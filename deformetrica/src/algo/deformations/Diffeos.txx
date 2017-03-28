/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Diffeos_txx
#define _Diffeos_txx

#include "Diffeos.h"

#include "GridFunctions.h"

#include "Landmark.h"
#include "LinearInterpImage.h"

#include "KernelFactory.h"

#include "AdjointEquationsIntegrator.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include <cassert>

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::Diffeos() : Superclass(), m_T0(0.0), m_TN(1.0), m_NumberOfTimePoints(10), m_KernelType(null), m_KernelWidth(1.0),
m_UseImprovedEuler(true), m_PaddingFactor(0.0), m_OutOfBox(true), m_ComputeTrueInverseFlow(true)
{
	this->SetDiffeosType();
}



template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::Diffeos(const Diffeos& other) : Superclass(other)
{
	this->SetDiffeosType();

	m_KernelType = other.m_KernelType;
	m_KernelWidth = other.m_KernelWidth;

	m_T0 = other.m_T0;
	m_TN = other.m_TN;
	m_NumberOfTimePoints = other.m_NumberOfTimePoints;

	m_StartPositions = other.m_StartPositions;
	m_StartMomentas = other.m_StartMomentas;

	m_DataDomain = other.m_DataDomain;
	m_PaddingFactor = other.m_PaddingFactor;
	m_OutOfBox = other.m_OutOfBox;

	/************************************ BEGIN CHANGE **********************************************/
	// m_CPDomain = other.m_CPDomain;
	/************************************ END CHANGE  **********************************************/

	m_ComputeTrueInverseFlow = other.m_ComputeTrueInverseFlow;
	m_UseImprovedEuler = other.m_UseImprovedEuler;
}


template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::~Diffeos() {}


// template <class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::CopyInformation(const Diffeos& other)
// {
// 	
// 	Superclass::CopyInformation(other);
// 	m_KernelType = other.m_KernelType;
// 	m_KernelWidth = other.m_KernelWidth;
// 
// 	m_T0 = other.m_T0;
// 	m_TN = other.m_TN;
// 	m_NumberOfTimePoints = other.m_NumberOfTimePoints;
// 	
// 
// 	m_DataDomain = other.m_DataDomain;
// 	m_PaddingFactor = other.m_PaddingFactor;
// 	m_OutOfBox = other.m_OutOfBox;
// 	
// 	m_ComputeTrueInverseFlow = other.m_ComputeTrueInverseFlow;
// 	m_UseImprovedEuler = other.m_UseImprovedEuler;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

// template<class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::ReverseFlow()
//  {
// 	int midPoint = static_cast<int>( 0.5*(this->GetNumberOfTimePoints() + 1) );
// 
// 	for (int t = 0; t < midPoint; t++)
// 	{
// 		MatrixType auxP = m_PositionsT[t];
// 		MatrixType auxM = - m_MomentasT[t];
// 
// 		m_PositionsT[t] = m_PositionsT[this->GetNumberOfTimePoints() - t - 1];
// 		m_MomentasT[t] = - m_MomentasT[this->GetNumberOfTimePoints() - t- 1];
// 
// 		m_PositionsT[this->GetNumberOfTimePoints() - t- 1] = auxP;
// 		m_MomentasT[this->GetNumberOfTimePoints() - t - 1] = auxM;
// 	}
// 
//  }


template<class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::Update()
{
	if (this->IsModified())
	{
		this->InitBoundingBox();
		this->Shoot();
		this->UnsetModified();
		Superclass::m_DeformableObjectModified = true;
	}
	
	if(Superclass::m_DeformableObjectModified)
	{
		if ( Superclass::m_IsLandmarkPoints )
			this->FlowLandmarkPointsTrajectory();
		if ( Superclass::m_IsImagePoints )
			this->FlowImagePointsTrajectory();
	}
	
	Superclass::m_DeformableObjectModified = false;
}


template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::InitBoundingBox()
{
	TScalar offset = m_PaddingFactor * m_KernelWidth /2;

	MatrixType boundingBox = m_DataDomain;
	boundingBox.set_column(0, boundingBox.get_column(0) - offset);
	boundingBox.set_column(1, boundingBox.get_column(1) + offset);

	m_BoundingBox = boundingBox;
	m_OutOfBox = false;
}



template <class TScalar, unsigned int Dimension>
bool
Diffeos<TScalar, Dimension>
::CheckBoundingBox(MatrixList& X, int t)
{
	MatrixType Xt = X[t];
	int outOfBox = 0;

	for (int d=0; d<Dimension; d++)
	{
		outOfBox += (Xt.get_column(d).min_value() < m_BoundingBox(d,0));
		outOfBox += (Xt.get_column(d).max_value() > m_BoundingBox(d,1));
	}

	return (outOfBox > 0);
}


template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>*
Diffeos<TScalar, Dimension>
::GetDeformedObjectAt(unsigned int t)
{
	if (t<0 || t>this->GetNumberOfTimePoints()-1)
		throw std::runtime_error("In Diffeos::GetDeformedObjectAt - Time t out of range");

	if (this->IsModified() || Superclass::m_DeformableObjectModified)
	{
		std::cerr << "In Diffeos::GetDeformedObjectAt() - The Diffeos was not updated or the objects not deformed" << std::endl;
		this->Update();
	}
		
	MatrixType LP, IP;
	if (Superclass::m_DeformableMultiObject->GetNumberOfLandmarkKindObjects())
		LP = m_LandmarkPointsT[t];

	if (Superclass::m_DeformableMultiObject->GetNumberOfImageKindObjects())
		IP = (m_ComputeTrueInverseFlow?m_InverseMapsT[t]:m_MapsT[m_NumberOfTimePoints - t - 1]);
	
	DeformableMultiObjectType* deformedObjects = Superclass::m_DeformableMultiObject->DeformedMultiObject(LP, IP);
	
	
	// DeformableMultiObjectType* deformedObjects = Superclass::m_DeformableMultiObject->Clone();
	// // std::cout << Superclass::m_DeformableMultiObject->GetNumberOfImagePoints() << "   " << deformedObjects->GetNumberOfImagePoints() << std::endl;
	// if (deformedObjects->GetNumberOfLandmarkKindObjects())
	// {
	// 	deformedObjects->SetLandmarkPoints(m_LandmarkPointsT[t]);
	// }
	// if (deformedObjects->GetNumberOfImageKindObjects())
	// {
	// 	if (m_ComputeTrueInverseFlow)
	// 		deformedObjects->SetImagePoints(m_InverseMapsT[t]);
	// 	else
	// 		deformedObjects->SetImagePoints(m_MapsT[m_NumberOfTimePoints - t - 1]);
	// }
	// deformedObjects->Update();
	
	// int numberOfObjects = Superclass::m_Template->GetTemplateObjects()->GetNumberOfObjects();
	// std::vector<int> numberOfPoints = Superclass::m_Template->GetTemplateObjects()->GetNumberOfPoints();
	// 
	// DeformableMultiObjectType* deformedObjects = Superclass::m_Template->GetTemplateObjects()->Clone();
	// 
	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// ImageTypePointer deformedImage;
	// 
	// MatrixList Y(numberOfObjects);
	// int counterLandmark = 0;
	// for (int i = 0; i < numberOfObjects; i++)
	// {
	// 	if ( Superclass::m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
	// 	{
	// 		Y[i] = m_LandmarkPointsT[t].extract(numberOfPoints[i], Dimension, counterLandmark);
	// 		counterLandmark += numberOfPoints[i];
	// 	}
	// 	else if ( Superclass::m_Template->GetObjectList()[i]->IsOfImageKind() )
	// 	{
	// 		// This is the flow
	// 		if (m_ComputeTrueInverseFlow)
	// 			Y[i] = m_InverseMapsT[t];
	// 			// BE CAREFUL: FOR t != numTimePoints-1 THIS IS NOT THE FLOW I_0\circ\phi_t^{-1}!!!! BUT I_0 \circ phi_t\circ\phi_1^{-1} instead
	// 			// the flow is backward in this case: needs to flip time, so that the deformed image is at time t = numTimePoints-1
	// 		else
	// 			Y[i] = m_MapsT[m_NumberOfTimePoints - t - 1];
	// 	}
	// 	else
	// 		throw std::runtime_error("In Diffeos::GetDeformedObjectAt(t) - Unknown object type");
	// }
	// 
	// deformedObjects->SetImagePointsAndPointData(Y);
	
	return deformedObjects;
}




template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::WriteFlow(const std::vector<std::string>& name, const std::vector<std::string>& extension)
{

	
	
	for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
	{
		DeformableMultiObjectType* deformedObjects = this->GetDeformedObjectAt(t);
		
		std::vector<std::string> fn(name.size());
		for (int i = 0; i < name.size(); i++)
		{
			std::ostringstream oss;
			oss << name[i] << "__t_" << t << extension[i] << std::ends;
			fn[i] = oss.str();
		}

		deformedObjects->WriteMultiObject(fn, m_LandmarkPointsVelocity);
	}
	
	
	// this->Update();
	// 
	// DeformableMultiObjectType* deformedObjects = Superclass::m_DeformableMultiObject->Clone();
	// 
	// for (int t = 0; t < m_NumberOfTimePoints; t++)
	// {
	// 	// remove to avoid to many copies of the template objects
	// 	// DeformableMultiObjectType* deformedObjects = this->GetDeformedObjectAt(t);
	// 	
	// 	if ( Superclass::m_IsLandmarkPoints )
	// 	{
	// 		deformedObjects->SetLandmarkPoints(m_LandmarkPointsT[t]);
	// 	}
	// 	if ( Superclass::m_IsImagePoints )
	// 	{
	// 		if (m_ComputeTrueInverseFlow)
	// 			deformedObjects->SetImagePoints(m_InverseMapsT[t]);
	// 		else
	// 			deformedObjects->SetImagePoints(m_MapsT[m_NumberOfTimePoints - t - 1]);
	// 	}
	// 	deformedObjects->Update(); // this is needed for images to update image intensities with new image point positions. Not really needed for landmark kind objects, though.
	// 	
	// 	for (int i = 0; i < deformedObjects->GetNumberOfObjects(); i++)
	// 	{
	// 		std::ostringstream oss;
	// 		oss << name[i] << "__t_" << t << extension[i] << std::ends;
	// 
	// 		deformedObjects->GetObjectList()[i]->WriteObject(oss.str());
	// 	}
	// }
}




template <class TScalar, unsigned int Dimension>
typename Diffeos<TScalar, Dimension>::MatrixType
Diffeos<TScalar, Dimension>
::SplatResidualImage(const DeformableMultiObjectType* target)
{
	this->Update();
	return Superclass::m_DeformableMultiObject->SplatDifferenceImage(target, m_ComputeTrueInverseFlow?m_InverseMapsT[m_NumberOfTimePoints-1]:m_MapsT[0]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::Shoot()
{
	unsigned int numCP = m_StartPositions.rows(); // WHY not creating a variable numberCP, they do not change and we could set it as const

	MatrixList& outPos = m_PositionsT;
	MatrixList& outMoms = m_MomentasT;

	outPos.resize(m_NumberOfTimePoints);
	outMoms.resize(m_NumberOfTimePoints);
	for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
	{
		outPos[t] = m_StartPositions;
		outMoms[t].set_size(numCP, Dimension);
		outMoms[t].fill(0);
	}
	outPos[0] = m_StartPositions;
	outMoms[0] = m_StartMomentas;

	// Special case: nearly zero momentas yield no motion
	if (outMoms[0].frobenius_norm() < 1e-20)
		return;

	//TScalar timeStep = 1.0 / (m_NumberOfTimePoints-1); // Why not T0 and Tm
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
	kernelObj->SetKernelWidth(this->GetKernelWidth());

	for (unsigned int t = 0; t < (m_NumberOfTimePoints-1); t++)
	{
		kernelObj->SetSources(outPos[t]);
		kernelObj->SetWeights(outMoms[t]);

		std::vector<typename KernelType::MatrixType> kGradMom =	kernelObj->ConvolveGradient(outPos[t]);

		typename KernelType::MatrixType dPos = kernelObj->Convolve(outPos[t]);

		typename KernelType::MatrixType dMom(numCP, Dimension, 0);
		for (unsigned int i = 0; i < numCP; i++)
			dMom.set_row(i,	(kGradMom[i].transpose() * outMoms[t].get_row(i)) );

		outPos[t+1] = outPos[t] + dPos * dt;
		outMoms[t+1] = outMoms[t] - dMom * dt;

		//		// Heun's method
		//		if (m_UseImprovedEuler)
		//		{
		//			kernelObj->SetSources(outPos[t+1]);
		//			kernelObj->SetWeights(outMoms[t+1]);
		//
		//			kGradMom = kernelObj->ConvolveGradient(outPos[t+1]);
		//
		//			typename KernelType::MatrixType dPos2 =
		//					kernelObj->Convolve(outPos[t+1]);
		//
		//			typename KernelType::MatrixType dMom2(numCP, Dimension, 0);
		//			for (unsigned int i = 0; i < numCP; i++)
		//				dMom2.set_row(i,
		//						(kGradMom[i].transpose() * outMoms[t+1].get_row(i)) );
		//
		//			outPos[t+1] = outPos[t] + (dPos + dPos2) * (timeStep * 0.5);
		//			outMoms[t+1] = outMoms[t] - (dMom + dMom2) * (timeStep * 0.5);
		//		}

		if (this->CheckBoundingBox(outPos,t+1))
		{
			std::cout << "Geodesic shooting: out of box at time t = " << t+1 << std::endl;
			break;
		}
	}

	delete kernelObj;
}


template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::FlowLandmarkPointsTrajectory()
{

	//TScalar dt = 1.0 / (m_NumberOfTimePoints - 1);
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
	kernelObj->SetKernelWidth(this->GetKernelWidth());

	MatrixList posT = this->GetTrajectoryPositions();
	MatrixList momT = this->GetTrajectoryMomentas();

	m_LandmarkPointsT.resize(m_NumberOfTimePoints);
	m_LandmarkPointsVelocity.resize(m_NumberOfTimePoints);
	for (unsigned int t = 0; t < m_NumberOfTimePoints ; t++)
	{
		m_LandmarkPointsT[t] = Superclass::m_LandmarkPoints;
		m_LandmarkPointsVelocity[t].fill(0.0);
	}

	if (momT[0].frobenius_norm() < 1e-20)
		return;

	for (unsigned int t = 0; t < m_NumberOfTimePoints-1 ; t++)
	{
		kernelObj->SetSources(posT[t]);
		kernelObj->SetWeights(momT[t]);

		m_LandmarkPointsVelocity[t] = kernelObj->Convolve(m_LandmarkPointsT[t]);

		m_LandmarkPointsT[t + 1] = m_LandmarkPointsT[t] + (m_LandmarkPointsVelocity[t] * dt);

		if (this->ImprovedEuler())
		{
			kernelObj->SetSources(posT[t + 1]);
			kernelObj->SetWeights(momT[t + 1]);

			m_LandmarkPointsVelocity[t + 1] = kernelObj->Convolve(m_LandmarkPointsT[t + 1]);

			m_LandmarkPointsT[t + 1] = m_LandmarkPointsT[t] + (m_LandmarkPointsVelocity[t] + m_LandmarkPointsVelocity[t + 1]) * (dt * 0.5);
		}

		if (this->CheckBoundingBox(m_LandmarkPointsT,t+1))
		{
			std::cout << "Landmark deformation: out of box at time t = " << t+1 << std::endl;
			break;
		}
	}

	// std::cout << "Done updating Point trajectory" << std::endl;

	delete kernelObj;

}


template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::FlowImagePointsTrajectory()
{

	if (m_ComputeTrueInverseFlow)
		this->IntegrateImagePointsWithTrueInverseFlow();
	else
		this->IntegrateImagePointsBackward();
}




// computes \phi_t\circ\phi_1^{-1}, which is equal to \phi_1^{-1} for t = 0
template <class TScalar, unsigned int Dimension>
void Diffeos<TScalar, Dimension>
::IntegrateImagePointsBackward()
{
	// typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	// LIImageType* tempLII =	dynamic_cast<LIImageType*>(
	// 	Superclass::m_Template->GetTemplateObjects()->GetIthObject(Superclass::m_Template->GetTemplateObjects()->GetIndexImage()));
	// ImageTypePointer downSampledWorkingImage = tempLII->GetDownSampledWorkingImage();
	// 
	// if (downSampledWorkingImage.IsNull())
	// 	throw std::runtime_error("DownSampledWorkingImage should have been set in LinearInterpImage object");
	// 
	// MatrixType downSampledY1 = Superclass::m_Template->GetTemplateObjects()->GetImagePoints();

	m_MapsT.resize(m_NumberOfTimePoints);
	for (long t = 0; t < m_NumberOfTimePoints; t++)
		m_MapsT[t] = Superclass::m_ImagePoints;

		// m_MapsT[t] = downSampledY1;

	// Special case: nearly zero momentas yield no motion
	if (m_MomentasT[0].frobenius_norm() < 1e-20)
		return;

	// Compute y(t)
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
	kernelObj->SetKernelWidth(m_KernelWidth);

	for (long t = this->GetNumberOfTimePoints() - 1; t >= 1; t--)
	{
		kernelObj->SetSources(m_PositionsT[t]);
		kernelObj->SetWeights(m_MomentasT[t]);

		MatrixType dY = kernelObj->Convolve(m_MapsT[t]);

		m_MapsT[t - 1] = m_MapsT[t] - dY * dt;

		// Heun's method
		if (m_UseImprovedEuler)
		{
			kernelObj->SetSources(m_PositionsT[t - 1]);
			kernelObj->SetWeights(m_MomentasT[t - 1]);

			MatrixType dY2 = kernelObj->Convolve(m_MapsT[t - 1]);

			m_MapsT[t - 1] = m_MapsT[t] - (dY + dY2) * (dt * 0.5);
		}

		if (this->CheckBoundingBox(m_MapsT, t - 1))
		{
			std::cout << "Image deformation: out of box at time t = " << t - 1 << std::endl;
			return;
		}
	}

	delete kernelObj;
}



// Computes the flow \phi_t^{-1}
template<class TScalar, unsigned int Dimension>
void Diffeos<TScalar, Dimension>
::IntegrateImagePointsWithTrueInverseFlow()
{
	// 
	// typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	// LIImageType* tempLII =	dynamic_cast<LIImageType*>(
	// 	Superclass::m_Template->GetTemplateObjects()->GetIthObject(Superclass::m_Template->GetTemplateObjects()->GetIndexImage()));
	// ImageTypePointer downSampledWorkingImage = tempLII->GetDownSampledWorkingImage();
	// 
	// // downSampledY1 is the template image
	// MatrixType downSampledY1 = Superclass::m_Template->GetTemplateObjects()->GetImagePoints();

	// ImageTypePointer downSampledWorkingImage = Superclass::m_Template->GetDownSampledWorkingImage();

	// Initialize all the time points to be the image at the first time point
	m_InverseMapsT.resize(m_NumberOfTimePoints);
	for (long t = 0; t < m_NumberOfTimePoints; t++)
		m_InverseMapsT[t] = Superclass::m_ImagePoints;
	
		// m_InverseMapsT[t] = downSampledY1;

	// Special case: nearly zero momentas yield no motion
	if (m_MomentasT[0].frobenius_norm() < 1e-20)
		return;

	// The time step
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

	// A list of Dimension images, each image stores a spatial coordinate representing the physical location of the image points
	std::vector<ImageTypePointer>  Yt;
	Yt.resize(Dimension);
	for (unsigned int d = 0; d < Dimension; d++)
	{
		// Create an image from each Dimension
		ImageTypePointer img = GridFunctionsType::VectorToImage(Superclass::m_DownSampledImage, Superclass::m_ImagePoints.get_column(d));
		// ImageTypePointer img = GridFunctionsType::VectorToImage(downSampledWorkingImage, downSampledY1.get_column(d));
		Yt[d] = img;
	}

	// The kernel is for computing v_t(y0)
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
	kernelObj->SetKernelWidth(this->GetKernelWidth());

	// Loop over time for the Euler scheme
	// This is a forward integration now ! (contrary to the backward integration for non true inverse flow)
	// JAMES: for (long t = numTimePoints-1; t>0; t--)
	// STANLEY:
	for (long t = 0; t < m_NumberOfTimePoints-1; t++)
	{
		kernelObj->SetSources(m_PositionsT[t]);
		kernelObj->SetWeights(m_MomentasT[t]);

		// The velocity is always v_t(y0)
		MatrixType VtY0 = kernelObj->Convolve(Superclass::m_ImagePoints);

		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
		std::vector<ImageTypePointer> gradImages;
		gradImages.resize(Dimension*Dimension);
		unsigned int indx = 0;
		// Compute the gradient in all directions
		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
		{
			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
			{
				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
				derivf->SetInput(Yt[dim1]);
				derivf->SetDirection(dim2);
				derivf->SetOrder(1);
				derivf->SetUseImageSpacingOn();
				derivf->Update();

				gradImages[indx] = derivf->GetOutput();
				indx++;
			}
		}

		// Now we have all the information we need to compute dY, but we need to loop over all the pixels
		// and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
		IteratorType it(Superclass::m_DownSampledImage, Superclass::m_DownSampledImage->GetLargestPossibleRegion());
		int i = 0;
		MatrixType dY(Superclass::m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension);

		// Loop over the grid to construct jacobian matrices and compute dY
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			typedef typename ImageType::IndexType ImageIndexType;
			ImageIndexType ind = it.GetIndex();
			MatrixType jacobian(Dimension, Dimension);

			unsigned int indx = 0;
			// Build the (DimensionxDimension) jacobian matrix
			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
			{
				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
				{
					ImageTypePointer curGradImage = gradImages[indx];
					jacobian(dim1, dim2) = curGradImage->GetPixel(ind);
					indx++;
				}
			}

			// Get the (Dimension) vector corresponding to this points v_t(y0)
			VectorType VtY0k = VtY0.get_row(i);

			// dY = d_y0 \phi_t * v_t(y0)
			dY.set_row(i, jacobian*VtY0k);
			i++;
		}

		// Update the inverse map by Euler scheme
		// JAMES
		// m_InverseMapsT[t-1] = m_InverseMapsT[t] - dY * dt;
		// STANLEY
		m_InverseMapsT[t+1] = m_InverseMapsT[t] - dY * dt;

		// Update Yt using the m_InverseMap matrix we just updated using Euler
		MatrixType inversem1 = m_InverseMapsT[t+1];
		for (unsigned int d = 0; d < Dimension; d++)
		{
			// Create an image from each Dimension
			ImageTypePointer img = GridFunctionsType::VectorToImage(Superclass::m_DownSampledImage, inversem1.get_column(d));
			Yt[d] = img;
		}

		if (this->CheckBoundingBox(m_InverseMapsT, t+1))
		{
			std::cout << "Image deformation: out of box at time t = " << t+1 << std::endl;
			return;
		}

	}

	delete kernelObj;
}


template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::IntegrateAdjointEquations(MatrixType& InitialConditionsLandmarkPoints, MatrixType& InitialConditionsImagePoints)
{
	MatrixList ListInitialConditionsLandmarkPoints;
	MatrixList ListInitialConditionsImagePoints;
	ListInitialConditionsLandmarkPoints.resize(1);
	ListInitialConditionsImagePoints.resize(1);

	ListInitialConditionsLandmarkPoints[0] = InitialConditionsLandmarkPoints;
	ListInitialConditionsImagePoints[0] = InitialConditionsImagePoints;

	std::vector<unsigned int> times;
	times.resize(0);

	this->IntegrateAdjointEquations(ListInitialConditionsLandmarkPoints, ListInitialConditionsImagePoints, times);
}


template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::IntegrateAdjointEquations(MatrixList& InitialConditionsLandmarkPoints, MatrixList& InitialConditionsImagePoints, std::vector<unsigned int> jumpTimes)
{
	// Upsample image maps, since initial condition of image objects are at full resolution
	// ImageTypePointer image = Superclass::m_Template->GetTemplateObjects()->GetImage();
	// The path of the image points over time
	MatrixList fullResMapsT(m_NumberOfTimePoints);
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

	if ( Superclass::m_IsImagePoints )
	{
		for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
			fullResMapsT[t] = GridFunctionsType::UpsampleImagePoints(Superclass::m_Image, Superclass::m_DownSampledImage, 
										 m_ComputeTrueInverseFlow?m_InverseMapsT[t]:m_MapsT[t]);	
	}	

	typedef AdjointEquationsIntegrator<TScalar, Dimension> AEIntegrator;
	AEIntegrator* integrator = new AEIntegrator();
	integrator->SetControlPointsTrajectory(m_PositionsT);
	integrator->SetMomentaTrajectory(m_MomentasT);
	if ( Superclass::m_IsLandmarkPoints )
	{
		integrator->SetLandmarkPointsTrajectory(m_LandmarkPointsT);
		integrator->SetInitialConditionsLandmarkPoints(InitialConditionsLandmarkPoints);
	}
	if ( Superclass::m_IsImagePoints )
	{
		integrator->SetImagePointsTrajectory(fullResMapsT);
		integrator->SetInitialConditionsImagePoints(InitialConditionsImagePoints);
	}
	integrator->SetKernelWidth(m_KernelWidth);
	integrator->SetKernelType(m_KernelType);
	integrator->SetNumberOfTimePoints(m_NumberOfTimePoints);
	integrator->SetFullResolutionImage(Superclass::m_Image);
	integrator->SetT0(m_T0);
	integrator->SetTN(m_TN);
	integrator->SetJumpTimes(jumpTimes);
	
	m_ComputeTrueInverseFlow?integrator->SetComputeTrueInverseFlow():integrator->UnsetComputeTrueInverseFlow();
	m_UseImprovedEuler?integrator->UseImprovedEuler():integrator->UseStandardEuler();
	
	integrator->Update();

	m_AdjointPosAt0 = integrator->GetAdjointPosAt(0);
	m_AdjointMomAt0 = integrator->GetAdjointMomAt(0);
	m_AdjointLandmarkPointsAt0 = integrator->GetAdjointLandmarkPointsAt(0);

	delete integrator;
}


// template <class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::TransportAlongGeodesic(MatrixList& InitialConditions, MatrixList& VectorsT, MatrixList& PointsT)
//  {
// 
// 	/*
// 	if (InitialConditions.size() != m_NumberOfObjects)
// 		throw std::runtime_error("Number of initial conditions should equal the number of objects");
// 
// 	VectorsT.resize(m_NumberOfTimePoints);
// 	PointsT.resize(m_NumberOfTimePoints);
// 
// 	int NumbPts = 0;
// 	for (int i = 0; i < m_NumberOfObjects; i++)
// 		NumbPts += m_NumberOfPoints[i];
// 
// 	MatrixType zeroM(NumbPts, Dimension, 0.0);
// 	for (int t = 0; t < m_NumberOfTimePoints; t++)
// 	{
// 		VectorsT[t] = zeroM;
// 		PointsT[t] = zeroM;
// 	}
// 
// 	int counter = 0;
// 	for (int i = 0; i < m_NumberOfObjects; i++)
// 	{
// 		MatrixList VecT;
// 		MatrixList PtsT;
// 		m_ObjectList[i]->TransportAlongGeodesic(InitialConditions[i], VecT, PtsT);
// 
// 		for (int t = 0; t < m_NumberOfTimePoints; t++)
// 			for (int r=0; r < m_NumberOfPoints[i]; r++)
// 			{
// 				PointsT[t].set_row(counter + r, PtsT[t].get_row(r));
// 				VectorsT[t].set_row(counter + r,VecT[t].get_row(r));
// 			}
// 
// 		counter += m_NumberOfPoints[i];
// 	}
// 	*/
// 
// 	int numberOfObjects = Superclass::m_Template->GetTemplateObjects()->GetNumberOfObjects();
// 	if (InitialConditions.size() != numberOfObjects)
// 		throw std::runtime_error(
// 				"In Diffeos::TransportAlongGeodesic() - Number of initial conditions should equal the number of objects");
// 
// 	int numberOfImagePoints = 0;
// 	int numberOfLandmarkPoints = 0;
// 	std::vector<int> numberOfPoints = Superclass::m_Template->GetTemplateObjects()->GetNumberOfPoints();
// 
// 	for (int i = 0; i < numberOfObjects; i++)
// 	{
// 		if ( Superclass::m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
// 			numberOfLandmarkPoints += numberOfPoints[i];
// 		else if ( Superclass::m_Template->GetObjectList()[i]->IsOfImageKind() )
// 			numberOfImagePoints += numberOfPoints[i];
// 		else
// 			throw std::runtime_error("In Diffeos::TransportAlongGeodesic() - Unknown object type");
// 	}
// 
// 	//
// 	// Split of InitialConditions into two matrices :
// 	//
// 	MatrixType initialConditionsImage(numberOfImagePoints, Dimension, 0.0);
// 	MatrixType initialConditionsLandmark(numberOfLandmarkPoints, Dimension, 0.0);
// 
// 	//std::cout << "initialConditionsImage" << initialConditionsImage << std::endl;
// 	//std::cout << "initialConditionsLandmark" << initialConditionsLandmark << std::endl;
// 
// 	int counterImage = 0;
// 	int counterLandmark = 0;
// 	for (int i = 0; i < numberOfObjects; i++)
// 	{
// 		if ( Superclass::m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
// 		{
// 			for (int r=0; r < InitialConditions[i].rows(); r++)
// 				initialConditionsLandmark.set_row(counterLandmark + r, InitialConditions[i].get_row(r));
// 
// 			counterLandmark += numberOfPoints[i];
// 		}
// 		else if ( Superclass::m_Template->GetObjectList()[i]->IsOfImageKind() )
// 		{
// 			for (int r=0; r < InitialConditions[i].rows(); r++)
// 				initialConditionsImage.set_row(counterImage + r, InitialConditions[i].get_row(r));
// 
// 			counterImage += numberOfPoints[i];
// 		}
// 	}
// 
// 	//
// 	// Transport along geodesic for image & landmark points :
// 	//
// 	MatrixList EtaT, YT, ThetaT, XT;
// 	if (m_ComputeTrueInverseFlow)
// 	{
// 		this->TransportAlongGeodesicBackwardForImagePoints(initialConditionsImage, EtaT, YT);
// 	}
// 	else
// 	{
// 		this->TransportAlongGeodesicForwardForImagePoints(initialConditionsImage, EtaT, YT);
// 	}
// 	TransportAlongGeodesicForLandmarkPoints(initialConditionsLandmark, ThetaT, XT);
// 
// 	//
// 	// Concatenation of the results in VectorsT & PointsT :
// 	//
// 	VectorsT.resize(m_NumberOfTimePoints);
// 	PointsT.resize(m_NumberOfTimePoints);
// 
// 	MatrixType zeroM(Superclass::m_Template->GetTemplateObjects()->GetTotalNumberOfPoints(), Dimension, 0.0);
// 	for (int t = 0; t < m_NumberOfTimePoints; t++)
// 	{
// 		VectorsT[t] = zeroM;
// 		PointsT[t] = zeroM;
// 	}
// 
// 	int counter = 0;
// 	counterImage = 0; counterLandmark = 0;
// 	for (int i = 0; i < numberOfObjects; i++)
// 	{
// 		if ( Superclass::m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
// 		{
// 			for (int t = 0; t < m_NumberOfTimePoints; t++)
// 				for (int r=0; r < numberOfPoints[i]; r++)
// 				{
// 					PointsT[t].set_row(counter + r, XT[t].get_row(counterLandmark + r));
// 					VectorsT[t].set_row(counter + r, ThetaT[t].get_row(counterLandmark + r));
// 				}
// 
// 			counterLandmark += numberOfPoints[i];
// 		}
// 		else if ( Superclass::m_Template->GetObjectList()[i]->IsOfImageKind() )
// 		{
// 			for (int t = 0; t < m_NumberOfTimePoints; t++)
// 				for (int r=0; r < numberOfPoints[i]; r++)
// 				{
// 					PointsT[t].set_row(counter + r, YT[t].get_row(counterImage + r));
// 					VectorsT[t].set_row(counter + r, EtaT[t].get_row(counterImage + r));
// 				}
// 
// 			counterImage += numberOfPoints[i];
// 		}
// 		counter += numberOfPoints[i];
// 	}
// 
// 	//
// 	// For debug (security) :
// 	//
// 	if (counter != Superclass::m_Template->GetTemplateObjects()->GetTotalNumberOfPoints())
// 		throw std::runtime_error("In Diffeos::TransportAlongGeodesic() - The number of points after concatenation is wrong.");
// 
//  }




// // This computes \dot{\eta}(t) = -\partial_1 G(t) \eta(t)
// template<class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::TransportAlongGeodesicBackwardForImagePoints(MatrixType& Eta1, MatrixList& EtaT, MatrixList& YT)
// {
// 	if (Superclass::m_Template->GetTemplateObjects()->GetNumberOfImageKindObjects() == 0)
// 		return;
// 
// 	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
// 	LIImageType* tempLII = dynamic_cast<LIImageType*>(
// 		Superclass::m_Template->GetTemplateObjects()->GetIthObject(
// 			Superclass::m_Template->GetTemplateObjects()->GetIndexImage()));
// 	ImageTypePointer image = tempLII->GetImage();
// 
// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 
// 	// The path of the image points over time
// 	YT.resize(m_NumberOfTimePoints);
// 	for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
// 	{
// 		//		YT[t] = GridFunctionsType::UpsampleImagePoints(tempLII->GetImage(), tempLII->GetDownSampledWorkingImage(), m_InverseMapsT[t]);
// 		YT[t] = tempLII->UpSampleImageMap(m_InverseMapsT[t]);
// 	}
// 
// 	MatrixType YT0 = YT[0];
// 
// 	// long numVoxels = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
// 	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);
// 
// 	// The aux variable Eta over time
// 	EtaT.resize(m_NumberOfTimePoints);
// 	// The source term for integration. There is a + sign here contrary to in TransportAlongGeodesicForward
// 	EtaT[m_NumberOfTimePoints-1] = Eta1;
// 
// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 
// 	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
// 	typedef typename KernelFactoryType::KernelBaseType KernelType;
// 	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();
// 
// 	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
// 	kernelObj->SetKernelWidth(this->GetKernelWidth());
// 
// 	// Integrate backwards in time
// 	for (long t = (m_NumberOfTimePoints - 1); t>=1; t--)
// 	{
// 		kernelObj->SetSources(m_PositionsT[t]);
// 		kernelObj->SetWeights(m_MomentasT[t]);
// 
// 		// The velocity is always v_t(y0)
// 		// MatrixType VtY0 = kernelObj->Convolve(m_DownSampledY1);
// 		MatrixType VtY0 = kernelObj->Convolve(YT0);
// 
// 		// This gives us a Dim*Dim matrix at every pixel of the original image
// 		// std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(m_DownSampledY1);
// 		std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(YT0);
// 
// 
// 		// We need to have eta(t) in the form of an image, so we can compute the jacobian
// 		std::vector<ImageTypePointer> etaTImage;
// 		etaTImage.resize(Dimension);
// 		for (unsigned int d = 0; d < Dimension; d++)
// 		{
// 			// Create an image from each Dimension
// 			ImageTypePointer img = GridFunctionsType::VectorToImage(image, EtaT[t].get_column(d));
// 			etaTImage[d] = img;
// 		}
// 
// 		// Compute jacobian of eta(t)
// 		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
// 		std::vector<ImageTypePointer> gradImages;
// 		gradImages.resize(Dimension*Dimension);
// 		unsigned int indx = 0;
// 		// Compute the gradient in all directions
// 		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 		{
// 			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 			{
// 				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
// 				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
// 				derivf->SetInput(etaTImage[dim1]);
// 				derivf->SetDirection(dim2);
// 				derivf->SetOrder(1);
// 				derivf->SetUseImageSpacingOn();
// 				derivf->Update();
// 
// 				gradImages[indx] = derivf->GetOutput();
// 				indx++;
// 			}
// 		}
// 
// 		// Now we have all the information we need to compute dEta, but we need to loop over all the pixels
// 		// and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
// 		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
// 		// IteratorType it(m_DownSampledReferenceImage, m_DownSampledReferenceImage->GetLargestPossibleRegion());
// 		// int matrixIndex = 0;
// 		// int numVoxelsDS = m_DownSampledReferenceImage->GetLargestPossibleRegion().GetNumberOfPixels();
// 		// MatrixType dEta(numVoxelsDS, Dimension, 0);
// 		IteratorType it(image, image->GetLargestPossibleRegion());
// 		int matrixIndex = 0;
// 		int numVoxels = image->GetLargestPossibleRegion().GetNumberOfPixels();
// 		MatrixType dEta(numVoxels, Dimension, 0);
// 
// 		// Loop over the grid to construct jacobian matrices and compute dY
// 		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
// 		{
// 			typedef typename ImageType::IndexType ImageIndexType;
// 			ImageIndexType imageIndex = it.GetIndex();
// 			MatrixType jacobian(Dimension, Dimension);
// 
// 			unsigned int tempMatrixIndex = 0;
// 			// Build the (DimensionxDimension) jacobian matrix
// 			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 			{
// 				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 				{
// 					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
// 					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
// 					tempMatrixIndex++;
// 				}
// 			}
// 
// 			// The trace of the first term
// 			TScalar traceValue = trace(firstTerm[matrixIndex]);
// 
// 			// Get the (Dimension) vector corresponding to this points v_t(y0)
// 			VectorType VtY0k = VtY0.get_row(matrixIndex);
// 
// 			// The final value for this pixel of dEta
// 			// JAMES
// 			// dEta.set_row(matrixIndex, traceValue*EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);
// 			// STANLEY
// 			dEta.set_row(matrixIndex, -traceValue*EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);
// 
// 
// 			matrixIndex++;
// 		}
// 
// 		// JAMES
// 		// EtaT[t-1] = EtaT[t] + dEta * dt;
// 		// STANLEY (backward integration: the differential changes sign)
// 		// also need to upsample image map
// 		// MatrixType dEta_HD = this->UpSampleImageMap(dEta);
// 		// EtaT[t-1] = EtaT[t] - dEta_HD * dt;
// 		EtaT[t-1] = EtaT[t] - dEta * dt;
// 
// 
// 	}
// 
// 	// We have computed Eta(t), the backwards propagator actually needs -(d_yp Y)^t Eta(t)
// 	for (long t = 0; t<m_NumberOfTimePoints; t++)
// 	{
// 		// We need to have Eta(t) in the form of an image, so we can compute the jacobian
// 		std::vector<ImageTypePointer> YTImage;
// 		YTImage.resize(Dimension);
// 		for (unsigned int d = 0; d < Dimension; d++)
// 		{
// 			// Create an image from each Dimension
// 			ImageTypePointer img = GridFunctionsType::VectorToImage(image, YT[t].get_column(d));
// 			YTImage[d] = img;
// 		}
// 
// 		// Compute jacobian of eta(t)
// 		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
// 		std::vector<ImageTypePointer> gradImages;
// 		gradImages.resize(Dimension*Dimension);
// 		unsigned int indx = 0;
// 		// Compute the gradient in all directions
// 		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 		{
// 			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 			{
// 				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
// 				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
// 				derivf->SetInput(YTImage[dim1]);
// 				derivf->SetDirection(dim2);
// 				derivf->SetOrder(1);
// 				derivf->SetUseImageSpacingOn();
// 				derivf->Update();
// 
// 				gradImages[indx] = derivf->GetOutput();
// 				indx++;
// 			}
// 		}
// 
// 		// Now we need to loop over all the pixels and construct the 3x3 jacobian matrix
// 		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
// 		// JAMES
// 		// IteratorType it(m_DownSampledReferenceImage, m_DownSampledReferenceImage->GetLargestPossibleRegion());
// 		// STANLEY: YT is high res image
// 		IteratorType it(image, image->GetLargestPossibleRegion());
// 		int matrixIndex = 0;
// 
// 		// Loop over the grid to construct jacobian matrices and compute dY
// 		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
// 		{
// 			typedef typename ImageType::IndexType ImageIndexType;
// 			ImageIndexType imageIndex = it.GetIndex();
// 			MatrixType jacobian(Dimension, Dimension);
// 
// 			unsigned int tempMatrixIndex = 0;
// 			// Build the (DimensionxDimension) jacobian matrix
// 			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 			{
// 				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 				{
// 					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
// 					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
// 					tempMatrixIndex++;
// 				}
// 			}
// 
// 			// The final value for this pixel of dEta
// 			//JAMES
// 			// EtaT[t].set_row(matrixIndex, jacobian.transpose()*EtaT[t].get_row(matrixIndex));
// 			//STANLEY
// 			EtaT[t].set_row(matrixIndex, -jacobian.transpose()*EtaT[t].get_row(matrixIndex));
// 			matrixIndex++;
// 		}
// 	}
// 
// 	// STANLEY
// 	// the backward propagator needs YT = Y0 for all t
// 	for (int t = 1; t < m_NumberOfTimePoints; t++)
// 	{
// 		YT[t] = YT[0];
// 	}
// 
// }



// template<class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::TransportAlongGeodesicForwardForImagePoints(MatrixType& Eta0, MatrixList& EtaT, MatrixList& YT)
// {
// 	if (Superclass::m_Template->GetTemplateObjects()->GetNumberOfImageKindObjects() == 0)
// 		return;
// 
// 	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
// 	LIImageType* tempLII = dynamic_cast<LIImageType*>(
// 		Superclass::m_Template->GetTemplateObjects()->GetIthObject(
// 			Superclass::m_Template->GetTemplateObjects()->GetIndexImage()));
// 	ImageTypePointer image = tempLII->GetImage();
// 
// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 
// 	YT.resize(m_NumberOfTimePoints);
// 	for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
// 	{
// 		//		YT[t] = GridFunctionsType::UpsampleImagePoints(tempLII->GetImage(), tempLII->GetDownSampledWorkingImage(), m_MapsT[t]);
// 		YT[t] = tempLII->UpSampleImageMap(m_MapsT[t]);
// 	}
// 
// 	long numVoxels = image->GetLargestPossibleRegion().GetNumberOfPixels();
// 
// 	// Propagate eta forward
// 	EtaT.resize(m_NumberOfTimePoints);
// 	// STANLEY: PREVIOUS
// 	// EtaT[0] = Eta0;
// 	// STANLEY: NOW (the sign is taken into account here)
// 	EtaT[0] = -Eta0;
// 
// 
// 	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);
// 
// 	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
// 	typedef typename KernelFactoryType::KernelBaseType KernelType;
// 	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();
// 
// 	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
// 	kernelObj->SetKernelWidth(this->GetKernelWidth());
// 
// 	for (long t = 0; t < (m_NumberOfTimePoints - 1); t++)
// 	{
// 		kernelObj->SetSources(m_PositionsT[t]);
// 		kernelObj->SetWeights(m_MomentasT[t]);
// 
// 		/*
// 		// Grad mom splatted at CP locations, with convolutions evaluated at y(t-1)
// 		// This computes alpha_p \nabla_1 K(y_k, c_p)
// 		std::vector<MatrixType> gradMom = kernelObj->ConvolveGradient(YT[t]);
// 
// 		// This computes \eta_k^t alpha_p /nabla_1 K(y_k, c_p)
// 		MatrixType dEta(numVoxels, Dimension, 0);
// 		for (unsigned int j = 0; j < numVoxels; j++)
// 		{
// 		dEta.set_row(j, gradMom[j].transpose() * EtaT[t].get_row(j));
// 		}
// 		*/
// 
// 		// Grad mom splatted at CP locations, with convolutions evaluated at y(t-1)
// 		// This computes \eta_k^t alpha_p /nabla_1 K(y_k, c_p)
// 		MatrixType dEta = kernelObj->ConvolveGradient(YT[t], EtaT[t]);
// 		assert(dEta.rows() == numVoxels);
// 		assert(dEta.columns() == Dimension);
// 
// 		EtaT[t + 1] = EtaT[t] - dEta * dt;
// 
// 		// Heun's method
// 		if (this->ImprovedEuler())
// 		{
// 			kernelObj->SetSources(m_PositionsT[t + 1]);
// 			kernelObj->SetWeights(m_MomentasT[t + 1]);
// 
// 			/*
// 			gradMom = kernelObj->ConvolveGradient(YT[t + 1]);
// 
// 			MatrixType dEta2(numVoxels, Dimension, 0);
// 			for (unsigned int j = 0; j < numVoxels; j++)
// 			{
// 			dEta2.set_row(j, gradMom[j].transpose() * EtaT[t + 1].get_row(j));
// 			}
// 			*/
// 			MatrixType dEta2 = kernelObj->ConvolveGradient(YT[t + 1], EtaT[t + 1]);
// 			assert(dEta.rows() == numVoxels);
// 			assert(dEta.columns() == Dimension);
// 
// 			EtaT[t + 1] = EtaT[t] - (dEta + dEta2) * (dt * 0.5);
// 		}
// 
// 	}
// 
// 	// std::cout << "EtaT[end] = " << std::endl;
// 	// std::cout << EtaT[numTimePoints-1].get_n_rows(0,20) << std::endl;
// 
// 
// 	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 	// ImageTypePointer testT = GridFunctionsType::VectorToImage(m_Image, YT[19].get_column(0));
// 	// typedef itk::ImageFileWriter<ImageType> WriterType;
// 	//   typename WriterType::Pointer writer = WriterType::New();
// 	//   writer->SetInput(testT);
// 	//   writer->SetFileName("YT.vtk");
// 	//   writer->Update();
// 
// 	delete kernelObj;
// }


// // UPdate of theta!
// template <class TScalar, unsigned int Dimension>
// void
// Diffeos<TScalar, Dimension>
// ::TransportAlongGeodesicForLandmarkPoints(MatrixType& Theta1, MatrixList& ThetaT, MatrixList& XT)
// {
// 	if (Superclass::m_Template->GetTemplateObjects()->GetNumberOfLandmarkKindObjects() == 0)
// 		return;
// 
// 	// if (this->isModified())
// 	// 	this->UpdatePointTrajectory();
// 
// 	XT = m_LandmarkPointsT;
// 	int numberOfPoints = m_LandmarkPointsT[0].rows();
// 
// 	// Propagate theta backward
// 	ThetaT.resize(m_NumberOfTimePoints); // <--- Initialize each matrix at 0 ?
// 	ThetaT[m_NumberOfTimePoints-1] = Theta1;
// 
// 	//TScalar dt = 1.0 / (m_NumberOfTimePoints - 1);
// 	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);
// 
// 	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
// 	typedef typename KernelFactoryType::KernelBaseType KernelType;
// 	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();
// 
// 	KernelType* kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
// 	kernelObj->SetKernelWidth(this->GetKernelWidth());
// 
// 	for (long t = m_NumberOfTimePoints-1; t > 0; t--)
// 	{
// 		kernelObj->SetSources(m_PositionsT[t]);
// 		kernelObj->SetWeights(m_MomentasT[t]);
// 
// 		// Grad mom splatted at CP locations, with convolutions evaluated at x(t-1)
// 		// SimpleTimer timer;
// 
// 		/*
// 		std::vector<MatrixType> gradMom = kernelObj->ConvolveGradient(m_LandmarkPointsT[t]);
// 
// 		MatrixType dTheta(numberOfPoints, Dimension, 0);
// 		for (unsigned int j = 0; j < numberOfPoints; j++)
// 		{
// 		dTheta.set_row(j, gradMom[j].transpose() * ThetaT[t].get_row(j) );
// 		}
// 		*/
// 		MatrixType dTheta = kernelObj->ConvolveGradient(m_LandmarkPointsT[t], ThetaT[t]);
// 		assert(dTheta.rows() == this->GetNumberOfTimePoints());
// 		assert(dTheta.columns() == Dimension);
// 
// 		// std::cout << "diff term: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 
// 
// 		ThetaT[t-1] = ThetaT[t] + dTheta * dt; // the plus is correct! dTheta should be negative.
// 
// 		// Heun's method
// 		if (m_UseImprovedEuler)
// 		{
// 			kernelObj->SetSources(m_PositionsT[t-1]);
// 			kernelObj->SetWeights(m_MomentasT[t-1]);
// 
// 			/*
// 			gradMom = kernelObj->ConvolveGradient(m_LandmarkPointsT[t-1]);
// 			MatrixType dTheta2(numberOfPoints, Dimension, 0);
// 			for (unsigned int j = 0; j < numberOfPoints; j++)
// 			{
// 			dTheta2.set_row(j, gradMom[j].transpose() * ThetaT[t-1].get_row(j) );
// 			}
// 			*/
// 			MatrixType dTheta2 = kernelObj->ConvolveGradient(m_LandmarkPointsT[t-1], ThetaT[t-1]);
// 			assert(dTheta2.rows() == this->GetNumberOfTimePoints());
// 			assert(dTheta2.columns() == Dimension);
// 
// 			ThetaT[t-1] = ThetaT[t] + (dTheta + dTheta2) * (dt * 0.5);
// 		}
// 
// 	}
// 
// 	delete kernelObj;
// }



#endif /* _Diffeos_txx */

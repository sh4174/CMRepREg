/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SSDImage_txx
#define _SSDImage_txx

#include "SSDImage.h"

#include "KernelFactory.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "GridFunctions.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
SSDImage<TScalar, Dimension>
::SSDImage() : Superclass()
{
	this->SetSSDImageType();
}



template <class TScalar, unsigned int Dimension>
SSDImage<TScalar, Dimension>
::SSDImage(const SSDImage& o) : Superclass(o)
{

}

template <class TScalar, unsigned int Dimension>
SSDImage<TScalar, Dimension>
::SSDImage(const SSDImage& ex, const MatrixType& IP) : Superclass(ex, IP)
{
	this->SetSSDImageType();
	this->Update();
}


template <class TScalar, unsigned int Dimension>
SSDImage<TScalar, Dimension>
::~SSDImage()
{
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class TScalar, unsigned int Dimension>
void
SSDImage<TScalar, Dimension>
::Update()
{
	if (this->IsModified())
	{
		Superclass::Update();
	}
	
	this->UnSetModified();
}

template<class TScalar, unsigned int Dimension>
TScalar
SSDImage<TScalar, Dimension>
::ComputeMatch(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	const SSDImage* targ = dynamic_cast<const SSDImage*>(target);
	
	this->Update();
	// target->Update();

	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	VectorType I0 = GridFunctionsType::VectorizeImage(this->GetImage());
	VectorType I1 = GridFunctionsType::VectorizeImage(targ->GetImage());

	// SSD norm between images  
	if (I0.size() != I1.size())
		throw std::runtime_error("image sizes mismatch");

	VectorType D = I0 - I1;
	TScalar match = D.squared_magnitude();

	return match;
}



// compute 2 (I_0(y(0)) - I_1) nabla_{y(0)}I_0
template<class TScalar, unsigned int Dimension>
typename SSDImage<TScalar, Dimension>::MatrixType
SSDImage<TScalar, Dimension>
::ComputeMatchGradient(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		std::cerr << "Deformable objects types mismatched: " << this->GetType() << " and " << target->GetType() << "\n";


	const SSDImage* targ = dynamic_cast<const SSDImage*>(target);

	this->Update();
	// target->Update();
	
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

//	MatrixType Yfinal = GridFunctionsType::UpsampleImagePoints(this->GetImage(), this->GetDownSampledWorkingImage(), Superclass::m_DownSampledY1);
	
	VectorType I0 = GridFunctionsType::VectorizeImage(this->GetImage());
	VectorType I1 = GridFunctionsType::VectorizeImage(targ->GetImage());
	VectorType D = I0 - I1;
	
	MatrixType gradMatch(Superclass::m_NumberOfVoxels, Dimension);
	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		VectorType gradI_d = GridFunctionsType::VectorizeImage(Superclass::m_GradientImages[dim]);
		for (unsigned int i = 0; i < Superclass::m_NumberOfVoxels; i++)
			gradI_d[i] *= D[i];
		
		gradI_d *= Superclass::m_FlipAxes[dim];
		gradMatch.set_column(Superclass::m_PermutationAxes[dim], gradI_d);
	}


	// VectorType m_image_vec = GridFunctionsType::VectorizeImage(m_Image);
	// std::cout << "m_Image = " << m_image_vec.max_value() << " " << m_image_vec.min_value() << std::endl;
	// 
	// std::cout << "Idef = " << Idef.max_value() << " " << Idef.min_value() << std::endl;
	// std::cout << "I1 = " << I1.max_value() << " " << I1.min_value() << std::endl;
	
		
	// ImageTypePointer DifferenceImage = GridFunctionsType::VectorToImage(m_Image, D);
	// typedef itk::ImageFileWriter<ImageType> WriterType;
	// typename WriterType::Pointer writer = WriterType::New();
	// writer->SetInput(DifferenceImage);
	// writer->SetFileName("ImageDifference.img");
	// writer->Update();

	// int numVoxels = Yfinal.rows();
	// int numVoxels = this->GetNumberOfPoints();
	// 
	// // Deform gradient image and interpolate
	// MatrixType Yup = this->UpSampleImageMap(Superclass::m_DownSampledY1);
	// MatrixType gradMatch(numVoxels, Dimension);
	// for (unsigned int dim = 0; dim < Dimension; dim++)
	// {
	// 	typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
	// 	typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
	// 	derivf->SetInput(this->GetReferenceImage());
	// 	derivf->SetDirection(dim);
	// 	derivf->SetOrder(1);
	// 	derivf->SetUseImageSpacingOn();
	// 	derivf->Update();
	// 	
	// 	ImageTypePointer grad_d = derivf->GetOutput();
	// 
	// 	VectorType gradI_d = GridFunctionsType::Interpolate(Yup, grad_d);
	// 
	// 	for (unsigned int i = 0; i < numVoxels; i++)
	// 		gradI_d[i] *= D[i];
	// 
	// 	// ImageTypePointer Gradient = GridFunctionsType::VectorToImage(m_Image, gradI_d);
	// 	// writer->SetInput(Gradient);
	// 	// if (dim==0)
	// 	// 	writer->SetFileName("gradientX.img");
	// 	// else if (dim==1)
	// 	// 	writer->SetFileName("gradientY.img");
	// 	// else
	// 	// 	writer->SetFileName("gradientZ.img");
	// 	// 	
	// 	// writer->Update();
	// 	
	// 	
	// 	// gradMatch.set_column(dim, gradI_d);
	// 	gradI_d *= Superclass::m_FlipAxes[dim];
	// 	gradMatch.set_column(Superclass::m_PermutationAxes[dim], gradI_d);
	// }

	gradMatch *= 2.0;

	return gradMatch;
}



#endif /* _SSDImage_txx */

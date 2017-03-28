/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _PointCloud_txx
#define _PointCloud_txx

#include "PointCloud.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include <cassert>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::PointCloud() : Superclass()
{
	this->SetPointCloudType();
	m_KernelWidth = 0;
	m_KernelType = null;
}



template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::~PointCloud()
{

}



template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::PointCloud(const PointCloud& other) : Superclass(other)
 {
	this->SetPointCloudType();

	m_PointWeights = other.m_PointWeights;

	m_KernelWidth = other.m_KernelWidth;
	m_KernelType = other.m_KernelType;
	m_NormSquared = other.m_NormSquared;

 }

template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::PointCloud(const PointCloud& ex, const MatrixType& LP) : Superclass(ex, LP)
 {
	this->SetPointCloudType();

	m_PointWeights = ex.m_PointWeights;

	m_KernelWidth = ex.m_KernelWidth;
	m_KernelType = ex.m_KernelType;
	
	this->UpdateSelfNorm();
	this->UnSetModified();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
::Update()
  {
	  if (this->IsModified())
	  {
		  Superclass::Update();
		  this->UpdatePointWeights();
		  this->UpdateSelfNorm();
	  }
	  this->UnSetModified();
 }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
TScalar
PointCloud<TScalar, Dimension>
::ComputeMatch(const PointCloud<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	const PointCloud* targetPointCloud = dynamic_cast<const PointCloud*>(target);

	if (m_KernelWidth != targetPointCloud->GetKernelWidth())
		throw std::runtime_error("Kernel width of point clouds mismatched");

	this->Update();
	// target->Update();

	MatrixType targPts = targetPointCloud->GetPointCoordinates();
	MatrixType targWts = targetPointCloud->GetPointWeights();
	MatrixType Pts = this->GetPointCoordinates();

	TScalar match = targetPointCloud->GetNormSquared() + this->GetNormSquared();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(Pts);
	kernelObject->SetWeights(m_PointWeights);

	MatrixType SdotT = kernelObject->Convolve(targPts);
	for (int i = 0; i < targetPointCloud->GetNumberOfPoints(); i++)
		match -= SdotT(i,0) * (2.0 * targWts(i,0));

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename PointCloud<TScalar, Dimension>::MatrixType
PointCloud<TScalar, Dimension>
::ComputeMatchGradient(const PointCloud<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	const PointCloud* targetPointCloud = dynamic_cast<const PointCloud*>(target);
	
	this->Update();
	// target->Update();

	if (m_KernelWidth != targetPointCloud->GetKernelWidth())
		throw std::runtime_error("Kernel width of point clouds mismatched");

	MatrixType targPts = targetPointCloud->GetPointCoordinates();
	MatrixType targWts = targetPointCloud->GetPointWeights();
	MatrixType Pts = this->GetPointCoordinates();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(Pts);
	kernelObject->SetWeights(m_PointWeights);
	MatrixType grad_SdotS = kernelObject->ConvolveGradient(Pts, m_PointWeights);

	kernelObject->SetSources(targPts);
	kernelObject->SetWeights(targWts);
	MatrixType grad_SdotT = kernelObject->ConvolveGradient(Pts, m_PointWeights);

	MatrixType gradmatch = (grad_SdotS + grad_SdotT) * 2.0;
	assert(gradmatch.rows() == this->GetNumberOfPoints());
	assert(gradmatch.columns() == Dimension);

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
::UpdatePointWeights()
 {
	m_PointWeights.set_size(Superclass::m_NumberOfPoints,1);
	m_PointWeights.fill(1.0);

	vtkSmartPointer<vtkPointData> pd = Superclass::m_PolyData->GetPointData();
	int numCmp = pd->GetNumberOfComponents();
	if (numCmp==0)
		std::cout << "Warning: No weights detected: use unit weight for each point" << std::endl;
	else if (numCmp==1)
	{
		vtkSmartPointer<vtkDataArray> scalars = pd->GetScalars();
		int numT = scalars->GetNumberOfTuples();
		if (numT != Superclass::m_NumberOfPoints)
			std::cout << "number of scalars and number of points mismatched! Defaulting to unit weights" << std::endl;
		else
			for (int i = 0; i < numT; i++)
				m_PointWeights.set_row(i, scalars->GetComponent(i,0)); // 0 because of scalar data
	}

 }



template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
::UpdateSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(this->GetPointCoordinates());
	kernelObject->SetWeights(m_PointWeights);

	MatrixType selfKW = kernelObject->Convolve(this->GetPointCoordinates());

	m_NormSquared = 0;

	for (unsigned int i = 0; i < this->GetNumberOfPoints(); i++)
		m_NormSquared += selfKW(i, 0) * m_PointWeights(i,0);

	delete kernelObject;
 }


#endif /* _PointCloud_txx */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedPolyLine_txx
#define _OrientedPolyLine_txx

#include "OrientedPolyLine.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include <cassert>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::OrientedPolyLine() : Superclass()
 {
	this->SetOrientedPolyLineType();
	m_KernelWidth = 0;
	m_KernelType = null;
 }



template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::~OrientedPolyLine()
{

}



template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::OrientedPolyLine(const OrientedPolyLine& other) : Superclass(other)
 {
	this->SetOrientedPolyLineType();

	m_Centers = other.m_Centers;
	m_Tangents = other.m_Tangents;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_KernelType = other.m_KernelType;
	m_NormSquared = other.m_NormSquared;

 }



template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::OrientedPolyLine(const OrientedPolyLine& ex, const MatrixType& LP) : Superclass(ex, LP)
 {
	this->SetOrientedPolyLineType();

	m_KernelWidth = ex.m_KernelWidth;
	m_KernelType = ex.m_KernelType;

	m_NumCells = ex.m_NumCells;

	// this is required since the call to Superclass(ex, LP) sets m_IsModified to false
	this->SetModified();
	this->Update();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	Superclass::SetPolyData(polyData);

	m_NumCells = polyData->GetNumberOfCells();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::Update()
 {
	 if (this->IsModified())
	 {
		 Superclass::Update();	
		 this->UpdateCentersTangents();
		 // replace the BoundingBox computed with vertices by the bounding box computed with centers of triangles 
		 this->UpdateBoundingBox();
		 this->UpdateSelfNorm();
	 }
	 this->UnSetModified();
 }



template <class TScalar, unsigned int Dimension>
TScalar
OrientedPolyLine<TScalar, Dimension>
::ComputeMatch(const OrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	const OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<const OrientedPolyLine*>(target);
	
	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
	{
		std::cerr << "DEBUG kw = "  << m_KernelWidth << " target kw = " << targetOrientedPolyLine->GetKernelWidth() << std::endl;
		throw std::runtime_error("Kernel width of curve currents mismatched");
	}

	this->Update();
	// target->Update();

	MatrixType targCenters = targetOrientedPolyLine->GetCenters();
	MatrixType targTangents = targetOrientedPolyLine->GetTangents();

	TScalar match = targetOrientedPolyLine->GetNormSquared() + this->GetNormSquared();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_Tangents);

	MatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetOrientedPolyLine->GetNumberOfCells(); i++)
		match -= 2.0 * dot_product( SdotT.get_row(i),  targTangents.get_row(i) );

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename OrientedPolyLine<TScalar, Dimension>::MatrixType
OrientedPolyLine<TScalar, Dimension>
::ComputeMatchGradient(const OrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	const OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<const OrientedPolyLine*>(target);
	
	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
		throw std::runtime_error("Kernel width of curve currents mismatched");

	this->Update();
	// target->Update();

	MatrixType targCenters = targetOrientedPolyLine->GetCenters();
	MatrixType targTangents = targetOrientedPolyLine->GetTangents();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_Tangents);
	MatrixType KtauS = kernelObject->Convolve(m_Centers);
	MatrixType gradKtauS = kernelObject->ConvolveGradient(m_Centers, m_Tangents);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targTangents);
	MatrixType KtauT = kernelObject->Convolve(m_Centers);
	MatrixType gradKtauT = kernelObject->ConvolveGradient(m_Centers, m_Tangents);

	MatrixType gradKtau = gradKtauS - gradKtauT;

	MatrixType gradmatch(this->GetNumberOfPoints(),Dimension,0.0);
	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int indM = ptIds->GetId(0);
		int indP = ptIds->GetId(1);

		VectorType Ktau = KtauS.get_row(f) - KtauT.get_row(f);
		Ktau *= 2.0;
		gradmatch.set_row(indM, gradmatch.get_row(indM) - Ktau);
		gradmatch.set_row(indP, gradmatch.get_row(indP) + Ktau);

		gradmatch.set_row(indM, gradmatch.get_row(indM) + gradKtau.get_row(f));
		gradmatch.set_row(indP, gradmatch.get_row(indP) + gradKtau.get_row(f));

	}

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::UpdateCentersTangents()
 {
	const MatrixType Pts = this->GetPointCoordinates();
	m_Centers.set_size(m_NumCells, Dimension);
	m_Tangents.set_size(m_NumCells, Dimension);

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		if (ptIds->GetNumberOfIds() != 2)
			throw std::runtime_error("Not a polygonal line!");

		VectorType p0 = Pts.get_row(ptIds->GetId(0));
		VectorType p1 = Pts.get_row(ptIds->GetId(1));
		m_Centers.set_row(i, (p0 + p1) / 2.0 );
		m_Tangents.set_row(i, p1 - p0);

	}
 }



template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::UpdateSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_Tangents);

	MatrixType selfKW = kernelObject->Convolve(m_Centers);

	m_NormSquared = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
		m_NormSquared += dot_product(selfKW.get_row(i), m_Tangents.get_row(i));

	delete kernelObject;
 }



template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::UpdateBoundingBox()
{
	// overload bounding box using the centers instead of the points (useful when using the results of matching pursuit as target)
	VectorType Min(Dimension, 1e20);
	VectorType Max(Dimension, -1e20);
	
	for (int i = 0; i < m_NumCells; i++)
	{
		VectorType p = m_Centers.get_row(i);
		for (unsigned int dim=0; dim<Dimension; dim++)
		{	
			Min[dim] = (p[dim] < Min[dim])?p[dim]:Min[dim];
			Max[dim] = (p[dim] > Max[dim])?p[dim]:Max[dim];
		}
	}
	
 	Superclass::Superclass::m_BoundingBox.set_column(0, Min);
 	Superclass::Superclass::m_BoundingBox.set_column(1, Max);
}


#endif /* _OrientedPolyLine_txx */

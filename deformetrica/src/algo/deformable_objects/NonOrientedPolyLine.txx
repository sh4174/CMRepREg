/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedPolyLine_txx
#define _NonOrientedPolyLine_txx

#include "NonOrientedPolyLine.h"

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



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::NonOrientedPolyLine() : Superclass()
 {
	this->SetNonOrientedPolyLineType();
	m_KernelWidth = 0;
	m_KernelType = null;
 }



template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::NonOrientedPolyLine(const NonOrientedPolyLine& other) : Superclass(other)
 {
	this->SetNonOrientedPolyLineType();

	m_Centers = other.m_Centers;
	m_Tangents = other.m_Tangents;
	m_MatrixTangents = other.m_MatrixTangents;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_KernelType = other.m_KernelType;
	m_NormSquared = other.m_NormSquared;
	
 }



template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::NonOrientedPolyLine(const NonOrientedPolyLine& ex, const MatrixType& LP) : Superclass(ex, LP)
 {
 	this->SetNonOrientedPolyLineType();

 	m_KernelWidth = ex.m_KernelWidth;
 	m_KernelType = ex.m_KernelType;

 	m_NumCells = ex.m_NumCells;
	
	// this is required since the call to Superclass(ex, LP) sets m_IsModified to false
	this->SetModified();
	this->Update();
	
 }



template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::~NonOrientedPolyLine()
{}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
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
NonOrientedPolyLine<TScalar, Dimension>
::Update()
{
 	if (this->IsModified())
	{
		Superclass::Update();
		this->UpdateCentersTangents();
        // replace the BoundingBox computed with vertices by the bounding box computed with centers of segments
		this->UpdateBoundingBox();
		this->UpdateSelfNorm();
	}
	this->UnSetModified();
}



template <class TScalar, unsigned int Dimension>
TScalar
NonOrientedPolyLine<TScalar, Dimension>
::ComputeMatch(const NonOrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");
	
	const NonOrientedPolyLine* targetNonOrientedPolyLine = dynamic_cast<const NonOrientedPolyLine*>(target);

	if (m_KernelWidth != targetNonOrientedPolyLine->GetKernelWidth())
	{
		std::cerr << "DEBUG kw = "  << m_KernelWidth << " target kw = " << targetNonOrientedPolyLine->GetKernelWidth() << std::endl;
		throw std::runtime_error("Kernel width of curve currents mismatched");
	}

	this->Update();
	// target->Update();

	TScalar match = targetNonOrientedPolyLine->GetNormSquared() + this->GetNormSquared();
	
	MatrixType targCenters = targetNonOrientedPolyLine->GetCenters();
	MatrixType targTangents = targetNonOrientedPolyLine->GetTangents();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixTangents);

	MatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetNonOrientedPolyLine->GetNumberOfCells(); i++)
	{
		VectorType Mi = special_product(SdotT.get_row(i), targTangents.get_row(i));
		match -= 2.0 * dot_product(targTangents.get_row(i), Mi);
	}

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedPolyLine<TScalar, Dimension>::MatrixType
NonOrientedPolyLine<TScalar, Dimension>
::ComputeMatchGradient(const NonOrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");
	
	const NonOrientedPolyLine* targetNonOrientedPolyLine = dynamic_cast<const NonOrientedPolyLine*>(target);

	if (m_KernelWidth != targetNonOrientedPolyLine->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	this->Update();
	// target->Update();

	MatrixType targCenters = targetNonOrientedPolyLine->GetCenters();
	MatrixType targMatrixTangents = targetNonOrientedPolyLine->GetMatrixTangents();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixTangents);
	MatrixType KtauS = kernelObject->Convolve(m_Centers);
	MatrixList gradKtauS = kernelObject->ConvolveGradient(m_Centers);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targMatrixTangents);
	MatrixType KtauT = kernelObject->Convolve(m_Centers);
	MatrixList gradKtauT = kernelObject->ConvolveGradient(m_Centers);

	MatrixType gradmatch(this->GetNumberOfPoints(),Dimension,0.0);

	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int ind0 = ptIds->GetId(0);
		int ind1 = ptIds->GetId(1);

		VectorType defN = m_Tangents.get_row(f);
		TScalar defN_mag2 = defN.squared_magnitude();

		if (defN_mag2 > 1e-20)
		{
			VectorType Ktau = special_product(KtauS.get_row(f) - KtauT.get_row(f), defN);
			Ktau *= 4.0;

			VectorType aux = special_product(m_MatrixTangents.get_row(f), Ktau);
			Ktau = (Ktau - aux / (2*defN_mag2) ) / sqrt(defN_mag2);

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) - Ktau );
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + Ktau );

			MatrixType delta = gradKtauS[f] - gradKtauT[f];
			VectorType gradKtau(Dimension);
			for (int p = 0; p < Dimension; p++)
			{
				VectorType Mf = special_product(delta.get_column(p), defN);
				gradKtau(p) = dot_product(Mf, defN);
			}

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau);
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau);
		}

	}

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
::UpdateCentersTangents()
 {
	const MatrixType Pts = this->GetPointCoordinates();
	m_Centers.set_size(m_NumCells, Dimension);
	m_Tangents.set_size(m_NumCells, Dimension);
	m_MatrixTangents.set_size(m_NumCells, Dimension*(Dimension+1)/2);

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		VectorType p0 = Pts.get_row(ptIds->GetId(0));
		VectorType p1 = Pts.get_row(ptIds->GetId(1));
		m_Centers.set_row(i, (p0 + p1) / 2.0 );

		VectorType Ti = p1-p0;
		Ti /= sqrt( Ti.magnitude() + 1e-20 ); // divided by norm^(1/2)
		m_Tangents.set_row(i, Ti);

		VectorType aux(Dimension*(Dimension+1)/2);
		int index = 0;
		for (int p = 0; p < Dimension; p++)
			for (int q = p; q < Dimension; q++)
				aux(index++) = Ti(p)*Ti(q);

		m_MatrixTangents.set_row(i, aux);
	}
 }



template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
::UpdateSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixTangents);

	MatrixType selfKW = kernelObject->Convolve(m_Centers);

	m_NormSquared = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VectorType Mi = special_product(selfKW.get_row(i), m_Tangents.get_row(i));
		m_NormSquared += dot_product(m_Tangents.get_row(i), Mi);
	}

	delete kernelObject;

 }



template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
::UpdateBoundingBox()
{
	// requires this->UpdateCentersTangents() to be called beforehand.

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



template <class TScalar, unsigned int Dimension>
typename NonOrientedPolyLine<TScalar, Dimension>::VectorType
NonOrientedPolyLine<TScalar, Dimension>
::special_product(VectorType M, VectorType X)
 {
	MatrixType Ms(Dimension,Dimension);
	int index = 0;
	for (int p = 0; p < Dimension; p++)
		for (int q = p; q < Dimension; q++)
		{
			Ms(p,q) = M(index++);
			Ms(q,p) = Ms(p,q);
		}

	return (Ms * X);
 }


#endif /* _NonOrientedPolyLine_txx */

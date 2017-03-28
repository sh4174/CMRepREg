/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedSurfaceMesh_txx
#define _NonOrientedSurfaceMesh_txx

#include "NonOrientedSurfaceMesh.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkTriangleFilter.h"
#include "myvtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include <cstring>
#include <iostream>
#include <sstream>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::NonOrientedSurfaceMesh() : Superclass()
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetNonOrientedSurfaceMeshType();
	m_KernelWidth = 0;
	m_KernelType = null;
 }



template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::~NonOrientedSurfaceMesh()
{}



template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::NonOrientedSurfaceMesh(const NonOrientedSurfaceMesh& other) : Superclass(other)
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetNonOrientedSurfaceMeshType();

	m_Centers = other.m_Centers;
	m_Normals = other.m_Normals;
	m_MatrixNormals = other.m_MatrixNormals;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_KernelType = other.m_KernelType;
	m_NormSquared = other.m_NormSquared;

 }



template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::NonOrientedSurfaceMesh(const NonOrientedSurfaceMesh& ex, const MatrixType& LP ) : Superclass(ex, LP)
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetNonOrientedSurfaceMeshType();

	m_KernelWidth = ex.m_KernelWidth;
	m_KernelType = ex.m_KernelType;
	
	// this is required since the call to Superclass(ex, LP) sets m_IsModified to false
	this->SetModified();
	this->Update();
	

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
 
template <class TScalar, unsigned int Dimension>
void
NonOrientedSurfaceMesh<TScalar, Dimension>
::Update()
 {
	 if (this->IsModified())
	 {
		 m_NumCells = Superclass::m_PolyData->GetNumberOfCells();
		 this->UpdateCentersNormals();
		 this->UpdateBoundingBox();
		 this->UpdateSelfNorm();
	 }
	this->UnSetModified();
 }


template <class TScalar, unsigned int Dimension>
TScalar
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatch(const NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	const NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<const NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface unsigned currents mismatched");

	this->Update();
	// target->Update();

	MatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	MatrixType targNormals = targetNonOrientedSurfaceMesh->GetNormals();

	TScalar match = targetNonOrientedSurfaceMesh->GetNormSquared() + this->GetNormSquared();
	
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixNormals);

	MatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetNonOrientedSurfaceMesh->GetNumberOfCells(); i++)
	{
		VectorType Mi = special_product(SdotT.get_row(i), targNormals.get_row(i));
		match -= 2.0 * dot_product(targNormals.get_row(i), Mi);		
	}

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedSurfaceMesh<TScalar, Dimension>::MatrixType
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatchGradient(const NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	const NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<const NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	this->Update();
	// target->Update();

	MatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	MatrixType targMatrixNormals = targetNonOrientedSurfaceMesh->GetMatrixNormals();

	MatrixType Pts = this->GetPointCoordinates();

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixNormals);
	MatrixType KtauS = kernelObject->Convolve(m_Centers);
	MatrixList gradKtauS = kernelObject->ConvolveGradient(m_Centers);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targMatrixNormals);
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
		int ind2 = ptIds->GetId(2);

		VectorType e0 = Pts.get_row(ind2) - Pts.get_row(ind1);
		VectorType e1 = Pts.get_row(ind0) - Pts.get_row(ind2);
		VectorType e2 = Pts.get_row(ind1) - Pts.get_row(ind0);

		VectorType defN = m_Normals.get_row(f);
		TScalar defN_mag2 = defN.squared_magnitude();

		if (defN_mag2 > 1e-20)
		{
			VectorType Ktau = special_product(KtauS.get_row(f) - KtauT.get_row(f), defN);
			Ktau *= 2.0;

			VectorType aux = special_product(m_MatrixNormals.get_row(f), Ktau);
			Ktau = (Ktau - aux / (2*defN_mag2) ) / sqrt(defN_mag2);

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + cross_3d(e0, Ktau) );
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + cross_3d(e1, Ktau) );
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + cross_3d(e2, Ktau) );

			MatrixType delta = gradKtauS[f] - gradKtauT[f];
			VectorType gradKtau(Dimension);
			for (int p = 0; p < Dimension; p++)
			{
				VectorType Mf = special_product(delta.get_column(p), defN);
				gradKtau(p) = dot_product(Mf, defN);
			}
			gradKtau *= 2.0/3.0;

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau);
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau);
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau);		
		}

	}

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

// make sure the ordering of point cells is consistent with the direction of the normal. Is that really need when orientation does not matter?
// template <class TScalar, unsigned int Dimension>
// void
// NonOrientedSurfaceMesh<TScalar, Dimension>
// ::CheckMeshAndNormals()
//  {
// 	vtkSmartPointer<vtkTriangleFilter> trif = vtkSmartPointer<vtkTriangleFilter>::New();
// 	trif->SetInput(Superclass::m_PolyData);
// 	trif->PassVertsOff();
// 	trif->PassLinesOff();
// 	trif->Update();
// 
// 	vtkSmartPointer<myvtkPolyDataNormals> normalf =	vtkSmartPointer<myvtkPolyDataNormals>::New();
// 	normalf->SetInput(trif->GetOutput());
// 	normalf->ComputePointNormalsOff();
// 	normalf->ComputeCellNormalsOn();
// 	normalf->SplittingOff();
// 	normalf->ConsistencyOff();
// 	normalf->AutoOrientNormalsOff();
// 
// 	normalf->Update();
// 
// 	Superclass::m_PolyData = normalf->GetOutput();
// 	Superclass::m_PolyData->BuildLinks();
// 
// 	m_NumCells = Superclass::m_PolyData->GetNumberOfCells();
// 
//  }


template <class TScalar, unsigned int Dimension>
void
NonOrientedSurfaceMesh<TScalar, Dimension>
::UpdateCentersNormals()
 {
	const MatrixType Pts = this->GetPointCoordinates();
	m_Centers.set_size(m_NumCells, Dimension);
	m_Normals.set_size(m_NumCells, Dimension);
	m_MatrixNormals.set_size(m_NumCells, Dimension*(Dimension+1)/2);

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		if (ptIds->GetNumberOfIds() != 3)
			throw std::runtime_error("Not a triangle cell!");

		VectorType p0 = Pts.get_row(ptIds->GetId(0));
		VectorType p1 = Pts.get_row(ptIds->GetId(1));
		VectorType p2 = Pts.get_row(ptIds->GetId(2));

		m_Centers.set_row(i, (p0 + p1 + p2) / 3.0 );

		VectorType Ni = cross_3d(p2 - p0, p1 - p0) / 2;
		Ni /= sqrt( Ni.magnitude() + 1e-20 ); // divided by norm^(1/2)
		m_Normals.set_row(i, Ni);

		VectorType aux(Dimension*(Dimension+1)/2);
		int index = 0;
		for (int p = 0; p < Dimension; p++)
			for (int q = p; q < Dimension; q++)
				aux(index++) = Ni(p)*Ni(q);

		m_MatrixNormals.set_row(i, aux);
	}
 }



template <class TScalar, unsigned int Dimension>
void
NonOrientedSurfaceMesh<TScalar, Dimension>
::UpdateSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixNormals);

	MatrixType selfKW = kernelObject->Convolve(m_Centers);

	m_NormSquared = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VectorType Mi = special_product(selfKW.get_row(i), m_Normals.get_row(i));
		m_NormSquared += dot_product(m_Normals.get_row(i), Mi);
	}

	delete kernelObject;
 }



template <class TScalar, unsigned int Dimension>
void
NonOrientedSurfaceMesh<TScalar, Dimension>
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



template <class TScalar, unsigned int Dimension>
typename NonOrientedSurfaceMesh<TScalar, Dimension>::VectorType
NonOrientedSurfaceMesh<TScalar, Dimension>
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



#endif

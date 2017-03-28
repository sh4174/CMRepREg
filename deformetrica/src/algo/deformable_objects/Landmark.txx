/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Landmark_txx
#define _Landmark_txx

#include "Landmark.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "SimpleTimer.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::Landmark() : Superclass()
 {
	this->SetLandmarkType();
	m_NumberOfPoints = 0;
 }



template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::~Landmark()
{

}



template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::Landmark(const Landmark& other) : Superclass(other)
 {
	this->SetLandmarkType();

	m_PolyData = vtkSmartPointer<vtkPolyData>::New();
	m_VTKMutex.Lock();
	m_PolyData->DeepCopy(other.m_PolyData);
	m_VTKMutex.Unlock();
	
	m_PointCoordinates = other.m_PointCoordinates;

	m_NumberOfPoints = other.m_NumberOfPoints;

 }



template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::Landmark(const Landmark& example, const MatrixType& LandmarkPoints) : Superclass(example, LandmarkPoints)
 {
	this->SetLandmarkType();

	m_NumberOfPoints = example.m_NumberOfPoints;
	
	if (m_NumberOfPoints != LandmarkPoints.rows())
		throw std::runtime_error("Number of LandmarkPoints mismatch in copy/update Landmark constructor");
	

	m_PolyData = vtkSmartPointer<vtkPolyData>::New();
	m_VTKMutex.Lock();
	m_PolyData->DeepCopy(example.m_PolyData);
	m_VTKMutex.Unlock();
	
	this->UpdatePolyDataWithPointCoordinates(LandmarkPoints);
	
	this->Update();	
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
{
	if (Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel().compare("LPS") != 0)
	{
		std::cout << "Reorienting polydata from " << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " to LPS" << std::endl;
		m_PolyData = ReorientPolyData(polyData, false);
	}
	else
		m_PolyData = polyData;

	m_NumberOfPoints = m_PolyData->GetNumberOfPoints();
	
	m_PointCoordinates.set_size(m_NumberOfPoints,Dimension);
	
	for (unsigned int i = 0; i < m_NumberOfPoints; i++)
	{
		// Here we used 3 since 2D points still have a z-coordinate that is equal to 0 in vtkPolyData.
	 	// This coordinate is removed in m_WorkingPointCoordinates
		 double p[3];
		 m_VTKMutex.Lock();
		 m_PolyData->GetPoint(i, p);
		 m_VTKMutex.Unlock();
		 for (int dim = 0; dim < Dimension; dim++)
		 	m_PointCoordinates(i, dim) = p[dim];
	}
	
	this->SetModified(); 
}


template <class TScalar, unsigned int Dimension>
vtkSmartPointer< vtkPolyData > 
Landmark<TScalar, Dimension>
::GetPolyData() const
{
	return m_PolyData;	
}


template <class TScalar, unsigned int Dimension>
vtkSmartPointer<vtkPolyData>
Landmark<TScalar, Dimension>
::ReorientPolyData(vtkPolyData* polyData, bool inverse) const
{
	double trfmtx[16];
	MatrixType orient_mtx;
	vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();

	if (inverse == false)
		orient_mtx = Superclass::m_AnatomicalOrientation.GetChangeOfBasisMatrix();
	else
		orient_mtx = Superclass::m_AnatomicalOrientation.GetInverseChangeOfBasisMatrix();

	int idx = 0;
	for (unsigned int i = 0; i < Dimension; i++)
	{
		for (unsigned int j = 0; j < Dimension; j++)
		{
			trfmtx[idx++] = orient_mtx(i, j);
			std::cout << trfmtx[idx-1] << " ";
		}
		trfmtx[idx++] = 0;
		std::cout << trfmtx[idx-1] << std::endl;
	}
	for (unsigned int j = 0; j < Dimension; j++)
	{
		trfmtx[idx++] = 0;
		std::cout << trfmtx[idx-1] << " ";
	}
	trfmtx[idx++] = 1;
	std::cout << trfmtx[idx-1] << std::endl;

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->SetMatrix(trfmtx);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyData =	vtkSmartPointer<vtkTransformPolyDataFilter>::New();

#if VTK_MAJOR_VERSION <= 5
	transformPolyData->SetInput(polyData);
    //transformPolyData->SetInputConnection(polyData->GetProducerPort());
#else
	transformPolyData->SetInputData(polyData);
#endif
	transformPolyData->SetTransform(transform);
	transformPolyData->SetOutput(newPolyData);
	transformPolyData->Update();

	return transformPolyData->GetOutput();
}



template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::UpdatePolyDataWithPointCoordinates(const MatrixType& Y)
 {
	// A polydata needs to be set before updating Point Coordinates
	if (m_PolyData == NULL)
		 throw std::runtime_error("a VTK PolyData should be set before setting new point coordinates");
	 
	if (Y.rows() != m_PolyData->GetNumberOfPoints())
		throw std::runtime_error("number of points mismatched");
	if (Y.columns() != Dimension)
		throw std::runtime_error("Dimension mismatched");

	
	for (unsigned int i = 0; i < m_NumberOfPoints; i++)
	{
		double p[3];
		p[2] = 0.0;
		for (int dim = 0; dim < Dimension; dim++)
			p[dim] = Y(i, dim);
	
		m_VTKMutex.Lock();
		m_PolyData->GetPoints()->SetPoint(i, p);
		m_VTKMutex.Unlock();
	}

	m_PointCoordinates = Y;

	this->SetModified();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::Update()
 {
	 if (m_PolyData == NULL)
		 throw std::runtime_error("A VTK PolyData should have been set in Landmark class or its children classes");
	 
	 
	 if (this->IsModified())
	 {	 
		 this->UpdateBoundingBox();
	 }
	 
 	this->UnSetModified(); 
 }


template <class TScalar, unsigned int Dimension>
TScalar
Landmark<TScalar, Dimension>
::ComputeMatch(const Landmark<TScalar, Dimension>::Superclass* target)
 {
	if (this->GetType() != target->GetType())
	{
		throw std::runtime_error("Deformable object types mismatched");
	}

	const Landmark* targetLandmark = dynamic_cast<const Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmark object should have the same number of points");

	MatrixType D = this->GetPointCoordinates() - targetLandmark->GetPointCoordinates();

	TScalar match = D.frobenius_norm();
	match *= match;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename Landmark<TScalar, Dimension>::MatrixType
Landmark<TScalar, Dimension>
::ComputeMatchGradient(const Landmark<TScalar, Dimension>::Superclass* target)
 {
	if ( this->GetType() != target->GetType())
	{
		throw std::runtime_error("Deformable object types mismatched");
	}


	const Landmark* targetLandmark = dynamic_cast<const Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmarks object should have the same number of points");

	MatrixType gradMatch = this->GetPointCoordinates() - targetLandmark->GetPointCoordinates();
	gradMatch *= 2.0;

	return gradMatch;
 }



template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::WriteObject(std::string str) const
{
	if (Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel().compare("LPS") != 0)
	{
		vtkSmartPointer<vtkPolyData> outData = vtkSmartPointer<vtkPolyData>::New();
		std::cout << "Reorienting polydata to " << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " before writing output file" << endl;
		outData = ReorientPolyData(m_PolyData, true);

		vtkSmartPointer<vtkPolyDataWriter> writer =
				vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
		writer->SetInput(outData);
#else
		writer->SetInputData(outData);
#endif
		writer->Update();
	}
	else
	{
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
		writer->SetInput(m_PolyData);
#else
		writer->SetInputData(m_PolyData); //
#endif
		writer->Update();
	}
}


template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::WriteObject(std::string str, const MatrixList& velocity) const
{
	m_VTKMutex.Lock();
	vtkSmartPointer<vtkPolyData> outData = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkDoubleArray> velocityVecs = vtkSmartPointer<vtkDoubleArray>::New();
	velocityVecs->SetNumberOfComponents(3);
	velocityVecs->SetNumberOfTuples(m_NumberOfPoints);
	velocityVecs->SetName("velocity");

	outData->DeepCopy(m_PolyData);
	m_VTKMutex.Unlock();

	for (int t=0; t<velocity.size(); t++)
	{
		MatrixType VT = velocity[t];

		for (int i=0; i<m_NumberOfPoints; i++)
		{
			double v[Dimension];

			for (int dim = 0; dim < Dimension; dim++)
			{
				v[dim] = VT(i, dim);
			}

			m_VTKMutex.Lock();
			velocityVecs->SetTuple(i, v);
			m_VTKMutex.Unlock();
		}
	}

	m_VTKMutex.Lock();
	//outData->GetFieldData()->AddArray(velocityVecs);
	// outData->GetCellData()->AddArray(velocityVecs);
	outData->GetPointData()->SetVectors(velocityVecs);
	m_VTKMutex.Unlock();

	m_PolyData->DeepCopy(outData);

	if (Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel().compare("LPS") != 0)
	{
		vtkSmartPointer<vtkPolyData> outData = vtkSmartPointer<vtkPolyData>::New();
		std::cout << "Reorienting polydata to " << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " before writing output file" << endl;
		outData = ReorientPolyData(m_PolyData, true);

		vtkSmartPointer<vtkPolyDataWriter> writer =
				vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
		writer->SetInput(outData);
#else
		writer->SetInputData(outData);
#endif
		writer->Update();
	}
	else
	{
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
		writer->SetInput(m_PolyData);
#else
		writer->SetInputData(m_PolyData); //
#endif
		writer->Update();
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::UpdateBoundingBox()
 {
 
 	VectorType Min(Dimension, 1e20);
 	VectorType Max(Dimension, -1e20);

 	for (int i = 0; i < m_NumberOfPoints; i++)
 	{
 		VectorType p = m_PointCoordinates.get_row(i);
 		for (unsigned int dim=0; dim<Dimension; dim++)
 		{	
 			Min[dim] = (p[dim] < Min[dim])?p[dim]:Min[dim];
 			Max[dim] = (p[dim] > Max[dim])?p[dim]:Max[dim];
 		}
 	}

	Superclass::m_BoundingBox.set_column(0, Min);
	Superclass::m_BoundingBox.set_column(1, Max);
 }



#endif /* _Landmark_txx */

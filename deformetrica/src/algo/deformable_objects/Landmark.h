/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Landmark_h
#define _Landmark_h

#include "DeformableObject.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkVersion.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief 		Landmarks (i.e. labelled point sets)
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The Landmark class inherited from DeformableObject represents a set of labelled points.
 *              This class assumes that the source and the target have the same number of points with a
 *              point-to-point correspondence.
 */
template <class TScalar, unsigned int Dimension>
class Landmark : public DeformableObject<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> Superclass;

	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;
	/// List of matrices type.
	typedef typename Superclass::MatrixList MatrixList;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Landmark();
	/// Copy constructor.
	Landmark(const Landmark& other);
	/// Constructor which copies the objects and update vertex coordinates.
	Landmark(const Landmark& example, const MatrixType& LandmarkPoints);
	
	/// Returns a Deformed version of the point set, where the deformation is given by the position of vertices in \e LandmarkPoints.
	virtual Landmark* DeformedObject(const MatrixType& LandmarkPoints) { return new Landmark(*this, LandmarkPoints); }

	virtual Landmark* Clone() { return new Landmark(*this); }

	virtual ~Landmark();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the pointer on a VTK object to a reference \e polyData.
	virtual void SetPolyData(vtkPolyData* polyData);
	virtual vtkSmartPointer< vtkPolyData > GetPolyData() const;



	/// Transforms polyData from native orientation system into ITK (LPS) default orientation,
	/// or inversely, if the flag "inverse" is set to true.
	vtkSmartPointer<vtkPolyData> ReorientPolyData(vtkPolyData* polyData, bool inverse) const;

	/// Update the PolyData with new coordinates of vertices. Need a call to Update() afterwards.
	void UpdatePolyDataWithPointCoordinates(const MatrixType& LandmarkPoints);

	/// Returns the vertex coordinates as a matrix.
	inline MatrixType GetPointCoordinates() const { return m_PointCoordinates; }

	/// Returns the number of points of the deformable object.
	inline int GetNumberOfPoints() const { return m_NumberOfPoints; }

	/// Returns the dimension of the discretized image, here the number of points times the dimension.
	virtual int GetDimensionOfDiscretizedObject() const { return (Dimension * m_NumberOfPoints); }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void Update();

	virtual TScalar ComputeMatch(const Superclass* target);
	virtual MatrixType ComputeMatchGradient(const Superclass* target);

	virtual void WriteObject(std::string filename) const;
	virtual void WriteObject(std::string filename, const MatrixList& velocity) const;


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the bounding box of the data.
	void UpdateBoundingBox();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Pointer on a geometric structure namely a VTK object.
	vtkSmartPointer<vtkPolyData> m_PolyData;

	///	Matrix coordinates of the points (Size : NumberOfPoints x Dimension).
	MatrixType m_PointCoordinates;

	///	Number of points of the deformable object.
	int m_NumberOfPoints;

	///	Object used to perform mutex (mutual exclusion) with m_ReferencePolyData (important for multithreaded programming).
	itk::SimpleFastMutexLock m_VTKMutex;


}; /* class Landmark */


#ifndef MU_MANUAL_INSTANTIATION
#include "Landmark.txx"
#endif


#endif /* _Landmark_h */

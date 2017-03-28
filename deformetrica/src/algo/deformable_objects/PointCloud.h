/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _PointCloud_h
#define _PointCloud_h

#include "DeformableObject.h"
#include "Landmark.h"

#include "Diffeos.h"

#include "KernelType.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      Point clouds.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The PointCloud class inherited from Landmark represents an unstructured unlabelled point set.
 *              This class uses the 0-current representation (aka measures) which does not assume point-to-point
 *              correspondence between source and target point clouds to be match.
 */
template <class TScalar, unsigned int Dimension>
class PointCloud : public Landmark<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Landmark type.
	typedef Landmark<TScalar, Dimension> Superclass;

	/// Deformable object type.
	typedef typename Superclass::Superclass DeformableObjectType;

	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;
	/// List of matrices type.
	typedef typename Superclass::MatrixList MatrixList;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	PointCloud();
	/// Copy constructor.
	PointCloud(const PointCloud& other);
	/// Constructor which copies the object and update vertex coordinates
	PointCloud(const PointCloud& example, const MatrixType& LandmarkPoints);
	
	/// Returns a Deformed version of the point cloud, where the deformation is given by the position of vertices in \e LandmarkPoints
	virtual PointCloud* DeformedObject(const MatrixType& LandmarkPoints) { return new PointCloud(*this, LandmarkPoints); }

	virtual PointCloud* Clone() { return new PointCloud(*this); }

	virtual ~PointCloud();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// void SetPolyData(vtkPolyData* polyData);

	/// Returns the weights associated to the points.
	inline MatrixType GetPointWeights() const { return m_PointWeights; }

	///	Returns the type of the kernel.
	inline KernelEnumType GetKernelType() const { return m_KernelType; }
	/// Sets the type of the kernel to \e kernelType.
	inline void SetKernelType(KernelEnumType kernelType) { m_KernelType = kernelType; this->SetModified(); }

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e kernelWidth.
	inline void SetKernelWidth(TScalar h) {	m_KernelWidth = h; this->SetModified(); }

	/// Returns the RKHS-norm of itself.
	inline TScalar GetNormSquared() const { return m_NormSquared; }

	// STANLEY
	/*
	virtual int GetDimensionOfDiscretizedObject() const
	{
		int d = 0;
		for (unsigned int dim = 0; dim < Dimension; dim ++)
		{
			d += floor( (Superclass::Superclass::m_BoundingBox(dim,1) - Superclass::Superclass::m_BoundingBox(dim,0)) / m_KernelWidth );
		}
		//d *= Dimension;
		return d;
	}
	*/
	// PIETRO
	virtual int GetDimensionOfDiscretizedObject() const
	{
		int d = 1;
		for (unsigned int dim = 0; dim < Dimension; dim ++)
		{
			d *= floor( (Superclass::Superclass::m_BoundingBox(dim,1) - Superclass::Superclass::m_BoundingBox(dim,0)) / m_KernelWidth + 1.0);
		}
		d *= Dimension;
		return d;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void Update();

	/// See DeformableObject::ComputeMatch(DeformableObject* target) for details.
	virtual TScalar ComputeMatch(const DeformableObjectType* target);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the weights associated to the points.
	void UpdatePointWeights();

	/// Computes the RKHS-norm of itself.
	void UpdateSelfNorm();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Weights associated to the points (Size : NumberOfPoints x 1).
	MatrixType m_PointWeights;

	///	Type of the kernel.
	KernelEnumType m_KernelType;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the point cloud.
	TScalar m_NormSquared;


}; /* class PointCloud */


#ifndef MU_MANUAL_INSTANTIATION
#include "PointCloud.txx"
#endif


#endif /* _PointCloud_h */

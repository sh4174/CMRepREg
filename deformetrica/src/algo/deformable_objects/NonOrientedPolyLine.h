/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedPolyLine_h
#define _NonOrientedPolyLine_h

#include "Landmark.h"

#include "KernelType.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      Non oriented curves.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The NonOrientedPolyLine class inherited from Landmark represents a set of polygonal lines.
 *	            This class uses the varifold representation of curves, which is insensitive to the local orientation of the curves.
 */

template <class TScalar, unsigned int Dimension>
class NonOrientedPolyLine : public Landmark<TScalar, Dimension>
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

	NonOrientedPolyLine();
	/// Copy constructor.
	NonOrientedPolyLine(const NonOrientedPolyLine& other);
	/// Constructor which copies the object and update vertex coordinates
	NonOrientedPolyLine(const NonOrientedPolyLine& example, const MatrixType& LandmarkPoints);
	
	/// Returns a Deformed version of the mesh, where the deformation is given by the position of vertices in \e LandmarkPoints
	virtual NonOrientedPolyLine* DeformedObject(const MatrixType& LandmarkPoints)
		{ return new NonOrientedPolyLine(*this, LandmarkPoints); }

	virtual NonOrientedPolyLine* Clone() { return new NonOrientedPolyLine(*this); }

	virtual ~NonOrientedPolyLine();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void SetPolyData(vtkPolyData* polyData);

	/// Returns the centers of the cells.
	inline MatrixType GetCenters() const { return m_Centers; }

	/// Returns the tangents of the cells.
	inline MatrixType GetTangents() const { return m_Tangents; }

	/// Returns the matrices of the type \f[ \frac{\tau_i \tau_j^T}{\left|\tau_i\right|^{1/2}\left|\tau_j\right|^{1/2}} \f],
	/// where \f$ \tau_i \f$ denotes the tangents.
	inline MatrixType GetMatrixTangents() const { return m_MatrixTangents; }

	/// Returns the number of cells.
	inline int GetNumberOfCells() const { return m_NumCells; }

	///	Returns the type of the kernel.
	inline KernelEnumType GetKernelType() const { return m_KernelType; }
	/// Sets the type of the kernel to \e kernelType.
	inline void SetKernelType(KernelEnumType kernelType) { m_KernelType = kernelType; this->SetModified(); }

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e h.
	inline void SetKernelWidth(TScalar h) {	m_KernelWidth = h; this->SetModified(); }

	/// Returns the squared RKHS-norm of itself.
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
		d *= Dimension;
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

	virtual void Update();

	/// See DeformableObject::ComputeMatch(DeformableObject* target) for details.
	virtual TScalar ComputeMatch(const DeformableObjectType* target);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the centers and the tangents from the points.
	void UpdateCentersTangents();

	/// Updates the bounding box.
	void UpdateBoundingBox();

	/*
	 *  \brief      Computes the centers and tangents from the points.
	 *
	 *  \details    Given a set of polygonal lines , this method computes from the vertices of the curve
	 *              the centers and the tangents of each cell. The tangents are normalised by the squared
	 *              root of the their norm.
	 *
	 *	\param[in]  Pts             The vertices of the curves.
	 *  \param[out] Centers         The centers of the cells (Size : NumCells x Dimension).
	 *  \param[out] Tangents        The tangents of the cells (Size : NumCells x Dimension).
	 *  \param[out] MatrixTangents  A matrix of vectors containing the products between all the couples
	 *                              of elements of the Tangents. Given Dimension=3 and being
	 *                              \f$ \alpha_i=[\alpha_{i1} \alpha_{i2} \alpha_{i3}] \f$ the tangent vector of
	 *                              cell i, MatrixTangents(i)= \f$ [\alpha_{i1}^2  \alpha_{i1}\alpha_{i2}
	 *                              \alpha_{i1}\alpha_{i3} \alpha_{i2}^2 \alpha_{i2}\alpha_{i3} \alpha_{i3}^2] \f$
	 *                              (Size : NumCells x (Dimension*(Dimension+1)/2)).
	 */

	/// Computes the RKHS-norm of itself in the framework of varifolds.
	void UpdateSelfNorm();

	/**
	 *  \brief      It transforms a vector of size (Dimension*(Dimension+1)/2)) in a symmetric matrix of
	 *              size Dimension x Dimension by considering the values of the vector as the upper
	 *              triangular part of the matrix. It then multiplies this matrix by the vector X.
	 *
	 *  \details    It transforms every row Mi of the MatrixTangents M in a symmetric matrix of size
	 *              Dimension x Dimension and it multiplies it by the tangent of a cell X.
	 *
	 *  \param[in]  Mi  Row of the MatrixTangents M  (Size : 1 x (Dimension*(Dimension+1)/2)).
	 *  \param[in]  X   Row of m_Tangents (Size : 1 x Dimension).
	 *
	 *  \return     Vector of size 1 x Dimension that can be used in the function dot_product as in OrientedPolyLine.
	 */
	VectorType special_product(VectorType Mi, VectorType X);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	MatrixType m_Centers;

	///	Matrix coordinates of the tangents of the cells  (Size : NumCells x Dimension).
	MatrixType m_Tangents;

	/// The matrix \f[ \frac{\tau_i \tau_j^T}{\left|\tau_i\right|^{1/2}\left|\tau_j\right|^{1/2}} \f], where \f$ \tau_i \f$
	/// denotes the tangents. Only the upper triangular part is stored (Size : NumCells x (Dimension*(Dimension+1)/2)).
	MatrixType m_MatrixTangents;

	/// Number of cells (i.e. polygonal lines).
	int m_NumCells;

	///	Type of the kernel.
	KernelEnumType m_KernelType;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented curve.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;



}; /* class NonOrientedPolyLine */


#ifndef MU_MANUAL_INSTANTIATION
#include "NonOrientedPolyLine.txx"
#endif

#endif /* _NonOrientedPolyLine_h */

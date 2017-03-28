/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _CMRep_h
#define _CMRep_h

#include "Landmark.h"
#include "NonOrientedSurfaceMesh.h"

#include "KernelType.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <cstring>
#include <iostream>
#include <sstream>


/**
 *  \brief      Continuous Medial Representation - Medial Surface and Reconstructed Shape Boundary
 *
 *  \copyright  New York University
 *  \version    
 *
 *  \details    CMRep
 *				This class uses the varifold representation of reconstructed surfaces, 
 *				which is insensitive to the local orientation of the mesh.
 */

template <class TScalar, unsigned int Dimension>
class CMRep : public Landmark<TScalar, Dimension>
{
public:
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Landmark type.
	typedef Landmark<TScalar, Dimension> Superclass;

	/// Target Deformable Object Type
	typedef NonOrientedSurfaceMesh<TScalar, Dimension> NonOrientedSurfaceMeshType;

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

	CMRep();
	/// Copy constructor.
	CMRep(const CMRep& other);
	/// Constructor which copies the object and update vertex coordinates
	CMRep(const CMRep& example, const MatrixType& LandmarkPoints);
	
	/// Returns a Deformed version of the mesh, where the deformation is given by the position of vertices in \e LandmarkPoints
	virtual CMRep* DeformedObject(const MatrixType& LandmarkPoints) {
		return new CMRep(*this, LandmarkPoints);
	}

	virtual CMRep* Clone() { return new CMRep(*this); }
	virtual ~CMRep();

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// void SetPolyData(vtkPolyData* polyData);

	/// Returns the centers of the cells.
	inline MatrixType GetCenters() const { return m_Centers; }

	/// Returns the normals of the cells.
	inline MatrixType GetNormals() const { return m_Normals; }

	/// Returns the matrices of the type \f[ \frac{n_i n_j^T}{\left|n_i\right|^{1/2}\left[n_j\right|^{1/2}} \f],
	/// where \f$ n_i \f$ denotes the normals.
	inline MatrixType GetMatrixNormals() const { return m_MatrixNormals; }

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
	// virtual void LoadRadiusFunction(const DeformableObjectType* target);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);

	/// Reconstruct Boundary from  CMRep
	void InitialReconstructBoundary();
	void ReconstructBoundary();
		

protected:
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the bounding box.
	virtual void UpdateBoundingBox();

	/// Updates the centers and the normals from the points.
	void UpdateCentersNormals();

	/*
	 *	\brief		Computes the centers and the normals from the points.
	 *
	 *	\details	Given a triangular mesh, this method computes from the vertices of the mesh
	 *				the centers and the normals of each face.
	 *
	 *	\param[in]	Pts				The vertices of the mesh.
	 *	\param[out]	Centers			The centers of the cells (Size : NumCells x Dimension).
	 *	\param[out]	Normals			The normals of the cells (Size : NumCells x Dimension).
	 *	\param[out]	MatrixNormals	A matrix of vectors containing the products between all the couples
	 *								of elements of the Tangents. Given Dimension=3 and being
	 *								\f$ \alpha_i=[\alpha_{i1} \alpha_{i2} \alpha_{i3}] \f$ the tangent vector of
	 *								cell i, MatrixTangents(i)= \f$ [\alpha_{i1}^2  \alpha_{i1}\alpha_{i2}
	 *								\alpha_{i1}\alpha_{i3} \alpha_{i2}^2 \alpha_{i2}\alpha_{i3} \alpha_{i3}^2] \f$
	 *								(Size : NumCells x (Dimension*(Dimension+1)/2))
	 */
	//void ComputeCentersNormals(MatrixType& Pts, MatrixType& Centers, MatrixType& Normals, MatrixType& MatrixNormals);

	/// Computes the RKHS-norm of itself.
	void UpdateSelfNorm();

	/**
	 *  \brief      It transforms a vector of size (Dimension*(Dimension+1)/2)) in a symmetric matrix of
	 *              size Dimension x Dimension by considering the values of the vector as the upper
	 *              triangular part of the matrix. It then multiplies this matrix by the vector X.
	 *
	 *  \details    It transforms every row Mi of the m_MatrixNormals M in a symmetric matrix of size
	 *	            Dimension x Dimension and it multiplies it by the normal of a cell X.
	 *
	 *  \param[in]  Mi  Row of the m_MatrixNormals M  (Size : 1 x (Dimension*(Dimension+1)/2)).
	 *  \param[in]  X   Row of m_Normals (Size : 1 x Dimension).
	 *
	 *  \return     Vector of size 1 x Dimension that can be used in the function dot_product as in OrientedSurfaceMesh.
	 */
	VectorType special_product(VectorType Mi, VectorType X);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	MatrixType m_Centers;

	/// Matrix coordinates of the normals of the cells  (Size : NumCells x Dimension).
	MatrixType m_Normals;

	/// The matrix \f[ \frac{n_i n_j^T}{\left|n_i\right|^{1/2}\left|n_j\right|^{1/2}} \f], where \f$ n_i \f$ denotes the normals.
	/// Only the upper triangular part is stored (Size : NumCells x (Dimension*(Dimension+1)/2)).
	MatrixType m_MatrixNormals;

	/// Number of cells (i.e. triangles) of the mesh.
	int m_NumCells;

	/// Type of the kernel.
	KernelEnumType m_KernelType;

	/// Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented surface.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;

	// Reconstructed Shape Boundary 
	vtkSmartPointer< vtkPolyData > m_ReconBndrPoly;
	NonOrientedSurfaceMeshType* m_ReconBndr;


}; /* class CMRep */


#ifndef MU_MANUAL_INSTANTIATION
#include "CMRep.txx"
#endif


#endif /* _CMRep_h */

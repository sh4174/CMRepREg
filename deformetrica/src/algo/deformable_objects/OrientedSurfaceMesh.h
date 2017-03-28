/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedSurfaceMesh_h
#define _OrientedSurfaceMesh_h

#include "Landmark.h"

#include "KernelType.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      Oriented surface meshes.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *  \details    The OrientedSurfaceMesh class inherited from Landmark represents a triangular mesh.
 *              This class uses the current representation of the surface, which is sensitive to the local orientation
 *              of the mesh (a flip in ordering of the vertex indices in the connectivity matrix changes the sign
 *              of the surface element (i.e. the triangle) in the space of currents).
 */
template <class TScalar, unsigned int Dimension>
class OrientedSurfaceMesh : public Landmark<TScalar, Dimension>
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

	OrientedSurfaceMesh();
	/// Copy constructor.
	OrientedSurfaceMesh(const OrientedSurfaceMesh& other);
	/// Constructor which copies the object and update vertex coordinates
	OrientedSurfaceMesh(const OrientedSurfaceMesh& example, const MatrixType& LandmarkPoints);
	
	/// Returns a Deformed version of the mesh, where the deformation is given by the position of vertices in \e LandmarkPoints
	virtual OrientedSurfaceMesh* DeformedObject(const MatrixType& LandmarkPoints) {
		return new OrientedSurfaceMesh(*this, LandmarkPoints);
	}

	virtual OrientedSurfaceMesh* Clone() { return new OrientedSurfaceMesh(*this); }

	virtual ~OrientedSurfaceMesh();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Uses VTK filter to consistently re-orient normals for genus 0 surfaces.
	inline void SetReorient() { m_Reorient = true; this->SetModified(); }

	/// Unsets the use of VTK filter to automatically re-orient surface normals.
	inline void UnSetReorient() { m_Reorient = false; this->SetModified(); }

	/// Returns the centers of the cells.
	inline MatrixType GetCenters() const { return m_Centers; }

	/// Returns the normals of the cells.
	inline MatrixType GetNormals() const { return m_Normals; }

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

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the bounding box.
	void UpdateBoundingBox();

	/// Possibly reorient normals according to vtk filters.
	/// \warning   Make sure the ordering of point cells is consistent with the direction of the normal.
	void CheckMeshAndNormals();

	/// Updates the centers and the normals from the points.
	void UpdateCentersNormals();

	/*
	 *	\brief		Computes the centers and the normals from the points.
	 *
	 *	\details	Given a triangular mesh, this method computes from the vertices of the mesh
	 *				the centers and the normals of each face.
	 *
	 *	\param[in]	Pts			The vertices of the mesh.
	 *	\param[out]	Centers		The centers of the cells (Size : NumCells x Dimension).
	 *	\param[out]	Normals		The normals of the cells (Size : NumCells x Dimension).
	 */
	//void ComputeCentersNormals(const MatrixType& Pts, MatrixType& Centers, MatrixType& Normals);

	/// Computes the RKHS-norm of itself.
	void UpdateSelfNorm();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	MatrixType m_Centers;

	///	Matrix coordinates of the normals of the cells  (Size : NumCells x Dimension).
	MatrixType m_Normals;

	/// Number of cells (i.e. triangles) of the mesh.
	int m_NumCells;

	/// true to use VTK filter to re-orient normals of genus-0 surfaces
	bool m_Reorient;

	///	Type of the kernel.
	KernelEnumType m_KernelType;
	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented surface.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;


}; /* class OrientedSurfaceMesh */


#ifndef MU_MANUAL_INSTANTIATION
#include "OrientedSurfaceMesh.txx"
#endif


#endif /* _OrientedSurfaceMesh_h */

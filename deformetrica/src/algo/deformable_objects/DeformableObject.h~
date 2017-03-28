/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObject_h
#define _DeformableObject_h

#include "LinearAlgebra.h"
#include <vector>

#include "AnatomicalCoordinateSystem.h"



/**
 *  \brief      Deformable objects.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The DeformableObject class encodes an object embedded in the current 2D or 3D space. The child classes contain
 *              the metric appearing in the fidelity-to-data term.\n \n
 *              See DeformableObject::DeformableObjectType for the list of available object types.
 */
template <class TScalar, unsigned int Dimension>
class DeformableObject
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Possible type of deformable object.
	typedef enum
	{
		null,			/*!< Null value. */
		SSDImage,		/*!< Image with linear interpolation of voxel values and Sum Of Squared Differences (SSD) Metric*/
		LCCImage,		/*!< Image with linear interpolation of voxel values and Local Correlation Coefficient (LCC) Metric*/
		EQLAImage,		/*!< Image with linear interpolation of voxel values and "Ecart Quadratique au modele Local Affine" (EQLA) Metric (variant of LCC metric)*/
		Landmark,		/*!< Landmark (see Landmark). */
		OrientedPolyLine,	/*!< Current representation of a curve (see OrientedPolyLine). */
		OrientedSurfaceMesh,	/*!< Current representation of a surface (see OrientedSurfaceMesh). */
		NonOrientedPolyLine,	/*!< Varifold representation of a curve (see NonOrientedPolyLine). */
		NonOrientedSurfaceMesh,	/*!< Varifold representation of a surface (see NonOrientedSurfaceMesh). */
		PointCloud		/*!< Point cloud (see PointCloud). */
	} DeformableObjectType;

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	/// List of matrices type.
	typedef std::vector<MatrixType> MatrixList;

	/// Anatomical Coordinate System type.
	typedef AnatomicalCoordinateSystem<TScalar, Dimension> AnatomicalCoordinateSystemType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObject();
	/// Copy constructor.
	DeformableObject(const DeformableObject& other);
	/// Copy constructor with an update of the landmark point coordinates or a resampling of the image
	/// according to the new voxel positions in \e LandmarkOrImagePoints
	DeformableObject(const DeformableObject& example, const MatrixType& LandmarkOrImagePoints);

	/// Returns a deformed version of the object
	virtual DeformableObject* DeformedObject(const MatrixType& LandmarkOrImagePoints) = 0;
	/// Makes a copy of the object.
	virtual DeformableObject* Clone() = 0;

	virtual ~DeformableObject();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

// inline TScalar GetTimePoint() { return m_Timepoint; }
// inline void SetTimePoint(TScalar t) { m_Timepoint = t; }

// inline unsigned int GetTimeIndex() { return m_TimeIndex; }
// inline void SetTimeIndex(unsigned int t) { m_TimeIndex = t; }

	/// Returns true if the parameters of the deformation have changed, false otherwise.
	inline bool IsModified() const { return m_Modified; }
	/// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
	inline void SetModified() { m_Modified = true; }
	/// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
	inline void UnSetModified() { m_Modified = false; }

	/// Return the bounding box of the data
	inline MatrixType GetBoundingBox() const { return m_BoundingBox; }
	
	/// Returns the type of the deformable object.
	inline DeformableObjectType GetType() const { return m_Type; }
	/// Returns true if the deformable object is of image kind, false otherwise.
	inline bool IsOfImageKind() const { return ( (m_Type == SSDImage) || (m_Type == LCCImage) || (m_Type == EQLAImage) ); }
	/// Returns true if the deformable object is of landmark kind, false otherwise.
	inline bool IsOfLandmarkKind() const { return ( (m_Type == Landmark) || (m_Type == PointCloud) ||
			(m_Type == OrientedPolyLine) || (m_Type == OrientedSurfaceMesh) ||
			(m_Type == NonOrientedPolyLine) || (m_Type == NonOrientedSurfaceMesh)  ); }
	/// Returns true if the deformable object is a mesh, false otherwise.
	inline bool IsMesh() const {return ( (m_Type == OrientedPolyLine) || (m_Type == NonOrientedPolyLine) ||
			(m_Type == OrientedSurfaceMesh) || (m_Type == NonOrientedSurfaceMesh) ); }
	/// Sets the type of the deformable object to SSDImage.
	inline void SetSSDImageType() { m_Type = SSDImage; }
	/// Sets the type of the deformable object to LCCImage.
	inline void SetLCCImageType() { m_Type = LCCImage; }
	/// Sets the type of the deformable object to EQLAImage.
	inline void SetEQLAImageType() { m_Type = EQLAImage; }
	/// Sets the type of the deformable object to Landmark.
	inline void SetLandmarkType() { m_Type = Landmark; }
	/// Sets the type of the deformable object to PointCloud.
	inline void SetPointCloudType() { m_Type = PointCloud; }
	/// Sets the type of the deformable object to OrientedPolyLine.
	inline void SetOrientedPolyLineType() { m_Type = OrientedPolyLine; }
	/// Sets the type of the deformable object to OrientedSurfaceMesh.
	inline void SetOrientedSurfaceMeshType() { m_Type = OrientedSurfaceMesh; }
	/// Sets the type of the deformable object to NonOrientedPolyLine.
	inline void SetNonOrientedPolyLineType() { m_Type = NonOrientedPolyLine; }
	/// Sets the type of the deformable object to NonOrientedSurfaceMesh.
	inline void SetNonOrientedSurfaceMeshType() { m_Type = NonOrientedSurfaceMesh;}

	/// Returns the number of points of the object.
	virtual int GetNumberOfPoints() const = 0;
	
	/// Returns the dimension of the discretization of the object, used in Bayesian atlas update.
	virtual int GetDimensionOfDiscretizedObject() const = 0;

	/// Returns the label of the anatomical orientation.
	inline std::string GetAnatomicalCoordinateSystemLabel() const { return m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel(); }
	/// Sets the label of the anatomical orientation to \e label.
	inline void SetAnatomicalCoordinateSystem(std::string label) { m_AnatomicalOrientation.SetAnatomicalCoordinateSystemLabel(label); }
	/// Returns the change of basis matrix of the anatomical orientation.
	inline MatrixType GetAnatomicalCoordinateSystemMatrix() const { return m_AnatomicalOrientation.GetChangeOfBasisMatrix(); }
	/// Sets the change of basis matrix of the anatomical orientation to \e m.
	inline void SetAnatomicalCoordinateSystem(MatrixType m) { m_AnatomicalOrientation.SetChangeOfBasisMatrix(m); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates several information concerning the deformable object (see details on child class).
	virtual void Update() = 0;

	/**
	 *  \brief      Returns the norm between itself and the target.
	 *
	 *  \details    Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
	 *              method computes the following norm :
	 *              \f[
	 *              \left\Vert S - T \right\Vert_{W}^2 ,
	 *              \f]
	 *              where the metric \f$ \left\Vert ~.~ \right\Vert_{W}^2 \f$ depends on the child class.
	 *
	 *  \param[in]  target	Target of same type as the deformable object.
	 *  \return     The norm of the difference.
	 */
	virtual TScalar ComputeMatch(const DeformableObject* target) = 0;

	/**
	 *  \brief      Returns the gradient of the norm between itself and the target.
	 *
	 *  \details    Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
	 *              method computes the following gradient :
	 *              \f[
	 *              \nabla\left\Vert S - T \right\Vert_{W}^2 ,
	 *              \f]
	 *              where the metric \f$ \left\Vert ~.~ \right\Vert_{W}^2 \f$ depends on the child class.
	 *
	 *  \param[in]  target	Target of same type as the deformable object.
	 *  \return	    The gradient of the norm of the difference.
	 */
	virtual MatrixType ComputeMatchGradient(const DeformableObject* target) = 0;

//	virtual void TransportAlongGeodesicWithJumps(MatrixList& gradDataTi,
//			std::vector<int>& jumpTimes, MatrixList& EtaT, MatrixList& YT) = 0;

// curvature-like measure for curves and surfaces only
// virtual TScalar ComputeSmoothness() = 0;
// virtual MatrixType ComputeSmoothnessGradient() = 0;
// virtual void Smooth() = 0;

	/// Saves in \e filename the deformable object (for Landmark type, it is saved in *.vtk format (VTK PolyData)).
	virtual void WriteObject(std::string filename) const = 0;
	virtual void WriteObject(std::string filename, const MatrixList& velocity) const
	{
		// No default implementation
	}


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Computes the union of two bounding box.
	MatrixType UnionBoundingBox(MatrixType BB1, MatrixType BB2);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Type of the deformable object.
	DeformableObjectType m_Type;

	/// Boolean which avoids re-updating computation of centers, tangents and self-norm
	bool m_Modified;
	
	/// Tight Bounding box, which contain the data (i.e. centers of mesh cells and centers of voxels)
	MatrixType m_BoundingBox;

	/// Anatomical coordinate system associated to the deformable object.
	AnatomicalCoordinateSystemType m_AnatomicalOrientation;


// TScalar m_Timepoint;

// unsigned int m_TimeIndex;


}; /* class DeformableObject */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableObject.txx"
#endif


#endif /* _DeformableObject_h */

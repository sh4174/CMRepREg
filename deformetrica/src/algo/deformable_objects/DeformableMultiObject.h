/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableMultiObject_h
#define _DeformableMultiObject_h

#include "DeformableObject.h"
#include "itkImage.h"

#include "LinearAlgebra.h"
#include <vector>

/**
 *  \brief      Collections of deformable objects (i.e. multi-object)
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The DeformableMultiObject class is used to deal with collections of deformable objects
 *              embedded in the current 2D or 3D space. It extends for such collections the methods
 *              of the deformable object class that are designed for a single object at a time.
 */
template <class TScalar, unsigned int Dimension>
class DeformableMultiObject
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	/// List of matrices type.
	typedef std::vector<MatrixType> MatrixList;

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	/// List of deformable objects type.
	typedef std::vector<DeformableObjectType*> DeformableObjectList;

	/// ITK image type.
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableMultiObject();
	/// Copy constructor.
	DeformableMultiObject(const DeformableMultiObject& other);
	/// Contructor which copies the \e Example and sets the coordinates of landmark points to \e LandmarkPoints
	/// and re-sample to image according to voxels positions in \e ImagePoints. It is used to get
	/// the deformed version of a DeformableMultiObject.
	DeformableMultiObject(const DeformableMultiObject& Example, const MatrixType& LandmarkPoints,
			const MatrixType& DownSampledImageMap);
	
	/// Returns a deformed version of the object.
	DeformableMultiObject* DeformedMultiObject(MatrixType& LandmarkPoints, MatrixType& ImagePoints) {
		return new DeformableMultiObject(*this, LandmarkPoints, ImagePoints);
	}

	/// Makes a copy of the object.
	DeformableMultiObject* Clone() { return new DeformableMultiObject(*this); }

	~DeformableMultiObject();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Return the list of the deformable objects.
	inline DeformableObjectList GetObjectList() const { return m_ObjectList; }
	inline DeformableObjectList GetCMRepObjectList() const { return m_CMRepObjectList; }
	
	/// Set the list of deformable objects as a copy of \e list.
	void SetObjectList(DeformableObjectList& list);
	void SetCMRepObjectList(DeformableObjectList& list);

	/// Returns the coordinates of the vertices for Landmark type (or children object).
	MatrixType GetLandmarkPoints() const;
	/// Returns the voxel positions of the downsample image in LinearInterpImage type objects (or children object).
	MatrixType GetDownSampledImageMap() const;
	/// Returns an example of downsampled image, whose voxel coordinates are given in this->GetDownSampledImageMap()
	ImageTypePointer GetDownSampledImage() const;
	/// Returns image at full resolution
	ImageTypePointer GetImage() const;
	
	// /// Set new coordinates of all landmark type objects (used mostly to compute a deformed object)
	// void SetLandmarkPoints(const MatrixType& Y);
	// /// Set new positions of the voxel in the image domain (used mostly to compute a deformed image)
	// void SetImagePoints(const MatrixType& Y);

	/// Return the coordinates of the vertices of the landmark-type objects and the intensities of the image object, if any (used mostly for updating template parameters)
	MatrixList GetImageIntensityAndLandmarkPointCoordinates() const;
	/// Set new vertex coordinates in landmark type objects and new image intensities in the image type object (used mostly for updating template parameters)
	void UpdateImageIntensityAndLandmarkPointCoordinates(const MatrixList & Y);

	/// Return the index of the image in the list of objects (returns -1 if there is no image).
	inline int GetImageIndex() const { return m_ImageIndex; }

	/// Return the number of objects of image kind.
	inline int GetNumberOfImageKindObjects() const { return m_NumberOfImageKindObjects; }
	/// Return the number of objects of landmark kind.
	inline int GetNumberOfLandmarkKindObjects() const { return m_NumberOfLandmarkKindObjects; } 
	/// Return the total number of objects.
	inline int GetNumberOfObjects() const { return m_NumberOfObjects; }

// TScalar GetTimePoint();

// unsigned int GetTimeIndex();

	/// Returns the list of each number of points associated to the deformable objects.
	inline std::vector<int> GetNumberOfPoints() const { return m_NumberOfPoints; }
	/// Returns the number of points of all deformable objects.
	inline int GetTotalNumberOfPoints() const { return (m_NumberOfLandmarkPoints + m_NumberOfImagePoints); }

	/// Returns the number of landmark points
	inline int GetNumberOfLandmarkPoints() const { return m_NumberOfLandmarkPoints; }
	
	/// Returns the number of image points
	inline int GetNumberOfImagePoints() const { return m_NumberOfImagePoints; }
	
	/// Returns a tight bounding box enclosing all objects
	inline MatrixType GetBoundingBox() const { return m_BoundingBox; };
	/// Computes a tight bounding box enclosing all objects.
	void UpdateBoundingBox();
	
	/// Get dimension of discretized objects for each object
	std::vector<int> GetDimensionOfDiscretizedObjects() const;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates information about the deformable objects and updates the bounding box.
	void Update();

	/// Computes DeformableObject::ComputeMatch(DeformableObject* target) for each deformable object.
	std::vector<TScalar> ComputeMatch(const DeformableMultiObject* target);

	/// Computes DeformableObject::ComputeMatchGradient(DeformableObject* target) for each deformable object.
	MatrixList ComputeMatchGradient(const DeformableMultiObject* target);

	/// Transforms concatenated data located at landmark points and image points into a list
	void ListToMatrices(const MatrixList& L, MatrixType& MLandmark, MatrixType& MImage) const;

	/// Transforms list of objects into concatenated data of landmark types and image types
	void MatricesToList(MatrixList& L, const MatrixType& MLandmark, const MatrixType& MImage) const;
	
	/// Splats the difference between image in this and image in \e targ
	MatrixType SplatDifferenceImage(const DeformableMultiObject* targ, const MatrixType& DownSampledImageMap) const;

	/// Calls the method DeformableObject::WriteObject() for each deformable object.
	void WriteMultiObject(std::vector<std::string>& str) const;
	void WriteMultiObject(std::vector<std::string>& str, const MatrixList& velocity) const;



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Collection of m_NumberOfObjects deformable objects.
	DeformableObjectList m_ObjectList;
	DeformableObjectList m_CMRepObjectList;

	/// Total number of objects.
	int m_NumberOfObjects;

	/// Number of objects of image kind.
	int m_NumberOfImageKindObjects;

	/// Number of objects of landmark kind.
	int m_NumberOfLandmarkKindObjects;

	/// Vector containing at each cell the number of points of the associated deformable object.
	std::vector<int> m_NumberOfPoints;
	
	/// Total number of vertices in the multi-object
	int m_NumberOfLandmarkPoints;
	
	/// Number of pixels/voxels in the image
	int m_NumberOfImagePoints;
	
	/// Index of the object of image type (at most one image allowed within a multi-object).
	int m_ImageIndex;
	
	/// Bounding Box tightly enclosing all objects
	MatrixType m_BoundingBox;


}; /* class DeformableMultiObject */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableMultiObject.txx"
#endif


#endif /* _DeformableMultiObject_h */

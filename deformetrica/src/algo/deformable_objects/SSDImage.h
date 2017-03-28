/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SSDImage_h
#define _SSDImage_h

#include "LinearInterpImage.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      An image using linear interpolation of voxel values and SSD Metric between images
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The SSDImage class inherited from LinearInterpImage ... TODO .
 */
template<class TScalar, unsigned int Dimension>
class SSDImage : public LinearInterpImage<TScalar, Dimension>
{
	public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Abstract image type.
	typedef LinearInterpImage<TScalar, Dimension> Superclass;
	
	/// Deformable Object type.
	typedef typename Superclass::Superclass DeformableObjectType;
	
	/// ITK image type.
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;
	/// ITK image index type.
	typedef typename ImageType::IndexType ImageIndexType;
	/// ITK image point type.
	typedef typename ImageType::PointType ImagePointType;
	/// ITK image region type.
	typedef typename ImageType::RegionType ImageRegionType;
	/// ITK image size type.
	typedef typename ImageType::SizeType ImageSizeType;
	/// ITK image spacing type.
	typedef typename ImageType::SpacingType ImageSpacingType;

	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;
	/// List of matrices type.
	typedef typename Superclass::MatrixList MatrixList;




	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SSDImage();
	virtual ~SSDImage();

	/// Copy constructor.
	SSDImage(const SSDImage& other);
	/// Constructor which copies the object and resample the image
	SSDImage(const SSDImage& example, const MatrixType& ImagePoints);


	/// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints
	virtual SSDImage* DeformedObject(const MatrixType& ImagePoints) { return new SSDImage(*this, ImagePoints); }
	
	virtual SSDImage* Clone() {	return new SSDImage(*this); }
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void Update();

	/// See DeformableObject::ComputeMatch(DeformableObject* target) for details.
	virtual TScalar ComputeMatch(const DeformableObjectType* target);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);

	/// Return the dimension of the discretized image, here the number of voxels of the original image
	virtual int GetDimensionOfDiscretizedObject() const { return Superclass::m_NumberOfVoxels; }


};


#ifndef MU_MANUAL_INSTANTIATION
#include "SSDImage.txx"
#endif


#endif /* _SSDImage_h */

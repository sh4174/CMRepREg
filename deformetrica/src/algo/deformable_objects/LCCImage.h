/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LCCImage_h
#define _LCCImage_h

#include "LinearInterpImage.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      An image using linear interpolation of voxel values and LCC Metric between images.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The LCCImage class inherited from LinearInterpImage ... TODO .
 */
template<class TScalar, unsigned int Dimension>
class LCCImage : public LinearInterpImage<TScalar, Dimension>
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

	LCCImage();
	virtual ~LCCImage();

	/// Copy constructor.
	LCCImage(const LCCImage& other);
	/// Constructor which copies the object and resample the image
	LCCImage(const LCCImage& example, const MatrixType& ImagePoints);

	/// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints
	virtual LCCImage* DeformedObject(const MatrixType& ImagePoints) { return new LCCImage(*this, ImagePoints); }
	
	virtual LCCImage* Clone() {	return new LCCImage(*this); }
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/// Set standard deviation of the Gaussian filter needed to compute local statistics
	void SetLCCKernelWidth(TScalar h) { m_LCCKernelWidth = h; }
	
	/// Get Local Mean Image
	ImageTypePointer GetLocalMeanImage() const {return m_LocalMeanImage; }
	/// Get Local Standard Deviation Image
	ImageTypePointer GetLocalVarianceImage() const { return m_LocalVarianceImage; }

	virtual int GetDimensionOfDiscretizedObject() const { throw std::runtime_error("No noise model for LCC Image!"); return 0; }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void Update();
	
	/// Compute LCC between this image (after deformation) and target image
	virtual TScalar ComputeMatch(const DeformableObjectType* target);

	/// Compute the gradient of the LCC metric between this image (once deformed) and target image
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Computes the local mean image from input image, namely the convolution between input image and a gaussian filter: W*I
	ImageTypePointer ComputeLocalMeanImage(const ImageType* img) const;

	/// Computes the local covariance image from 2 input images (requires to have Local Mean Image computed): W*(I.J) - (W*I)(W*J)
	ImageTypePointer ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
			const ImageType* img2, const ImageType* LocalMeanImg2, bool sameImage) const;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Image of local mean (convolution of the image I with kernel W: W*I)
	ImageTypePointer m_LocalMeanImage;
	
	/// Image of local standard deviation: sqrt( (W*(I^2)) )
	ImageTypePointer m_LocalVarianceImage;


	/// Standard deviation of the Gaussian function used to compute the local mean and correlation
	TScalar m_LCCKernelWidth;
	
	/// Minimum image size
	// int m_SmallestImageDimension;
	
	
};


#ifndef MU_MANUAL_INSTANTIATION
#include "LCCImage.txx"
#endif


#endif /* _LCCImage_h */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _EQLAImage_h
#define _EQLAImage_h

#include "LinearInterpImage.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      An image using linear interpolation of voxel values and EQLA Metric between images.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The EQLAImage class inherited from LinearInterpImage... TODO .
 */
template<class TScalar, unsigned int Dimension>
class EQLAImage : public LinearInterpImage<TScalar, Dimension>
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

	EQLAImage();
	virtual ~EQLAImage();

	/// Copy constructor.
	EQLAImage(const EQLAImage& other);
	/// Constructor which copies the object and resample the image
	EQLAImage(const EQLAImage& example, const MatrixType& ImagePoints);

	/// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints
	virtual EQLAImage* DeformedObject(const MatrixType& ImagePoints) { return new EQLAImage(*this, ImagePoints); }

	virtual EQLAImage* Clone() { return new EQLAImage(*this); }
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/// Sets input image and compute local mean image and local standard deviation image.
	// void SetImage(ImageType* img);
	
	/// Sets standard deviation of the Gaussian filter needed to compute local statistics.
	void SetEQLAKernelWidth(TScalar h) { m_EQLAKernelWidth = h; }
	
	/// Returns local mean image.
	ImageTypePointer GetLocalMeanImage() const { return m_LocalMeanImage; }
	/// Returns the local standard deviation image.
	ImageTypePointer GetLocalVarianceImage() const { return m_LocalVarianceImage; }

	virtual int GetDimensionOfDiscretizedObject() const { throw std::runtime_error("No noise model for EQLA Image!"); return 0; }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
	virtual void Update();

	/// Computes EQLA between this image (after deformation) and target image.
	virtual TScalar ComputeMatch(const DeformableObjectType* target);

	/// Computes the gradient of the EQLA metric between this image (once deformed) and target image.
	virtual MatrixType ComputeMatchGradient(const DeformableObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Computes the local mean image from input image, namely the convolution between input image and a gaussian filter: W*I.
	ImageTypePointer ComputeLocalMeanImage(const ImageType* img) const;

	/// Computes the local covariance image from 2 input images (requires to have Local Mean Image computed): W*(I.J) - (W*I)(W*J)
	ImageTypePointer ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
			const ImageType* img2, const ImageType* LocalMeanImg2) const;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Image of local mean (convolution of the image I with kernel W: W*I).
	ImageTypePointer m_LocalMeanImage;
	
	/// Image of local standard deviation: sqrt( (W*(I^2)) ).
	ImageTypePointer m_LocalVarianceImage;

	/// Standard deviation of the Gaussian function used to compute the local mean and correlation.
	TScalar m_EQLAKernelWidth;
	
	/// Minimum image size.
	// int m_SmallestImageDimension;
	

};


#ifndef MU_MANUAL_INSTANTIATION
#include "EQLAImage.txx"
#endif


#endif /* _EQLAImage_h */

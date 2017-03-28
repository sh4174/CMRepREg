/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _GridFunctions_h
#define _GridFunctions_h

#include "itkImage.h"

#include "LinearAlgebra.h"
#include <vector>


/**
 *  \brief      TODO .
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The GridFunctions class TODO .
 */
template <class TScalar, unsigned int Dimension>
class GridFunctions
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// VectorType type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;

	/// TODO .
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// TODO .
	typedef typename ImageType::Pointer ImagePointer;
	/// TODO .
	typedef typename ImageType::IndexType ImageIndexType;
	/// TODO .
	typedef typename ImageType::PointType ImagePointType;
	/// TODO .
	typedef typename ImageType::RegionType ImageRegionType;
	/// TODO .
	typedef typename ImageType::SizeType ImageSizeType;
	/// TODO .
	typedef typename ImageType::SpacingType ImageSpacingType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Convert image to set of values, ordered using standard ITK traversal order
	static VectorType VectorizeImage(const ImageType* img);

	/// Generate matrix of point locations from image grid.
	static MatrixType ImageToPoints(const ImageType* img);

	/// Returns a pointer on a (itk) image based on \e img whose the intensity of the voxels are given by \e values.
	static ImagePointer VectorToImage(const ImageType* img, const VectorType& values);

	/// Interpolate values from image at location \e pos.
	static VectorType Interpolate(const MatrixType& pos, const ImageType* values);

	/// TODO .
	static ImagePointer SplatToImage(const ImageType* example, const MatrixType& pos, const VectorType& values, long gridPadding=2);

	/// Downsample image to a size roughly old size / factor
	static ImagePointer DownsampleImage(const ImageType* img, TScalar factor);

	/// Upsample \e img to the space defined by \e exImg.
	static ImagePointer UpsampleImage(const ImageType* exImg, const ImageType* img);

	/// TODO .
	static MatrixType UpsampleImagePoints(const ImageType* img, const ImageType* downSampledImg, const MatrixType& downSampledPos);

	/// Virtual method to avoid instanciation of the class GridFunctions.
	virtual void Abstract() = 0;



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	static inline void _getInterpolationWeightsAndGridPoints(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x,
			const ImageType* img);

	/// Specialized linear interpolation weights for 2D.
	static inline void _getInterpolationWeightsAndGridPoints_2(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x,
			const ImageType* img);

	/// Specialized linear interpolation weights for 3D.
	static inline void _getInterpolationWeightsAndGridPoints_3(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x,
			const ImageType* img);


}; /* class GridFunctions */


#ifndef MU_MANUAL_INSTANTIATION
#include "GridFunctions.txx"
#endif


#endif /* _GridFunctions_h */

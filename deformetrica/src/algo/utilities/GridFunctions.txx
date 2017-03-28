/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _GridFunctions_txx
#define _GridFunctions_txx

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkResampleImageFilter.h"



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::MatrixType
GridFunctions<TScalar, Dimension>
::ImageToPoints(const ImageType* img)
 {
	MatrixType pointsM(
			img->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension, 0);

	long r = 0;

	typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
	IteratorType it(img, img->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		ImageIndexType ind = it.GetIndex();

		ImagePointType p;
		img->TransformIndexToPhysicalPoint(ind, p);

		for (long dim = 0; dim < Dimension; dim++)
			pointsM(r, dim) = p[dim];
		r++;
	}

	return pointsM;
 }



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::VectorType
GridFunctions<TScalar, Dimension>
::VectorizeImage(const ImageType* img)
 {
	VectorType values(img->GetLargestPossibleRegion().GetNumberOfPixels(), 0);

	long r = 0;

	typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
	IteratorType it(img, img->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		values[r] = it.Get();

		r++;
	}

	return values;
 }



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::ImagePointer
GridFunctions<TScalar, Dimension>
::VectorToImage(const ImageType* imgEx, const VectorType& values)
 {
	if (imgEx->GetLargestPossibleRegion().GetNumberOfPixels() != values.size())
		throw std::runtime_error("Cannot set image voxels values: vector dimension mismatch");

	ImagePointer img = ImageType::New();
	img->SetRegions(imgEx->GetLargestPossibleRegion());
	img->CopyInformation(imgEx);
	img->Allocate();
	img->FillBuffer(0);

	typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
	IteratorType it(img, img->GetLargestPossibleRegion());

	unsigned int r = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		it.Set(values[r]);

		r++;
	}

	return img;
 }



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::VectorType
GridFunctions<TScalar, Dimension>
::Interpolate(
		const MatrixType& X, const ImageType* img)
		{
	long numPoints = X.rows();

	VectorType values(numPoints, 0.0);

	std::vector<TScalar> weights;
	std::vector<ImageIndexType> gridIndices;

	for (unsigned int i = 0; i < numPoints; i++)
	{
		_getInterpolationWeightsAndGridPoints(
				weights, gridIndices, X.get_row(i), img);

		TScalar v = 0;

		for (unsigned int j = 0; j < weights.size(); j++)
		{
			TScalar w = weights[j];
			ImageIndexType ind = gridIndices[j];
			v += w*img->GetPixel(ind);
		}

		values[i] = static_cast<TScalar>(v);
	}

	return values;
		}



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::ImagePointer
GridFunctions<TScalar, Dimension>
::SplatToImage(const ImageType* exImg, const MatrixType& X, const VectorType& values, long gridPadding)
 {
	long numPoints = X.rows();

	ImagePointer img = ImageType::New();
	img->SetRegions(exImg->GetLargestPossibleRegion());
	img->CopyInformation(exImg);
	img->Allocate();
	img->FillBuffer(0);

	ImageSizeType size =  exImg->GetLargestPossibleRegion().GetSize();

	for (unsigned int i = 0; i < numPoints; i++)
	{
		std::vector<TScalar> weights;
		std::vector<ImageIndexType> gridIndices;

		_getInterpolationWeightsAndGridPoints(weights, gridIndices, X.get_row(i), img);

		for (unsigned int j = 0; j < weights.size(); j++)
		{
			TScalar w = weights[j];
			ImageIndexType ind = gridIndices[j];

			bool isborder = false;
			for (unsigned int d = 0; d < Dimension; d++)
				if (ind[d] <= gridPadding || ind[d] >= ((long)size[d]-gridPadding))
					isborder = true;
			if (isborder)
				continue;

			img->SetPixel(ind, img->GetPixel(ind) + w*values[i]);
		}
	}

	/*
	ImageSpacingType spacing = img->GetSpacing();

	TScalar minSpacing = spacing[0];
	for (unsigned int dim = 1; dim < Dimension; dim++)
	if (spacing[dim] < minSpacing)
	  minSpacing = spacing[dim];

	typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
	GaussFilterType;
	typename GaussFilterType::Pointer gaussf = GaussFilterType::New();
	gaussf->SetInput(img);
	gaussf->SetVariance(minSpacing*minSpacing * 0.5);
	gaussf->Update();

	img = gaussf->GetOutput();
	 */

	return img;
 }



template <class TScalar, unsigned int Dimension>
void
GridFunctions<TScalar, Dimension>
::_getInterpolationWeightsAndGridPoints(
		std::vector<TScalar>& weights,
		std::vector<ImageIndexType>& gridIndices,
		const VectorType& x,
		const ImageType* img)
		{
	if (Dimension == 2)
	{
		// Bilinear
		_getInterpolationWeightsAndGridPoints_2(weights, gridIndices, x, img);
	}
	else if (Dimension == 3)
	{
		// Trilinear
		_getInterpolationWeightsAndGridPoints_3(weights, gridIndices, x, img);
	}
	else
	{
		// Nearest neighbor
		ImagePointType p;
		for (unsigned int k = 0; k < Dimension; k++)
			p[k] = x[k];

		ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

		// we need to make sure that this method well behaves when point coordinates are negative (floor instead of int conversion)
		ImageIndexType ind;
		img->TransformPhysicalPointToIndex(p, ind);

		weights.clear();
		gridIndices.clear();

		bool isOut = false;
		for (unsigned int k = 0; k < Dimension; k++)
			if (ind[k] < 0 || ind[k] >= (long)size[k])
			{
				isOut = true;
				break;
			}

		if (!isOut)
		{
			weights.push_back(1.0);
			gridIndices.push_back(ind);
		}
	}
		}

// template <class TScalar, unsigned int Dimension>
// void
// GridFunctions<TScalar, Dimension>
// ::_getInterpolationWeightsAndGridPoints_2(
//   std::vector<TScalar>& weights,
//   std::vector<ImageIndexType>& gridIndices,
//   const VectorType& p,
//   const ImageType* img)
// {
//   ImagePointType origin = img->GetOrigin();
//   ImageSpacingType spacing = img->GetSpacing();
// 
//   ImageSizeType size = img->GetLargestPossibleRegion().GetSize();
// 
//   std::vector<TScalar> contInd(Dimension);
//   for (unsigned int d = 0; d < Dimension; d++)
//     contInd[d] = (p[d]-origin[d]) / spacing[d];
// 
//   //
//   // Bilinear interpolation
//   //
// 
//   TScalar x = p[0];
//   TScalar y = p[1];
// 
//   // Get the 4 grid positions
// 
// 	// MARCEL
//   // int ix1 = (int)contInd[0];
//   // int iy1 = (int)contInd[1];
// 
// 	// STANLEY (problem when point coordinates become negative)
//   int ix1 = floor(contInd[0]);
//   int iy1 = floor(contInd[1]);
// 
// 
//   int ix2 = ix1 + 1;
//   int iy2 = iy1 + 1;
// 
//   TScalar x1 = ix1 * spacing[0] + origin[0];
//   TScalar y1 = iy1 * spacing[1] + origin[1];
// 
//   TScalar x2 = ix2 * spacing[0] + origin[0];
//   TScalar y2 = iy2 * spacing[1] + origin[1];
// 
//   TScalar fx = x - x1;
//   TScalar fy = y - y1;
// 
//   TScalar gx = x2 - x;
//   TScalar gy = y2 - y;
// 
//   TScalar A = spacing[0] * spacing[1];
// 
//   // Add valid grid positions and corresponding weights
//   weights.clear();
//   gridIndices.clear();
// 
// #define interpWeightMacro2(ix, iy, w) \
//   if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
//     (0 <= (iy)) && ((iy) < (int)size[1])) \
//   { \
//     ImageIndexType ind; \
//     ind[0] = (ix); ind[1] = (iy); \
//     if (w > 0) \
//     { \
//       gridIndices.push_back(ind); \
//       weights.push_back(w); \
//     } \
//   }
// 
//   interpWeightMacro2(ix1, iy1, gx*gy / A);
//   interpWeightMacro2(ix1, iy2, gx*fy / A);
//   interpWeightMacro2(ix2, iy1, fx*gy / A);
//   interpWeightMacro2(ix2, iy2, fx*fy / A);
// 
// #undef interpWeightMacro2
// 
// }



template <class TScalar, unsigned int Dimension>
void
GridFunctions<TScalar, Dimension>
::_getInterpolationWeightsAndGridPoints_2(
		std::vector<TScalar>& weights,
		std::vector<ImageIndexType>& gridIndices,
		const VectorType& p,
		const ImageType* img)
		{
	//  ImagePointType origin = img->GetOrigin();
	//  ImageSpacingType spacing = img->GetSpacing();

	ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

	// MARCEL
	// std::vector<TScalar> contInd(Dimension);
	// for (unsigned int d = 0; d < Dimension; d++)
	//   contInd[d] = (p[d]-origin[d]) / spacing[d];

	//STANLEY: solve problem when axes are flipped and transform contInd[d] = (p[d]-origin[d]) / spacing[d] is inaccurate
	typedef itk::Point<TScalar,Dimension> PointType;
	typedef itk::ContinuousIndex<TScalar,Dimension> ContIndexType;
	PointType pt;
	for (unsigned int d = 0; d < Dimension; d++)
		pt[d] = p[d];

	ContIndexType contInd; 
	img->TransformPhysicalPointToContinuousIndex(pt, contInd);


	//
	// Bilinear interpolation
	//

	// Get the 4 grid positions
	// STANLEY: floor and not int conversion, since contIndex may be negative (in case deformation field ends up outside image boundary for instance)
	// but since it may also have a neighborhing voxel *on* the boundary, one wants to keep it
	int ix1 = floor(contInd[0]);
	int iy1 = floor(contInd[1]);

	int ix2 = ix1 + 1;
	int iy2 = iy1 + 1;

	TScalar fx = contInd[0] - ix1;
	TScalar fy = contInd[1] - iy1;

	TScalar gx = ix2 - contInd[0];
	TScalar gy = iy2 - contInd[1];

	// Add valid grid positions and corresponding weights
	weights.clear();
	gridIndices.clear();

#define interpWeightMacro2(ix, iy, w) \
		if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
				(0 <= (iy)) && ((iy) < (int)size[1])) \
				{ \
			ImageIndexType ind; \
			ind[0] = (ix); ind[1] = (iy); \
			if (w > 0) \
			{ \
				gridIndices.push_back(ind); \
				weights.push_back(w); \
			} \
				}

	interpWeightMacro2(ix1, iy1, gx*gy);
	interpWeightMacro2(ix1, iy2, gx*fy);
	interpWeightMacro2(ix2, iy1, fx*gy);
	interpWeightMacro2(ix2, iy2, fx*fy);

#undef interpWeightMacro2

		}


// template <class TScalar, unsigned int Dimension>
// void
// GridFunctions<TScalar, Dimension>
// ::_getInterpolationWeightsAndGridPoints_3(
//   std::vector<TScalar>& weights,
//   std::vector<ImageIndexType>& gridIndices,
//   const VectorType& p,
//   const ImageType* img)
// {
//   ImagePointType origin = img->GetOrigin();
//   ImageSpacingType spacing = img->GetSpacing();
// 
//   ImageSizeType size = img->GetLargestPossibleRegion().GetSize();
// 
//   std::vector<TScalar> contInd(Dimension);
//   for (unsigned int d = 0; d < Dimension; d++)
//     contInd[d] = (p[d]-origin[d]) / spacing[d];
// 
//   //
//   // Trilinear interpolation
//   //
// 
//   TScalar x = p[0];
//   TScalar y = p[1];
//   TScalar z = p[2];
// 
//   // Get the 8 grid positions
// 	// MARCEL
//   // int ix1 = (int)contInd[0];
//   // int iy1 = (int)contInd[1];
//   // int iz1 = (int)contInd[2];
// 
// 	// STANLEY (problem when point coordinates become negative...)
//   int ix1 = floor(contInd[0]);
//   int iy1 = floor(contInd[1]);
//   int iz1 = floor(contInd[2]);
// 	
// 	
// 	
// 	
//   int ix2 = ix1 + 1;
//   int iy2 = iy1 + 1;
//   int iz2 = iz1 + 1;
// 
//   TScalar x1 = ix1 * spacing[0] + origin[0];
//   TScalar y1 = iy1 * spacing[1] + origin[1];
//   TScalar z1 = iz1 * spacing[2] + origin[2];
// 
//   TScalar x2 = ix2 * spacing[0] + origin[0];
//   TScalar y2 = iy2 * spacing[1] + origin[1];
//   TScalar z2 = iz2 * spacing[2] + origin[2];
// 
//   TScalar V = spacing[0] * spacing[1] * spacing[2];
// 
//   // Get distances to the image grid
//   TScalar fx = x - x1;
//   TScalar fy = y - y1;
//   TScalar fz = z - z1;
// 
//   TScalar gx = x2 - x;
//   TScalar gy = y2 - y;
//   TScalar gz = z2 - z;
// 
//   // Add valid grid positions and corresponding weights
//   weights.clear();
//   gridIndices.clear();
// 
// #define interpWeightMacro3(ix, iy, iz, w) \
//   if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
//     (0 <= (iy)) && ((iy) < (int)size[1]) && \
//     (0 <= (iz)) && ((iz) < (int)size[2])) \
//   { \
//     ImageIndexType ind; \
//     ind[0] = (ix); ind[1] = (iy); ind[2] = (iz); \
//     if (w > 0) \
//     { \
//       gridIndices.push_back(ind); \
//       weights.push_back(w); \
//     } \
//   }
// 
// 
// 	std::cout << ix1 << " " << " " << iy1 << " " << iz1 << " " << gx*gy*gz / V << std::endl;
// 
//   interpWeightMacro3(ix1, iy1, iz1, gx*gy*gz / V);
//   interpWeightMacro3(ix1, iy1, iz2, gx*gy*fz / V);
//   interpWeightMacro3(ix1, iy2, iz1, gx*fy*gz / V);
//   interpWeightMacro3(ix1, iy2, iz2, gx*fy*fz / V);
//   interpWeightMacro3(ix2, iy1, iz1, fx*gy*gz / V);
//   interpWeightMacro3(ix2, iy1, iz2, fx*gy*fz / V);
//   interpWeightMacro3(ix2, iy2, iz1, fx*fy*gz / V);
//   interpWeightMacro3(ix2, iy2, iz2, fx*fy*fz / V);
// 
// #undef interpWeightMacro3
// 
// 
// }


template <class TScalar, unsigned int Dimension>
void
GridFunctions<TScalar, Dimension>
::_getInterpolationWeightsAndGridPoints_3(
		std::vector<TScalar>& weights,
		std::vector<ImageIndexType>& gridIndices,
		const VectorType& p,
		const ImageType* img)
		{
	//  ImagePointType origin = img->GetOrigin();
	//  ImageSpacingType spacing = img->GetSpacing();

	ImageSizeType size = img->GetLargestPossibleRegion().GetSize();


	// MARCEL 
	//   std::vector<TScalar> contInd(Dimension);
	//   for (unsigned int d = 0; d < Dimension; d++)
	//   {
	// 	contInd[d] = (p[d]-origin[d]) / spacing[d];
	// }  

	//STANLEY: solve problem when axes are flipped and transform contInd[d] = (p[d]-origin[d]) / spacing[d] is inaccurate
	typedef itk::Point<TScalar,Dimension> PointType;
	typedef itk::ContinuousIndex<TScalar,Dimension> ContIndexType;
	PointType pt;
	for (unsigned int d = 0; d < Dimension; d++)
		pt[d] = p[d];

	ContIndexType contInd; 
	img->TransformPhysicalPointToContinuousIndex(pt, contInd);

	//
	// Trilinear interpolation
	//

	// Get the 8 grid positions
	// STANLEY: floor and not int conversion, since contIndex may be negative (in case deformation field ends up outside image boundary for instance)
	// but since it may also have a neighborhing voxel *on* the boundary, one wants to keep it
	int ix1 = floor(contInd[0]);
	int iy1 = floor(contInd[1]);
	int iz1 = floor(contInd[2]);

	int ix2 = ix1 + 1;
	int iy2 = iy1 + 1;
	int iz2 = iz1 + 1;

	// Get distances to the image grid
	TScalar fx = contInd[0] - ix1;
	TScalar fy = contInd[1] - iy1;
	TScalar fz = contInd[2] - iz1;

	TScalar gx = ix2 - contInd[0];
	TScalar gy = iy2 - contInd[1];
	TScalar gz = iz2 - contInd[2];

	// Add valid grid positions and corresponding weights
	weights.clear();
	gridIndices.clear();

#define interpWeightMacro3(ix, iy, iz, w) \
		if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
				(0 <= (iy)) && ((iy) < (int)size[1]) && \
				(0 <= (iz)) && ((iz) < (int)size[2])) \
				{ \
			ImageIndexType ind; \
			ind[0] = (ix); ind[1] = (iy); ind[2] = (iz); \
			if (w > 0) \
			{ \
				gridIndices.push_back(ind); \
				weights.push_back(w); \
			} \
				}

	interpWeightMacro3(ix1, iy1, iz1, gx*gy*gz);
	interpWeightMacro3(ix1, iy1, iz2, gx*gy*fz);
	interpWeightMacro3(ix1, iy2, iz1, gx*fy*gz);
	interpWeightMacro3(ix1, iy2, iz2, gx*fy*fz);
	interpWeightMacro3(ix2, iy1, iz1, fx*gy*gz);
	interpWeightMacro3(ix2, iy1, iz2, fx*gy*fz);
	interpWeightMacro3(ix2, iy2, iz1, fx*fy*gz);
	interpWeightMacro3(ix2, iy2, iz2, fx*fy*fz);

#undef interpWeightMacro3


		}



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::ImagePointer
GridFunctions<TScalar, Dimension>
::DownsampleImage(const ImageType* img, TScalar factor)
 {
	if (factor <= 1.0)
	{
		typedef itk::ImageDuplicator<ImageType> DuperType;
		typename DuperType::Pointer dupef = DuperType::New();
		dupef->SetInputImage(img);
		dupef->Update();
		return dupef->GetOutput();
	}

	ImageSpacingType spacing = img->GetSpacing();
	ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

	TScalar minSpacing = spacing[0];
	for (unsigned int dim = 1; dim < Dimension; dim++)
		if (spacing[dim] < minSpacing)
			minSpacing = spacing[dim];

	ImageSpacingType downSpacing;
	ImageSizeType downSize;

	// MARCEL
	// for (unsigned int dim = 0; dim < Dimension; dim++)
	// {
	//   TScalar length = size[dim]*spacing[dim];
	//
	//   //downSize[dim] = static_cast<long>(length / factor);
	//   downSize[dim] = static_cast<long>(size[dim] / factor);
	//
	//   downSpacing[dim] = length / downSize[dim];
	// }

	// STANLEY: this function is used to define a coarser voxel grid, and then to interpolate it to give the coordinate of voxels at a finest resolution
	// in this case, the downsampled image needs to cover a *larger* domain than the fine image. Otherwise, voxels in the fine grid outside the domain covered by the coarse grid will be interpolated with zero values in the neighborhood located outside the domain
	// If this is desirable to interpolate intensities (assuming black background), this has dramatic consequences when interpolating voxel coordinates


	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		downSize[dim] = ceil(size[dim] / factor) + 1;
		downSpacing[dim] = spacing[dim] * factor;
	}

	typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
	GaussFilterType;
	typename GaussFilterType::Pointer gaussf = GaussFilterType::New();
	gaussf->SetInput(img);
	gaussf->SetVariance(minSpacing*minSpacing);
	//gaussf->SetMaximumKernelWidth(5);
	//gaussf->SetMaximumError(0.1);
	gaussf->Update();


	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	typename ResamplerType::Pointer resf = ResamplerType::New();
	//resf->SetInput(img);
	resf->SetInput(gaussf->GetOutput());
	resf->SetOutputParametersFromImage(img);
	resf->SetSize(downSize);
	resf->SetOutputSpacing(downSpacing);
	resf->Update();

	return resf->GetOutput();
 }



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::ImagePointer
GridFunctions<TScalar, Dimension>
::UpsampleImage(const ImageType* exImg, const ImageType* img)
 {
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	typename ResamplerType::Pointer resf = ResamplerType::New();
	resf->SetInput(img);
	resf->SetOutputParametersFromImage(exImg);
	resf->Update();

	return resf->GetOutput();
 }



template <class TScalar, unsigned int Dimension>
typename GridFunctions<TScalar, Dimension>::MatrixType
GridFunctions<TScalar, Dimension>
::UpsampleImagePoints(const ImageType* img, const ImageType* downSampledImg, const MatrixType& downSampledPos)
 {
	unsigned int numPixels = img->GetLargestPossibleRegion().GetNumberOfPixels();

	MatrixType Yup(numPixels, Dimension, 0.0);

	// STANLEY
	MatrixType PointImage = GridFunctions<TScalar, Dimension>::ImageToPoints(img);
	for (unsigned int d = 0; d < Dimension; d++)
	{
		ImagePointer Hd = GridFunctions<TScalar, Dimension>::VectorToImage(downSampledImg, downSampledPos.get_column(d));
		// MARCEL
		// Hd = GridFunctionsType::UpsampleImage(m_Image, Hd);
		// Yup.set_column(d, GridFunctionsType::VectorizeImage(Hd));

		// STANLEY
		Yup.set_column(d, GridFunctions<TScalar, Dimension>::Interpolate(PointImage, Hd));

	}

	return Yup;
 }


#endif /* _GridFunctions_txx */

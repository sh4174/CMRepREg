/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _P3MKernel_txx
#define _P3MKernel_txx

#include "P3MKernel.h"

#include "itkConfigure.h"
#if defined(USE_FFTWF)
  #if ITK_VERSION_MAJOR < 4
    #include "itkFFTWComplexConjugateToRealImageFilter.h"
    #include "itkFFTWRealToComplexConjugateImageFilter.h"
  #else
    #include "itkFFTWForwardFFTImageFilter.h"
    #include "itkFFTWInverseFFTImageFilter.h"
  #endif
#else
  #if ITK_VERSION_MAJOR < 4
    #include "itkVnlFFTRealToComplexConjugateImageFilter.h"
    #include "itkVnlFFTComplexConjugateToRealImageFilter.h"
  #else
    #include "itkVnlForwardFFTImageFilter.h"
    #include "itkVnlInverseFFTImageFilter.h"
  #endif
#endif

#include "itkVersion.h"

#include <itkFFTShiftImageFilter.h>

#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"

#include <exception>
#include <stdexcept>


// #include "itkImageFileWriter.h"
// #include "itkMinimumMaximumImageFilter.h"
// #include "itkRescaleIntensityImageFilter.h"

#define DO_NEAR_FIELD 0
#define DO_SORT_ID 0

//
// Support class
//

class TagIndex
{
public:
	unsigned int index;
	unsigned int j;
	inline TagIndex& operator=(const TagIndex& o)
		{ this->index = o.index; this->j = o.j; return (*this); }
	inline bool operator<(const TagIndex& o) const
		{ return this->index < o.index; }
};

//
// P3MKernel class definition
//

// Static variables
template <class TScalar, unsigned int PointDim>
itk::SimpleFastMutexLock
P3MKernel<TScalar, PointDim>::m_FFTMutex;

template <class TScalar, unsigned int PointDim>
itk::SimpleFastMutexLock
P3MKernel<TScalar, PointDim>::m_CacheMutex;

template<class TScalar, unsigned int PointDim>
std::vector<TScalar>
P3MKernel<TScalar, PointDim>::m_CacheFFTKernelWidths;

template<class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ImageSizeType>
P3MKernel<TScalar, PointDim>::m_CacheFFTKernelSizes;

template<class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ComplexImagePointer>
P3MKernel<TScalar, PointDim>::m_CacheFFTKernelImages;

template<class TScalar, unsigned int PointDim>
std::vector<TScalar>
P3MKernel<TScalar, PointDim>::m_CacheFFTGradientKernelWidths;

template<class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ImageSizeType>
P3MKernel<TScalar, PointDim>::m_CacheFFTGradientKernelSizes;

template<class TScalar, unsigned int PointDim>
std::vector< std::vector<typename P3MKernel<TScalar, PointDim>::ComplexImagePointer> >
P3MKernel<TScalar, PointDim>::m_CacheFFTGradientKernelImages;

template<class TScalar, unsigned int PointDim>
std::vector<TScalar>
P3MKernel<TScalar, PointDim>::m_CacheFFTHessianKernelWidths;

template<class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ImageSizeType>
P3MKernel<TScalar, PointDim>::m_CacheFFTHessianKernelSizes;

template<class TScalar, unsigned int PointDim>
std::vector< std::vector<typename P3MKernel<TScalar, PointDim>::ComplexImagePointer> >
P3MKernel<TScalar, PointDim>::m_CacheFFTHessianKernelImages;


// Methods

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>
::P3MKernel(): Superclass()
 {
	this->Init();
	this->SetKernelWidth(1.0);
 }

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>
::P3MKernel(const P3MKernel& o)
 {
	//this->SetSources(o.GetSources());
	//this->SetWeights(o.GetWeights());
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();

	m_GridOrigin = o.m_GridOrigin;
	m_GridSpacing = o.m_GridSpacing;
	m_GridSize = o.m_GridSize;
	m_GridPadding = o.m_GridPadding;

	m_FFTKernel = o.m_FFTKernel;
	m_FFTGradientKernels = o.m_FFTGradientKernels;
	m_FFTHessianKernels = o.m_FFTHessianKernels;

	// m_MeshPreConvList  = o.m_MeshPreConvList;

	m_MeshList  = o.m_MeshList;
	// m_MeshWYList  = o.m_MeshWYList;
	// m_MeshWYYtList  = o.m_MeshWYYtList;

	m_MeshGradientList = o.m_MeshGradientList;
	m_MeshHessianList = o.m_MeshHessianList;

	m_HessianUpdated = o.m_HessianUpdated;

	m_SourceIDs = o.m_SourceIDs;

	m_IDLookupTable = o.m_IDLookupTable;

	m_NearThresholdScale = o.m_NearThresholdScale;
 }

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>
::P3MKernel(const MatrixType& X, double h): Superclass(X, h)
 {
	this->Init();
	this->SetKernelWidth(h);
	this->SetModified();
 }

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>
::P3MKernel(
		const MatrixType& X, const MatrixType& W,
		double h): Superclass(X, W, h)
		{
	this->Init();
	this->SetKernelWidth(h);
	this->SetModified();
		}

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>
::~P3MKernel()
{
	this->ClearGrids();
	m_FFTKernel = 0;
}

template <class TScalar, unsigned int PointDim>
P3MKernel<TScalar, PointDim>*
P3MKernel<TScalar, PointDim>
::Clone() const
 {
	return new P3MKernel(*this);
 }

template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::ClearGrids()
{
/*
  for (unsigned int k = 0; k < m_MeshList.size(); k++)
    m_MeshList[k] = 0;
  m_MeshList.clear();

  for (unsigned int i = 0; i < m_MeshGradientList.size(); i++)
  {
    ImageListType list = m_MeshGradientList[i];
    for (unsigned int j = 0; j < list.size(); j++)
      list[j] = 0;
    list.clear();
  }
  m_MeshGradientList.clear();

  for (unsigned int k = 0; k < m_MeshHessianList.size(); k++)
  {
    ImageMatrixType mat = m_MeshHessianList[k];
    for (unsigned int i = 0; i < mat.size(); i++)
    {
      ImageListType list = mat[i];
      for (unsigned int j = 0; j < list.size(); j++)
        list[j] = 0;
      list.clear();
    }
    mat.clear();
  }
  m_MeshHessianList.clear();
 */
	m_MeshList.clear();
	m_MeshGradientList.clear();
	m_MeshHessianList.clear();
}

template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::Init()
 {
	m_GridOrigin.Fill(0.0);
	m_GridSpacing.Fill(1.0);
	m_GridSize.Fill(128);
	m_GridPadding = 0;
	m_NearThresholdScale = 0.5;
	m_FFTKernel = 0;

	m_HessianUpdated = true;

	m_WorkingSpacingRatio = 1.0;
	m_PaddingFactor = 0.0;

 }

template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::UpdateGrids()
 {

	// Update grid size and spacing
	this->DetermineGrids();

	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	unsigned int sourceDim = Y.columns();

	if (sourceDim != PointDim)
		throw std::runtime_error("Can only handle certain dimension");
	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	ImageRegionType region;
	region.SetSize(m_GridSize);

	this->ClearGrids();

	for (unsigned int k = 0; k < weightDim; k++)
	{
		ImagePointer img = ImageType::New();
		img->SetRegions(region);
		img->Allocate();
		img->SetOrigin(m_GridOrigin);
		img->SetSpacing(m_GridSpacing);
		img->FillBuffer(0);

		m_MeshList.push_back(img);
	}

	// Splat weight values to meshes
	for (unsigned int i = 0; i < Superclass::m_Sources.rows(); i++)
	{
		VectorType yi = Superclass::m_Sources.get_row(i);

		std::vector<TScalar> weights;
		std::vector<ImageIndexType> gridIndices;

		this->_getInterpolationWeightsAndGridPoints(
				weights, gridIndices, yi);

		for (unsigned int j = 0; j < weights.size(); j++)
		{
			TScalar w = weights[j];
			ImageIndexType ind = gridIndices[j];

			bool isborder = false;
			for (unsigned int d = 0; d < PointDim; d++)
				if (ind[d] <= m_GridPadding || ind[d] >= ((long)m_GridSize[d]-m_GridPadding))
					isborder = true;
			if (isborder)
				continue;

			for (unsigned int k = 0; k < weightDim; k++)
				m_MeshList[k]->SetPixel(ind,
						m_MeshList[k]->GetPixel(ind) + w*Superclass::m_Weights(i, k));
		}
	}

	// apply FFT to the weight meshes
	#if defined(USE_FFTWF)

	#if ITK_VERSION_MAJOR < 4
	  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
	  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
	#else
	  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
	  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
	#endif

	#else

	#if ITK_VERSION_MAJOR < 4
	  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
	  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
	#else
	  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
	  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
	#endif

	#endif

	  // TODO:
	  // fftw execute is thread safe, but planner is not
	  // for thread safety and speed may need to write own ITK fftw wrapper

	  m_MeshListFFT.resize(weightDim);
	  for (unsigned int k = 0; k < weightDim; k++)
	  {
		  m_MeshList[k]->ReleaseDataFlagOff();

		  //P3MKernel::m_FFTMutex.Lock();

		  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
		  forwardFFT->SetInput(m_MeshList[k]);
		  forwardFFT->SetNumberOfThreads(1);
		  forwardFFT->Update();

		  //P3MKernel::m_FFTMutex.Unlock();

		  m_MeshListFFT[k] = forwardFFT->GetOutput();
	  }

	  m_FFTKernel = this->BuildFFTKernel();
	  m_FFTGradientKernels = this->BuildFFTGradientKernels();

//   // Convolve the splatted values in mesh
// 
//   m_MeshPreConvList.resize(weightDim);
//   for (unsigned int k = 0; k < weightDim; k++)
//   {
//     ImagePointer oldmesh = m_MeshList[k]; // Avoid cyclic smart pointer
//     m_MeshList[k] = this->ApplyKernelFFT(m_FFTKernel, oldmesh);
//     m_MeshPreConvList[k] = oldmesh;
//   }
// 
//   for (unsigned int k = 0; k < m_MeshWYList.size(); k++)
//   {
//     ImagePointer oldmesh = m_MeshWYList[k]; // Avoid cyclic smart pointer
//     m_MeshWYList[k] = this->ApplyKernelFFT(m_FFTKernel, oldmesh);
//   }
// 
// #if DO_NEAR_FIELD
//   // Label all source points with nearest grid id
//   unsigned int numSources = Y.rows();
// #if DO_SORT_ID
//   // Sort IDs for faster search
//   TaggedIndex* temp = new TaggedIndex[numSources];
//   for (unsigned int i = 0; i < numSources; i++)
//   {
//     VectorType yi = m_Sources.get_row(i);
//     TaggedIndex[i].index = this->_getPointID(yi);
//     TaggedIndex[i].j = i;
//   }
// 
//   sort(temp, temp+numSources);
// 
//   MatrixType movedSources(numSources, PointDim);
// 
//   m_SourceIDs.resize(numSources);
//   for (unsigned int i = 0; i < numSources; i++)
//   {
//     unsigned int j = temp[i].j;
//     movedSources.set_row(i, m_Sources.get_row(j));
//     m_SourceIDs[i] = temp[i].index;
//   }
// 
//   Superclass::m_Sources = movedSources;
// 
//   delete [] temp;
// #else
//   m_SourceIDs.resize(numSources);
//   for (unsigned int i = 0; i < numSources; i++)
//   {
//     VectorType yi = Superclass::m_Sources.get_row(i);
//     m_SourceIDs[i] = this->_getPointID(yi);
//   }
// #endif
// #endif

  m_HessianUpdated = false;
}

template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::UpdateHessianGrids()
 {
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	unsigned int sourceDim = Y.columns();
//	unsigned int weightDim = W.columns();

	if (sourceDim != PointDim)
		throw std::runtime_error("Can only handle certain dimension");
	if (Y.rows() != W.rows())
		throw std::runtime_error("Source and weights count mismatch");

	// ImageRegionType region;
	// region.SetSize(m_GridSize);
	//
	// m_MeshWYYtList.clear();
	// for (unsigned int k = 0; k < weightDim*(PointDim+1)*PointDim/2; k++)
	// {
	//   ImagePointer img = ImageType::New();
	//   img->SetRegions(region);
	//   img->Allocate();
	//   img->SetOrigin(m_GridOrigin);
	//   img->SetSpacing(m_GridSpacing);
	//   img->FillBuffer(0);
	//
	//   m_MeshWYYtList.push_back(img);
	// }
	//
	// for (unsigned int k = 0; k < weightDim; k++)
	// {
	//   ImageIteratorType it(m_MeshPreConvList[k], region);
	//   for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	//   {
	//     TScalar w = it.Get();
	//
	//     ImageIndexType ind = it.GetIndex();
	//
	//     ImagePointType p;
	//     m_MeshPreConvList[k]->TransformIndexToPhysicalPoint(ind, p);
	//
	//     unsigned int which = k*PointDim*(PointDim+1)/2;
	//     for (unsigned int r = 0; r < PointDim; r++)
	//       for (unsigned int c = r; c < PointDim; c++)
	//         m_MeshWYYtList[which++]->SetPixel(ind, w*p[r]*p[c]);
	//   }
	// }
	//
	// for (unsigned int k = 0; k < m_MeshWYYtList.size(); k++)
	// {
	//   ImagePointer oldmesh = m_MeshWYYtList[k]; // Avoid cyclic smart pointer
	//   m_MeshWYYtList[k] = this->ApplyKernelFFT(m_FFTKernel, oldmesh);
	// }

	m_FFTHessianKernels = this->BuildFFTHessianKernels();

	m_HessianUpdated = true;

}



template<class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::DetermineGrids()
{
	// parameters of the working image
	//   ImageSpacingType spacing = m_WorkingImage->GetSpacing();
	//   ImageSizeType size = m_WorkingImage->GetLargestPossibleRegion().GetSize();
	// ImagePointType origin = m_WorkingImage->GetOrigin();	
	VectorType Xmin = m_DataDomain.get_column(0);
	VectorType Xmax = m_DataDomain.get_column(1);

	TScalar h = this->GetKernelWidth();

	if (h <= 1e-20)
	{
		//     this->SetGridOrigin(origin);
		//     this->SetGridSize(size);
		//     this->SetGridSpacing(spacing);
		// this->SetGridPadding(0.0);
		std::cout << "kernel width too small to set up a p3m grid!" << std::endl;
		return;
	}


	ImageSpacingType grid_spacing;
	ImageSizeType grid_size;
	ImagePointType grid_origin;

	grid_spacing.Fill(m_WorkingSpacingRatio * h);

	TScalar length;
	TScalar padded_length;
	for (unsigned int d = 0; d < PointDim; d++)
	{
		length = Xmax[d] - Xmin[d]; //size[d] * spacing[d];
		padded_length = length + 2*m_PaddingFactor * h;

		grid_size[d] = (long) (padded_length / grid_spacing[d]);
	}

	std::vector<long> pads(PointDim, 0);
	for (unsigned int d = 0; d < PointDim; d++)
	{
		// log_2(size)
		TScalar l2 = log((TScalar) grid_size[d]) / log(2.0);

		// Round up
		long trunc_l2 = (long) l2;
		if ((l2 - trunc_l2) > 0)
			l2 = trunc_l2 + 1;
		else
			l2 = trunc_l2;

		// Nearest power of 2
		grid_size[d] = (long) pow(2, l2);

		TScalar shift = ( grid_size[d]*grid_spacing[d] - (Xmax[d] - Xmin[d]) ) / 2; //(grid_size[d]*grid_spacing[d] - size[d]*spacing[d]) / 2;
		grid_origin[d] = Xmin[d] - shift; //origin[d] - shift;

	}

	// if (grid_size != m_GridSize)
	// std::cout << "Grid changed: origin = " << grid_origin << " size = " << grid_size << " spacing = " << grid_spacing << std::endl;


	this->SetGridOrigin(grid_origin);
	this->SetGridSize(grid_size);
	this->SetGridSpacing(grid_spacing);

}



template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::_getInterpolationWeightsAndGridPoints(
  std::vector<TScalar>& weights,
  std::vector<ImageIndexType>& gridIndices,
  const VectorType& x)
  {

	//TODO sync with the one in GridFunctions
	//typedef GridFunctions<TScalar, PointDim> GridFunctionsType;
	//GridFunctionsType::_getInterpolationWeightsAndGridPoints(weights, gridIndices, x, m_MeshList[0]);

	if (PointDim == 2)
	{
		// Bilinear
		_getInterpolationWeightsAndGridPoints_2(weights, gridIndices, x);
	}
	else if (PointDim == 3)
	{
		// Trilinear
		_getInterpolationWeightsAndGridPoints_3(weights, gridIndices, x);
	}
	else
	{
		// Nearest neighbor
		ImageIndexType ind;
		for (unsigned int d = 0; d < PointDim; d++)
			ind[d] = (long)((x[d]-m_GridOrigin[d]) / m_GridSpacing[d]);
		gridIndices.clear();

		bool isOut = false;
		for (unsigned int k = 0; k < PointDim; k++)
			if (ind[k] < 0 || ind[k] >= (long)m_GridSize[k])
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

template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::_getInterpolationWeightsAndGridPoints_2(
  std::vector<TScalar>& weights,
  std::vector<ImageIndexType>& gridIndices,
  const VectorType& p)
{
	// it works because origin is always the point with the smallest coordinate in each direction (i.e. axes are not flipped) -- different from same function in GridFunctions.txx
	std::vector<TScalar> contInd(PointDim);
	for (unsigned int d = 0; d < PointDim; d++)
		contInd[d] = (p[d]-m_GridOrigin[d]) / m_GridSpacing[d];

	//
	// Bilinear interpolation
	//

	// Get the 4 grid positions

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
  if ((0 <= (ix)) && ((ix) < (int)m_GridSize[0]) && \
    (0 <= (iy)) && ((iy) < (int)m_GridSize[1])) \
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



template <class TScalar, unsigned int PointDim>
void
P3MKernel<TScalar, PointDim>
::_getInterpolationWeightsAndGridPoints_3(
  std::vector<TScalar>& weights,
  std::vector<ImageIndexType>& gridIndices,
  const VectorType& p)
  {
	// it works because origin is always the point with the smallest coordinate in each direction (i.e. axes are not flipped) -- different from same function in GridFunctions.txx
	std::vector<TScalar> contInd(PointDim);
	for (unsigned int d = 0; d < PointDim; d++)
		contInd[d] = (p[d]-m_GridOrigin[d]) / m_GridSpacing[d];

	//
	// Trilinear interpolation
	//

	// Get the 8 grid positions
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
  if ((0 <= (ix)) && ((ix) < (int)m_GridSize[0]) && \
    (0 <= (iy)) && ((iy) < (int)m_GridSize[1]) && \
    (0 <= (iz)) && ((iz) < (int)m_GridSize[2])) \
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


template <class TScalar, unsigned int PointDim>
long
P3MKernel<TScalar, PointDim>
::_getPointID(const VectorType& x)
{
	std::vector<TScalar> contInd(PointDim);
	for (unsigned int d = 0; d < PointDim; d++)
		contInd[d] = (x[d]-m_GridOrigin[d]) / m_GridSpacing[d];

	// Check if out of bounds
	bool isOut = false;
	for (unsigned int d = 0; d < PointDim; d++)
		if (contInd[d] < 0 || contInd[d] > (m_GridSize[d]-1))
		{
			isOut = true;
			break;
		}

	if (isOut)
	{
		long outid = 1;
		for (unsigned int d = 0; d < PointDim; d++)
			outid *= m_GridSize[d];
		outid += 1;
		return outid;
	}

	// Round up continiuous index to nearest point
	ImageIndexType ind;
	for (unsigned int d = 0; d < PointDim; d++)
	{
		ind[d] = (long)(contInd[d]+0.5);
	}

	return this->_getPointID(ind);
}

template <class TScalar, unsigned int PointDim>
long
P3MKernel<TScalar, PointDim>
::_getPointID(const ImageIndexType& ind)
{
	long id = ind[PointDim-1];
	long mult = m_GridSize[PointDim-1];
	for (long d = PointDim-2; d >= 0; d--)
	{
		id += ind[d]*mult;
		mult *= m_GridSize[d];
	}

	return id;
}

template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::VectorType
P3MKernel<TScalar, PointDim>
::Interpolate(
  const VectorType& x, const std::vector<ImagePointer>& images)
{
	unsigned int numInterps = images.size();

	VectorType values(numInterps, 0.0);

	std::vector<TScalar> weights;
	std::vector<ImageIndexType> gridIndices;

	this->_getInterpolationWeightsAndGridPoints(
			weights, gridIndices, x);

	for (unsigned int i = 0; i < weights.size(); i++)
	{
		TScalar w = weights[i];
		ImageIndexType ind = gridIndices[i];
		for (unsigned int j = 0; j < numInterps; j++)
			values[j] += w*images[j]->GetPixel(ind);
	}

	return values;
}

template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::ComplexImagePointer
P3MKernel<TScalar, PointDim>
::BuildFFTKernel()
{
	TScalar kernelWidth = this->GetKernelWidth();

	m_CacheMutex.Lock();

	unsigned int numCache = m_CacheFFTKernelWidths.size();

	unsigned int cacheIndex = numCache+1;
	for (unsigned int i = 0; i < numCache; i++)
		if (m_CacheFFTKernelWidths[i] == kernelWidth)
		{
			bool sameSize = true;
			ImageSizeType size_i = m_CacheFFTKernelSizes[i];
			for (unsigned int d = 0; d < PointDim; d++)
				if (size_i[d] != m_GridSize[d])
				{
					sameSize = false;
					break;
				}
			if (sameSize)
			{
				cacheIndex = i;
				break;
			}
		}

	m_CacheMutex.Unlock();

	if (cacheIndex < numCache)
	{
		m_CacheMutex.Lock();
		ComplexImagePointer kImg = m_CacheFFTKernelImages[cacheIndex];
		m_CacheMutex.Unlock();

		return kImg;
	}
	else
	{
		std::cout << "New grid size: " << m_GridSize << " with spacing = " << m_GridSpacing << std::endl;	
	}


//std::cout << "Building kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

  ImagePointer kernelImg = ImageType::New();
  kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
  kernelImg->Allocate();
  kernelImg->CopyInformation(m_MeshList[0]);
  kernelImg->FillBuffer(0);

  TScalar nearThres = 0;
  for (unsigned int d = 0; d < PointDim; d++)
	  nearThres += m_GridSpacing[d]*m_GridSpacing[d];
  nearThres *= m_NearThresholdScale*m_NearThresholdScale;

  VectorType center(PointDim, 0.0);

#if 1
  for (unsigned int d = 0; d < PointDim; d++)
  {
	  // MARCEL
	  // center[d] =
	  //   m_GridOrigin[d] + 0.5 * (m_GridSize[d]-1) * m_GridSpacing[d];
	  // STANLEY
	  center[d] =
			  m_GridOrigin[d] + 0.5 * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
		  kernelImg, kernelImg->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
	  ImageIndexType ind = it.GetIndex();

	  bool isborder = false;
	  for (unsigned int d = 0; d < PointDim; d++)
		  if (ind[d] <= m_GridPadding || ind[d] >= ((long)m_GridSize[d]-m_GridPadding))
			  isborder = true;
	  if (isborder)
		  continue;

	  ImagePointType p;
	  kernelImg->TransformIndexToPhysicalPoint(ind, p);

	  VectorType x(PointDim, 0.0);
	  for (unsigned int d = 0; d < PointDim; d++)
		  x[d] = p[d];

	  /*
    VectorType dvec = x - center;
    if (dvec.squared_magnitude() < nearThres)
      continue;
	   */

	  TScalar k = this->EvaluateKernel(x, center);

	  it.Set(k);
  }


  // Shift and flip halfplanes
  typedef itk::FFTShiftImageFilter<ImageType, ImageType>
  ShifterType;
  typename ShifterType::Pointer shiftf = ShifterType::New();
  shiftf->SetInput(kernelImg);
  shiftf->SetInverse(true);
  shiftf->SetNumberOfThreads(1);
  shiftf->Update();
  kernelImg = shiftf->GetOutput();
#else
  ImageIteratorType it(
    kernelImg, kernelImg->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageIndexType ind = it.GetIndex();

    for (unsigned int d = 0; d < PointDim; d++)
    { 
      long half = m_GridSize[d]/2 + 2;
      if (ind[d] > half)
        ind[d] = half - (ind[d]%half);
    }

    ImagePointType p;
    kernelImg->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    TScalar k = this->EvaluateKernel(x, center);
  
    it.Set(k);
  }
#endif

	// STANLEY: I remove kernels's normalization, so that the integral of the kernel does not sum to 1, but to sqrt(2*pi*kernelwidth^d)
  // Make sure the discrete version is still a unit kernel
	//   TScalar deltaVol = m_GridSpacing[0];
	//   for (unsigned int d = 1; d < PointDim; d++)
	//     deltaVol *= m_GridSpacing[d];
	//   
	//   TScalar sumK = 1e-20;
	// ImageIteratorType it2(
	//     kernelImg, kernelImg->GetLargestPossibleRegion());
	// 
	//   for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2)
	//     sumK += it2.Get();
	//   sumK *= deltaVol;
	// 
	//   for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2)
	// 	it2.Set( it2.Get() / sumK );


#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

  //P3MKernel::m_FFTMutex.Lock();

  kernelImg->ReleaseDataFlagOff();

  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
  forwardFFT->SetInput(kernelImg);
  forwardFFT->SetNumberOfThreads(1);
  forwardFFT->Update();


  //P3MKernel::m_FFTMutex.Unlock();

  ComplexImagePointer kernelFFTImage = forwardFFT->GetOutput();

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTKernelWidths.push_back(kernelWidth);
  m_CacheFFTKernelSizes.push_back(m_GridSize);
  m_CacheFFTKernelImages.push_back(kernelFFTImage);
  P3MKernel::m_CacheMutex.Unlock();

  return kernelFFTImage;
}

template <class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ComplexImagePointer>
P3MKernel<TScalar, PointDim>
::BuildFFTGradientKernels()
{
	TScalar kernelWidth = this->GetKernelWidth();

	P3MKernel::m_CacheMutex.Lock();

	unsigned int numCache = m_CacheFFTGradientKernelWidths.size();

	unsigned int cacheIndex = numCache+1;
	for (unsigned int i = 0; i < numCache; i++)
		if (m_CacheFFTGradientKernelWidths[i] == kernelWidth)
		{
			bool sameSize = true;
			ImageSizeType size_i = m_CacheFFTGradientKernelSizes[i];
			for (unsigned int d = 0; d < PointDim; d++)
				if (size_i[d] != m_GridSize[d])
				{
					sameSize = false;
					break;
				}
			if (sameSize)
			{
				cacheIndex = i;
				break;
			}
		}

	P3MKernel::m_CacheMutex.Unlock();

	if (cacheIndex < numCache)
	{
		return m_CacheFFTGradientKernelImages[cacheIndex];
	}

	//std::cout << "Building grad kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

	std::vector<ImagePointer> kernList;

	for (unsigned int dim = 0; dim < PointDim; dim++)
	{
		ImagePointer kernelImg = ImageType::New();
		kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
		kernelImg->Allocate();
		kernelImg->CopyInformation(m_MeshList[0]);
		kernelImg->FillBuffer(0);

		kernList.push_back(kernelImg);
	}

	VectorType center(PointDim, 0.0);

//	TScalar sumK = 1e-20;

#if 1
  for (unsigned int d = 0; d < PointDim; d++)
  {
	  center[d] =
			  m_GridOrigin[d] + 0.5 * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
		  kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
	  ImageIndexType ind = it.GetIndex();

	  bool isborder = false;
	  for (unsigned int d = 0; d < PointDim; d++)
		  if (ind[d] <= m_GridPadding || ind[d] >= ((long)m_GridSize[d]-m_GridPadding))
			  isborder = true;
	  if (isborder)
		  continue;

	  ImagePointType p;
	  kernList[0]->TransformIndexToPhysicalPoint(ind, p);

	  VectorType x(PointDim, 0.0);
	  for (unsigned int d = 0; d < PointDim; d++)
		  x[d] = p[d];

	  VectorType g = this->EvaluateKernelGradient(x, center);

	  for (unsigned int dim = 0; dim < PointDim; dim++)
		  kernList[dim]->SetPixel(ind, g[dim]);

	  // TScalar k = this->EvaluateKernel(x, center);
	  //
	  // sumK += k;
  }

  // Shift and flip halfplanes
  for (unsigned int dim = 0; dim < PointDim; dim++)
  {
	  typedef itk::FFTShiftImageFilter<ImageType, ImageType>
	  ShifterType;
	  typename ShifterType::Pointer shiftf = ShifterType::New();
	  shiftf->SetInput(kernList[dim]);
	  shiftf->SetInverse(true);
	  shiftf->Update();
	  kernList[dim] = shiftf->GetOutput();
  }
#else
  ImageIteratorType it(
    kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageIndexType ind = it.GetIndex();

    for (unsigned int d = 0; d < PointDim; d++)
    { 
      long half = m_GridSize[d]/2 + 1;
      if (ind[d] > half)
        ind[d] = half - (ind[d]%half);
    }
  
    ImagePointType p;
    kernList[0]->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    VectorType g = this->EvaluateKernelGradient(x, center);

    for (unsigned int dim = 0; dim < PointDim; dim++)
      kernList[dim]->SetPixel(ind, g[dim]);

    TScalar k = this->EvaluateKernel(x, center);

    sumK += k;
  }
#endif

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

	// STANLEY: does the gradient of the kernel need to sum to 1?? Anyway: NEED TO REDEFINE ITERATOR AFTER HAVING SHIFTED THE HALF-PLANES!!
  // // Make sure the discrete version is still a unit kernel
  // TScalar deltaVol = m_GridSpacing[0];
  // for (unsigned int d = 1; d < PointDim; d++)
  //   deltaVol *= m_GridSpacing[d];
  // 
  // sumK *= deltaVol;
  // 
  // for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  // {
  //   ImageIndexType ind = it.GetIndex();
  //   for (unsigned int dim = 0; dim < PointDim; dim++)
  //     kernList[dim]->SetPixel(ind, kernList[dim]->GetPixel(ind) / sumK);
  // }

  std::vector<ComplexImagePointer> gradFFTKernels;

  for (unsigned int dim = 0; dim < PointDim; dim++)
  {
	  //P3MKernel::m_FFTMutex.Lock();

	  kernList[dim]->ReleaseDataFlagOff();

	  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
	  forwardFFT->SetInput(kernList[dim]);
	  forwardFFT->SetNumberOfThreads(1);
	  forwardFFT->Update();

	  //P3MKernel::m_FFTMutex.Unlock();

	  gradFFTKernels.push_back(forwardFFT->GetOutput());
  } // for dim

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTGradientKernelWidths.push_back(kernelWidth);
  m_CacheFFTGradientKernelSizes.push_back(m_GridSize);
  m_CacheFFTGradientKernelImages.push_back(gradFFTKernels);
  P3MKernel::m_CacheMutex.Unlock();

  return gradFFTKernels;
}

template <class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::ComplexImagePointer>
P3MKernel<TScalar, PointDim>
::BuildFFTHessianKernels()
{
	TScalar kernelWidth = this->GetKernelWidth();

	P3MKernel::m_CacheMutex.Lock();

	unsigned int numCache = m_CacheFFTHessianKernelWidths.size();

	unsigned int cacheIndex = numCache+1;
	for (unsigned int i = 0; i < numCache; i++)
		if (m_CacheFFTHessianKernelWidths[i] == kernelWidth)
		{
			bool sameSize = true;
			ImageSizeType size_i = m_CacheFFTHessianKernelSizes[i];
			for (unsigned int d = 0; d < PointDim; d++)
				if (size_i[d] != m_GridSize[d])
				{
					sameSize = false;
					break;
				}
			if (sameSize)
			{
				cacheIndex = i;
				break;
			}
		}

	P3MKernel::m_CacheMutex.Unlock();

	if (cacheIndex < numCache)
	{
		return m_CacheFFTHessianKernelImages[cacheIndex];
	}

	//std::cout << "Building hess kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

	std::vector<ImagePointer> kernList;

	for (unsigned int r = 0; r < PointDim; r++)
		for (unsigned int c = r; c < PointDim; c++)
		{
			ImagePointer kernelImg = ImageType::New();
			kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
			kernelImg->Allocate();
			kernelImg->CopyInformation(m_MeshList[0]);
			kernelImg->FillBuffer(0);

			kernList.push_back(kernelImg);
		}

	VectorType center(PointDim, 0.0);

//	TScalar sumK = 1e-20;

#if 1
  for (unsigned int d = 0; d < PointDim; d++)
  {
	  center[d] =
			  m_GridOrigin[d] + 0.5 * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
		  kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
	  ImageIndexType ind = it.GetIndex();

	  bool isborder = false;
	  for (unsigned int d = 0; d < PointDim; d++)
		  if (ind[d] <= m_GridPadding || ind[d] >= ((long)m_GridSize[d]-m_GridPadding))
			  isborder = true;
	  if (isborder)
		  continue;

	  ImagePointType p;
	  kernList[0]->TransformIndexToPhysicalPoint(ind, p);

	  VectorType x(PointDim, 0.0);
	  for (unsigned int d = 0; d < PointDim; d++)
		  x[d] = p[d];

	  MatrixType H = this->EvaluateKernelHessian(x, center);

	  unsigned int ikernel = 0;
	  for (unsigned int r = 0; r < PointDim; r++)
		  for (unsigned int c = r; c < PointDim; c++)
			  kernList[ikernel++]->SetPixel(ind, H(r, c));

	  // TScalar k = this->EvaluateKernel(x, center);
	  //
	  // sumK += k;
  }

  // Shift and flip halfplanes
  for (unsigned int j = 0; j < kernList.size(); j++)
  {
	  typedef itk::FFTShiftImageFilter<ImageType, ImageType>
	  ShifterType;
	  typename ShifterType::Pointer shiftf = ShifterType::New();
	  shiftf->SetInput(kernList[j]);
	  shiftf->SetInverse(true);
	  shiftf->Update();
	  kernList[j] = shiftf->GetOutput();
  }
#else
  ImageIteratorType it(
    kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageIndexType ind = it.GetIndex();

    for (unsigned int d = 0; d < PointDim; d++)
    { 
      long half = m_GridSize[d]/2 + 1;
      if (ind[d] > half)
        ind[d] = half - (ind[d]%half);
    }
  
    ImagePointType p;
    kernList[0]->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    MatrixType H = this->EvaluateKernelHessian(x, center);
  
    unsigned int ikernel = 0;
    for (unsigned int r = 0; r < PointDim; r++)
      for (unsigned int c = r; c < PointDim; c++)
        kernList[ikernel++]->SetPixel(ind, H(r, c));

    TScalar k = this->EvaluateKernel(x, center);

    sumK += k;
  }
#endif

	// STANLEY: does the gradient of the kernel need to sum to 1?? Anyway: NEED TO REDEFINE ITERATOR AFTER HAVING SHIFTED THE HALF-PLANES!!
  // Make sure the discrete version is still a unit kernel
  // TScalar deltaVol = m_GridSpacing[0];
  // for (unsigned int d = 1; d < PointDim; d++)
  //   deltaVol *= m_GridSpacing[d];
  // 
  // sumK *= deltaVol;
  // 
  // for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  // {
  //   ImageIndexType ind = it.GetIndex();
  // 
  //   unsigned int ikernel = 0;
  //   for (unsigned int r = 0; r < PointDim; r++)
  //     for (unsigned int c = r; c < PointDim; c++)
  //     {
  //       kernList[ikernel]->SetPixel(ind,
  //         kernList[ikernel]->GetPixel(ind) / sumK);
  //       ikernel++;
  //     }
  // }

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

  std::vector<ComplexImagePointer> hessFFTKernels;

  for (unsigned int j = 0; j < kernList.size(); j++)
  {
	  //P3MKernel::m_FFTMutex.Lock();

	  kernList[j]->ReleaseDataFlagOff();

	  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
	  forwardFFT->SetInput(kernList[j]);
	  forwardFFT->SetNumberOfThreads(1);
	  forwardFFT->Update();

	  //P3MKernel::m_FFTMutex.Unlock();

	  hessFFTKernels.push_back(forwardFFT->GetOutput());
  } // for dim

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTHessianKernelWidths.push_back(kernelWidth);
  m_CacheFFTHessianKernelSizes.push_back(m_GridSize);
  m_CacheFFTHessianKernelImages.push_back(hessFFTKernels);
  P3MKernel::m_CacheMutex.Unlock();

  return hessFFTKernels;
}

template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::ImagePointer
P3MKernel<TScalar, PointDim>
::ApplyKernelFFT(ComplexImageType* kernelImg, ImageType* img)
{
	if (kernelImg == 0)
	{

		std::cout << "no kernelImg!... implementation doubtful..." << std::endl;

		TScalar h = this->GetKernelWidth();

		//typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
		typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>
		GaussFilterType;
		typename GaussFilterType::Pointer gaussf = GaussFilterType::New();
		gaussf->SetInput(img);
		//gaussf->SetVariance(h*h);
		gaussf->SetSigma(h);
		gaussf->SetNormalizeAcrossScale(false);
		gaussf->Update();

		ImagePointer outimg = gaussf->GetOutput();

		TScalar Z = pow(2.0 * vnl_math::pi * h, (TScalar)PointDim / 2);

		ImageIteratorType it(outimg, outimg->GetLargestPossibleRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
			it.Set(it.Get() * Z);

		return outimg;
	}

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#endif

  // TODO:
  // fftw execute is thread safe, but planner is not
  // for thread safety and speed may need to write own ITK fftw wrapper

  //P3MKernel::m_FFTMutex.Lock();

  img->ReleaseDataFlagOff();

  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
  forwardFFT->SetInput(img);
  forwardFFT->SetNumberOfThreads(1);
  forwardFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ComplexImagePointer fftImg = forwardFFT->GetOutput();
  fftImg->ReleaseDataFlagOff();

  typedef itk::ImageRegionIteratorWithIndex<ComplexImageType>
  ComplexImageIteratorType;

  ComplexImageIteratorType f_it(fftImg, fftImg->GetLargestPossibleRegion());
  for (f_it.GoToBegin(); !f_it.IsAtEnd(); ++f_it)
  {
	  ImageIndexType ind = f_it.GetIndex();
	  f_it.Set( f_it.Get()*kernelImg->GetPixel(ind) );
  }

  //P3MKernel::m_FFTMutex.Lock();

  typename InverseFFTType::Pointer invFFT = InverseFFTType::New();
  invFFT->SetInput(fftImg);
  invFFT->SetNumberOfThreads(1);
  invFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ImagePointer outimg = invFFT->GetOutput();

  return outimg;
}


template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::ImagePointer
P3MKernel<TScalar, PointDim>
::ApplyKernelFFT(ComplexImageType* kernelImg, ComplexImageType* img)
{
	if (kernelImg == 0)
		throw std::runtime_error("should give kernelImg");


#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#endif

  typedef itk::ImageDuplicator< ComplexImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(img);
  duplicator->Update();
  typename ComplexImageType::Pointer AuxImg = duplicator->GetOutput();
  AuxImg->ReleaseDataFlagOff();

  typedef itk::ImageRegionIteratorWithIndex<ComplexImageType>
  ComplexImageIteratorType;

  ComplexImageIteratorType f_it(AuxImg, AuxImg->GetLargestPossibleRegion());
  for (f_it.GoToBegin(); !f_it.IsAtEnd(); ++f_it)
  {
	  ImageIndexType ind = f_it.GetIndex();
	  f_it.Set( f_it.Get()*kernelImg->GetPixel(ind) );
  }

  // TODO:
  // fftw execute is thread safe, but planner is not
  // for thread safety and speed may need to write own ITK fftw wrapper

  //P3MKernel::m_FFTMutex.Lock();

  typename InverseFFTType::Pointer invFFT = InverseFFTType::New();
  invFFT->SetInput(AuxImg);
  invFFT->SetNumberOfThreads(1);
  invFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ImagePointer outimg = invFFT->GetOutput();

  return outimg;
}



template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::MatrixType
P3MKernel<TScalar, PointDim>
::Convolve(const MatrixType& X)
{
	
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  MatrixType V(X.rows(), weightDim, 0.0);

  // TScalar nearThres = 0;
  // for (unsigned int d = 0; d < PointDim; d++)
  //   nearThres += m_GridSpacing[d]*m_GridSpacing[d];
  // nearThres *= m_NearThresholdScale*m_NearThresholdScale;


	// Convolve weights splatted on mesh with kernel
	std::vector<ImagePointer> img(weightDim);
  for (unsigned int k = 0; k < weightDim; k++)
    img[k] = this->ApplyKernelFFT(m_FFTKernel, m_MeshListFFT[k]);
		
  for (unsigned int i = 0; i < X.rows(); i++)
  {
    VectorType xi = X.get_row(i);

    VectorType vi = this->Interpolate(xi, img);

#if DO_NEAR_FIELD
    //long id = this->_getPointID(xi);

    std::vector<TScalar> weights;
    std::vector<ImageIndexType> gridIndices;
    this->_getInterpolationWeightsAndGridPoints(
      weights, gridIndices, xi);

    std::vector<long> idlist;
    for (unsigned int j = 0; j < gridIndices.size(); j++)
      idlist.push_back(this->_getPointID(gridIndices[j]));

    for (unsigned int j = 0; j < m_SourceIDs.size(); j++)
    {
#if DO_SORT_ID
      if (m_SourceIDs[j] > id)
        break;
#endif

      //if (m_SourceIDs[j] != id)
      //  continue;
      bool isneigh = false;
      for (unsigned int k = 0; k < idlist.size(); k++)
        if (idlist[k] == m_SourceIDs[j])
        {
          isneigh = true;
          break;
        }
      if (!isneigh)
        continue;

      VectorType yj = Superclass::m_Sources.get_row(j);

      VectorType dvec = xi - yj;

      if (dvec.squared_magnitude() < nearThres)
      {
        TScalar sij = this->EvaluateKernel(xi, yj);
        for (unsigned int k = 0; k < Superclass::m_Weights.columns(); k++)
          vi[k] += sij * Superclass::m_Weights(j, k);
      }
    }

#endif

    V.set_row(i, vi);
  }

  return V;
}



template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::MatrixType
P3MKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, const MatrixType& alpha)
 {
	std::vector<MatrixType> gradMom = this->ConvolveGradient(X);

	MatrixType result(X.rows(), PointDim, 0);
	for (unsigned int j = 0; j < X.rows(); j++)
		result.set_row(j, gradMom[j].transpose() * alpha.get_row(j) );

	return result;
 }



template <class TScalar, unsigned int PointDim>
std::vector<typename P3MKernel<TScalar, PointDim>::MatrixType>
P3MKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X)
 {
	if (this->IsModified())
	{
		this->UpdateGrids();
		this->UnsetModified();
	}

	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	unsigned int sourceDim = Y.columns();
	unsigned int weightDim = W.columns();

	if (sourceDim != PointDim)
		throw std::runtime_error("Can only handle certain dimension");
	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	// TScalar nearThres = 0;
	// for (unsigned int d = 0; d < PointDim; d++)
	//   nearThres += m_GridSpacing[d]*m_GridSpacing[d];
	// nearThres *= m_NearThresholdScale*m_NearThresholdScale;

	std::vector<MatrixType> gradK(
			X.rows(), MatrixType(weightDim, PointDim, 0.0));

	std::vector<ImagePointer> img(weightDim*PointDim);
	for (unsigned int k = 0; k < weightDim; k++)
		for (unsigned int p = 0; p < PointDim; p++)
			img[p + PointDim*k] = this->ApplyKernelFFT(m_FFTGradientKernels[p], m_MeshListFFT[k]);


	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		VectorType wi = this->Interpolate(xi, img);

		MatrixType& G = gradK[i];
		for (unsigned int k = 0; k < weightDim; k++)
			for (unsigned int p = 0; p < PointDim; p++)
				G(k,p) = wi[p + PointDim*k];

	}
	return gradK;
	

	


//   for (unsigned int k = 0; k < weightDim; k++)
//   {
//     ImageListType interpImages;
//     interpImages.push_back(m_MeshList[k]);
// 
//     for (unsigned int d = 0; d < PointDim; d++)
//       interpImages.push_back(m_MeshWYList[k*PointDim+d]);
// 
//     for (unsigned int i = 0; i < X.rows(); i++)
//     {
//       VectorType xi = X.get_row(i);
//       VectorType wi = this->Interpolate(xi, interpImages);
// 
//       VectorType gi(PointDim, 0.0);
//       for (unsigned int d = 0; d < PointDim; d++)
//         gi[d] = wi[0]*xi[d] - wi[d+1];
//       gi *= (-2.0 / Superclass::m_KernelWidthSquared);
// 
//       MatrixType& G = gradK[i];
//       G.set_row(k, G.get_row(k) + gi);
//     }
//   }
// 
// #if DO_NEAR_FIELD
//   for (unsigned int i = 0; i < X.rows(); i++)
//   {
//     VectorType xi = X.get_row(i);
// 
//     //long id = this->_getPointID(xi);
// 
//     std::vector<TScalar> weights;
//     std::vector<ImageIndexType> gridIndices;
//     this->_getInterpolationWeightsAndGridPoints(
//       weights, gridIndices, xi);
// 
//     std::vector<long> idlist;
//     for (unsigned int j = 0; j < gridIndices.size(); j++)
//       idlist.push_back(this->_getPointID(gridIndices[j]));
// 
//     for (unsigned int j = 0; j < m_SourceIDs.size(); j++)
//     {
// #if DO_SORT_ID
//       if (m_SourceIDs[j] > id)
//         break;
// #endif
// 
//       //if (m_SourceIDs[j] != id)
//       //  continue;
//       bool isneigh = false;
//       for (unsigned int k = 0; k < idlist.size(); k++)
//         if (idlist[k] == m_SourceIDs[j])
//         {
//           isneigh = true;
//           break;
//         }
//       if (!isneigh)
//         continue;
// 
//       VectorType yj = Superclass::m_Sources.get_row(j);
// 
//       VectorType dvec = xi - yj;
// 
//       if (dvec.squared_magnitude() < nearThres)
//       {
//         VectorType gij = this->EvaluateKernelGradient(xi, yj);
//         MatrixType& G = gradK[i];
//         for (unsigned int k = 0; k < Superclass::m_Weights.columns(); k++)
//           G.set_row(k, G.get_row(k) + (gij * Superclass::m_Weights(j, k)) );
//       }
//     }
// 
//   }
// #endif
  // 
  // return gradK;
}


template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::MatrixType
P3MKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, unsigned int dim)
{
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
	if ( (dim < 0)||(dim >= PointDim) )
		throw std::runtime_error("dimension parameter out of bounds");

  MatrixType gradK(X.rows(), weightDim, 0.0);

	std::vector<ImagePointer> img(weightDim);
  for (unsigned int k = 0; k < weightDim; k++)
    	img[k] = this->ApplyKernelFFT(m_FFTGradientKernels[dim], m_MeshListFFT[k]);
	
  for (unsigned int i = 0; i < X.rows(); i++)
  {
    VectorType xi = X.get_row(i);

    VectorType wi = this->Interpolate(xi, img);

		gradK.set_row(i, wi);
	}
	
	return gradK;
}







template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::VectorType
P3MKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, unsigned int k, unsigned int dp)
{
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");
  if (dp >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  VectorType gradK(X.rows(), 0.0);

  // ImageListType interpImages;
  // interpImages.push_back(m_MeshList[k]);
  // interpImages.push_back(m_MeshWYList[k*PointDim+dp]);
  // 
  // for (unsigned int i = 0; i < X.rows(); i++)
  // {
  //   VectorType xi = X.get_row(i);
  //   VectorType wi = this->Interpolate(xi, interpImages);
  // 
  //   gradK[i] = wi[0]*xi[dp] - wi[1];
  // }
  // 
  // gradK *= (-2.0 / Superclass::m_KernelWidthSquared);

  return gradK;
}




template <class TScalar, unsigned int PointDim>
std::vector< std::vector<typename P3MKernel<TScalar, PointDim>::MatrixType> >
P3MKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X)
{
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  std::vector< std::vector<MatrixType> > hessK;
  for (unsigned int i = 0; i < X.rows(); i++)
  {
    MatrixType M(PointDim, PointDim, 0.0);
    std::vector<MatrixType> H(weightDim, M);
    hessK.push_back(H);
  }


	int ptsD2 = PointDim*(PointDim + 1)/2;
	std::vector<ImagePointer> img(weightDim*ptsD2);
	for (unsigned int p = 0; p < ptsD2; p++)
		for (unsigned int k = 0; k < weightDim; k++)
			img[k + weightDim*p] = this->ApplyKernelFFT(m_FFTHessianKernels[p], m_MeshListFFT[k]);


  for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);
		VectorType w = this->Interpolate(xi, img);
		
		for (unsigned int k = 0; k < weightDim; k++)
		{
			MatrixType Hik(PointDim,PointDim, 0.0);
	    unsigned int idx = 0;
	    for (unsigned int r = 0; r < PointDim; r++)
				for (unsigned int c = r; c < PointDim; c++)
	      {
	      	Hik(r, c) = w[k + weightDim * idx];
	        Hik(c, r) = Hik(r, c);
					idx ++;
	      }
	
			hessK[i][k] = Hik;
		}
	}


      
	

//   for (unsigned int k = 0; k < weightDim; k++)
//   {
//     ImageListType wy_k;
//     wy_k.push_back(m_MeshList[k]);
//     for (unsigned int d = 0; d < PointDim; d++)
//       wy_k.push_back(m_MeshWYList[k*PointDim + d]);
//     ImageListType wyyt_k;
//     for (unsigned int c = 0; c < PointDim*(PointDim+1)/2; c++)
//       wyyt_k.push_back(m_MeshWYYtList[k*PointDim*(PointDim+1)/2 + c]);
// 
//     for (unsigned int i = 0; i < X.rows(); i++)
//     {
//       VectorType xi = X.get_row(i);
// 
//       VectorType wy = this->Interpolate(xi, wy_k);
//       VectorType wyyt = this->Interpolate(xi, wyyt_k);
// 
//       MatrixType Hik(PointDim, PointDim, 0.0);
// 
//       unsigned int idx = 0;
//       for (unsigned int r = 0; r < PointDim; r++)
//         for (unsigned int c = r; c < PointDim; c++)
//         {
//           Hik(r, c) = wyyt[idx++];
//           Hik(c, r) = Hik(r, c);
//         }
// 
//       TScalar w = wy[0];
// 
//       for (unsigned int r = 0; r < PointDim; r++)
//         for (unsigned int c = 0; c < PointDim; c++)
//           Hik(r, c) += xi[r]*xi[c]*w - wy[r+1]*xi[c] - xi[r]*wy[c+1];
// 
//       Hik *= (4.0 / (Superclass::m_KernelWidthSquared*Superclass::m_KernelWidthSquared));
// 
//       for (unsigned int r = 0; r < PointDim; r++)
//         Hik(r, r) -= 2.0 / Superclass::m_KernelWidthSquared * w;
// 
//       hessK[i][k] = Hik;
//     } // for i
// 
//   } // for k
// 
// #if DO_NEAR_FIELD
//   for (unsigned int i = 0; i < X.rows(); i++)
//   {
//     VectorType xi = X.get_row(i);
// 
//     //long id = this->_getPointID(xi);
// 
//     std::vector<TScalar> weights;
//     std::vector<ImageIndexType> gridIndices;
//     this->_getInterpolationWeightsAndGridPoints(
//       weights, gridIndices, xi);
// 
//     std::vector<long> idlist;
//     for (unsigned int j = 0; j < gridIndices.size(); j++)
//       idlist.push_back(this->_getPointID(gridIndices[j]));
// 
//     for (unsigned int j = 0; j < m_SourceIDs.size(); j++)
//     {
// #if DO_SORT_ID
//       if (m_SourceIDs[j] > id)
//         break;
// #endif
// 
//       //if (m_SourceIDs[j] != id)
//       //  continue;
//       bool isneigh = false;
//       for (unsigned int k = 0; k < idlist.size(); k++)
//         if (idlist[k] == m_SourceIDs[j])
//         {
//           isneigh = true;
//           break;
//         }
//       if (!isneigh)
//         continue;
// 
//       VectorType yj = Superclass::m_Sources.get_row(j);
// 
//       VectorType dvec = xi - yj;
// 
//       if (dvec.squared_magnitude() < nearThres)
//       {
//         MatrixType Hij = this->EvaluateKernelHessian(xi, yj);
//         for (unsigned int k = 0; k < Superclass::m_Weights.columns(); k++)
//           hessK[i][k] += Hij * Superclass::m_Weights(j, k);
//       }
//     }
// 
//   }
// #endif

  return hessK;
}




template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::MatrixType
P3MKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X, unsigned int row, unsigned int col)
{
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
	if (row >= PointDim || col >= PointDim)
		throw std::runtime_error("dimension index out of bounds");

  MatrixType hessK(X.rows(), weightDim, 0.0);
 
	int index = (row <= col)?(col + PointDim*row - row*(row+1)/2):(row + PointDim*col - col*(col+1)/2);
	std::vector<ImagePointer> img(weightDim);
	for (unsigned int k = 0; k < weightDim; k++)
		img[k] = this->ApplyKernelFFT(m_FFTHessianKernels[index], m_MeshListFFT[k]);
	
	for (int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);
		VectorType wi = this->Interpolate(xi,img);
		
		hessK.set_row(i, wi);
	}



//   for (unsigned int k = 0; k < weightDim; k++)
//   {
//     ImageListType wy_k;
//     wy_k.push_back(m_MeshList[k]);
//     for (unsigned int d = 0; d < PointDim; d++)
//       wy_k.push_back(m_MeshWYList[k*PointDim + d]);
//     ImageListType wyyt_k;
//     for (unsigned int c = 0; c < PointDim*(PointDim+1)/2; c++)
//       wyyt_k.push_back(m_MeshWYYtList[k*PointDim*(PointDim+1)/2 + c]);
// 
//     for (unsigned int i = 0; i < X.rows(); i++)
//     {
//       VectorType xi = X.get_row(i);
// 
//       VectorType wy = this->Interpolate(xi, wy_k);
//       VectorType wyyt = this->Interpolate(xi, wyyt_k);
// 
//       MatrixType Hik(PointDim, PointDim, 0.0);
// 
//       unsigned int idx = 0;
//       for (unsigned int r = 0; r < PointDim; r++)
//         for (unsigned int c = r; c < PointDim; c++)
//         {
//           Hik(r, c) = wyyt[idx++];
//           Hik(c, r) = Hik(r, c);
//         }
// 
//       TScalar w = wy[0];
// 
//       for (unsigned int r = 0; r < PointDim; r++)
//         for (unsigned int c = 0; c < PointDim; c++)
//           Hik(r, c) += xi[r]*xi[c]*w - wy[r+1]*xi[c] - xi[r]*wy[c+1];
// 
//       Hik *= (4.0 / (Superclass::m_KernelWidthSquared*Superclass::m_KernelWidthSquared));
// 
//       for (unsigned int r = 0; r < PointDim; r++)
//         Hik(r, r) -= 2.0 / Superclass::m_KernelWidthSquared * w;
// 
//       hessK[i][k] = Hik;
//     } // for i
// 
//   } // for k
// 
// #if DO_NEAR_FIELD
//   for (unsigned int i = 0; i < X.rows(); i++)
//   {
//     VectorType xi = X.get_row(i);
// 
//     //long id = this->_getPointID(xi);
// 
//     std::vector<TScalar> weights;
//     std::vector<ImageIndexType> gridIndices;
//     this->_getInterpolationWeightsAndGridPoints(
//       weights, gridIndices, xi);
// 
//     std::vector<long> idlist;
//     for (unsigned int j = 0; j < gridIndices.size(); j++)
//       idlist.push_back(this->_getPointID(gridIndices[j]));
// 
//     for (unsigned int j = 0; j < m_SourceIDs.size(); j++)
//     {
// #if DO_SORT_ID
//       if (m_SourceIDs[j] > id)
//         break;
// #endif
// 
//       //if (m_SourceIDs[j] != id)
//       //  continue;
//       bool isneigh = false;
//       for (unsigned int k = 0; k < idlist.size(); k++)
//         if (idlist[k] == m_SourceIDs[j])
//         {
//           isneigh = true;
//           break;
//         }
//       if (!isneigh)
//         continue;
// 
//       VectorType yj = Superclass::m_Sources.get_row(j);
// 
//       VectorType dvec = xi - yj;
// 
//       if (dvec.squared_magnitude() < nearThres)
//       {
//         MatrixType Hij = this->EvaluateKernelHessian(xi, yj);
//         for (unsigned int k = 0; k < Superclass::m_Weights.columns(); k++)
//           hessK[i][k] += Hij * Superclass::m_Weights(j, k);
//       }
//     }
// 
//   }
// #endif

  return hessK;
}



template <class TScalar, unsigned int PointDim>
typename P3MKernel<TScalar, PointDim>::VectorType
P3MKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X, unsigned int k,
  unsigned int dp, unsigned int dq)
{
  if (this->IsModified())
  {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

  MatrixType& Y = this->GetSources();
  MatrixType& W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");
  if (dp >= PointDim || dq >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  VectorType hessK(X.rows(), 0.0);

  // unsigned int wyyt_idx = 0;
  // for (unsigned int r = 0; r < PointDim; r++)
  // {
  //   bool found = false;
  //   for (unsigned int c = r; c < PointDim; c++)
  //   {
  //     if ((r == dp && c == dq) || (r == dq && c == dp))
  //     {
  //       found = true;
  //       break;
  //     }
  //     wyyt_idx++;
  //   }
  //   if (found)
  //     break;
  // }
  // 
  // ImageListType interpImages;
  // interpImages.push_back(m_MeshList[k]);
  // interpImages.push_back(m_MeshWYList[k*PointDim + dp]);
  // interpImages.push_back(m_MeshWYList[k*PointDim + dq]);
  // interpImages.push_back(m_MeshWYYtList[k*PointDim*(PointDim+1)/2 + wyyt_idx]);
  // 
  // for (unsigned int i = 0; i < X.rows(); i++)
  // {
  //   VectorType xi = X.get_row(i);
  // 
  //   VectorType iVals = this->Interpolate(xi, interpImages);
  // 
  //   TScalar h_pq = iVals[3];
  // 
  //   h_pq += xi[dp]*xi[dq]*iVals[0] - iVals[1]*xi[dq] - xi[dp]*iVals[2];
  //   h_pq *= (4.0 / (Superclass::m_KernelWidthSquared*Superclass::m_KernelWidthSquared));
  // 
  //   if (dp == dq)
  //     h_pq -= 2.0 / Superclass::m_KernelWidthSquared * iVals[0];
  // 
  //   hessK[i] = h_pq;
  // } // for i

  return hessK;
}

#endif

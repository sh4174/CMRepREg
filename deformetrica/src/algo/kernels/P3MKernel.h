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

#ifndef _P3MKernel_h
#define _P3MKernel_h

#include "ExactKernel.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkSimpleFastMutexLock.h"

/**
 *	\brief      P3M kernels.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The P3MKernel class inherited from AbstractKernel implements operations with kernels
 *              using approximations based on FFTs and projections/interpolation on regular lattices.
 */
template <class TScalar, unsigned int PointDim>
class P3MKernel: public ExactKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> Superclass;

	/// Vector type
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;

	/// Image type (itk).
	typedef itk::Image<TScalar, PointDim> ImageType;
	/// Image pointer type (itk).
	typedef typename ImageType::Pointer ImagePointer;
	/// Image index type (itk).
	typedef typename ImageType::IndexType ImageIndexType;
	/// Image point type (itk).
	typedef typename ImageType::PointType ImagePointType;
	/// Image region type (itk).
	typedef typename ImageType::RegionType ImageRegionType;
	/// Image size type (itk).
	typedef typename ImageType::SizeType ImageSizeType;
	/// Image spacing type (itk).
	typedef typename ImageType::SpacingType ImageSpacingType;

	/// Complex image type (itk).
	typedef itk::Image< std::complex<TScalar>, PointDim> ComplexImageType;
	/// Complex image pointer type (itk).
	typedef typename ComplexImageType::Pointer ComplexImagePointer;

	/// Image iterator type (itk).
	typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;

	/// List of images type (itk).
	typedef std::vector<ImagePointer> ImageListType;
	/// Image matrix type (itk).
	typedef std::vector<ImageListType> ImageMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	P3MKernel();
	/// Copy constructor.
	P3MKernel(const P3MKernel& o);
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, double h) for details.
	P3MKernel(const MatrixType& X, double h);
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, const MatrixType& W, double h) for details.
	P3MKernel(const MatrixType& X, const MatrixType& W, double h);

	virtual P3MKernel* Clone() const;

	virtual ~P3MKernel();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the origin of the grid to \e d.
	void SetGridOrigin(const ImagePointType& p) { m_GridOrigin = p; }

	/// Sets the grid spacing to \e d.
	void SetGridSpacing(const ImageSpacingType& s) { m_GridSpacing = s; }

	/// Sets the size of the grid to \e size.
	void SetGridSize(const ImageSizeType& size) { m_GridSize = size; }

	/// TODO .
	void SetGridPadding(long n) { m_GridPadding = n; }

// void SetWorkingImage(ImageType* img) { m_WorkingImage = img; }

	/// TODO .
	void SetDataDomain(MatrixType DD) { m_DataDomain = DD; }

	/// TODO .
	void SetWorkingSpacingRatio(const TScalar d) { m_WorkingSpacingRatio = d; }

	/// Sets the padding factor to \e d.
	void SetPaddingFactor(const TScalar d) { m_PaddingFactor = d; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Automatically called before convolution calls
	/// TODO .
	virtual void UpdateGrids();
	/// TODO .
	virtual void UpdateHessianGrids();

	/// TODO .
	void SetNearThresholdScale(double s) { m_NearThresholdScale = s; }

	virtual MatrixType Convolve(const MatrixType& X);

	virtual MatrixType ConvolveGradient(const MatrixType& X, const MatrixType& alpha);
	virtual std::vector<MatrixType> ConvolveGradient(const MatrixType& X);

	virtual VectorType ConvolveGradient(const MatrixType& X, unsigned int k, unsigned int dp);
	virtual MatrixType ConvolveGradient(const MatrixType& X, unsigned int dim);

	virtual std::vector< std::vector<MatrixType> > ConvolveHessian(const MatrixType & X);

	virtual VectorType ConvolveHessian(const MatrixType & X, unsigned int k, unsigned int dp, unsigned int dq);
	virtual MatrixType ConvolveHessian(const MatrixType& X, unsigned int row, unsigned int col);

	/// TODO .
	ComplexImagePointer BuildFFTKernel();
	/// TODO .
	std::vector<ComplexImagePointer> BuildFFTGradientKernels();
	/// TODO .
	std::vector<ComplexImagePointer> BuildFFTHessianKernels();

	// Methods to override the FFT kernel used
	/// TODO .
	void UseFFTKernel(ComplexImagePointer kernImg) { m_FFTKernel = kernImg; }
	/// TODO .
	void UseFFTGradientKernels(std::vector<ComplexImagePointer>& klist) { m_FFTGradientKernels = klist; }
	/// TODO .
	void UseFFTHessianKernels(std::vector<ComplexImagePointer>& klist) { m_FFTHessianKernels = klist; }



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	void Init();

	/// TODO .
	void ClearGrids();
	/// TODO .
	void DetermineGrids();

	/// TODO .
	VectorType Interpolate(const VectorType& x, const std::vector<ImagePointer>& imgs);

	/// TODO .
	void SplatToGrid(ImageType* mesh, const MatrixType& X, const VectorType& values);

	/// TODO .
	ImagePointer ApplyKernelFFT(ComplexImageType* kernelImg, ImageType* img);
	/// TODO .
	ImagePointer ApplyKernelFFT(ComplexImageType* kernelImg, ComplexImageType* img);


	/// TODO .
	void inline _getInterpolationWeightsAndGridPoints(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x);
	/// Specialized linear interpolation weights for 2D.
	void inline _getInterpolationWeightsAndGridPoints_2(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x);
	/// Specialized linear interpolation weights for 3D.
	void inline _getInterpolationWeightsAndGridPoints_3(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VectorType& x);

	/// TODO .
	long inline _getPointID(const VectorType& x);
	/// TODO .
	long inline _getPointID(const ImageIndexType& x);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

// ImagePointer m_WorkingImage;
	/// TODO .
	MatrixType m_DataDomain;

	/// TODO .
	ImagePointType m_GridOrigin;
	/// TODO .
	ImageSpacingType m_GridSpacing;
	/// TODO .
	ImageSizeType m_GridSize;
	/// TODO .
	long m_GridPadding;
	//STANLEY
	/// TODO .
	// 1/5 gives a relative approximation error of about 5%, use 0.3 for increased speed
	TScalar m_WorkingSpacingRatio;
	/// Padding factor. It will enlarge the grid by m_PaddingFactor x m_KernelWidth to avoid side effects
	/// (FFTs have circular boundary conditions). It is also used to define a bounding box.
	TScalar m_PaddingFactor;
	//STANLEY


	/// \cond HIDE_FOR_DOXYGEN

	ComplexImagePointer m_FFTKernel;
	std::vector<ComplexImagePointer> m_FFTGradientKernels;
	std::vector<ComplexImagePointer> m_FFTHessianKernels;

// Splatted values before convolution
// std::vector<ImagePointer> m_MeshPreConvList;

	// Splatted values after convolution
	std::vector<ImagePointer> m_MeshList;
	std::vector<ComplexImagePointer> m_MeshListFFT;
// std::vector<ImagePointer> m_MeshWYList;
// std::vector<ImagePointer> m_MeshWYYtList;

	ImageMatrixType m_MeshGradientList;
	std::vector< ImageMatrixType > m_MeshHessianList;

	bool m_HessianUpdated;

	std::vector<unsigned int> m_SourceIDs;

	std::vector< std::vector<unsigned int> > m_IDLookupTable;

	TScalar m_NearThresholdScale;

	static itk::SimpleFastMutexLock m_FFTMutex;

	static itk::SimpleFastMutexLock m_CacheMutex;

	// Cache for BuildFFT...Kernel()
	static std::vector<TScalar> m_CacheFFTKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTKernelSizes;
	static std::vector<ComplexImagePointer> m_CacheFFTKernelImages;

	static std::vector<TScalar> m_CacheFFTGradientKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTGradientKernelSizes;
	static std::vector< std::vector<ComplexImagePointer> >  m_CacheFFTGradientKernelImages;

	static std::vector<TScalar> m_CacheFFTHessianKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTHessianKernelSizes;
	static std::vector< std::vector<ComplexImagePointer> > m_CacheFFTHessianKernelImages;

	/// \endcond

}; /* class P3MKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "P3MKernel.txx"
#endif


#endif /* _P3MKernel_h */

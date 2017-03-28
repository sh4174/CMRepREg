/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _EQLAImage_txx
#define _EQLAImage_txx

#include "EQLAImage.h"

#include "KernelFactory.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkSqrtImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkAddImageFilter.h"

#include "GridFunctions.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
EQLAImage<TScalar, Dimension>
::EQLAImage() : Superclass()
{
	this->SetEQLAImageType();
	m_EQLAKernelWidth = 5.0;
	// m_SmallestImageDimension = 0;
}



template <class TScalar, unsigned int Dimension>
EQLAImage<TScalar, Dimension>
::EQLAImage(const EQLAImage& o) : Superclass(o)
{

	m_EQLAKernelWidth = o.m_EQLAKernelWidth;
	m_LocalMeanImage = o.m_LocalMeanImage;
	m_LocalVarianceImage = o.m_LocalVarianceImage;
	// m_SmallestImageDimension = o.m_SmallestImageDimension;

}


template <class TScalar, unsigned int Dimension>
EQLAImage<TScalar, Dimension>
::EQLAImage(const EQLAImage& ex, const MatrixType& IP) : Superclass(ex, IP)
{
	this->SetEQLAImageType();
	m_EQLAKernelWidth = ex.m_EQLAKernelWidth;
	this->Update();
}



template <class TScalar, unsigned int Dimension>
EQLAImage<TScalar, Dimension>
::~EQLAImage()
{
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation methods
////////////////////////////////////////////////////////////////////////////////////////////////////

// template <class TScalar, unsigned int Dimension>
// void
// EQLAImage<TScalar, Dimension>
// ::SetImage(ImageType* img)
// {
// 	Superclass::m_Image = img;
// 	
// 	ImageSizeType size = img->GetLargestPossibleRegion().GetSize();
// 	m_SmallestImageDimension = size[0];
// 	for (int dim = 1; dim < Dimension; dim++)
// 	{
// 		if (m_SmallestImageDimension>size[dim])
// 			m_SmallestImageDimension = size[dim];
// 	}
// 	
// 	m_LocalMeanImage = this->ComputeLocalMeanImage(img);
// 	ImageTypePointer LocalVarianceImage = this->ComputeLocalCovarianceImage(img, m_LocalMeanImage, img, m_LocalMeanImage);
// 
// 	// small perturbation to avoid "division by zero" issue in image regions of constant intensity	
// 	typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddImageFilterType;
// 	typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
// 	addFilter->SetInput1(LocalVarianceImage);
// 	addFilter->SetConstant2(pow(10,-12));
// 	addFilter->Update();
// 
// 	m_LocalVarianceImage = addFilter->GetOutput();
// 	
// 	this->SetModified();
// }



template <class TScalar, unsigned int Dimension>
void
EQLAImage<TScalar, Dimension>
::Update()
{
	if (this->IsModified())
	{
		Superclass::Update();
		
		m_LocalMeanImage = this->ComputeLocalMeanImage(this->GetImage());
		ImageTypePointer LocalVarianceImage = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, this->GetImage(), m_LocalMeanImage);

		// small perturbation to avoid "division by zero" issue in image regions of constant intensity	
		typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddImageFilterType;
		typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
		addFilter->SetInput1(LocalVarianceImage);
		addFilter->SetConstant2(pow(10,-12));
		addFilter->Update();

		m_LocalVarianceImage = addFilter->GetOutput();
	}
	
	this->UnSetModified();
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int Dimension>
TScalar EQLAImage<TScalar, Dimension>
::ComputeMatch(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	// EQLAImage* Idef = dynamic_cast<EQLAImage*>(this->GetDeformedObjectAt(nbTimePoints-1));	
	// ImageTypePointer Img = this->GetImageFromDownSampledImagePoints();
	const EQLAImage* targ = dynamic_cast<const EQLAImage*>(target);
	
	this->Update();
	// target->Update();

	// Get the Local Mean Image and Standard Deviation Image of the target image
	// ImageTypePointer targLocalMean = targ->GetLocalMeanImage();
	// ImageTypePointer targLocalVariance = targ->GetLocalVarianceImage();
	// 
	// // Compute the Local Mean Image and Standard Deviation Image of the deformed Image
	// ImageTypePointer ImgLocalMean = this->ComputeLocalMeanImage(Img);
	// ImageTypePointer ImgLocalVariance = this->ComputeLocalCovarianceImage(Img, ImgLocalMean, Img, ImgLocalMean);
	
	// Compute the Local Covariance between deformed source and target
	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage());

	// typename ImageType::IndexType pixelIndex;
	// pixelIndex[0] = 60;      // x position of the pixel
	// pixelIndex[1] = 25;      // y position of the pixel
	// 
	// std::cout << " covariance = " << LocalCovariance->GetPixel(pixelIndex) << std::endl;
	// std::cout << " VAR 1 = " << defImgLocalVariance->GetPixel(pixelIndex) << std::endl;
	// std::cout << " VAR 2 = " << targLocalVariance->GetPixel(pixelIndex) << std::endl;
	// 	
	// typedef itk::Image<unsigned char, Dimension> OutImageType;
	// typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
	// typedef itk::ImageFileWriter<OutImageType> WriterType;
	// 
	// typename CasterType::Pointer castf = CasterType::New();
	// castf->SetInput(defImgLocalMean); //targLocalVariance
	// castf->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// castf->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// castf->Update();
	// 
	// typename WriterType::Pointer writer = WriterType::New();		
	// writer->SetInput(castf->GetOutput());
	// writer->SetFileName("defImgLocalMean.png");
	// writer->Update();
	// 
	// typename CasterType::Pointer castf2 = CasterType::New();
	// castf2->SetInput(defImgLocalVariance);
	// castf2->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// castf2->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// castf2->Update();
	// 
	// typename WriterType::Pointer writer2 = WriterType::New();		
	// writer2->SetInput(castf2->GetOutput());
	// writer2->SetFileName("defImgLocalVariance.png");
	// writer2->Update();
	// 
	// typename CasterType::Pointer castf3 = CasterType::New();
	// castf3->SetInput(LocalCovariance);
	// castf3->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// castf3->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// castf3->Update();
	// 
	// typename WriterType::Pointer writer3 = WriterType::New();		
	// writer3->SetInput(castf3->GetOutput());
	// writer3->SetFileName("LocalCovariance.png");
	// writer3->Update();
	// 
	// typename CasterType::Pointer castf4 = CasterType::New();
	// castf4->SetInput(targLocalMean);
	// castf4->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// castf4->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// castf4->Update();
	// 
	// typename WriterType::Pointer writer4 = WriterType::New();		
	// writer4->SetInput(castf4->GetOutput());
	// writer4->SetFileName("targLocalMean.png");
	// writer4->Update();
	
	
	// Vectorize images
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	VectorType LC = GridFunctionsType::VectorizeImage(LocalCovariance);
	VectorType Var1 = GridFunctionsType::VectorizeImage(m_LocalVarianceImage);
	VectorType Var2 = GridFunctionsType::VectorizeImage(targ->GetLocalVarianceImage());
	
	TScalar match = 0.0;
	for (int i = 0; i < LC.size(); i++)
	{
		match += Var1.get(i) - ( LC.get(i) * LC.get(i) / Var2.get(i) );
	}
	
	
	// ImageTypePointer T1 = GridFunctionsType::VectorToImage(Superclass::m_Image, Var1);
	// typedef itk::Image<TScalar, Dimension> OutImageType2;
	// typedef itk::ImageFileWriter<OutImageType2> WriterType2;
	// typename WriterType2::Pointer writer2 = WriterType2::New();		
	// writer2->SetInput(T1);
	// writer2->SetFileName("T1.img");
	// writer2->Update();
	// 
	// ImageTypePointer T2 = GridFunctionsType::VectorToImage(Superclass::m_Image, Term2);
	// typedef itk::Image<TScalar, Dimension> OutImageType3;
	// typedef itk::ImageFileWriter<OutImageType2> WriterType3;
	// typename WriterType3::Pointer writer3 = WriterType3::New();		
	// writer3->SetInput(T2);
	// writer3->SetFileName("T2.img");
	// writer3->Update();

	// ImageTypePointer VI = GridFunctionsType::VectorToImage(Superclass::m_Image, V);
	// typedef itk::Image<TScalar, Dimension> OutImageType4;
	// typedef itk::ImageFileWriter<OutImageType4> WriterType4;
	// typename WriterType4::Pointer writer4 = WriterType4::New();		
	// writer4->SetInput(VI);
	// writer4->SetFileName("VI.img");
	// writer4->Update();
	
	
	
	// std::cout << EQLA << std::endl;
	// std::cout << STD1.extract(10) << std::endl;
	// std::cout << STD2.extract(10) << std::endl;
	// std::cout << LC.get(1) / (STD1.get(1) * STD1.get(1)) << std::endl;
	
	// typedef itk::Image<unsigned short, Dimension> OutImageType;
	// typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
	// typedef itk::ImageFileWriter<OutImageType> WriterType;
	// 
	// typename CasterType::Pointer castf5 = CasterType::New();
	// castf5->SetInput(GridFunctionsType::VectorToImage(Superclass::m_Image, EQLA));
	// castf5->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// castf5->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// castf5->Update();
	// 
	// typename WriterType::Pointer writer5 = WriterType::New();		
	// writer5->SetInput(castf5->GetOutput());
	// writer5->SetFileName("EQLA.png");
	// writer5->Update();
		
//	std::cout << "match = " << match << std::endl;
	
	return match;
}



template<class TScalar, unsigned int Dimension>
typename EQLAImage<TScalar, Dimension>::MatrixType EQLAImage<TScalar, Dimension>
::ComputeMatchGradient(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		std::cerr << "Deformable objects types mismatched: " << this->GetType() << " and " << target->GetType() << "\n";

	/*
	// Need Yfinal to interpolate the deformed gradient image \nabla I
	MatrixType Yfinal;
	if (Superclass::m_ComputeTrueInverseFlow)
	{
		Yfinal = this->UpSampleImageMap(Superclass::m_InverseMapsT[nbTimePoints-1]);
	}
	else
	{
		Yfinal = this->UpSampleImageMap(Superclass::m_MapsT[0]);
	}
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	ImageTypePointer Img = GridFunctionsType::VectorToImage(Superclass::m_Image,
		GridFunctionsType::Interpolate(Yfinal, Superclass::m_Image));
	*/
	
	const EQLAImage* targ = dynamic_cast<const EQLAImage*>(target);
	
	this->Update();
	// target->Update();
	
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

// //	MatrixType Yfinal = GridFunctionsType::UpsampleImagePoints(this->GetImage(), this->GetDownSampledWorkingImage(), Superclass::m_DownSampledY1);
// 	MatrixType Yfinal = this->UpSampleImageMap(Superclass::m_DownSampledY1);
// 	ImageTypePointer Img = this->GetImageFromDownSampledImagePoints();
// 	
// 	// Get the Local Mean Image and Standard Deviation Image of the target image
// 	ImageTypePointer targLocalMean = targ->GetLocalMeanImage();
// 	ImageTypePointer targLocalVariance = targ->GetLocalVarianceImage();
// 	
// 	// Compute the Local Mean Image and Standard Deviation Image of the deformed Image
// 	ImageTypePointer ImgLocalMean = this->ComputeLocalMeanImage(Img);
	// ImageTypePointer ImgLocalVariance = this->ComputeLocalCovarianceImage(this->GetWorkingImage(), m_LocalMeanImage, this->GetWorkingImage(), m_LocalMeanImage);
	
	// Compute the Local Covariance between deformed source and target
	// ImageTypePointer targImg =  targ->GetImage();
	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage());

	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyImageFilterType;
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractImageFilterType;
	typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DivideImageFilterType;
	
	// Compute Cov(Idef,targ)/Var(targ)
	typename DivideImageFilterType::Pointer divideFilter = DivideImageFilterType::New();
	divideFilter->SetInput1(LocalCovariance);
	divideFilter->SetInput2(targ->GetLocalVarianceImage());
	divideFilter->Update();
	
	// Compute W*( Cov(Idef,targ)/Var(targ) )
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->SetInput(divideFilter->GetOutput());
	gaussf->Update();
	
	ImageTypePointer Convolution1 = gaussf->GetOutput();
	
	// Compute W*( targMean * Cov(Idef,targ)/Var(targ) )
	typename MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
	multiplyFilter->SetInput1(targ->GetLocalMeanImage());
	multiplyFilter->SetInput2(divideFilter->GetOutput());
	multiplyFilter->Update();

	typename GaussianFilterType::Pointer gaussf2 = GaussianFilterType::New();
	gaussf2->SetSigma( m_EQLAKernelWidth );
	gaussf2->SetInput(multiplyFilter->GetOutput());
	gaussf2->Update();
	
	ImageTypePointer Convolution2 = gaussf2->GetOutput();
	
	// Compute W*( IdefMean)
	typename GaussianFilterType::Pointer gaussf3 = GaussianFilterType::New();
	gaussf3->SetSigma( m_EQLAKernelWidth );
	gaussf3->SetInput(m_LocalMeanImage);
	gaussf3->Update();

	ImageTypePointer IdefMeanMean = gaussf3->GetOutput();
	
	
	// Vectorize required images to compute the gradient
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	VectorType VtargImg = GridFunctionsType::VectorizeImage(targ->GetImage());
	VectorType VdefImg = GridFunctionsType::VectorizeImage(this->GetImage());
	VectorType VdefImgMeanMean = GridFunctionsType::VectorizeImage(IdefMeanMean);
	VectorType VConv1 = GridFunctionsType::VectorizeImage(Convolution1);
	VectorType VConv2 = GridFunctionsType::VectorizeImage(Convolution2);
	
	// Compute the gradient of the deformed source image
	MatrixType gradMatch(Superclass::m_NumberOfVoxels, Dimension, 0.0);
	MatrixType gradI(Superclass::m_NumberOfVoxels, Dimension);

	// ImageTypePointer DifferenceImage = GridFunctionsType::VectorToImage(m_Image, D);
	// typedef itk::ImageFileWriter<ImageType> WriterType;
	// typename WriterType::Pointer writer = WriterType::New();
	// writer->SetInput(DifferenceImage);
	// writer->SetFileName("ImageDifference.img");
	// writer->Update();

	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		VectorType gradI_d = GridFunctionsType::VectorizeImage(Superclass::m_GradientImages[dim]);
		gradI_d *= Superclass::m_FlipAxes[dim];
		gradI.set_column(Superclass::m_PermutationAxes[dim], gradI_d);
	}


	// MatrixType Yup = this->UpSampleImageMap(Superclass::m_DownSampledY1);
	// for (unsigned int dim = 0; dim < Dimension; dim++)
	// {
	// 	typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
	// 	typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
	// 	derivf->SetInput(this->GetReferenceImage());
	// 	derivf->SetDirection(dim);
	// 	derivf->SetOrder(1);
	// 	derivf->SetUseImageSpacingOn();
	// 	derivf->Update();
	// 	
	// 	ImageTypePointer grad_d = derivf->GetOutput();
	// 	
	// 	VectorType gradI_d = GridFunctionsType::Interpolate(Yup, grad_d);
	// 	
	// 	gradI_d *= Superclass::m_FlipAxes[dim];
	// 	gradI.set_column(Superclass::m_PermutationAxes[dim], gradI_d);		
	// }
	
	// VectorType V(numVoxels);
	// VectorType V1(numVoxels);
	// VectorType V2(numVoxels);
	
	for (int i = 0; i < Superclass::m_NumberOfVoxels; i++)
	{
		TScalar val = VdefImg.get(i) - VdefImgMeanMean.get(i) - VtargImg.get(i) * VConv1.get(i) + VConv2.get(i);
		gradMatch.set_row(i, val * gradI.get_row(i) );
		// V1(i) = VdefImg.get(i) - VdefImgMeanMean.get(i);
		// V2(i) = VtargImg.get(i) * VConv1.get(i) - VConv2.get(i);
		// V(i) = V1(i) - V2(i);
	}

	gradMatch *= 2.0;

	// std::cout << gradMatch.get_column(0) << std::endl << std::endl;
	// 
	// std::cout << gradMatch.get_column(1) << std::endl << std::endl;
	
	// ImageTypePointer VI1 = GridFunctionsType::VectorToImage(Superclass::m_Image, V1);
	// typedef itk::Image<TScalar, Dimension> OutImageType2;
	// typedef itk::ImageFileWriter<OutImageType2> WriterType2;
	// typename WriterType2::Pointer writer2 = WriterType2::New();		
	// writer2->SetInput(VI1);
	// writer2->SetFileName("VI1.img");
	// writer2->Update();
	// 
	// ImageTypePointer VI2 = GridFunctionsType::VectorToImage(Superclass::m_Image, V2);
	// typedef itk::Image<TScalar, Dimension> OutImageType3;
	// typedef itk::ImageFileWriter<OutImageType2> WriterType3;
	// typename WriterType3::Pointer writer3 = WriterType3::New();		
	// writer3->SetInput(VI2);
	// writer3->SetFileName("VI2.img");
	// writer3->Update();
	// 
	// ImageTypePointer VI = GridFunctionsType::VectorToImage(Superclass::m_Image, V);
	// typedef itk::Image<TScalar, Dimension> OutImageType4;
	// typedef itk::ImageFileWriter<OutImageType4> WriterType4;
	// typename WriterType4::Pointer writer4 = WriterType4::New();		
	// writer4->SetInput(VI);
	// writer4->SetFileName("VI.img");
	// writer4->Update();


	// for (int dim = 0; dim < Dimension; dim++)
	// {
	// 	ImageTypePointer Gradient = GridFunctionsType::VectorToImage(Superclass::m_Image, gradMatch.get_column(dim));
	// 	typedef itk::Image<TScalar, Dimension> OutImageType;
	// 	// typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
	// 	// typename CasterType::Pointer castf = CasterType::New();
	// 	// castf->SetInput(Gradient);
	// 	// castf->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// 	// castf->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// 	// castf->Update();
	// 
	// 	typedef itk::ImageFileWriter<OutImageType> WriterType;
	// 	typename WriterType::Pointer writer = WriterType::New();		
	// 	// writer->SetInput(castf->GetOutput());
	// 	writer->SetInput(Gradient);
	// 	if (dim==0)
	// 		writer->SetFileName("gradientX.img");
	// 	else if (dim==1)
	// 		writer->SetFileName("gradientY.img");
	// 	else
	// 		writer->SetFileName("gradientZ.img");
	// 
	// 	writer->Update();
	// }


	// std::cout << "save m_LocalMeanImage" << std::endl;
	// typedef itk::Image<TScalar, Dimension> OutImageType;
	// // typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
	// // typename CasterType::Pointer castf = CasterType::New();
	// // castf->SetInput(m_LocalMeanImage);
	// // castf->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// // castf->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// // castf->Update();
	// 
	// // write output image
	// typedef itk::ImageFileWriter<OutImageType> WriterType;
	// typename WriterType::Pointer writer = WriterType::New();
	// // writer->SetInput(castf->GetOutput());
	// writer->SetInput(m_LocalMeanImage);
	// writer->SetFileName("LocalMeanImage.img");
	// writer->Update();


	return gradMatch;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TScalar, unsigned int Dimension>
typename EQLAImage<TScalar, Dimension>::ImageTypePointer
EQLAImage<TScalar, Dimension>
::ComputeLocalMeanImage(const ImageType* img) const
{
	// convolution between input image and a gaussian filter: W*I
	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetInput(img);
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->Update();

	return gaussf->GetOutput();
}



template<class TScalar, unsigned int Dimension>
typename EQLAImage<TScalar, Dimension>::ImageTypePointer
EQLAImage<TScalar, Dimension>
::ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
		const ImageType* img2, const ImageType* LocalMeanImg2) const
{
	// multiply img1 and img2 voxel-wise in multiplyFilter->GetOutput()
	typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
	multiplyFilter->SetInput1(img1);
	multiplyFilter->SetInput2(img2);
	multiplyFilter->Update();

	// convolve the multiplied image: W*(I.J) in gaussf->GetOutput()
	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->SetInput(multiplyFilter->GetOutput());
	gaussf->Update();

	// multiply the local mean images (W*I).(W*J) in multiplyFilter->GetOutput()
	typename MultiplyImageFilterType::Pointer multiplyFilter2 = MultiplyImageFilterType::New();
	multiplyFilter2->SetInput1(LocalMeanImg1);
	multiplyFilter2->SetInput2(LocalMeanImg2);
	multiplyFilter2->Update();
	
	// subtract the last two images: W*(I.J) - (W*I).(W*J)
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractImageFilterType;
	typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
	subtractFilter->SetInput1(gaussf->GetOutput());
	subtractFilter->SetInput2(multiplyFilter2->GetOutput());
	subtractFilter->Update();
	
	return subtractFilter->GetOutput();

	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// VectorType V = GridFunctionsType::VectorizeImage(subtractFilter->GetOutput());
	// std::cout << V << std::endl;
}



#endif /* _EQLAImage_txx */

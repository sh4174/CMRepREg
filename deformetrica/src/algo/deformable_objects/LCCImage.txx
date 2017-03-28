/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LCCImage_txx
#define _LCCImage_txx

#include "LCCImage.h"

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
LCCImage<TScalar, Dimension>
::LCCImage() : Superclass()
{
	this->SetLCCImageType();
	m_LCCKernelWidth = 5.0;
	// m_SmallestImageDimension = 0;
}



template <class TScalar, unsigned int Dimension>
LCCImage<TScalar, Dimension>
::LCCImage(const LCCImage& o) : Superclass(o)
{
	m_LCCKernelWidth = o.m_LCCKernelWidth;
	m_LocalMeanImage = o.m_LocalMeanImage;
	m_LocalVarianceImage = o.m_LocalVarianceImage;
	// m_SmallestImageDimension = o.m_SmallestImageDimension;
}


template <class TScalar, unsigned int Dimension>
LCCImage<TScalar, Dimension>
::LCCImage(const LCCImage& ex, const MatrixType& IP) : Superclass(ex, IP)
{
	this->SetLCCImageType();
	m_LCCKernelWidth = ex.m_LCCKernelWidth;
	this->Update();
}



template <class TScalar, unsigned int Dimension>
LCCImage<TScalar, Dimension>
::~LCCImage()
{
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation methods
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
LCCImage<TScalar, Dimension>
::Update()
{
	if (this->IsModified())
	{
		Superclass::Update();
		
		m_LocalMeanImage = this->ComputeLocalMeanImage( this->GetImage() );
		m_LocalVarianceImage = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, this->GetImage(), m_LocalMeanImage, 1);
	}
	this->UnSetModified();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int Dimension>
TScalar LCCImage<TScalar, Dimension>
::ComputeMatch(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	// /*
	//   // Deform this image
	// int nbTimePoints = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	// // LCCImage* Idef;
	// // Idef =  dynamic_cast<LCCImage*>(this->GetDeformedObjectAt(nbTimePoints-1));
	//  */
	// ImageTypePointer Img = this->GetImageFromDownSampledImagePoints();

	const LCCImage* targ = dynamic_cast<const LCCImage*>(target);
	
	this->Update();
	// target->Update();
	
	// Get the Local Mean Image and Standard Deviation Image of the target image
	// ImageTypePointer targLocalMean = targ->GetLocalMeanImage();
	// ImageTypePointer targLocalVariance = targ->GetLocalVarianceImage();
	
	// // Compute the Local Mean Image and Standard Deviation Image of the deformed Image
	// ImageTypePointer defImgLocalMean = Idef->GetLocalMeanImage(); // A VERIFIER A VERIFIER
	// ImageTypePointer defImgLocalVariance = Idef->GetLocalVarianceImage();
	// ImageTypePointer ImgLocalMean = this->ComputeLocalMeanImage(Img);
	// ImageTypePointer ImgLocalVariance = this->ComputeLocalCovarianceImage(Img, ImgLocalMean, Img, ImgLocalMean, 1);
	
	// Compute the Local Covariance between deformed source and target
	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage(), 0);

	// Vectorize images
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	VectorType LC = GridFunctionsType::VectorizeImage(LocalCovariance);
	VectorType Var1 = GridFunctionsType::VectorizeImage(m_LocalVarianceImage);
	VectorType Var2 = GridFunctionsType::VectorizeImage(targ->GetLocalVarianceImage());
	
	// for (int i = 0; i < 100; i++)
	// {
	// 	std::cout << "LC = " << LC.get(i) << " Var1 = " << Var1.get(i) << " Var2 = " << Var2.get(i) << " match = " << ( LC.get(i) * LC.get(i) / ( Var2.get(i) * Var1.get(i) ) ) << std::endl;
	// }
	
	TScalar match = 0.0;
	for (int i = 0; i < LC.size(); i++)
	{
		match += ( LC.get(i) * LC.get(i) / ( Var2.get(i) * Var1.get(i) ) );
	}
	match /= Superclass::m_NumberOfVoxels; // normalize by the number of voxels
	match = 1 - match;

	std::cout << "match = " << match << std::endl;
	
	return match;
}



template<class TScalar, unsigned int Dimension>
typename LCCImage<TScalar, Dimension>::MatrixType LCCImage<TScalar, Dimension>
::ComputeMatchGradient(const DeformableObjectType* target)
{
	if (this->GetType() != target->GetType())
		std::cerr << "Deformable objects types mismatched: " << this->GetType() << " and " << target->GetType() << "\n";

	/*
	int nbTimePoints = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	// LCCImage* Idef;
	// Idef = dynamic_cast<LCCImage*>(this->GetDeformedObjectAt(nbTimePoints-1));
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
	*/
	const LCCImage* targ = dynamic_cast<const LCCImage*>(target);

	this->Update();
	// target->Update();
	
	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// ImageTypePointer Img = this->GetImageFromDownSampledImagePoints(); //GridFunctionsType::VectorToImage(Superclass::m_Image,GridFunctionsType::Interpolate(Yfinal, Superclass::m_Image));
	
//	MatrixType Yfinal = GridFunctionsType::UpsampleImagePoints(this->GetImage(), this->GetDownSampledWorkingImage(), Superclass::m_DownSampledY1);
	// MatrixType Yfinal = this->UpSampleImageMap(Superclass::m_DownSampledY1);
	// 
	// // Get the Local Mean Image and Standard Deviation Image of the target image
	// ImageTypePointer targLocalMean = targ->GetLocalMeanImage();
	// ImageTypePointer targLocalVariance = targ->GetLocalVarianceImage();
	
	// Compute the Local Mean Image and Standard Deviation Image of the deformed Image
	// ImageTypePointer ImgLocalMean = this->ComputeLocalMeanImage(Img);
	// ImageTypePointer ImgLocalVariance = this->ComputeLocalCovarianceImage(Img, ImgLocalMean, Img, ImgLocalMean, 1);
	
	// Compute the Local Covariance between deformed source and target
	// ImageTypePointer targImg =  targ->GetWorkingImage();
	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage(), 0);

	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyImageFilterType;
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractImageFilterType;
	typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DivideImageFilterType;
	
	// Compute Cov(Idef,targ)/Var(Idef)
	typename DivideImageFilterType::Pointer divideFilter1 = DivideImageFilterType::New();
	divideFilter1->SetInput1(LocalCovariance);
	divideFilter1->SetInput2(m_LocalVarianceImage);
	divideFilter1->Update();
	
	// Compute Cov(Idef,targ) / (Var(Idef) * Var(targ))
	typename DivideImageFilterType::Pointer divideFilter2 = DivideImageFilterType::New();
	divideFilter2->SetInput1(divideFilter1->GetOutput());
	divideFilter2->SetInput2(targ->GetLocalVarianceImage());
	divideFilter2->Update();
	
	
	// Compute W*( Cov(Idef,targ)/ ( Var(Idef) * Var(targ) )
	typename GaussianFilterType::Pointer gaussf1 = GaussianFilterType::New();
	gaussf1->SetSigma( m_LCCKernelWidth );
	gaussf1->SetInput(divideFilter2->GetOutput());
	gaussf1->Update();
	
	ImageTypePointer Convolution1 = gaussf1->GetOutput();
	
	// Compute W*( Cov(Idef,targ)^2 / (Var(Idef)^2 * Var(targ) )
	typename MultiplyImageFilterType::Pointer multiplyFilter1 = MultiplyImageFilterType::New();
	multiplyFilter1->SetInput1(divideFilter1->GetOutput());
	multiplyFilter1->SetInput2(divideFilter2->GetOutput());
	multiplyFilter1->Update();

	typename GaussianFilterType::Pointer gaussf2 = GaussianFilterType::New();
	gaussf2->SetSigma( m_LCCKernelWidth );
	gaussf2->SetInput(multiplyFilter1->GetOutput());
	gaussf2->Update();
	
	ImageTypePointer Convolution2 = gaussf2->GetOutput();
	
	// Compute Cov(Idef,targ)^2 * defImgLocalMean / ( Var(Idef)^2 * Var(targ) )
	typename MultiplyImageFilterType::Pointer multiplyFilter2 = MultiplyImageFilterType::New();
	multiplyFilter2->SetInput1(multiplyFilter1->GetOutput());
	multiplyFilter2->SetInput2(m_LocalMeanImage);
	multiplyFilter2->Update();

	// Compute targLocalMean * Cov(Idef,targ) / ( Var(Idef) * Var(targ) )
	typename MultiplyImageFilterType::Pointer multiplyFilter3 = MultiplyImageFilterType::New();
	multiplyFilter3->SetInput1(divideFilter2->GetOutput());
	multiplyFilter3->SetInput2(targ->GetLocalMeanImage());
	multiplyFilter3->Update();
	
	// Subtract the last two and convolve the output
	typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
	subtractFilter->SetInput1(multiplyFilter2->GetOutput());
	subtractFilter->SetInput2(multiplyFilter3->GetOutput());
	subtractFilter->Update();
	
	typename GaussianFilterType::Pointer gaussf3 = GaussianFilterType::New();
	gaussf3->SetSigma( m_LCCKernelWidth );
	gaussf3->SetInput(subtractFilter->GetOutput());
	gaussf3->Update();

	ImageTypePointer Convolution3 = gaussf3->GetOutput();
	
	
	// Vectorize required images to compute the gradient
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	VectorType VtargImg = GridFunctionsType::VectorizeImage(targ->GetImage());
	VectorType VdefImg = GridFunctionsType::VectorizeImage(this->GetImage());
	VectorType VConv1 = GridFunctionsType::VectorizeImage(Convolution1);
	VectorType VConv2 = GridFunctionsType::VectorizeImage(Convolution2);
	VectorType VConv3 = GridFunctionsType::VectorizeImage(Convolution3);
	
	
	// for (int i = 0; i < 100; i++)
	// {
	// 	std::cout << VtargImg.get(i) << " " << VConv1.get(i) << " " <<  VdefImg.get(i) << " " << VConv2.get(i) << " " << VConv3.get(i) << "total = " << -1.0 * VtargImg.get(i) * VConv1.get(i) + VdefImg.get(i) * VConv2.get(i) - VConv3.get(i) << std::endl;
	// }
	
	// Compute the gradient of the deformed source image
	MatrixType gradMatch(Superclass::m_NumberOfVoxels, Dimension, 0.0);
	MatrixType gradI(Superclass::m_NumberOfVoxels, Dimension);

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
	
	for (int i = 0; i < Superclass::m_NumberOfVoxels; i++)
	{
		TScalar val = -1.0 * VtargImg.get(i) * VConv1.get(i) + VdefImg.get(i) * VConv2.get(i) - VConv3.get(i);
		gradMatch.set_row(i, val * gradI.get_row(i) );
	}

	gradMatch *= ( 2.0 / Superclass::m_NumberOfVoxels );

	// std::cout << gradMatch.get_column(0) << std::endl << std::endl;
	// 
	// std::cout << gradMatch.get_column(1) << std::endl << std::endl;
	

	// for (int dim = 0; dim < Dimension; dim++)
	// {
	// 	ImageTypePointer Gradient = GridFunctionsType::VectorToImage(Superclass::m_Image, gradMatch.get_column(dim));
	// 	typedef itk::Image<TScalar, Dimension> OutImageType;
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

	// for (int dim = 0; dim < Dimension; dim++)
	// {
	// 	ImageTypePointer Gradient = GridFunctionsType::VectorToImage(Superclass::m_Image, gradMatch.get_column(dim));
	// 	typedef itk::Image<unsigned char, Dimension> OutImageType;	
	// 	typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
	// 	typename CasterType::Pointer castf = CasterType::New();
	// 	castf->SetInput(Gradient);
	// 	castf->SetOutputMinimum(Superclass::m_MinIntensityOutput);
	// 	castf->SetOutputMaximum(Superclass::m_MaxIntensityOutput);
	// 	castf->Update();
	// 
	// 	typedef itk::ImageFileWriter<OutImageType> WriterType;
	// 	typename WriterType::Pointer writer = WriterType::New();		
	// 	writer->SetInput(castf->GetOutput());
	// 	if (dim==0)
	// 		writer->SetFileName("gradientX.png");
	// 	else if (dim==1)
	// 		writer->SetFileName("gradientY.png");
	// 	else
	// 		writer->SetFileName("gradientZ.png");
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
typename LCCImage<TScalar, Dimension>::ImageTypePointer
LCCImage<TScalar, Dimension>
::ComputeLocalMeanImage(const ImageType* img) const
{
	// convolution between input image and a gaussian filter: W*I
	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetInput(img);
	gaussf->SetSigma( m_LCCKernelWidth );
	gaussf->Update();

	return gaussf->GetOutput();
}



template<class TScalar, unsigned int Dimension>
typename LCCImage<TScalar, Dimension>::ImageTypePointer
LCCImage<TScalar, Dimension>
::ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
		const ImageType* img2, const ImageType* LocalMeanImg2, bool sameImage) const
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
	gaussf->SetSigma( m_LCCKernelWidth );
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
	
	if (sameImage)
	{
		// We want to avoid division-by-zero in areas of constant intensities. For that, we consider each image corrupted by a white noise with tiny variance \epsilon.
		// This amounts to add \epsilon to the local variance of each image (sameImage = true) and do nothing for the local correlation between distinct images (sameImage = false) since white noise are independent for each image
		typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddImageFilterType;
		typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
		addFilter->SetInput1(subtractFilter->GetOutput());
		addFilter->SetConstant2(pow(10,-7));
		addFilter->Update();
		return addFilter->GetOutput();
	}
	else
	{
		return subtractFilter->GetOutput();
	}
	
	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// VectorType V = GridFunctionsType::VectorizeImage(subtractFilter->GetOutput());
	// std::cout << V << std::endl;

}


#endif /* _LCCImage_txx */

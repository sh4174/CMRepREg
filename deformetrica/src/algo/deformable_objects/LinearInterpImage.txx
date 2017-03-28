/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LinearInterpImage_txx
#define _LinearInterpImage_txx

#include "LinearInterpImage.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "GridFunctions.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkOrientImageFilter.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
LinearInterpImage<TScalar, Dimension>
::LinearInterpImage() : Superclass()
{
	m_MinIntensityOutput = 0;
	m_MaxIntensityOutput = (Dimension==2)?255:1.0;

	std::vector<int> p(Dimension);
	std::vector<int> f(Dimension);
	for (unsigned int d=0; d<Dimension; d++)
	{
		p[d] = d;
		f[d] = 1;
	}
	m_PermutationAxes = p;
	m_FlipAxes = f;
	
}



template <class TScalar, unsigned int Dimension>
LinearInterpImage<TScalar, Dimension>
::LinearInterpImage(const LinearInterpImage& o) : Superclass(o)
{
    typedef typename itk::ImageDuplicator< ImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

	if (!o.m_Image.IsNull())
	{
	    duplicator->SetInputImage(o.m_Image);
	    duplicator->Update();
	    m_Image = duplicator->GetOutput();		
	}

	if (!o.m_DownSampledImage.IsNull())
	{
		duplicator->SetInputImage(o.m_DownSampledImage);
		duplicator->Update();
		m_DownSampledImage = duplicator->GetOutput();
	}
	
	m_GradientImages.resize(Dimension);
	for (int dim = 0; dim < Dimension; dim++)
	{
		if (!o.m_GradientImages[dim].IsNull())
		{
			duplicator->SetInputImage(o.m_GradientImages[dim]);
			duplicator->Update();
			m_GradientImages[dim] = duplicator->GetOutput();
		}
	}
	
	m_MinIntensityOutput = o.m_MinIntensityOutput;
	m_MaxIntensityOutput = o.m_MaxIntensityOutput;

	m_PermutationAxes = o.m_PermutationAxes;
	m_FlipAxes = o.m_FlipAxes;
	
	m_NumberOfVoxels = o.m_NumberOfVoxels;
	
}


template <class TScalar, unsigned int Dimension>
LinearInterpImage<TScalar, Dimension>
::LinearInterpImage(const LinearInterpImage& example, const MatrixType& DownSampledImageMap) : Superclass(example, DownSampledImageMap)
{
	if (example.m_Image.IsNull() || example.m_DownSampledImage.IsNull())
		throw std::runtime_error("Could not resample image if no image has been set!");
	
    typedef typename itk::ImageDuplicator< ImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
	
    duplicator->SetInputImage(example.m_Image);
    duplicator->Update();
    m_Image = duplicator->GetOutput();		
	
	duplicator->SetInputImage(example.m_DownSampledImage);
	duplicator->Update();
	m_DownSampledImage = duplicator->GetOutput();
	
	m_GradientImages.resize(Dimension);
	for (int dim = 0; dim < Dimension; dim++)
	{
		duplicator->SetInputImage(example.m_GradientImages[dim]);
		duplicator->Update();
		m_GradientImages[dim] = duplicator->GetOutput();
	}
		
	m_MinIntensityOutput = example.m_MinIntensityOutput;
	m_MaxIntensityOutput = example.m_MaxIntensityOutput;

	m_PermutationAxes = example.m_PermutationAxes;
	m_FlipAxes = example.m_FlipAxes;
	
	m_NumberOfVoxels = example.m_NumberOfVoxels;
	
	// Upsample ImagePoints to be at the same resolution as the reference image
	if (example.m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels() != DownSampledImageMap.rows())
		throw std::runtime_error("Number of voxels in downsampled image maps mismatch in LinearInterpImage constructor");
	
	MatrixType Yup = this->UpSampleImageMap(DownSampledImageMap);
	VectorType I1 = GridFunctionsType::Interpolate(Yup, m_Image);
	
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType it(m_Image, m_Image->GetLargestPossibleRegion());

    unsigned int r = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      it.Set(I1[r++]);
    }

	for (int dim = 0; dim < Dimension; dim++)
	{
		VectorType I1g = GridFunctionsType::Interpolate(Yup, m_GradientImages[dim]);
	
	    IteratorType itg(m_GradientImages[dim], m_GradientImages[dim]->GetLargestPossibleRegion());

	    unsigned int r = 0;
	    for (itg.GoToBegin(); !itg.IsAtEnd(); ++itg)
	    {
	      itg.Set(I1g[r++]);
	    }
		
	}

	// DownsampledReferenceImage is not re-sampled since only the positions of its voxels is used.
	// Call to Update() is done in children classes
}




template <class TScalar, unsigned int Dimension>
LinearInterpImage<TScalar, Dimension>
::~LinearInterpImage()
{
	// for (int dim = 0; dim < Dimension; dim++)
	// 	delete m_GradientImages[dim];
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
template <class TScalar, unsigned int Dimension>
void LinearInterpImage<TScalar, Dimension>
::SetOriginalImage(ImageType *img)
{
	ImageTypePointer oimg;
	void *poimg = &oimg;

	//timg = ReorientImage(img);

	if (img != NULL)
	{
		m_OrigImage = img;

		std::cout << "Original image " << img << std::endl;
		std::cout << "Original direction "<< " " << img->GetDirection() << std::endl;

		ReorientImage(img, poimg);

	}
	// Added by ABGF@03mars14

	if (oimg.IsNotNull())
	{
		std::cout << "m_Image (LPS) is a copy of oriented image " << oimg.GetPointer() << std::endl;
		// Set working image
		m_Image = oimg;
		//std::cout << "m_Image (LPS): " << wimg << std::endl;
		//std::cout << "m_Image (LPS): " << wimg->GetDirection() << std::endl;
	}
	else
	{
		std::cout << "m_Image is a copy of original image " << img << std::endl;
		m_Image = img;
		//std::cout << "m_Image (original): " << 		Superclass::m_Image->GetDirection() << std::endl;
	}

	this->SetImage();
	//this->SetModified();
}
*/

template <class TScalar, unsigned int Dimension>
void
LinearInterpImage<TScalar, Dimension>
::SetImageAndDownSamplingFactor(ImageType* img, int ds)
{
	if (ds<0)
		throw std::runtime_error("Negative downsampling factor not allowed");
	
	m_Image = img;
	m_NumberOfVoxels = img->GetLargestPossibleRegion().GetNumberOfPixels();
	
	if (ds==0 || ds == 1)
	{
		m_DownSampledImage = m_Image; // I don't think we need to duplicate
	}
	else
	{
		m_DownSampledImage = GridFunctionsType::DownsampleImage(m_Image, ds);
	}
	
	
	m_GradientImages.resize(Dimension);
	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
		typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
		derivf->SetInput(m_Image);
		derivf->SetDirection(dim);
		derivf->SetOrder(1);
		derivf->SetUseImageSpacingOn();
		derivf->Update();
	
		m_GradientImages[dim] = derivf->GetOutput();
	}
	
	this->SetModified();
}


template <class TScalar, unsigned int Dimension>
void
LinearInterpImage<TScalar, Dimension>
::UpdateImageIntensity(const MatrixType& I)
{
	if (m_Image.IsNull())
		throw std::runtime_error("ITK image should have been set before setting new intensities");

	if (I.rows() != m_NumberOfVoxels)
		throw std::runtime_error("image size and number of pixels mismatched");
	
	if (I.cols() != 1)
		throw std::runtime_error("only scalar images allowed");
		
	int i=0;
	
	typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
	ImageIterator it(m_Image, m_Image->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		it.Set( I(i++,0) );

	this->SetModified(); 
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
template <class TScalar, unsigned int Dimension>
void LinearInterpImage<TScalar, Dimension>
::ReorientImage(ImageType* img, void* oimg)
{
   if ((Dimension > 2) && (img != NULL)) // orientation does not make sense with 2D data ?
   {
      /// If image has a non-canonical anatomical orientation, save original orientation and then reorient it to LPS referential.
      itk::Matrix<double,Dimension,Dimension> dir = img->GetDirection();
      itk::Matrix<double,Dimension,Dimension> id;

      id.SetIdentity();

      // not LPS orientation ?
      if (dir != id)
      {
    	  MatrixType orient_mtx(Dimension, Dimension, 0.0);
    	  for (unsigned int i = 0; i < Dimension; i++)
    		  for (unsigned int j = 0; j < Dimension; j++)
    			  	  orient_mtx[i][j] = dir(i,j);


    	  //m_anatOrient.LabelToChangeOfBasisMatrix("LPS"))GetChangeOfBasisMatrix()
	      m_anatOrient.SetChangeOfBasisMatrix(orient_mtx);
	      std::cout << "Deformable Object is an image originally oriented in the " << m_anatOrient.GetAnatomicalCoordinateSystemLabel() << " coordinate system." << std::endl;
	      std::cout << "Reorienting image to match the LPS standard." << std::endl;

	      // reorient image (LPS orientation / RAI code for ITK filter)
	      typename ImageType::Pointer *outimg = (typename ImageType::Pointer *)oimg;
	      typename itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
	      orienter->UseImageDirectionOn();
	      orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
	      orienter->SetInput(img);
	      orienter->Update();
	      *outimg = orienter->GetOutput();
	      std::cout << "Reoriented image matrix (LPS): " << outimg << std::endl; // << " " << (*outimg)->GetDirection() << std::endl;
      }
   }
}
*/

template <class TScalar, unsigned int Dimension>
void 
LinearInterpImage<TScalar, Dimension>
::Update()
{
	if (m_Image.IsNull() || m_DownSampledImage.IsNull())
		throw std::runtime_error("Reference image and downsampled version of it should be set in LinearInterpImage");

	if (this->IsModified())
	{
			
		this->UpdateBoundingBox();

		// // If m_DownSampledY1 has been set, it is used to re-sample the reference image, otherwise one assumes m_DownSampledY1 to be a regular sampling of m_DownSampledReferenceImage (i.e. no deformation)	
		// if (m_DownSampledY1.rows() == 0)
		// {
		// 	m_DownSampledY1.set_size(m_DownSampledReferenceImage->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension);
		// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
		// 	m_DownSampledY1 = GridFunctionsType::ImageToPoints(m_DownSampledReferenceImage);
		// }
		// 
		// this->UpdateWorkingImageFromDownSampledImagePoints();
	}
	// this->UnSetModified() not necessary, since this method could only be called by children classes, which will take care of this.
}



template<class TScalar, unsigned int Dimension>
typename LinearInterpImage<TScalar, Dimension>::MatrixType
LinearInterpImage<TScalar, Dimension>
::SplatDifferenceImage(const LinearInterpImage<TScalar, Dimension>* target, const MatrixType& DownSampledImageMap )
{
	MatrixType Y0 = this->UpSampleImageMap(DownSampledImageMap);
	VectorType I0 = GridFunctionsType::Interpolate(Y0, m_Image);
	VectorType I1 = GridFunctionsType::VectorizeImage(target->GetImage());
	VectorType Residual = I0 - I1;

	ImageTypePointer splat = GridFunctionsType::SplatToImage(m_Image, Y0, Residual, 0);
	
	MatrixType out(m_NumberOfVoxels,1);
	
	out.set_column(0 ,GridFunctionsType::VectorizeImage(splat));
	
	return out;
}



/*
template <class TScalar, unsigned int Dimension>
void LinearInterpImage<TScalar, Dimension>
::WriteObject(std::string str)
{
	/// rescale and convert intensities of the working image, then save it.
	/// If initial intensity range is not included in [0,1], intensity values are converted to unsigned short. Otherwise, Tscalar (float or double) is used
	if ( (m_MinIntensityOutput >= 0.0) && (m_MaxIntensityOutput > 1.0) )
	{
		typedef itk::Image<unsigned char, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// write output image
		typedef itk::ImageFileWriter<OutImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput(castf->GetOutput());
		writer->SetFileName(str.c_str());
		writer->Update();
	}
	else
	{
		typedef itk::Image<TScalar, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// write output image
		typedef itk::ImageFileWriter<OutImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput(castf->GetOutput());
		writer->SetFileName(str.c_str());
		writer->Update();
	}
	
}*/



template <class TScalar, unsigned int Dimension>
void LinearInterpImage<TScalar, Dimension>
::WriteObject(std::string str) const
{
	// Rescale and convert intensities of the working image, then save it.
	// If initial intensity range is not included in [0,1], intensity values are converted to unsigned short.
	// Otherwise, Tscalar (float or double) is used
	if ( (m_MinIntensityOutput >= 0.0) && (m_MaxIntensityOutput > 1.0) )
	{
		typedef itk::Image<unsigned char, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// std::cout << "Defined cast filter" << std::endl;
		// If original image (> 2D) was not oriented according to LPS reference, switch back to original coordinate system
		std::string lps_label("LPS");

		if ((Dimension ==3) && (lps_label.compare(Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel()) != 0))
		{
			typedef itk::Image<unsigned char, 3> OutImageType_dim_3;
			typename itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::Pointer orienter = itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::New();
			orienter->UseImageDirectionOn();
			orienter->SetDesiredCoordinateOrientation(Superclass::m_AnatomicalOrientation.GetITKAnatomicalCoordinateSystemLabel());
			orienter->SetInput(dynamic_cast<OutImageType_dim_3*>( castf->GetOutput() ));
			orienter->Update();

			// std::cout << "Orientation filter" << std::endl;
			typedef itk::ImageFileWriter<OutImageType_dim_3> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(orienter->GetOutput()); //
			writer->SetFileName(str.c_str());
			writer->Update();
			// std::cout << "Writing filter for oriented image" << std::endl;
		}
		else
		{
			// write output image
			typedef itk::ImageFileWriter<OutImageType> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(castf->GetOutput());
			writer->SetFileName(str.c_str());
			writer->Update();
			// std::cout << "Writing filter" << std::endl;
		}
	}
	else
	{
		typedef itk::Image<TScalar, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// If original image was not oriented according to LPS reference, switch back to original coordinate system
		std::string lps_label("LPS");

		if ((Dimension ==3) && (lps_label.compare(Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel()) != 0))
		{
			std::cout << "Reorienting output image from LPS orientation to the original " << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " orientation." << std::endl;

			typedef itk::Image<TScalar, 3> OutImageType_dim_3;
			typename itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::Pointer orienter = itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::New();
			orienter->UseImageDirectionOn();
			orienter->SetDesiredCoordinateOrientation(Superclass::m_AnatomicalOrientation.GetITKAnatomicalCoordinateSystemLabel());
			orienter->SetInput(dynamic_cast<OutImageType_dim_3*>( castf->GetOutput() ));
			orienter->Update();

			typedef itk::ImageFileWriter<OutImageType_dim_3> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(orienter->GetOutput()); //
			writer->SetFileName(str.c_str());
			writer->Update();
		}
		else
		{
	    // write output image
		typedef itk::ImageFileWriter<OutImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput(castf->GetOutput()); //
		writer->SetFileName(str.c_str());
		writer->Update();
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void 
LinearInterpImage<TScalar, Dimension>
::UpdateBoundingBox()
{
	// Not really efficient... but the safiest given possible weird coordinate systems...
	MatrixType Yex = GridFunctionsType::ImageToPoints(m_Image);
	
	for (unsigned int d=0; d<Dimension; d++)
	{
		Superclass::m_BoundingBox(d,0) = Yex.get_column(d).min_value();
		Superclass::m_BoundingBox(d,1) = Yex.get_column(d).max_value();
	}
	
}



template <class TScalar, unsigned int Dimension>
typename LinearInterpImage<TScalar, Dimension>::MatrixType
LinearInterpImage<TScalar, Dimension>
::UpSampleImageMap(const MatrixType& Y)
{
	if (Y.rows() != m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels())
		throw std::runtime_error("Number points in image map and size of downsampled image mismatch");
	
	return GridFunctionsType::UpsampleImagePoints(m_Image, m_DownSampledImage, Y);
	
	// MatrixType Yup(m_NumberOfVoxels, Dimension, 0.0);
	// 
	// MatrixType ImageMap = GridFunctionsType::ImageToPoints(m_Image);
	// for (unsigned int d = 0; d < Dimension; d++)
	// {
	// 	ImageTypePointer Hd = GridFunctionsType::VectorToImage(m_DownSampledImage, Y.get_column(d));
	// 	Yup.set_column(d, GridFunctionsType::Interpolate(ImageMap, Hd));
	// }

	// 
	// /*
	// MatrixType otherYup = GridFunctionsType::UpsampleImagePoints(m_Image, m_DownSampledWorkingImage, m_DownSampledY1);
	// TScalar l2_error = 0.0;
	// for(int i = 0; i < otherYup.rows(); i++)
	// 	l2_error += (otherYup.get_row(i) - Yup.get_row(i)).squared_magnitude();
	// if ( l2_error > 0 )
	// 	std::cerr << "In LinearInterpImage::UpSampleImageMap - l2_error = " << l2_error << std::endl;
	//  */
	// 
	// return Yup;
}




// template <class TScalar, unsigned int Dimension>
// typename LinearInterpImage<TScalar, Dimension>::ImageTypePointer
// LinearInterpImage<TScalar, Dimension>
// ::RescaleImage(TScalar minVal, TScalar maxVal)
// {
// 	typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
// 	
// 	typename RescalerType::Pointer resf = RescalerType::New();
// 	resf->SetInput(m_Image);
// 	resf->SetOutputMinimum(minVal);
// 	resf->SetOutputMaximum(maxVal);
// 	resf->Update();
//   
// 	// this->SetImage(resf->GetOutput()); // this also sets m_Modified = true;
// 	return resf->GetOuptput();
// }



// template <class TScalar, unsigned int Dimension>
// void
// LinearInterpImage<TScalar, Dimension>
// ::UpdateWorkingImageFromDownSampledImagePoints()
//  {
// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 
// 	MatrixType Y1 = this->UpSampleImageMap(m_DownSampledY1);
// 	VectorType I1 = GridFunctionsType::Interpolate(Y1, m_ReferenceImage);
// 	m_WorkingImage = GridFunctionsType::VectorToImage(m_ReferenceImage, I1);
// 
// }



// template<class TScalar, unsigned int Dimension>
// void LinearInterpImage<TScalar, Dimension>
// ::TransportAlongGeodesicWithJumps(MatrixList& gradDataTi, std::vector<int>& jumpTimes, MatrixList& EtaT, MatrixList& YT)
// {
// 	// The path of the image points over time
// 	YT = this->GetFlow();
// 
// 	long numTimePoints = Superclass::m_Def->GetNumberOfTimePoints();
// 	long numVoxels = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
// 	TScalar dt = (Superclass::m_Def->GetTN() - Superclass::m_Def->GetT0()) / (numTimePoints-1);
// 
// 	// The aux variable Eta over time
// 	EtaT.resize(numTimePoints);
// 	// The source term for integration
// 	MatrixType zeroMat(this->GetNumberOfPoints(), Dimension, 0.0);
// 	EtaT[numTimePoints - 1] = zeroMat;
// 
// 	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
// 
// 	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
// 	typedef typename KernelFactoryType::KernelBaseType KernelType;
// 	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();
// 
// 	KernelType* kernelObj = kFactory->CreateKernelObject();
// 	kernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());
// 
// 	// The path of control points and momenta over time
// 	MatrixList PositionsT = Superclass::m_Def->GetTrajectoryPositions();
// 	MatrixList MomentasT = Superclass::m_Def->GetTrajectoryMomentas();
// 
// 	int subjIndex = jumpTimes.size()-1;
// 
// 	// Integrate backwards in time
// 	for (long t = (numTimePoints - 1); t>=1; t--)
// 	{
// 		kernelObj->SetSources(PositionsT[t]);
// 		kernelObj->SetWeights(MomentasT[t]);
// 
// 		// The velocity is always v_t(y0)
// 		MatrixType VtY0 = kernelObj->Convolve(m_DownSampledY1);
// 
// 		// This gives us a Dim*Dim matrix at every pixel of the original image
// 		std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(m_DownSampledY1);
// 
// 		// We need to have eta(t) in the form of an image, so we can compute the jacobian
// 		std::vector<ImageTypePointer> etaTImage;
// 		etaTImage.resize(Dimension);
// 		for (unsigned int d = 0; d < Dimension; d++)
// 		{
// 			// Create an image from each Dimension
// 			ImageTypePointer img = GridFunctionsType::VectorToImage(m_Image, EtaT[t].get_column(d));
// 			etaTImage[d] = img;
// 		}
// 
// 		// Compute jacobian of eta(t)
// 		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
// 		std::vector<ImageTypePointer> gradImages;
// 		gradImages.resize(Dimension*Dimension);
// 		unsigned int indx = 0;
// 		// Compute the gradient in all directions
// 		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 		{
// 			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 			{
// 				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
// 				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
// 				derivf->SetInput(etaTImage[dim1]);
// 				derivf->SetDirection(dim2);
// 				derivf->SetOrder(1);
// 				derivf->SetUseImageSpacingOn();
// 				derivf->Update();
// 
// 				gradImages[indx] = derivf->GetOutput();
// 				indx++;
// 			}
// 		}
// 
// 		// Now we have all the information we need to compute dEta, but we need to loop over all the pixels
// 		// and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
// 		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
// 		IteratorType it(m_DownSampledWorkingImage, m_DownSampledWorkingImage->GetLargestPossibleRegion());
// 		int matrixIndex = 0;
// 		MatrixType dEta(numVoxels, Dimension, 0);
// 
// 		// Loop over the grid to construct jacobian matrices and compute dY
// 		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
// 		{
// 			ImageIndexType imageIndex = it.GetIndex();
// 			MatrixType jacobian(Dimension, Dimension);
// 
// 			unsigned int tempMatrixIndex = 0;
// 			// Build the (DimensionxDimension) jacobian matrix
// 			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 			{
// 				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 				{
// 					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
// 					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
// 					tempMatrixIndex++;
// 				}
// 			}
// 
// 			// The trace of the first term
// 			TScalar traceValue = trace(firstTerm[matrixIndex]);
// 
// 			// Get the (Dimension) vector corresponding to this points v_t(y0)
// 			VectorType VtY0k = VtY0.get_row(matrixIndex);
// 
// 			// The final value for this pixel of dEta
// 			dEta.set_row(matrixIndex, traceValue*EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);
// 
// 			matrixIndex++;
// 		}
// 
// 		// Add the jump condition
// 		if (t == jumpTimes[subjIndex])
// 		{
// 			dEta += gradDataTi[subjIndex];
// 			subjIndex--;
// 			if (subjIndex < 0)
// 				subjIndex = 0;
// 		}
// 
// 		EtaT[t-1] = EtaT[t] + dEta * dt;
// 	}
// 
// 	// We have computed Eta(t), the backwards propagator actually needs (d_yp Y)^t Eta(t)
// 	for (long t = 0; t<numTimePoints; t++)
// 	{
// 		// We need to have Eta(t) in the form of an image, so we can compute the jacobian
// 		std::vector<ImageTypePointer> YTImage;
// 		YTImage.resize(Dimension);
// 		for (unsigned int d = 0; d < Dimension; d++)
// 		{
// 			// Create an image from each Dimension
// 			ImageTypePointer img = GridFunctionsType::VectorToImage(m_Image, YT[t].get_column(d));
// 			YTImage[d] = img;
// 		}
// 
// 		// Compute jacobian of eta(t)
// 		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
// 		std::vector<ImageTypePointer> gradImages;
// 		gradImages.resize(Dimension*Dimension);
// 		unsigned int indx = 0;
// 		// Compute the gradient in all directions
// 		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 		{
// 			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 			{
// 				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
// 				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
// 				derivf->SetInput(YTImage[dim1]);
// 				derivf->SetDirection(dim2);
// 				derivf->SetOrder(1);
// 				derivf->SetUseImageSpacingOn();
// 				derivf->Update();
// 
// 				gradImages[indx] = derivf->GetOutput();
// 				indx++;
// 			}
// 		}
// 
// 		// Now we need to loop over all the pixels and construct the 3x3 jacobian matrix
// 		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
// 		IteratorType it(m_DownSampledWorkingImage, m_DownSampledWorkingImage->GetLargestPossibleRegion());
// 		int matrixIndex = 0;
// 
// 		// Loop over the grid to construct jacobian matrices and compute dY
// 		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
// 		{
// 			ImageIndexType imageIndex = it.GetIndex();
// 			MatrixType jacobian(Dimension, Dimension);
// 
// 			unsigned int tempMatrixIndex = 0;
// 			// Build the (DimensionxDimension) jacobian matrix
// 			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
// 			{
// 				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
// 				{
// 					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
// 					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
// 					tempMatrixIndex++;
// 				}
// 			}
// 
// 			// The final value for this pixel of dEta
// 			EtaT[t].set_row(matrixIndex, jacobian.transpose()*EtaT[t].get_row(matrixIndex));
// 
// 			matrixIndex++;
// 		}
// 
// 	}
// }



#endif /* _LinearInterpImage_txx */

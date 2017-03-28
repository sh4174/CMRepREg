/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectReader_txx
#define _DeformableObjectReader_txx

#include "LinearInterpImage.h"
#include "SSDImage.h"
#include "LCCImage.h"
#include "EQLAImage.h"
#include "Landmark.h"
#include "PointCloud.h"
#include "OrientedPolyLine.h"
#include "NonOrientedPolyLine.h"
#include "OrientedSurfaceMesh.h"
#include "NonOrientedSurfaceMesh.h"
#include "CMRep.h"

#include "KernelFactory.h"

#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNiftiImageIO.h"

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "itksys/SystemTools.hxx"


#include "itkOrientImageFilter.h"
//#include "itkMatrix.h"


#define instantiateCurrentOrVarifold(name) \
	{ \
		vtkSmartPointer<vtkPolyDataReader> reader = \
				vtkSmartPointer<vtkPolyDataReader>::New(); \
		reader->SetFileName(m_FileName); \
		reader->Update(); \
		vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput(); \
\
		object##name->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem()); \
		object##name->SetPolyData(PolyData); \
\
		m_Min = object##name->GetMin(); \
		m_Max = object##name->GetMax(); \
\
		m_Bbox.set_column(0,m_Min); \
		m_Bbox.set_column(1,m_Max); \
		typedef KernelFactory<TScalar, Dimension> KernelFactoryType; \
		KernelFactoryType* kfac = KernelFactoryType::Instantiate(); \
		kfac->SetDataDomain(m_Bbox); \
\
		object##name->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str())); \
		object##name->SetKernelWidth(m_ParamObject->GetKernelWidth()); \
	}




////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformableObjectReader<TScalar, Dimension>
::DeformableObjectReader()
 {
	//m_IsTemplate = false;
 }



template <class TScalar, unsigned int Dimension>
DeformableObjectReader<TScalar, Dimension>
::~DeformableObjectReader()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableObjectReader<TScalar, Dimension>
::Update()
 {
	typedef Landmark<TScalar, Dimension> LandmarkType;
	typedef PointCloud<TScalar, Dimension> PointCloudType;
	typedef OrientedPolyLine<TScalar, Dimension> OrientedPolyLineType;
	typedef NonOrientedPolyLine<TScalar, Dimension> NonOrientedPolyLineType;
	typedef OrientedSurfaceMesh<TScalar, Dimension> OrientedSurfaceMeshType;
	typedef NonOrientedSurfaceMesh<TScalar, Dimension> NonOrientedSurfaceMeshType;
	typedef CMRep<TScalar, Dimension> CMRepType;
	
	typedef LinearInterpImage<TScalar, Dimension> LinearInterpImageType;
	typedef SSDImage<TScalar, Dimension> SSDImageType;
	typedef LCCImage<TScalar, Dimension> LCCImageType;
	typedef EQLAImage<TScalar, Dimension> EQLAImageType;
	
	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

	const char* ObjectType = m_ParamObject->GetDeformableObjectType().c_str();

	if (itksys::SystemTools::Strucmp(ObjectType,"Landmark") == 0)
	{
		LandmarkType* objectLandmark = new LandmarkType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectLandmark->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectLandmark->SetPolyData(PolyData);
		}
		m_Object = objectLandmark;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"PointCloud") == 0)
	{
		PointCloudType* objectPointCloud = new PointCloudType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectPointCloud->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectPointCloud->SetPolyData(PolyData);

			objectPointCloud->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectPointCloud->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectPointCloud;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"OrientedPolyLine") == 0)
	{
		OrientedPolyLineType* objectOrientedPolyLine = new OrientedPolyLineType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectOrientedPolyLine->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectOrientedPolyLine->SetPolyData(PolyData);

			objectOrientedPolyLine->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectOrientedPolyLine;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"NonOrientedPolyLine") == 0)
	{
		NonOrientedPolyLineType* objectNonOrientedPolyLine = new NonOrientedPolyLineType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectNonOrientedPolyLine->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectNonOrientedPolyLine->SetPolyData(PolyData);

			objectNonOrientedPolyLine->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectNonOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectNonOrientedPolyLine;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"OrientedSurfaceMesh") == 0)
	{
		bool ReOrientSurface = m_ParamObject->ReOrient();

		OrientedSurfaceMeshType* objectOrientedSurfaceMesh = new OrientedSurfaceMeshType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectOrientedSurfaceMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectOrientedSurfaceMesh->SetPolyData(PolyData);
			(ReOrientSurface)?objectOrientedSurfaceMesh->SetReorient():objectOrientedSurfaceMesh->UnSetReorient();

			objectOrientedSurfaceMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectOrientedSurfaceMesh;
	}

	else if ( m_ParamObject->GetDeformableObjectType() == "NonOrientedSurfaceMesh" )
	{
		NonOrientedSurfaceMeshType* objectNonOrientedSurfaceMesh = new NonOrientedSurfaceMeshType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectNonOrientedSurfaceMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectNonOrientedSurfaceMesh->SetPolyData(PolyData);

			objectNonOrientedSurfaceMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectNonOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());			
		}
		m_Object = objectNonOrientedSurfaceMesh;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"CMRep") == 0)
	{
		CMRepType* objectCMRep= new CMRepType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectCMRep->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
			objectCMRep->SetPolyData(PolyData);

			objectCMRep->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
			objectCMRep->SetKernelWidth(m_ParamObject->GetKernelWidth());			
		}
		m_Object = objectCMRep;
	}

	else if ( (itksys::SystemTools::Strucmp(ObjectType,"SSDImage") == 0) || (itksys::SystemTools::Strucmp(ObjectType,"LCCImage") == 0) || (itksys::SystemTools::Strucmp(ObjectType,"EQLAImage") == 0) )
	{
		typedef itk::ImageFileReader<ImageType> ReaderType;
		ImageTypePointer ExampleImg;

		LinearInterpImageType* objectImg;
		if (itksys::SystemTools::Strucmp(ObjectType,"SSDImage") == 0)
		{
			SSDImageType* objectImgAux = new SSDImageType();
			objectImg = objectImgAux;
		}
		else if (itksys::SystemTools::Strucmp(ObjectType,"EQLAImage") == 0)	
		{
			EQLAImageType* objectImgAux = new EQLAImageType();
			objectImgAux->SetEQLAKernelWidth( m_ParamObject->GetKernelWidth() );
			objectImg = objectImgAux;
		}
		else
		{
			LCCImageType* objectImgAux = new LCCImageType();
			objectImgAux->SetLCCKernelWidth( m_ParamObject->GetKernelWidth() );
			objectImg = objectImgAux;
		}
/*
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(m_FileName);

		typedef itk::NiftiImageIO NiftiIOType;
		typename NiftiIOType::Pointer NiftiIOPointer = NiftiIOType::New();

		ImageTypePointer img;
		if ( NiftiIOPointer->CanReadFile(m_FileName) )
		{
			NiftiIOPointer->SetFileName(m_FileName);
			reader->SetImageIO(NiftiIOPointer);
			reader->Update();
			img = reader->GetOutput();
			// std::cout << "Read " << m_FileName << " as Nifti image" << std::endl;
		}
		else
		{
			reader->Update();
			img = reader->GetOutput();
		}

		typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
		typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		minMaxFilter->SetInput(img);
		minMaxFilter->Update();
		double minI = minMaxFilter->GetMinimum();
		double maxI = minMaxFilter->GetMaximum();
		
		typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
		typename RescalerType::Pointer resf = RescalerType::New();
		resf->SetInput(img);
		resf->SetOutputMinimum(0.0);
		resf->SetOutputMaximum(1.0);
		resf->Update();
		ImageTypePointer rescaledImage = resf->GetOutput();
		
		TScalar ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();
		
		objectImg->SetImageAndDownSamplingFactor(rescaledImage, ImageGridDownsampling);
		objectImg->SetMinIntensityOutput(minI);
		objectImg->SetMaxIntensityOutput(maxI);

		typename ImageType::PointType minAux = img->GetOrigin();
		typename ImageType::PointType maxAux = minAux;
		typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
		typename ImageType::SpacingType spacing = img->GetSpacing();

		std::cout << "Image read: origin = " << img->GetOrigin() << " size = " << size << " spacing = " << spacing << " min value = " << minI << " max value = " << maxI << ", rescaled between 0 and 1" << std::endl;
		std::cout << "orientation =\n" << img->GetDirection() << std::endl;
		// minMaxFilter->SetInput(objectImg->GetImage());
		// minI = minMaxFilter->GetMinimum();
		// maxI = minMaxFilter->GetMaximum();
		// std::cout << "after renormalization " << minI << " " << maxI << std::endl;

		itk::Matrix<double,Dimension,Dimension> dir = img->GetDirection();
		std::vector<int> permutation(Dimension);
		std::vector<int> flip_axes(Dimension);
		for (unsigned int d=0; d < Dimension; d++)
		{
			int ind = 0;
			double val = fabs(dir(0,d));
			for (unsigned int k=1; k < Dimension; k++)
			{
				if (fabs(dir(k,d))>val)
				{
					ind = k;
					val = fabs(dir(k,d));
				}
			}
			permutation[d] = ind;
			flip_axes[d] = (dir(ind,d)>0.0)?1:-1;
		}
		objectImg->SetPermutationAxes(permutation);
		objectImg->SetFlipAxes(flip_axes);

		// std::cout << "permutation\n" << permutation[0] << permutation[1] << permutation[2] << std::endl;
		// std::cout << "flip axes\n" << flip_axes[0] << flip_axes[1] << flip_axes[2] << std::endl;

		// MatrixType Yex = GridFunctionsType::ImageToPoints(img);
		// for (unsigned int d=0; d<Dimension; d++)
		// {
		// 	m_Max[d] = Yex.get_column(d).max_value();
		// 	m_Min[d] = Yex.get_column(d).min_value();
		// 	// wrong if axes are flipped
		// 	// m_Max[d] = maxAux[d] + size[d]*spacing[d];
		// 	// m_Min[d] = minAux[d];
		// }

		// if (m_IsTemplate)
		// {
		// TScalar ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();
		// ImageTypePointer ExampleDownsampledImage;
		// if (ImageGridDownsampling <= 1.0)
		// 	ExampleDownsampledImage = img;
		// else
		// 	ExampleDownsampledImage = GridFunctionsType::DownsampleImage(img, ImageGridDownsampling);
		// 
		// 	// MatrixType DownsampledY1 = GridFunctionsType::ImageToPoints(ExampleDownsampledImage);
		// 
		// 	// std::cout << "DownsampledImage = origin " <<  ExampleDownsampledImage->GetOrigin() << " size = " << ExampleDownsampledImage->GetLargestPossibleRegion().GetSize() << " spacing = " << ExampleDownsampledImage->GetSpacing() << std::endl;
		// objectImg->SetDownSampledReferenceImage(ExampleDownsampledImage);
		// // }


		m_Object = objectImg;*/

		{

			typename ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(m_FileName);

			typedef itk::NiftiImageIO NiftiIOType;
			typename NiftiIOType::Pointer NiftiIOPointer = NiftiIOType::New();

			ImageTypePointer img, inputimg;
			// img = working image, which can be the same as the original image, if no orientation issues.
			if ( NiftiIOPointer->CanReadFile(m_FileName) )
			{
				//ImageTypePointer inputimg;
				NiftiIOPointer->SetFileName(m_FileName);
				reader->SetImageIO(NiftiIOPointer);
				reader->Update();
				inputimg = reader->GetOutput();

				// Reorient image if necessary, saving original coordinate system info in the final corresponding deformable object
				itk::Matrix<double,Dimension,Dimension> dir = inputimg->GetDirection();
				itk::Matrix<double,Dimension,Dimension> id;
				id.SetIdentity();

				// If not LPS orientation
				if (dir != id)
				{
					if (Dimension != 3)
					{
						std::cout << "Warning: cannot reorient images that are not 3D" << std::endl;
					}
					else
					{
						MatrixType orient_mtx(Dimension, Dimension, 0.0);
						for (unsigned int i = 0; i < Dimension; i++)
							for (unsigned int j = 0; j < Dimension; j++)
								orient_mtx(i, j) = dir(i,j);

						//std::cout << "Orientation matrix " << orient_mtx << std::endl;
						objectImg->SetAnatomicalCoordinateSystem(orient_mtx);
						std::cout << "Reorienting input image from " <<  objectImg->GetAnatomicalCoordinateSystemLabel() << " orientation to the internal LPS standard." << std::endl;

						//std::cout << "Deformable Object is an image originally oriented in the " << objectImg->GetAnatomicalCoordinateSystemLabel() << " coordinate system." << std::endl;
						//std::cout << "Reorienting image to match the LPS standard." << std::endl;

						// reorient image (LPS orientation / RAI code for ITK filter)
						typedef itk::Image<TScalar, 3> ImageType_dim_3;
						typename itk::OrientImageFilter<ImageType_dim_3,ImageType_dim_3>::Pointer orienter = itk::OrientImageFilter<ImageType_dim_3,ImageType_dim_3>::New();
						orienter->UseImageDirectionOn();
						orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
						orienter->SetInput( dynamic_cast<ImageType_dim_3*>( inputimg.GetPointer() ));
						orienter->Update();
						img = dynamic_cast<ImageType*>( orienter->GetOutput() );

						/* typedef itk::ImageFileWriter<ImageType> WriterType;
					    			typename WriterType::Pointer writer = WriterType::New();
					    			writer->SetInput(img); //
					    			writer->SetFileName(string(m_FileName).append("LPS.nii"));
					    			writer->Update(); */

						//std::cout << "Reoriented image matrix (LPS): " << img << std::endl; // << " " << (*outimg)->GetDirection() << std::endl;
					}
				}
				//std::cout << "Read " << m_FileName << " as Nifti image" << std::endl;
				else
					img = inputimg;
			}
			else
			{
				reader->Update();
				img = reader->GetOutput();
			}

			//std::cout << "Working image = " << img << std::endl;

			typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
			typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
			minMaxFilter->SetInput(img);
			minMaxFilter->Update();
			double minI = minMaxFilter->GetMinimum();
			double maxI = minMaxFilter->GetMaximum();


			typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
			typename RescalerType::Pointer resf = RescalerType::New();
			resf->SetInput(img);
			resf->SetOutputMinimum(0.0);
			resf->SetOutputMaximum(1.0);
			resf->Update();
			ImageTypePointer rescaledImage = resf->GetOutput();

			TScalar ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();

			objectImg->SetImageAndDownSamplingFactor(rescaledImage, ImageGridDownsampling);
			objectImg->SetMinIntensityOutput(minI);
			objectImg->SetMaxIntensityOutput(maxI);

//			typename ImageType::PointType minAux = img->GetOrigin();
//			typename ImageType::PointType maxAux = minAux;
			typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
			typename ImageType::SpacingType spacing = img->GetSpacing();

			std::cout << "Image read: origin = " << img->GetOrigin() << " size = " << size << " spacing = " << spacing << " min value = " << minI << " max value = " << maxI << ", rescaled between 0 and 1" << std::endl;
			
			// At this stage, img->GetDirection should be the identity matrix for a 3D Nifti image. The following could be useful for other types of images.
			itk::Matrix<double,Dimension,Dimension> dir = img->GetDirection();
			std::vector<int> permutation(Dimension);
			std::vector<int> flip_axes(Dimension);
			for (unsigned int d=0; d < Dimension; d++)
			{
				int ind = 0;
				double val = fabs(dir(0,d));
				for (unsigned int k=1; k < Dimension; k++)
				{
					if (fabs(dir(k,d))>val)
					{
						ind = k;
						val = fabs(dir(k,d));
					}
				}
				permutation[d] = ind;
				flip_axes[d] = (dir(ind,d)>0.0)?1:-1;
			}
			objectImg->SetPermutationAxes(permutation);
			objectImg->SetFlipAxes(flip_axes);

			// std::cout << "permutation\n" << permutation[0] << permutation[1] << permutation[2] << std::endl;
			// std::cout << "flip axes\n" << flip_axes[0] << flip_axes[1] << flip_axes[2] << std::endl;

/*			MatrixType Yex = GridFunctionsType::ImageToPoints(img);
			for (unsigned int d=0; d<Dimension; d++)
			{
				m_Max[d] = Yex.get_column(d).max_value();
				m_Min[d] = Yex.get_column(d).min_value();
				// wrong if axes are flipped
				// m_Max[d] = maxAux[d] + size[d]*spacing[d];
				// m_Min[d] = minAux[d];
			}
			*/

		}
		m_Object = objectImg;
	}
	
	else
	{
		std::cerr << "Unknown object type: " << m_ParamObject->GetDeformableObjectType() << 
				"\nCurrently available types are: SSDImage, LCCImage, EQLAImage, Landmark, PointCloud, OrientedPolyLine, "
				"NonOrientedPolyLine, OrientedSurfaceMesh, NonOrientedSurfaceMesh" << std::endl;
		return;
	}
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
KernelEnumType
DeformableObjectReader<TScalar, Dimension>
::StringToKernelEnumType(const char* kernelType)
 {
	KernelEnumType result = null;
	if (itksys::SystemTools::Strucmp(kernelType, "fgt") == 0)
		result = FGT;
	else if (itksys::SystemTools::Strucmp(kernelType, "p3m") == 0)
		result = P3M;
#ifdef USE_CUDA
	else if (itksys::SystemTools::Strucmp(kernelType, "cudaexact") == 0)
		result = CUDAExact;
#endif
	else
	{
		if (itksys::SystemTools::Strucmp(kernelType, "exact") != 0)
			std::cerr << "Unknown kernel type for the deformable object : defaulting to exact" << std::endl;
		result = Exact;
	}
	return result;
 }

#endif

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformationFieldIO_txx
#define _DeformationFieldIO_txx


#include "GridFunctions.h"

//#include "Diffeos.h"
#include "DeformationFieldIO.h"

//#include "writeMatrixDLM.txx"

#include <iostream>
#include <iomanip>
#include <sstream>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformationFieldIO<TScalar, Dimension>
::DeformationFieldIO()
{
    m_Factor = 5;
	m_Def = NULL;
	m_CoordSystem.SetAnatomicalCoordinateSystemLabel("LPS");
}



template <class TScalar, unsigned int Dimension>
DeformationFieldIO<TScalar, Dimension>
::~DeformationFieldIO()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void DeformationFieldIO<TScalar, Dimension>::Update()
{
	this->WriteDeformationField("output", true);
}

template <class TScalar, unsigned int Dimension>
void DeformationFieldIO<TScalar, Dimension>::WriteDeformationField(const std::string& fn, bool finalField)
{
    if ( m_Def != NULL )
	{
    	std::string prefix(fn);

   		prefix.append("_defFieldLPS");

   		WriteDeformationField3DImage(prefix, finalField);
/*
		switch (Dimension)
		{
		case 2:
			WriteDeformationField2DImage(prefix, finalField);
		case 3:
			WriteDeformationField3DImage(prefix, finalField);
		default:
			cout << "DeformationFieldIO::WriteDeformationField: unsupported field dimension. No output produced." << endl;
								return;
		}*/
	}
	else
	{
		std::cerr << "DeformationFieldIO::control diffeomorphism (abstract deformation) object should have been set before calling WriteDeformationField!" << std::endl;
	}


}




template <class TScalar, unsigned int Dimension>
bool DeformationFieldIO<TScalar, Dimension>::WriteDeformationField3DImage(const std::string& fn, bool finalField)
{
	MatrixType DD = m_Def->GetDataDomain();

	TScalar img_spacing = m_Def->GetKernelWidth();

	// Create image
	DisplacementField3DTypePointer fieldImg = DisplacementField3DType::New();

	// Image properties
	typename DisplacementField3DType::RegionType region;
	typename DisplacementField3DType::IndexType startIdx;
	typename DisplacementField3DType::SizeType size;
	typename DisplacementField3DType::SpacingType spacing;
	typename DisplacementField3DType::PointType origin;
	typename DisplacementField3DType::DirectionType direction;

	VectorType minCoord = DD.get_column(0);
	VectorType maxCoord = DD.get_column(1);

	// Compute information for building the deformation field image
	direction.Fill(0);

	for (unsigned int i = 0; i < Dimension; i++)
	{

		startIdx[i] = 0;
		spacing[i] = img_spacing / m_Factor; // en fonction de la taille du noyau des deformations
		size[i] = ceil(((maxCoord[i] - minCoord[i]) + 1) / spacing[i]);
		origin[i] = minCoord[i];
		direction(i,i) = 1;

	}

	// complete info about z axis if dealing with 2D image
	if (Dimension == 2)
	{
		startIdx[Dimension] = 0;
		spacing[Dimension] = img_spacing / m_Factor; // en fonction de la taille du noyau des deformations
		size[Dimension] = 1;
		origin[Dimension] = 0;
		direction(Dimension,Dimension) = 1;
	}

	region.SetIndex(startIdx);
	region.SetSize(size);

	fieldImg->SetRegions(region);
	fieldImg->SetSpacing(spacing);
	fieldImg->SetOrigin(origin);
	fieldImg->SetDirection(direction);

	fieldImg->Allocate();

	// Fill image
	Vector3DType nullvector;
	nullvector.Fill(0.0);

	fieldImg->FillBuffer(nullvector);

	// Convert image points (coordinates) into vtkPoints
	long npoints = fieldImg->GetLargestPossibleRegion().GetNumberOfPixels();

	vtkSmartPointer<vtkPoints> fieldPoints = vtkSmartPointer<vtkPoints>::New();

	long r = 0;

	typedef itk::ImageRegionConstIteratorWithIndex<DisplacementField3DType> IteratorType;
	IteratorType it(fieldImg, fieldImg->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		typename DisplacementField3DType::IndexType ind = it.GetIndex();

		typename DisplacementField3DType::PointType p;
		fieldImg->TransformIndexToPhysicalPoint(ind, p);

		fieldPoints->InsertNextPoint(p[0], p[1], p[2]);

		r++;
	}

	// save it to vtkPolyData container required by the Landmark object
	vtkSmartPointer<vtkPolyData> fieldPointData = vtkSmartPointer<vtkPolyData>::New();
	fieldPointData->SetPoints(fieldPoints);

	// create a Landmark object and set points/polyData accordingly
	Landmark<TScalar, Dimension> ldmkObj;
	ldmkObj.SetAnatomicalCoordinateSystem("LPS");
	ldmkObj.SetPolyData(fieldPointData);


	// set the Landmark object to the corresponding Diffeo object
	DeformableObjectList list;
	list.push_back(&ldmkObj);

	DeformableMultiObjectType* DMO = new DeformableMultiObjectType();
	DMO->SetObjectList(list);
    DMO->Update();

	m_Def->SetDeformableMultiObject(DMO);

	// update points according to the computed diffeomorphism
	m_Def->Update();


	// get object flow and compute displacement field
	DeformableMultiObjectType *UDMO, *UDMO2;
	DeformableObjectList ulist, ulist2;

//	int NumbTimePoints = m_Def->GetNumberOfTimePoints();
	MatrixList deformedFieldPoints(2); // for two consecutive timepoints

	Point3DType pt0, ptN, or_pt0, or_ptN;
	Vector3DType displacement, orient_displacement;
	typename DisplacementField3DType::IndexType idx;

	int T0 = m_Def->GetT0();
	int TN = m_Def->GetNumberOfTimePoints();
	int step;

	MatrixType baseMtx(3,3);

	if (Dimension == 2)
	{
		MatrixType baseMtx2D = m_CoordSystem.GetInverseChangeOfBasisMatrix();
		for (int i = 0; i < Dimension; i++)
		{
			for (int j = 0; j < Dimension; j++)
			{
				baseMtx(i,j) = baseMtx2D(i,j);
			}
		}
		baseMtx(Dimension,Dimension) = 1;
	}
	else
		baseMtx = m_CoordSystem.GetInverseChangeOfBasisMatrix();

	typedef itk::ImageFileWriter< DisplacementField3DType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	typename itk::OrientImageFilter<DisplacementField3DType,DisplacementField3DType>::Pointer orienter = itk::OrientImageFilter<DisplacementField3DType,DisplacementField3DType>::New();
	orienter->UseImageDirectionOn();
	orienter->SetDesiredCoordinateOrientation(m_CoordSystem.GetITKAnatomicalCoordinateSystemLabel());

	std::stringstream str; //output string stream
	std::string fnMetaData;

	if ((finalField == true) && (TN > 1))
		step = TN - 1;
	else
	{
		step = 1;

		// save deformation field at time 0 (null field)
		fnMetaData.assign(fn + "_t0.mha");
		writer->SetFileName(fnMetaData);
		orienter->SetInput(fieldImg);
		orienter->Update();
		writer->SetInput(orienter->GetOutput());

		// Write MetaImage
		try
		{
		writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}

	for (int t = T0 + step; t < TN; t += step)
	{
		// compute deformation fields between timestamps
		// get data from T0
		UDMO = m_Def->GetDeformedObjectAt(t - step);
		ulist = UDMO->GetObjectList();

		// for each timestep, get the deformed objects list, traverse it, get updated point coordinates into deformedFieldPoints
		UDMO2 = m_Def->GetDeformedObjectAt(t);
		ulist2 = UDMO2->GetObjectList();

		for (DeformableObjectListIterator it = ulist.begin() ; it != ulist.end(); ++it)
		{
			Landmark<TScalar, Dimension> *uldmkObj;
			uldmkObj = dynamic_cast<Landmark<TScalar, Dimension> *> (*it);
			deformedFieldPoints[0] = uldmkObj->GetPointCoordinates();
		}

		for (DeformableObjectListIterator it = ulist2.begin() ; it != ulist2.end(); ++it)
		{
			Landmark<TScalar, Dimension> *uldmkObj2;
			uldmkObj2 = dynamic_cast<Landmark<TScalar, Dimension> *> (*it);
			deformedFieldPoints[1] = uldmkObj2->GetPointCoordinates();
		}

		// compute displacements from flow and save them in a proper vector image format
		for (long i = 0; i < npoints; i++)
		{

			for (unsigned int j = 0; j < Dimension; j++)
			{
				pt0[j] = (deformedFieldPoints[0])(i,j);
				ptN[j] = (deformedFieldPoints[1])(i,j);
				displacement[j] = ptN[j] - pt0[j];
			}

			if (Dimension == 2)
			{
				displacement[Dimension] = 0;
			}

			bool isInside = fieldImg->TransformPhysicalPointToIndex(pt0, idx);

			for (int m = 0; m < 3; m++)
			{
				orient_displacement[m] = (baseMtx(m, 0) * displacement[0]) + (baseMtx(m, 1) * displacement[1]) + (baseMtx(m, 2) * displacement[2]);
				//displacement[m] = (baseMtx(m, 0) * displacement[0]) + (baseMtx(m, 1) * displacement[1]) + (baseMtx(m, 2) * displacement[2]);
				//std::cout << displacement[m] << " " << orient_displacement[m] << std::endl;
			}

			if (isInside)
			{
				fieldImg->SetPixel(idx, orient_displacement);
			}
		}

		// save deformation field at time t (null field)
		str << t;

		fnMetaData.assign(fn + "_t" + str.str() + ".mha");
		writer->SetFileName(fnMetaData);
		orienter->SetInput(fieldImg);
		orienter->Update();
		writer->SetInput(orienter->GetOutput());

		// Write MetaImage
		try
		{
		writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}


	return true;
}



#endif /* _DeformationFieldIO_txx */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AnatomicalCoordinateSystem_txx
#define _AnatomicalCoordinateSystem_txx

#include "AnatomicalCoordinateSystem.h"

#include <iostream>
#include <cmath>

#include "itksys/SystemTools.hxx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
AnatomicalCoordinateSystem<TScalar, Dimension>
::AnatomicalCoordinateSystem()
 {
	m_ACSLabel = "LPS";
	m_ChangeOfBasisMatrix.set_size(Dimension, Dimension);
	m_ChangeOfBasisMatrix.set_identity();
 }



template <class TScalar, unsigned int Dimension>
AnatomicalCoordinateSystem<TScalar, Dimension>
::AnatomicalCoordinateSystem(std::string label)
 {
	this->SetAnatomicalCoordinateSystemLabel(label);
 }



template <class TScalar, unsigned int Dimension>
AnatomicalCoordinateSystem<TScalar, Dimension>
::AnatomicalCoordinateSystem(MatrixType changeOfBasisMatrix)
 {
	this->SetChangeOfBasisMatrix(changeOfBasisMatrix);
 }



template <class TScalar, unsigned int Dimension>
AnatomicalCoordinateSystem<TScalar, Dimension>
::~AnatomicalCoordinateSystem()
{}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
bool
AnatomicalCoordinateSystem<TScalar, Dimension>
::SetAnatomicalCoordinateSystemLabel(const std::string label)
 {
	bool axial_dir = false, sag_dir = false, cor_dir = false;

	/*if(label.size() != Dimension) {
		std::cerr << "The dimension and the number of characters of the label do not match." << std::endl;
	}*/
	m_ACSLabel = itksys::SystemTools::UpperCase(label);

	// Check if label contains valid characters : only one from R/L; only one from A/P; only one from I/S
	for (int i = 0; i < label.size(); i++)
	{
		switch (m_ACSLabel.at(i))
		{
		case 'R':
		case 'L':
			if (sag_dir == true)
			{
				std::cerr << "Invalid sequence of label characters." << std::endl;
				return false;
			}
			else
				sag_dir = true;
			break;
		case 'A':
		case 'P':
			if (cor_dir == true)
			{
				std::cerr << "Invalid sequence of label characters." << std::endl;
				return false;
			}
			else
				cor_dir = true;
			break;
		case 'S':
		case 'I':
			if (axial_dir == true)
			{
				std::cerr << "Invalid sequence of label characters." << std::endl;
				return false;
			}
			else
				axial_dir = true;
			break;
		default: // unknown orientation
			std::cerr << "Invalid label character." << std::endl;
			return false;
		}
	}

	if (!(sag_dir && cor_dir && axial_dir))
	{
		std::cerr << "Invalid sequence of label characters." << std::endl;
		return false;
	}

	m_ChangeOfBasisMatrix = this->LabelToChangeOfBasisMatrix(m_ACSLabel);

	return true;
 }



template <class TScalar, unsigned int Dimension>
void
AnatomicalCoordinateSystem<TScalar, Dimension>
::SetChangeOfBasisMatrix(const MatrixType& changeOfBasisMatrix)
 {
	if( (changeOfBasisMatrix.rows() != Dimension) || (changeOfBasisMatrix.columns() != Dimension) ) {
		throw std::runtime_error("In AnatomicalCoordinateSystem : The dimensions of the change of basis matrix are not correct");
	}
	m_ChangeOfBasisMatrix = changeOfBasisMatrix;
	m_ACSLabel = this->ChangeOfBasisMatrixToLabel(m_ChangeOfBasisMatrix);
 }



template <class TScalar, unsigned int Dimension>
itk::SpatialOrientation::ValidCoordinateOrientationFlags AnatomicalCoordinateSystem<TScalar, Dimension>
::GetITKAnatomicalCoordinateSystemLabel() const
 {
	std::string codestr(m_ACSLabel);
	std::size_t codestr_sz = codestr.length();

//	itk::SpatialOrientation::ValidCoordinateOrientationFlags label;

	// invert code characters for consistency with itk spatial orientation flags
	for (int i = 0; i < codestr_sz; i++)
	{
		switch (codestr.at(i))
		{
		case 'R':
			codestr.at(i) = 'L';
			break;
		case 'L':
			codestr.at(i) = 'R';
			break;
		case 'A':
			codestr.at(i) = 'P';
			break;
		case 'P':
			codestr.at(i) = 'A';
			break;
		case 'S':
			codestr.at(i) = 'I';
			break;
		case 'I':
			codestr.at(i) = 'S';
			break;
		default: // unknown orientation
			return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;
		}
	}

	if ( codestr.compare("RAI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
	}
	else if ( codestr.compare("IRA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA;
	}
	else if ( codestr.compare("AIR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR;
	}
	else if ( codestr.compare("ARI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI;
	}
	else if ( codestr.compare("IAR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR;
	}
	else if ( codestr.compare("RIA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA;
	}
	else if ( codestr.compare("RAS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
	}
	else if ( codestr.compare("SRA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA;
	}
	else if ( codestr.compare("ASR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR;
	}
	else if ( codestr.compare("ARS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS;
	}
	else if ( codestr.compare("SAR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR;
	}
	else if ( codestr.compare("RSA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
	}
	else if ( codestr.compare("RPI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI;
	}
	else if ( codestr.compare("IRP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP;
	}
	else if ( codestr.compare("PIR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR;
	}
	else if ( codestr.compare("PRI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI;
	}
	else if ( codestr.compare("IPR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR;
	}
	else if ( codestr.compare("RIP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP;
	}
	else if ( codestr.compare("RPS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS;
	}
	else if ( codestr.compare("SRP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP;
	}
	else if ( codestr.compare("PSR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR;
	}
	else if ( codestr.compare("PRS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS;
	}
	else if ( codestr.compare("SPR") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR;
	}
	else if ( codestr.compare("RSP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP;
	}
	else if ( codestr.compare("LAI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI;
	}
	else if ( codestr.compare("ILA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA;
	}
	else if ( codestr.compare("AIL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL;
	}
	else if ( codestr.compare("ALI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI;
	}
	else if ( codestr.compare("IAL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL;
	}
	else if ( codestr.compare("LIA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA;
	}
	else if ( codestr.compare("LAS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS;
	}
	else if ( codestr.compare("SLA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA;
	}
	else if ( codestr.compare("ASL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;
	}
	else if ( codestr.compare("ALS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS;
	}
	else if ( codestr.compare("SAL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL;
	}
	else if ( codestr.compare("LSA") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA;
	}
	else if ( codestr.compare("LPI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
	}
	else if ( codestr.compare("ILP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP;
	}
	else if ( codestr.compare("PIL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL;
	}
	else if ( codestr.compare("PLI") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI;
	}
	else if ( codestr.compare("IPL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL;
	}
	else if ( codestr.compare("LIP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP;
	}
	else if ( codestr.compare("LPS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
	}
	else if ( codestr.compare("SLP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP;
	}
	else if ( codestr.compare("PSL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL;
	}
	else if ( codestr.compare("PLS") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS;
	}
	else if ( codestr.compare("SPL") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL;
	}
	else if ( codestr.compare("LSP") == 0 )
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP;
	}
	else
	{
		return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;
	}
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
AnatomicalCoordinateSystem<TScalar, Dimension>
::PrintSelf() const
 {
	std::cout << "The label is " << this->GetAnatomicalCoordinateSystemLabel();
	std::cout << " and the associated matrix is :" << std::endl << this ->GetChangeOfBasisMatrix();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
std::string
AnatomicalCoordinateSystem<TScalar, Dimension>
::ChangeOfBasisMatrixToLabel(MatrixType changeOfBasisMatrix)
 {
	std::string result = "___";
	int index;
	double value;
	char letter;

	for (int j = 0; j < Dimension; j++)
	{
		index = 0;
		value = fabs(changeOfBasisMatrix(0,j));
		for (int i=1; i < Dimension; i++)
		{
			if (fabs(changeOfBasisMatrix(i,j)) > value)
			{
				index = i;
				value = fabs(changeOfBasisMatrix(i,j));
			}
		}
		switch (index)
		{
		case 0 : letter = (changeOfBasisMatrix(index,j)>0.0)?'L':'R';
		break;
		case 1 : letter = (changeOfBasisMatrix(index,j)>0.0)?'P':'A';
		break;
		case 2 : letter = (changeOfBasisMatrix(index,j)>0.0)?'S':'I';
		break;
		default :
			throw std::runtime_error("In AnatomicalCoordinateSystem : Invalid change of basis matrix");
		}
		result[j] = letter;
	}

	return result;
 }



template <class TScalar, unsigned int Dimension>
typename LinearAlgebra<TScalar>::Matrix
AnatomicalCoordinateSystem<TScalar, Dimension>
::LabelToChangeOfBasisMatrix(std::string label)
 {
	MatrixType result(Dimension, Dimension, 0.0);
	int index, value;

	for(int j = 0; j < Dimension; j++) {
		switch (label[j])
		{
		case 'L' : index = 0; value =  1;
		break;
		case 'P' : index = 1; value =  1;
		break;
		case 'S' : index = 2; value =  1;
		break;
		case 'R' : index = 0; value = -1;
		break;
		case 'A' : index = 1; value = -1;
		break;
		case 'I' : index = 2; value = -1;
		break;
		default :
			throw std::runtime_error("In AnatomicalCoordinateSystem::::LabelToChangeOfBasisMatrix(label) : Invalid label");
		}
		result(index, j) = value;
	}
	return result;
 }



#endif /* _AnatomicalCoordinateSystem_txx */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AnatomicalCoordinateSystem_h
#define _AnatomicalCoordinateSystem_h

#include "LinearAlgebra.h"

#include "itkSpatialOrientation.h"

#include <string>



/**
 *  \brief      Anatomical Coordinate System.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The AnatomicalCoordinateSystem class contains all the operations to convert the anatomical coordinate system
 *              of a given object to our frame of reference (which is LPS in Deformetrica) and vice versa.
 *              When an instance of AnatomicalCoordinateSystem is created on default, we assume that the anatomical coordinate system
 *              and the frame of reference are the same (i.e. LPS).
 */
template <class TScalar, unsigned int Dimension>
class AnatomicalCoordinateSystem {

public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;

    // /// ITK Coordinate System enumeration type.
	// typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags ITKCoordinateLabel;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	AnatomicalCoordinateSystem();
	/// Constructor which initializes the anatomical coordinate system with the label \e label.
	AnatomicalCoordinateSystem(std::string label);
	/// Constructor which initializes the anatomical coordinate system with the matrix \e changeOfBasisMatrix.
	AnatomicalCoordinateSystem(MatrixType changeOfBasisMatrix);

	~AnatomicalCoordinateSystem();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the label of the coordinates system.
	inline std::string GetAnatomicalCoordinateSystemLabel() const { return m_ACSLabel; }
	/// Sets the label of the coordinates system to \e label.
	bool SetAnatomicalCoordinateSystemLabel(const std::string label);

	/// Returns the (ITK) label of the coordinates system.
	itk::SpatialOrientation::ValidCoordinateOrientationFlags GetITKAnatomicalCoordinateSystemLabel() const;
	// void SetITKAnatomicalCoordinateSystemLabel(const std::string label);

	/// Returns the change-of-basis matrix (towards the frame of reference).
	inline MatrixType GetChangeOfBasisMatrix() const { return m_ChangeOfBasisMatrix; }
	/// Sets the change-of-basis matrix (towards the frame of reference) to \e changeOfBasisMatrix.
	void SetChangeOfBasisMatrix(const MatrixType& changeOfBasisMatrix);


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the change-of-basis matrix (towards the anatomical coordinate system).
	MatrixType GetInverseChangeOfBasisMatrix() const { return m_ChangeOfBasisMatrix.transpose(); }

	/// Returns the change-of-basis matrix (towards the anatomical coordinate system).
	void PrintSelf() const;



protected :

	///	Matrix containing the change-of-basis from our anatomical coordinate system to our frame of reference.
	MatrixType m_ChangeOfBasisMatrix;

	///	Label containing the anatomical coordinate system (e.g. LPS or RAS).
	std::string m_ACSLabel;

	///	(ITK) Label containing the anatomical coordinate system (e.g. LPS or RAS).
	itk::SpatialOrientation::ValidCoordinateOrientationFlags m_ItkCoordinateLabel;


private :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Converts a label to a change-of-basis matrix.
	MatrixType LabelToChangeOfBasisMatrix(std::string label);

	/// Converts a change-of-basis matrix to label.
	std::string ChangeOfBasisMatrixToLabel(MatrixType changeOfBasisMatrix);


}; /* class AnatomicalCoordinateSystem */


#ifndef MU_MANUAL_INSTANTIATION
#include "AnatomicalCoordinateSystem.txx"
#endif


#endif /* _AnatomicalCoordinateSystem_h */

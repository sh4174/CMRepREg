/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectReader_h
#define _DeformableObjectReader_h

#include "LinearAlgebra.h"

#include "KernelType.h"

#include "DeformableObject.h"
#include "DeformableObjectParametersXMLFile.h"

#include <cstring>
#include <iostream>
#include <sstream>

template <class TScalar, unsigned int Dimension>
class DeformableObjectReader
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type.
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;

	/// ITK Image type
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK Image Pointer type
	typedef typename ImageType::Pointer ImageTypePointer;

	/// Deformable object type.
	typedef DeformableObject<TScalar,Dimension> DeformableObjectType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectReader();

	~DeformableObjectReader();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the file name to \e fn.
	void SetFileName(char* fn) { m_FileName = fn; }

	void SetObjectParameters(DeformableObjectParameters::Pointer param) { m_ParamObject = param; }

	void SetTemplateType() { m_IsTemplate = true; }

	/// Returns the ouput i.e. the deformable object.
	DeformableObjectType* GetOutput() { return m_Object; }

	// /// Returns the bounding box.
	// MatrixType GetBoundingBox() { return m_Bbox; }

	// /// Returns data-sigma squared.
	// double GetDataSigmaSquared() const { return (m_ParamObject->GetDataSigma() * m_ParamObject->GetDataSigma()); }
	// /************************************** CHANGE BEGIN ******************************************/
	// double GetWeightData() const { return m_ParamObject->GetWeightData(); }
	// double GetPriorData() const { return m_ParamObject->GetPriorData(); }
	// double GetKernelWidth() const { return m_ParamObject->GetKernelWidth(); }
	// /************************************** CHANGE END ******************************************/


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void Update();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	KernelEnumType StringToKernelEnumType(const char* kernelType);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Name of the file where the deformable object will be read.
	char* m_FileName;	

	DeformableObjectParameters::Pointer m_ParamObject;

	/// Our deformable object which will be extracted from the file.
	DeformableObjectType* m_Object;


	// MatrixType m_Bbox;

	bool m_IsTemplate;

	// VectorType m_Min;
	// VectorType m_Max;


};

#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableObjectReader.txx"
#endif


#endif

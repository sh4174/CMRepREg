/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AbstractDeformations_txx
#define _AbstractDeformations_txx

#include "KernelFactory.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::AbstractDeformations() :
m_Type(null), m_DeformableMultiObject(NULL), m_Modified(true), m_DeformableObjectModified(true)
 {
 }



template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::~AbstractDeformations()
{}



template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::AbstractDeformations(const AbstractDeformations& other)
 {
	m_Type = other.m_Type;
	m_Modified = true;
	m_DeformableObjectModified = true;

	if (other.m_DeformableMultiObject != NULL)
    {
		m_DeformableMultiObject = other.m_DeformableMultiObject->Clone();
    }
    else
    {
        m_DeformableMultiObject = NULL;
    }
         
	
	m_LandmarkPoints = other.m_LandmarkPoints;
	m_ImagePoints = other.m_ImagePoints;
	m_IsLandmarkPoints = other.m_IsLandmarkPoints;
	m_IsImagePoints = other.m_IsImagePoints;
	
 }


// template <class TScalar, unsigned int Dimension>
// void
// AbstractDeformations<TScalar, Dimension>
// ::CopyInformation(const AbstractDeformations& other)
// {
//  	m_Type = other.m_Type;
//  	m_Modified = true;
//  	m_DeformableObjectModified = true;
// 	
// }



template <class TScalar, unsigned int Dimension>
void
AbstractDeformations<TScalar, Dimension>
::SetDeformableMultiObject(DeformableMultiObjectType* DMO)
{
	m_DeformableMultiObject = DMO;
	m_LandmarkPoints = DMO->GetLandmarkPoints();
	m_ImagePoints = DMO->GetDownSampledImageMap();
	m_DownSampledImage = DMO->GetDownSampledImage();
	m_Image = DMO->GetImage();
	
	m_IsLandmarkPoints = (DMO->GetNumberOfLandmarkKindObjects() != 0);
	m_IsImagePoints = (DMO->GetNumberOfImageKindObjects() == 1);
	
	m_DeformableObjectModified = true; 
	
}



#endif /* _AbstractDeformations_txx */

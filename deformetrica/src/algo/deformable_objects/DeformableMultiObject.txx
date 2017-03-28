/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableMultiObject_txx
#define _DeformableMultiObject_txx

#include "DeformableMultiObject.h"
#include "DeformableObject.h"
#include "Landmark.h"
#include "LinearInterpImage.h"

#include "GridFunctions.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#include "KernelFactory.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::DeformableMultiObject()
 {
	m_ImageIndex = -1;

	m_NumberOfImageKindObjects = 0;
	m_NumberOfLandmarkKindObjects = 0;
	m_NumberOfObjects = 0;
	
	m_NumberOfLandmarkPoints = 0;
	m_NumberOfImagePoints = 0;
	
	m_ObjectList.resize(0);
	m_NumberOfPoints.resize(0);
	
 }



template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::~DeformableMultiObject()
{ 
	for (unsigned int i = 0; i < m_ObjectList.size(); i++)
		delete m_ObjectList[i];
}



template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::DeformableMultiObject(const DeformableMultiObject& other)
 {
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(other.m_ObjectList.size());
	for (unsigned int i = 0; i < other.m_ObjectList.size(); i++)
		m_ObjectList[i] = other.m_ObjectList[i]->Clone();
		
	
	m_NumberOfObjects = other.m_NumberOfObjects;
	m_NumberOfImageKindObjects = other.m_NumberOfImageKindObjects;
	m_NumberOfLandmarkKindObjects = other.m_NumberOfLandmarkKindObjects;

	m_ImageIndex = other.m_ImageIndex;

	m_NumberOfLandmarkPoints = other.m_NumberOfLandmarkPoints;
	m_NumberOfImagePoints = other.m_NumberOfImagePoints;
	
	m_NumberOfPoints.resize(m_ObjectList.size());
	for (unsigned int i = 0; i < m_ObjectList.size(); i++)
		m_NumberOfPoints[i] = other.m_NumberOfPoints[i];
	
	m_BoundingBox = other.m_BoundingBox;
 }


template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::DeformableMultiObject(const DeformableMultiObject& Example, const MatrixType& LandmarkPoints, const MatrixType& DownSampledImageMap)
{
	if (Example.m_NumberOfLandmarkPoints != LandmarkPoints.rows() || Example.m_NumberOfImagePoints != DownSampledImageMap.rows() )
		throw std::runtime_error("Number of Landmark and/or Image points mismatch in copy/update DeformableMultiObject constructor ");

	MatrixList List(Example.m_NumberOfObjects);
	int counter = 0;
	for (int i = 0; i < Example.m_NumberOfObjects; i++)
	{
		if (Example.m_ObjectList[i]->IsOfLandmarkKind())
		{
			List[i] = LandmarkPoints.get_n_rows(counter, Example.m_NumberOfPoints[i]);
			counter += Example.m_NumberOfPoints[i];
		}
		else if (Example.m_ObjectList[i]->IsOfImageKind())
		{
			List[i] = DownSampledImageMap;
		}
		else
			throw std::runtime_error("unknown object type");
	}
			
	m_ObjectList.resize(Example.m_ObjectList.size());
	for (unsigned int i = 0; i < Example.m_ObjectList.size(); i++)
		m_ObjectList[i] = Example.m_ObjectList[i]->DeformedObject(List[i]);
		
	this->Update();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::SetObjectList(DeformableObjectList& objList)
 {
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(objList.size());
	for (unsigned int i = 0; i < objList.size(); i++)
		m_ObjectList[i] = objList[i]->Clone();
 }


template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::SetCMRepObjectList(DeformableObjectList& objList)
 {
	// Deep copy so we can safely modify the objects
	m_CMRepObjectList.resize(objList.size());
	for (unsigned int i = 0; i < objList.size(); i++)
		m_CMRepObjectList[i] = objList[i]->Clone();
 }


template<class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::ImageTypePointer
DeformableMultiObject<TScalar, Dimension>
::GetImage() const
{
	if (m_NumberOfImageKindObjects == 0)
		return NULL;
		
	if(m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetImage() - There are more than one image (it should not be).");

	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	LIImageType* tempLII = dynamic_cast<LIImageType*>( m_ObjectList[m_ImageIndex] );
	return tempLII->GetImage();
	
}

template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::ImageTypePointer
DeformableMultiObject<TScalar, Dimension>
::GetDownSampledImage() const
 {
	if(m_NumberOfImageKindObjects == 0)
		return NULL;

	if(m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetDownSampledImage() - There are more than one image (it should not be).");

	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	LIImageType* tempLII = dynamic_cast<LIImageType*>(m_ObjectList[m_ImageIndex]);

	return tempLII->GetDownSampledImage();
 }

//template<class TScalar, unsigned int Dimension>
//TScalar DeformableMultiObject<TScalar, Dimension>
//::GetTimePoint()
// {
//	return m_ObjectList[0]->GetTimePoint();
// }



//template<class TScalar, unsigned int Dimension>
//unsigned int DeformableMultiObject<TScalar, Dimension>
//::GetTimeIndex()
// {
//	return m_ObjectList[0]->GetTimeIndex();
// }

// 
// template <class TScalar, unsigned int Dimension>
// void
// DeformableMultiObject<TScalar, Dimension>
// ::SetImagePoints(const MatrixType& Y)
// {
// 	if (Y.rows() != m_NumberOfImagePoints)
// 		throw std::runtime_error("number of image points mismatch");
// 	
// 	typedef LinearInterpImage<TScalar, Dimension> LIImageType;	
// 	LIImageType* obj = dynamic_cast<LIImageType*>(m_ObjectList[m_ImageIndex]);
// 	
// 	obj->SetImagePoints(Y);
// }
// 
// 

// template <class TScalar, unsigned int Dimension>
// void
// DeformableMultiObject<TScalar, Dimension>
// ::SetLandmarkPoints(const MatrixType& Y)
// {
// 	if (Y.rows() != m_NumberOfLandmarkPoints)
// 		throw std::runtime_error("number of landmark points mismatch");
// 	
// 	typedef Landmark<TScalar, Dimension> LandmarkType;
// 	
// 	int counter = 0;
// 	for (int i = 0; i < m_NumberOfObjects; i++)
// 		if ( m_ObjectList[i]->IsOfLandmarkKind() )
// 		{
// 			LandmarkType* obj = dynamic_cast<LandmarkType*>(m_ObjectList[i]);
// 			obj->SetPointCoordinates(
// 				Y.get_n_rows(counter, m_NumberOfPoints[i]) );
// 		
// 			counter += m_NumberOfPoints[i];
// 		}
// }


template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::MatrixType
DeformableMultiObject<TScalar, Dimension>
::GetLandmarkPoints() const
 {
 	if(m_NumberOfLandmarkKindObjects == 0)
 		return MatrixType(0,0);
	
	typedef Landmark<TScalar, Dimension> LandmarkType;

	MatrixType X0(m_NumberOfLandmarkPoints,Dimension);
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			MatrixType Xaux = dynamic_cast<LandmarkType*>(m_ObjectList[i])->GetPointCoordinates();
			for (int r = 0; r < m_NumberOfPoints[i] ; r++)
				X0.set_row(counter + r, Xaux.get_row(r));
			counter += m_NumberOfPoints[i];
		}
	return X0;
}


template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::MatrixType
DeformableMultiObject<TScalar, Dimension>
::GetDownSampledImageMap() const
 {
	if(m_NumberOfImageKindObjects == 0)
		return MatrixType(0,0);

	if(m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetDownSampledImageMap() - There are more than one image (it should not be).");

	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	LIImageType* tempLII = dynamic_cast<LIImageType*>(m_ObjectList[m_ImageIndex]);

	return tempLII->GetDownSampledImageMap();
 }


template <class TScalar, unsigned int Dimension>
typename std::vector<int>
DeformableMultiObject<TScalar, Dimension>
::GetDimensionOfDiscretizedObjects() const
 {
	 std::vector<int> out(m_NumberOfObjects);
	 for (int i = 0; i < m_NumberOfObjects; i++)
		 out[i] = m_ObjectList[i]->GetDimensionOfDiscretizedObject();
	 
	 return out;
 }



template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::MatrixList
DeformableMultiObject<TScalar, Dimension>
::GetImageIntensityAndLandmarkPointCoordinates() const
 {
	MatrixList Y(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			typedef Landmark<TScalar, Dimension> LandmarkType;
			LandmarkType* obj = dynamic_cast<LandmarkType*>(m_ObjectList[i]);
			Y[i] = obj->GetPointCoordinates();
		}
		else if ( m_ObjectList[i]->IsOfImageKind() )
		{
			typedef LinearInterpImage<TScalar, Dimension> LIImageType;
			typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
			LIImageType* obj = dynamic_cast<LIImageType*>(m_ObjectList[i]);
			MatrixType Yaux(obj->GetImage()->GetLargestPossibleRegion().GetNumberOfPixels(), 1);
			Yaux.set_column(0,
					GridFunctionsType::VectorizeImage(obj->GetImage()) );
			Y[i] = Yaux;
		}
		else
			throw std::runtime_error("Unknown object type");
	}

	return Y;
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::UpdateImageIntensityAndLandmarkPointCoordinates(const MatrixList& Y)
 {
	if (Y.size() != m_NumberOfObjects)
		throw std::runtime_error("size of list and number of objects mismatch");

	typedef Landmark<TScalar, Dimension> LandmarkType;
	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			LandmarkType* obj = dynamic_cast<LandmarkType*>(m_ObjectList[i]);
			obj->UpdatePolyDataWithPointCoordinates(Y[i]);
		}
		else if ( m_ObjectList[i]->IsOfImageKind() )
		{
			LIImageType* obj = dynamic_cast<LIImageType*>(m_ObjectList[i]);
			obj->UpdateImageIntensity(Y[i]);
		}
		else
			throw std::runtime_error("Unknown object type");
	}
 }


template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar,Dimension>::MatrixType
DeformableMultiObject<TScalar, Dimension>
::SplatDifferenceImage(const DeformableMultiObject* targ, const MatrixType& DownSampledImageMap) const
{

	if (m_NumberOfImageKindObjects)
	{
	 	typedef LinearInterpImage<TScalar, Dimension> LIImageType;
	 	LIImageType* tempImage = dynamic_cast<LIImageType*>(this->GetObjectList()[ m_ImageIndex ]);
	 	LIImageType* targImage = dynamic_cast<LIImageType*>(targ->GetObjectList()[ targ->GetImageIndex() ]);
	 	MatrixType splat = tempImage->SplatDifferenceImage(targImage, DownSampledImageMap);
		return splat;		
	}

	return MatrixType(0,0);
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::Update()
{

	m_NumberOfObjects = m_ObjectList.size();

	m_NumberOfImageKindObjects = 0;
	m_NumberOfLandmarkKindObjects = 0;
	m_NumberOfLandmarkPoints = 0;
	m_NumberOfImagePoints = 0;

	m_NumberOfPoints.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		// Update each object in the list
		if (m_ObjectList[i] == NULL)
			throw std::runtime_error("All deformable objects have not been set");
		
		m_ObjectList[i]->Update();
		
		int nb = m_ObjectList[i]->GetNumberOfPoints();
		m_NumberOfPoints[i] = nb;

		if ( m_ObjectList[i]->IsOfImageKind() )
		{
			m_NumberOfImageKindObjects += 1;
			m_NumberOfImagePoints = nb;
			m_ImageIndex = i;
		}

		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{	
			m_NumberOfLandmarkPoints += nb;
			m_NumberOfLandmarkKindObjects += 1;
		}
	}

	if ( m_NumberOfImageKindObjects > 1 )
		throw std::runtime_error("In DeformableMultiObject::SetObjectList() - There are more than one image (it should not be).");
	// in case of multi-modalities, it would be better to add a new class multi-image, since the same image domain should be shared by all images
	// and the eta maps in TransportAlongGeodesic could be added, as well as the initial eta0 maps (provided they are weighted by their respective dataSigma)
	
	this->UpdateBoundingBox();
	
}

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::UpdateBoundingBox()
{
	MatrixType BB(Dimension,2);
	
	BB = m_ObjectList[0]->GetBoundingBox();
	
	for (int i = 1; i < m_NumberOfObjects; i++)
	{
		for (int d = 0; d < Dimension; d++)
		{
			MatrixType BBaux = m_ObjectList[i]->GetBoundingBox();
			BB(d,0) = (BB(d,0)<BBaux(d,0)?BB(d,0):BBaux(d,0));
			BB(d,1) = (BB(d,1)>BBaux(d,1)?BB(d,1):BBaux(d,1));
		}
	}

	m_BoundingBox = BB;

}



template <class TScalar, unsigned int Dimension>
typename std::vector<TScalar>
DeformableMultiObject<TScalar, Dimension>
::ComputeMatch(const DeformableMultiObject* target)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

	DeformableObjectList targetList = target->GetObjectList();
	DeformableObjectList cmRepList = target->GetCMRepObjectList();

	std::vector<TScalar> match(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		match[i] = m_ObjectList[i]->ComputeMatch(targetList[i]);
	}

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::MatrixList
DeformableMultiObject<TScalar, Dimension>
::ComputeMatchGradient(const DeformableMultiObject* target)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

	DeformableObjectList targetList = target->GetObjectList();
	MatrixList gradmatch(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradmatch[i] = m_ObjectList[i]->ComputeMatchGradient(targetList[i]);

	return gradmatch;
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::ListToMatrices(const MatrixList& L, MatrixType& MLandmark, MatrixType& MImage) const
{
	if (L.size() != m_NumberOfObjects)
		throw std::runtime_error("number of objects mismatch in DeformableMultiObject::ListToMatrices");
	
	int dim = L[0].columns();
	
	MLandmark.set_size(m_NumberOfLandmarkPoints, dim);
	MImage.set_size(m_NumberOfImagePoints, dim);
	
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (m_ObjectList[i]->IsOfLandmarkKind())
		{
			for (int r = 0; r < m_NumberOfPoints[i]; r++)
			{
				MLandmark.set_row(counter + r, L[i].get_row(r));
			}
			counter += m_NumberOfPoints[i];
		}
		else if (m_ObjectList[i]->IsOfImageKind())
		{
			MImage = L[i];
		}
		else
			throw std::runtime_error("unknown object type");
	}
}


template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::MatricesToList(MatrixList& L, const MatrixType& MLandmark, const MatrixType& MImage) const
{
	L.resize(m_NumberOfObjects);
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (m_ObjectList[i]->IsOfLandmarkKind())
		{
			L[i] = MLandmark.get_n_rows(counter, m_NumberOfPoints[i]);
			counter += m_NumberOfPoints[i];
		}
		else if (m_ObjectList[i]->IsOfImageKind())
		{
			L[i] = MImage;
		}
		else
			throw std::runtime_error("unknown object type");
	}
}



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::WriteMultiObject(std::vector<std::string>& outfn) const
 {
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteObject(outfn[i]);
 }


template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::WriteMultiObject(std::vector<std::string>& outfn, const MatrixList& velocity) const
 {
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteObject(outfn[i], velocity);
 }

#endif /* _DeformableMultiObject_txx */

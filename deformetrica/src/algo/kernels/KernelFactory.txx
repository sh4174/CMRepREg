/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _KernelFactory_txx
#define _KernelFactory_txx

#ifdef USE_CUDA
	#include "CUDAExactKernel.h"
#endif
#include "ExactKernel.h"
#include "FGTKernel.h"
#include "P3MKernel.h"


#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>*
KernelFactory<TScalar, PointDim>::m_SingletonInstance = 0;


template<class TScalar, unsigned int PointDim>
itk::SimpleFastMutexLock
KernelFactory<TScalar, PointDim>::m_Mutex;


template<class TScalar, unsigned int PointDim>
bool
KernelFactory<TScalar, PointDim>::m_IsInstantiated = false;


template<class TScalar, unsigned int PointDim>
typename KernelFactory<TScalar, PointDim>::MatrixType
KernelFactory<TScalar, PointDim>::m_DataDomain;


template<class TScalar, unsigned int PointDim>
TScalar
KernelFactory<TScalar, PointDim>::m_PaddingFactor = 0;


template<class TScalar, unsigned int PointDim>
TScalar
KernelFactory<TScalar, PointDim>::m_WorkingSpacingRatio = 0;



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>
::KernelFactory() {}



template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>
::~KernelFactory() {}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>*
KernelFactory<TScalar, PointDim>
::Instantiate()
 {
	if (!m_IsInstantiated)
	{
		m_Mutex.Lock();
		if (m_SingletonInstance == 0)
			m_SingletonInstance = new KernelFactory<TScalar, PointDim>();
		m_Mutex.Unlock();
		m_IsInstantiated = true;
	}

	return m_SingletonInstance;
 }



template<class TScalar, unsigned int PointDim>
typename KernelFactory<TScalar, PointDim>::KernelBaseType*
KernelFactory<TScalar, PointDim>
::CreateKernelObject(KernelEnumType kernelType)
 {
	// if (kernelType == null)
	// 	throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");

	switch (kernelType)
	{
	case Exact:
		return new ExactKernel<TScalar, PointDim>();
#ifdef USE_CUDA
	case CUDAExact:
		return new CUDAExactKernel<TScalar, PointDim>();
#endif
	case FGT:
		return new FGTKernel<TScalar, PointDim>();
	case P3M:
	{
		typedef P3MKernel<TScalar, PointDim> P3MKernelType;
		P3MKernelType* obj = new P3MKernelType();
		if (m_DataDomain.size())
			obj->SetDataDomain(this->GetDataDomain());

		obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
		obj->SetPaddingFactor(this->GetPaddingFactor());
		return obj;
	}
	default:
		throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");
	}
	return 0;
 }

 template<class TScalar, unsigned int PointDim>
 typename KernelFactory<TScalar, PointDim>::KernelBaseType*
 KernelFactory<TScalar, PointDim>
 ::CreateKernelObject(KernelEnumType kernelType, MatrixType DataDomain)
  {
 	// if (kernelType == null)
 	// 	throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");

 	switch (kernelType)
 	{
 	case Exact:
 		return new ExactKernel<TScalar, PointDim>();
 #ifdef USE_CUDA
 	case CUDAExact:
 		return new CUDAExactKernel<TScalar, PointDim>();
 #endif
 	case FGT:
 		return new FGTKernel<TScalar, PointDim>();
 	case P3M:
 	{
 		typedef P3MKernel<TScalar, PointDim> P3MKernelType;
 		P3MKernelType* obj = new P3MKernelType();
		obj->SetDataDomain(DataDomain);
 		obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
 		obj->SetPaddingFactor(this->GetPaddingFactor());
 		return obj;
 	}
 	default:
 		throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");
 	}
 	return 0;
  }


template<class TScalar, unsigned int PointDim>
typename KernelFactory<TScalar, PointDim>::KernelBaseType*
KernelFactory<TScalar, PointDim>
::CreateKernelObject(
		const MatrixType& X,
		const MatrixType& W,
		TScalar h,
		KernelEnumType kernelType)
{
	switch (kernelType)
	{
	case Exact:
		return new ExactKernel<TScalar, PointDim>(X, W, h);
#ifdef USE_CUDA
		return new CUDAExactKernel<TScalar, PointDim>(X, W, h);
#endif
	case FGT:
		return new FGTKernel<TScalar, PointDim>(X, W, h);
	case P3M:
	{
		typedef P3MKernel<TScalar, PointDim> P3MKernelType;
		P3MKernelType* obj = new P3MKernelType(X, W, h);
		if (m_DataDomain.size())
			obj->SetDataDomain(this->GetDataDomain());
		obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
		obj->SetPaddingFactor(this->GetPaddingFactor());
		return obj;
	}
	default:
		throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");	
	}
	return 0;
}



#endif /* _KernelFactory_txx */

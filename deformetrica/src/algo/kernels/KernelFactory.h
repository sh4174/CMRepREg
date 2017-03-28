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

#ifndef _KernelFactory_h
#define _KernelFactory_h

// #include "itkImage.h"
#include "itkSimpleFastMutexLock.h"

#include "KernelType.h"

#include "ExactKernel.h"


/**
 *	\brief      A kernel factory.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The KernelFactory class enables to instantiate objects whose type is derived from an
 *              abstract type (it is the principle of the factory method pattern). On top of that, this
 *              class implements the singleton pattern, that is to say that you can only instantiate
 *              one factory.
 */
template <class TScalar, unsigned int PointDim>
class KernelFactory
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> KernelBaseType;

	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;




	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

/*	/// Returns the type of the kernel.
	KernelEnumType GetKernelCode() const { return m_WhichKernel; }
	/// Sets the type of the kernel to ExactKernel.
	void UseExactKernel() { m_Mutex.Lock(); m_WhichKernel = Exact; m_Mutex.Unlock();}
#ifdef USE_CUDA
	/// Sets the type of the kernel to CUDAExactKernel.
	void UseCUDAExactKernel() { m_Mutex.Lock(); m_WhichKernel = CUDAExact; m_Mutex.Unlock();}
#endif
	/// Sets the type of the kernel to FGTKernel.
	void UseFGTKernel() { m_Mutex.Lock(); m_WhichKernel = FGT; m_Mutex.Unlock();}
	/// Sets the type of the kernel to P3MKernel.
	void UseP3MKernel() { m_Mutex.Lock(); m_WhichKernel = P3M; m_Mutex.Unlock();}*/

	/// See AbstractDeformations::GetDataDomain() for details.
	static MatrixType GetDataDomain() { return m_DataDomain; }
	/// See AbstractDeformations::SetDataDomain() for details.
	static void SetDataDomain(MatrixType DD) { m_Mutex.Lock(); m_DataDomain = DD; m_Mutex.Unlock(); }

	/// See AbstractDeformations::GetWorkingSpacingRatio() for details.
	static TScalar GetWorkingSpacingRatio() { return m_WorkingSpacingRatio; }
	/// See AbstractDeformations::SetWorkingSpacingRatio() for details.
	static void SetWorkingSpacingRatio(TScalar d) { m_Mutex.Lock(); m_WorkingSpacingRatio = d; m_Mutex.Unlock(); }

	/// See AbstractDeformations::GetPaddingFactor() for details.
	static TScalar GetPaddingFactor() { return m_PaddingFactor; }
	/// See AbstractDeformations::SetPaddingFactor() for details.
	static void SetPaddingFactor(TScalar d) { m_Mutex.Lock(); m_PaddingFactor = d; m_Mutex.Unlock(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Instantiates an object of KernelFactory type with the Singleton strategy.
	static KernelFactory<TScalar, PointDim>* Instantiate();

	/// Returns the instance of the object, NULL in case of error.
	KernelBaseType* CreateKernelObject(KernelEnumType kernelType);

	/// Returns the instance of the object, NULL in case of error.
	KernelBaseType* CreateKernelObject(KernelEnumType kernelType, MatrixType DataDomain);

	/// Returns the instance of the object, NULL in case of error.
	KernelBaseType* CreateKernelObject(
			const MatrixType& X,
			const MatrixType& W,
			TScalar h, KernelEnumType kernelType);

// NOTE: make sure to assign the target image when using P3M
//   void SetWorkingImage(ImageType* img);
// ImagePointer GetWorkingImage() const { return m_WorkingImage; }



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	KernelFactory();

	~KernelFactory();

// void DetermineP3MGrid(
//   ImagePointType& outOrigin,
//   ImageSpacingType& outSpacing,
//   ImageSizeType& outSize,
//   ImageType* img, TScalar h);

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

//	/// Type of the kernel.
//	KernelEnumType m_WhichKernel;

	// Information of image domain for heuristically determining parameters
	// ImagePointer m_WorkingImage;



private:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// See AbstractDeformations::m_DataDomain for details.
	static MatrixType m_DataDomain;

	/// See P3MKernel::m_WorkingSpacingRatio for details.
	static TScalar m_WorkingSpacingRatio;

	/// See P3MKernel::m_PaddingFactor for details.
	static TScalar m_PaddingFactor;

	///	Object used to perform mutex (important for multithreaded programming).
	static itk::SimpleFastMutexLock m_Mutex;

	/// Boolean which enables to instantiate once a kernel factory.
	static bool m_IsInstantiated;

	/// The unique kernel factory.
	static KernelFactory* m_SingletonInstance;


}; /* class KernelFactory */


#ifndef MU_MANUAL_INSTANTIATION
#include "KernelFactory.txx"
#endif


#endif /* _KernelFactory_h */

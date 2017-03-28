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

#ifndef _CUDAExactKernel_h
#define _CUDAExactKernel_h

#include "ExactKernel.h"


/**
 *	\brief      An exact kernel with CUDA implementation.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.1
 *
 *	\details    The CUDAExactKernel class inherited from ExactKernel is a reimplementation of some methods
 *              with Cuda in order to take advantage of the parallization with the GPU.
 */
template <class TScalar, unsigned int PointDim>
class CUDAExactKernel : public ExactKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> Superclass;
	/// Abstract kernel type.
	typedef typename Superclass::Superclass AbstractKernel;

	/// Vector type
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	CUDAExactKernel(): Superclass() { }
	/// Copy constructor.
	CUDAExactKernel(const CUDAExactKernel& o);
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, double h).
	CUDAExactKernel(const MatrixType& X, double h): Superclass(X, h) {}
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, const MatrixType& W, double h).
	CUDAExactKernel(const MatrixType& X, const MatrixType& W, double h): Superclass(X, W, h) {}

	virtual CUDAExactKernel* Clone() const {
		return new CUDAExactKernel(*this);
	}

	virtual ~CUDAExactKernel() { }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual MatrixType Convolve(const MatrixType& X);

	virtual MatrixType ConvolveGradient(const MatrixType& X, const MatrixType& alpha);

	virtual MatrixType ConvolveSpecialHessian(const MatrixType& xi);


protected:


}; /* class CUDAExactKernel */



#ifndef MU_MANUAL_INSTANTIATION
#include "CUDAExactKernel.txx"
#endif



#endif /* _CUDAExactKernel_h */

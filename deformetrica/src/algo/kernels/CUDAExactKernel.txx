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

#ifndef _CUDAExactKernel_txx
#define _CUDAExactKernel_txx

#include "CUDAExactKernel.h"

#include "../../lib/cuda_convolutions/GpuConv1D.h"
#include "../../lib/cuda_convolutions/GpuConv2D.h"

#include "writeMatrixDLM.txx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
CUDAExactKernel<TScalar, PointDim>
::CUDAExactKernel(const CUDAExactKernel& o)
 {
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::MatrixType
CUDAExactKernel<TScalar, PointDim>
::Convolve(const MatrixType& X) {

	if (this->GetSources().rows() != this->GetWeights().rows())
		throw std::runtime_error("Sources and weights count mismatch");

	int DimVect = this->GetWeights().columns();
	TScalar* gammap = new TScalar [X.rows()*DimVect];

	if (DimVect == PointDim)
		GaussGpuEvalConv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				X.vectorise_row_wise().memptr(), this->GetSources().vectorise_row_wise().memptr(), this->GetWeights().vectorise_row_wise().memptr(),
				gammap, X.rows(), this->GetSources().rows());
	else if (DimVect == 2*PointDim)
		GaussGpuEvalConv1D<TScalar, PointDim, (2*PointDim)>(this->GetKernelWidth(),
				X.vectorise_row_wise().memptr(), this->GetSources().vectorise_row_wise().memptr(), this->GetWeights().vectorise_row_wise().memptr(),
				gammap, X.rows(), this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::Convolve(X) - Invalid number of rows of beta !");


	MatrixType gamma(gammap, DimVect, X.rows());
	return gamma.transpose();

}



template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::MatrixType
CUDAExactKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, const MatrixType& alpha)
 {
	int DimVect = this->GetWeights().columns();
	TScalar* gammap = new TScalar [X.rows()*DimVect];

	if (DimVect == PointDim)
		GaussGpuGrad1Conv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				alpha.vectorise_row_wise().memptr(), X.vectorise_row_wise().memptr(),
				this->GetSources().vectorise_row_wise().memptr(), this->GetWeights().vectorise_row_wise().memptr(),
				gammap, X.rows(), this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::ConvolveGradient(X, alpha) - Problem with the number of rows of beta !");

	MatrixType gamma(gammap, DimVect, X.rows());

	return gamma.transpose();
 }



template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::MatrixType
CUDAExactKernel<TScalar, PointDim>
::ConvolveSpecialHessian(const MatrixType& xi)
{
	if (this->GetSources().rows() != this->GetWeights().rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (this->GetSources().rows() != xi.rows())
		throw std::runtime_error("Y and Xi count mismatch");

	/*
	int DimVect = this->GetWeights().columns();
	TScalar *xi_p, *y_p, *beta_p, *gamma_p;
//	xi_p    = (TScalar*)xi.memory_pointer_row_wise();
	xi_p    = xi.vectorise_row_wise().memptr();
	y_p     = this->GetSources().vectorise_row_wise().memptr();
	beta_p  = this->GetWeights().vectorise_row_wise().memptr();
	gamma_p = new TScalar[xi.rows()*DimVect];

	if (DimVect == PointDim)
		GaussGpuGradDiffConv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				y_p, beta_p, xi_p, gamma_p, this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::ConvolveSpecialHessian(xi) - Invalid number of rows of beta !");

	MatrixType gamma(gamma_p, DimVect, xi.rows());
	return (gamma.transpose() * (- 0.5));
	*/
	int DimVect = this->GetWeights().columns();
	TScalar* gamma_p = new TScalar[xi.rows()*DimVect];

	if (DimVect == PointDim)
		GaussGpuGradDiffConv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				this->GetSources().vectorise_row_wise().memptr(), this->GetWeights().vectorise_row_wise().memptr(),
				xi.vectorise_row_wise().memptr(), gamma_p, this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::ConvolveSpecialHessian(xi) - Invalid number of rows of beta !");

	MatrixType gamma(gamma_p, DimVect, xi.rows());
	return (- 0.5 * gamma.transpose());
}



#endif /* _CUDAExactKernel_txx */

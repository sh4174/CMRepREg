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

#ifndef _ExactKernel_txx

#include "ExactKernel.h"

#include <exception>
#include <stdexcept>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
ExactKernel<TScalar, PointDim>
::ExactKernel(const ExactKernel& o)
 {
	//this->SetSources(o.GetSources());
	//this->SetWeights(o.GetWeights());
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();

 }



template <class TScalar, unsigned int PointDim>
ExactKernel<TScalar, PointDim>*
ExactKernel<TScalar, PointDim>
::Clone() const
 {
	return new ExactKernel(*this);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::MatrixType
ExactKernel<TScalar, PointDim>
::Convolve(const MatrixType& X)
 {
#ifdef DEFORMETRICA_WITH_OLD_ALGORITHM
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	MatrixType V(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		VectorType vi(weightDim, 0.0);
		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			vi += wj * this->EvaluateKernel(xi, yj);
		}

		V.set_row(i, vi);
	}

	return V;
#else
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	MatrixType V(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			TScalar Kij = this->EvaluateKernel(X, Y, i, j);

			for(unsigned int k=0; k< weightDim; k++)
				V(i,k) += W(j,k) * Kij;
		}
	}

	return V;
#endif
 }



template <class TScalar, unsigned int PointDim>
std::vector<typename ExactKernel<TScalar, PointDim>::MatrixType>
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X)
 {
#ifdef DEFORMETRICA_WITH_OLD_ALGORITHM
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector<MatrixType> gradK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		MatrixType Gi(weightDim, PointDim, 0.0);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			VectorType g = this->EvaluateKernelGradient(xi, yj);
			for (unsigned int k = 0; k < weightDim; k++)
				Gi.set_row(k, Gi.get_row(k) + g*wj[k]);
		}

		gradK.push_back(Gi);
	}

	return gradK;
#else

	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector<MatrixType> gradK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
	// 	VectorType xi = X.get_row(i);

		MatrixType Gi(weightDim, PointDim, 0.0);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			// VectorType yj = Y.get_row(j);

			VectorType g = this->EvaluateKernelGradient(X, Y, i, j);
			for (unsigned int k = 0; k < weightDim; k++)
			{
				TScalar Wjk = W(j,k);
				for (unsigned int l = 0; l < PointDim; l++)
				{
					Gi(k,l) +=  g[l]*Wjk;
				}
			}
		}

		gradK.push_back(Gi);
	}

	return gradK;
#endif
}



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::MatrixType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, const MatrixType& alpha)
 {
	std::vector<MatrixType> convolveGradient = this->ConvolveGradient(X);

	MatrixType result(X.rows(), PointDim, 0);
	for (unsigned int j = 0; j < X.rows(); j++)
		result.set_row(j, convolveGradient[j].transpose() * alpha.get_row(j) );

	return result;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::MatrixType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, unsigned int dim)
 {
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (dim >= PointDim || dim < 0)
		throw std::runtime_error("dimension index out of bounds");

	unsigned int weightDim = W.columns();

	MatrixType gradK(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			VectorType g = this->EvaluateKernelGradient(xi, yj);
			gradK.set_row(i, gradK.get_row(i) + g[dim]*wj);
		}
	}

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VectorType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const MatrixType& X, unsigned int k, unsigned int dp)
 {
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	if (k >= weightDim)
		throw std::runtime_error("Invalid weight index");

	if (dp >= PointDim)
		throw std::runtime_error("Invalid derivative direction");

	unsigned int numPoints = Y.rows();

	VectorType gradK(numPoints, 0);

	for (unsigned int i = 0; i < numPoints; i++)
	{
		VectorType xi = X.get_row(i);

		gradK[i] = 0;

		for (unsigned int j = 0; j < numPoints; j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			VectorType g = this->EvaluateKernelGradient(xi, yj);

			gradK[i] += wj[k] * g[dp];
		}
	}

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
std::vector< std::vector<typename ExactKernel<TScalar, PointDim>::MatrixType> >
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X)
 {
#ifdef DEFORMETRICA_WITH_OLD_ALGORITHM
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector< std::vector<MatrixType> > hessK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		std::vector<MatrixType> Hi;
		for (unsigned int k = 0; k < weightDim; k++)
			Hi.push_back(MatrixType(PointDim, PointDim, 0.0));

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			MatrixType H = this->EvaluateKernelHessian(xi, yj);
			for (unsigned int k = 0; k < weightDim; k++)
				Hi[k] = Hi[k] + H*wj[k];
		}

		hessK.push_back(Hi);
	}

	return hessK;
#else
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector< std::vector<MatrixType> > hessK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
//		VectorType xi = X.get_row(i);

		std::vector<MatrixType> Hi;
		for (unsigned int k = 0; k < weightDim; k++)
			Hi.push_back(MatrixType(PointDim, PointDim, 0.0));

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
//			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			MatrixType H = this->EvaluateKernelHessian(X, Y, i, j);
			for (unsigned int k = 0; k < weightDim; k++)
				Hi[k] = Hi[k] + H*wj[k];
		}

		hessK.push_back(Hi);
	}

	return hessK;
#endif
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::MatrixType
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X, unsigned int row, unsigned int col)
 {
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (row >= PointDim || col >= PointDim)
		throw std::runtime_error("Dimension index out of bounds");

	unsigned int weightDim = W.columns();

	MatrixType hessK(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VectorType xi = X.get_row(i);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			MatrixType H = this->EvaluateKernelHessian(xi, yj);
			hessK.set_row(i, hessK.get_row(i) + wj*H(row,col) );
		}
	}

	return hessK;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VectorType
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const MatrixType& X, unsigned int k,
		unsigned int dp, unsigned int dq)
		{
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	if (k >= weightDim)
		throw std::runtime_error("Invalid weight index");

	if (dp >= PointDim || dq >= PointDim)
		throw std::runtime_error("Invalid derivative direction");

	unsigned int numPoints = X.rows();

	VectorType hessK(numPoints, 0);

	for (unsigned int i = 0; i < numPoints; i++)
	{
		VectorType xi = X.get_row(i);

		hessK[i] = 0;

		// MARCEL
		// for (unsigned int j = 0; j < numPoints; j++)
		// STANLEY
		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VectorType yj = Y.get_row(j);
			VectorType wj = W.get_row(j);

			MatrixType H = this->EvaluateKernelHessian(xi, yj);

			hessK[i] += wj[k] * H(dp, dq);

			// DEBUG
			// if (k==0 && dp==0 && dq==1)
			// {
			// 	std::cout << "H(dp,dq) = " << H(dp,dq) << std::endl;
			// 	std::cout << "hessK[i] = " << hessK[i] << std::endl;
			// }

		}
	}

	return hessK;
		}


#endif /* _ExactKernel_txx */

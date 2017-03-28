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

#ifndef _AbstractKernel_txx

#include "AbstractKernel.h"

#include <exception>
#include <stdexcept>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename AbstractKernel<TScalar, PointDim>::MatrixType
AbstractKernel<TScalar, PointDim>
::ConvolveSpecialHessian(const MatrixType& xi)
{
	MatrixType& Y = this->GetSources();
	MatrixType& W = this->GetWeights();

	MatrixType beta = W;  // Deep copy (because of this->SetWeights() below)

	const unsigned int N = Y.rows();

	MatrixType dXi4(N, PointDim, 0);

	// dXi4 first term
	// timer.Start();
	std::vector< std::vector<MatrixType> > HessA = this->ConvolveHessian(Y);
	for (unsigned int i = 0; i < N; i++)
	{
		VectorType XiMomi = xi.get_row(i);
		MatrixType Auxi(PointDim, PointDim, 0.0);
		for (unsigned int dim = 0; dim < PointDim; dim++)
			Auxi -= HessA[i][dim] * beta(i,dim);

		dXi4.set_row(i, Auxi * XiMomi);
	}
	// dXi4 second term
	std::vector< MatrixType > H(N);
	MatrixType Zeros(PointDim, PointDim, 0.0);
	for (int i = 0; i < N; i++)
		H[i] = Zeros;

	for (unsigned r = 0; r < PointDim; r++)
	{
		MatrixType Wd(N,PointDim,0.0);
		for (int i = 0; i < N; i++)
			Wd.set_row(i, beta.get_row(i) * xi(i,r) );

		this->SetWeights(Wd);
		MatrixType Hrr = this->ConvolveHessian(Y,r,r);
		for (int i = 0; i < N; i ++)
			for (unsigned int q = 0; q < PointDim; q++)
				H[i](r,q) += Hrr(i,q);

		for (unsigned q = r+1; q < PointDim; q++)
		{
			MatrixType Wtmp(N, 2*PointDim, 0.0);
			for (int i = 0; i < N; i++)
				for (unsigned int d = 0; d < PointDim; d++)
				{
					Wtmp(i,d) = Wd(i,d);
					Wtmp(i, d + PointDim) = beta(i,d) * xi(i,q);
				}

			this->SetWeights(Wtmp);
			MatrixType Hrq = this->ConvolveHessian(Y,r,q);
			for (int i = 0; i < N; i ++)
				for (unsigned int d = 0; d < PointDim; d++)
				{
					H[i](q,d) += Hrq(i,d);
					H[i](r,d) += Hrq(i,d + PointDim);
				}
		}
	}

	for (int i = 0; i < N; i++)
		dXi4.set_row(i, dXi4.get_row(i) + H[i] * beta.get_row(i));

	this->SetWeights(beta);


	return dXi4;

}


#endif /* _AbstractKernel_txx */

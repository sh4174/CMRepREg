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

#ifndef _FGTKernel_h
#define _FGTKernel_h

#include "ExactKernel.h"
#include "FastGauss.h"

/**
 *	\brief      Kernels with Fast Gauss Transform approximation.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.1
 *
 *	\details    The ExactKernel class inherited from AbstractKernel implements the operations
 *              with Gaussian kernels using the Fast Gauss Transforms (i.e. multipole approximation).
 */
template <class TScalar, unsigned int PointDim>
class FGTKernel: public ExactKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> Superclass;

	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;

	/// Fast Gauss type.
	typedef FastGauss<TScalar> FGCalculatorType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//TODO: constructor builds fgt object with weights w and source*w?
	// for conv and convgrad -> grad is x*wconv - source_w_conv
	// ???
	// for hessian needs w and (x-y)*(x-y)' = (x-y)*x' - (x-y)*y'
	// = x*x' - y*x' - x*y' + y*y'
	// so can get away with source*w as well

	FGTKernel();
	/// Copy constructor.
	FGTKernel(const FGTKernel& o);
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, double h) for details.
	FGTKernel(const MatrixType& X, double h);
	/// See AbstractKernel::AbstractKernel(const MatrixType& X, const MatrixType& W, double h) for details.
	FGTKernel(const MatrixType& X, const MatrixType& W, double h);

	virtual FGTKernel* Clone() const;

	virtual ~FGTKernel();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the number of clusters to \e n.
	void SetNumberOfClusters(unsigned int n) { m_NumberOfClusters = n; }

	/// TODO .
	void SetFarRatio(double f) { m_FarRatio = f; }

	/// Sets the truncation order to \e n.
	void SetTruncationOrder(unsigned int n) { m_TruncationOrder = n; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual MatrixType Convolve(const MatrixType& X);

	virtual MatrixType ConvolveGradient(const MatrixType& X, const MatrixType& alpha);
	virtual std::vector<MatrixType> ConvolveGradient(const MatrixType& X);
	virtual VectorType ConvolveGradient(const MatrixType& X, unsigned int k, unsigned int dp);

	virtual std::vector< std::vector<MatrixType> > ConvolveHessian(const MatrixType & X);
	virtual VectorType ConvolveHessian(const MatrixType& X, unsigned int k, unsigned int dp, unsigned int dq);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	void Init();

	/// TODO .
	FGCalculatorType* BuildFGT(const MatrixType& sources, const MatrixType& weights, double h);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	FGCalculatorType* m_FGTObject;

	/// Number of clusters where the point set will be gathered.
	unsigned int m_NumberOfClusters;

	/// Truncation order of the Taylor expansion.
	unsigned int m_TruncationOrder;

	/// TODO .
	double m_FarRatio;


}; /* class FGTKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "FGTKernel.txx"
#endif


#endif /* _FGTKernel_h */

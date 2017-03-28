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

#ifndef _AbstractKernel_h
#define _AbstractKernel_h

#include "LinearAlgebra.h"
#include <vector>

#include <cassert>
#include <cmath>

/**
 *	\brief 		Abstract kernels.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The AbstractKernel class defines standard operations with kernels used in Deformetrica.\n \n
 *	            See children classes (ExactKernel, P3MKernel, FGTKernel) for implementation of these operations.
 */
template <class TScalar, unsigned int PointDim>
class AbstractKernel
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	AbstractKernel()
	{
		this->SetKernelWidth(1.0);
		m_Modified = false;
	}

	/// Constructor with sources initialized to \e X, kernel width to \e h and weights to 1.
	AbstractKernel(const MatrixType& X, double h)
	{
		MatrixType W(X.rows(), 1, 1.0);
		this->SetSources(X); this->SetWeights(W);
		this->SetKernelWidth(h);
		m_Modified = true;
	}

	/// Constructor with sources initialized to \e X, width to \e W and kernel width to \e h.
	AbstractKernel(const MatrixType& X, const MatrixType& W, double h)
	{
		this->SetSources(X); this->SetWeights(W);
		this->SetKernelWidth(h);
		m_Modified = true;
	}

	/// Makes a copy of the object.
	virtual AbstractKernel* Clone() const =  0;

	virtual ~AbstractKernel() { }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the size of the kernel.
	inline double GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e h.
	void SetKernelWidth(double h)
	{
		if (h < 1e-10)
			h = 1e-10;
		m_KernelWidth = h;
		m_KernelWidthSquared = h*h;
		m_Modified = true;
	}

	/// Returns the sources.
	MatrixType& GetSources() { return m_Sources; }
	/// Sets the sources to \e X.
	void SetSources(const MatrixType& X)	{ assert(X.columns() == PointDim); m_Sources = X; m_Modified = true; }

	/// Returns the weights.
	MatrixType& GetWeights() { return m_Weights; }
	/// Sets the weights to \e W.
	void SetWeights(const MatrixType& W) { m_Weights = W; m_Modified = true; }

	/// Returns true if one of the parameters of the kernel has changed, false otherwise.
	inline bool IsModified() const { return m_Modified; }
	/// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
	inline void SetModified() { m_Modified = true; }
	/// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
	inline void UnsetModified() { m_Modified = false; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Evaluates \f$ K(x,y) \f$.
	virtual TScalar EvaluateKernel(const VectorType& x, const VectorType& y) = 0;

	/// Evaluates the gradient of \f$ K(x,y) \f$ at \e x.
	virtual VectorType EvaluateKernelGradient(const VectorType& x, const VectorType& y) = 0;

	/// Evaluates the hessian of \f$ K(x,y) \f$ at \e x.
	virtual MatrixType EvaluateKernelHessian(const VectorType& x, const VectorType& y) = 0;

	/// Computes convolution of the "weights" located at the "source" points with the kernel and provides results at output point \e X .
	virtual MatrixType Convolve(const MatrixType& X) = 0;

	/// Derivative of convolved weight w in rows in direction dp in columns.
	virtual std::vector<MatrixType> ConvolveGradient(const MatrixType& X) = 0;
	/// TODO .
	virtual MatrixType ConvolveGradient(const MatrixType& X, const MatrixType& alpha) = 0;
	/// Derivative of convolved weight \e k at direction \e dp.
	virtual VectorType ConvolveGradient(const MatrixType& X, unsigned int k, unsigned int dp) = 0;
	/// TODO .
	virtual MatrixType ConvolveGradient(const MatrixType& X, unsigned int dim) = 0;

	/// List of Hessian of convolved weight \e k H[point_index][weight_index](dir1, dir2).
	virtual std::vector< std::vector<MatrixType> > ConvolveHessian(const MatrixType & X) = 0;
	/// Second derivative of weight \e k at directions \e dp and \e dq.
	virtual VectorType ConvolveHessian(const MatrixType & X, unsigned int k, unsigned int dp, unsigned int dq) = 0;
	/// TODO .
	virtual MatrixType ConvolveHessian(const MatrixType& X, unsigned int row, unsigned int col) = 0;

	/// Computes the matrix of \f$\gamma_i\f$ where \f$\gamma_i = \sum_j \beta_i^T\beta_j Hess(y_i, y_j)^T(\xi_j-\xi_i)\f$.
	virtual MatrixType ConvolveSpecialHessian(const MatrixType& xi);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix containing coordinates of the sources i.e. the positions (Size : NumberOfPoints x Dimension).
	MatrixType m_Sources;

	///	Matrix containing coordinates of the weights i.e. the momentas  (Size : NumberOfPoints x Dimension).
	MatrixType m_Weights;

	/// Size of the kernel.
	TScalar m_KernelWidth;

	/// Size of the kernel squared.
	TScalar m_KernelWidthSquared;

	///	Boolean which avoids computing the trajectory (via Update()) if no parameter has changed.
	bool m_Modified;


}; /* class AbstractKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "AbstractKernel.txx"
#endif


#endif /* _AbstractKernel_h */

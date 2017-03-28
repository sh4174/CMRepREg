/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _EigenMatrixWrapper_h
#define _EigenMatrixWrapper_h

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Dense>

#define EIGEN_NO_DEBUG

template<class TScalar>
class EigenMatrixWrapper;

template<class TScalar>
class EigenVectorWrapper;



template<class TScalar>
EigenMatrixWrapper<TScalar> operator-(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);
template<class TScalar>
EigenMatrixWrapper<TScalar> operator+(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);
template<class TScalar>
EigenMatrixWrapper<TScalar> operator*(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);

template<class TScalar>
EigenMatrixWrapper<TScalar> operator/(const EigenMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
template<class TScalar>
EigenMatrixWrapper<TScalar> operator*(const EigenMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
template<class TScalar>
EigenMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const EigenMatrixWrapper<TScalar> & rightMatrix);

template<class TScalar>
EigenVectorWrapper<TScalar> operator*(EigenMatrixWrapper<TScalar> const &leftMatrix, EigenVectorWrapper<TScalar> const &rightVector);

template<class TScalar>
TScalar vnl_trace(EigenMatrixWrapper<TScalar> const& M);
template<class TScalar>
EigenMatrixWrapper<TScalar> vnl_inverse(EigenMatrixWrapper<TScalar> const& M);
template<class TScalar>
EigenMatrixWrapper<TScalar> inverse_sympd(EigenMatrixWrapper<TScalar> const& M);

template <class TScalar>
std::ostream & operator<<(std::ostream& os, EigenMatrixWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical matrix class (Eigen).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The EigenMatrixWrapper class is a wrapping for mathematical matrix class, templated by type of element
 *              and using the Eigen library.
 */
template<class TScalar>
class EigenMatrixWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend EigenMatrixWrapper<TScalar> operator-<>(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);
	friend EigenMatrixWrapper<TScalar> operator+<>(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);
	friend EigenMatrixWrapper<TScalar> operator*<>(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right);

	friend EigenMatrixWrapper<TScalar> operator/<>(const EigenMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
	friend EigenMatrixWrapper<TScalar> operator*<>(const EigenMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
	friend EigenMatrixWrapper<TScalar> operator*<>(TScalar const& leftScalar, const EigenMatrixWrapper<TScalar> & rightMatrix);

	friend EigenVectorWrapper<TScalar> operator*<>(EigenMatrixWrapper<TScalar> const &leftMatrix, EigenVectorWrapper<TScalar> const &rightVector);

	friend TScalar vnl_trace<>(EigenMatrixWrapper<TScalar> const& M);
	friend EigenMatrixWrapper<TScalar> vnl_inverse<>(EigenMatrixWrapper<TScalar> const& M);
	friend EigenMatrixWrapper<TScalar> inverse_sympd<>(EigenMatrixWrapper<TScalar> const& M);

	friend std::ostream & operator<<<>(std::ostream& os, EigenMatrixWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename arma::Col<TScalar> EigenVectorType;
	/// Matrix type.
	typedef typename Eigen::Matrix<TScalar, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;
//	typedef typename Eigen::Matrix<TScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline EigenMatrixWrapper() : m_Matrix() {}

	explicit inline EigenMatrixWrapper(unsigned r, unsigned c, TScalar const& v0) { m_Matrix = EigenMatrixType::Constant(r, c, v0) ; }

	explicit inline EigenMatrixWrapper(unsigned r, unsigned c) : m_Matrix(r, c) {}

	inline EigenMatrixWrapper(const Eigen::Matrix<TScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & m) : m_Matrix(m) {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return number of rows
	inline unsigned rows()    const { return m_Matrix.rows(); }

	//: Return number of columns
	inline unsigned columns()  const { return m_Matrix.cols(); }

	//: Return number of columns
	inline unsigned cols()    const { return m_Matrix.cols(); }

	//: Return number of elements
	inline unsigned size()    const { return m_Matrix.innerSize(); }

	//: Access the contiguous block storing the elements in the matrix row-wise. O(1).
	inline TScalar const* data_block() const { exit(42); return 0; }

	//: Access the contiguous block storing the elements in the matrix row-wise. O(1).
	inline TScalar      * data_block() { exit(42); return 0; }

	//: Resize to r rows by c columns. Old data lost.
	// Returns true if size changed.
	inline bool set_size(unsigned r, unsigned c) { m_Matrix.resize(r, c); return true; }

	//: Set all elements of matrix to specified value.
	// Complexity $O(r.c)$
	inline EigenMatrixWrapper<TScalar>& fill(TScalar const& scalar) { m_Matrix.fill(scalar); return *this; }

	//: Get n rows beginning at rowstart
	inline EigenMatrixWrapper<TScalar> get_n_rows(unsigned rowstart, unsigned n) const {
		return EigenMatrixWrapper<TScalar>(m_Matrix.block(rowstart, 0, n, m_Matrix.cols()));
	}

	//: Get n columns beginning at colstart
	inline EigenMatrixWrapper<TScalar> get_n_columns(unsigned colstart, unsigned n) const {
		return EigenMatrixWrapper<TScalar>(m_Matrix.block(0, colstart, m_Matrix.rows(), n));
	}

	//: Set columns to those in M, starting at starting_column, then return *this.
	inline EigenMatrixWrapper<TScalar> & set_columns(unsigned starting_column, EigenMatrixWrapper<TScalar> const& M) {
		m_Matrix.block(0, starting_column, M.rows(), M.cols()) = M.m_Matrix;
		return *this;
	}

	//: Get a vector equal to the given row
	inline EigenVectorWrapper<TScalar> get_row(unsigned r) const { return EigenVectorWrapper<TScalar>(m_Matrix.row(r).transpose()); }

	//: Get a vector equal to the given column
	inline EigenVectorWrapper<TScalar> get_column(unsigned c) const { return EigenVectorWrapper<TScalar>(m_Matrix.col(c)); }

	//: Set the i-th row
	inline void set_row(unsigned i, EigenVectorWrapper<TScalar> const& v) { m_Matrix.row(i) = v.toEigen(); }

	//: Set the elements of the i'th row to value, then return *this.
	inline EigenMatrixWrapper<TScalar> & set_row(unsigned i, TScalar value ) { m_Matrix.row(i).fill(value); return *this; }

	//: Set j-th column to v
	inline void set_column(unsigned i, EigenVectorWrapper<TScalar> const& v) { m_Matrix.col(i) = v.toEigen(); }

	//: Set this matrix to an identity matrix
	//  Abort if the matrix is not square
	inline void set_identity() { m_Matrix.setIdentity(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return Frobenius norm of matrix (sqrt of sum of squares of its elements)
	inline TScalar frobenius_norm() const { return m_Matrix.norm(); }

	//: Return transpose
	inline EigenMatrixWrapper<TScalar> transpose() const { return EigenMatrixWrapper<TScalar>(m_Matrix.transpose()); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar const & operator()(long r, long c) const { return m_Matrix(r, c); }
	inline TScalar       & operator()(long r, long c)       { return m_Matrix(r, c); }

	inline EigenMatrixWrapper<TScalar>& operator/=(TScalar rhsScalar)
							{ m_Matrix /= rhsScalar; return *this; }
	inline EigenMatrixWrapper<TScalar>& operator*=(TScalar rhsScalar)
							{ m_Matrix *= rhsScalar; return *this; }
	inline EigenMatrixWrapper<TScalar>& operator+=(TScalar rhsScalar)
							{ m_Matrix += rhsScalar; return *this; }

	inline EigenMatrixWrapper<TScalar>& operator*=(EigenMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix *= rhsMatrix.m_Matrix; return *this; }
	inline EigenMatrixWrapper<TScalar>& operator+=(EigenMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix += rhsMatrix.m_Matrix; return *this; }
	inline EigenMatrixWrapper<TScalar>& operator-=(EigenMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix -= rhsMatrix.m_Matrix; return *this; }

	inline EigenMatrixWrapper<TScalar> operator-() const { return EigenMatrixWrapper<TScalar>(- m_Matrix); }



private :

	EigenMatrixType m_Matrix;

};



#include "EigenMatrixWrapper_friend.h"


#endif /* _EigenMatrixWrapper_h */

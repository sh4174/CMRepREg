/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _BlazeMatrixWrapper_h
#define _BlazeMatrixWrapper_h

#include <blaze/Math.h>
#include <blaze/math/DenseSubmatrix.h>
#include <blaze/math/DenseColumn.h>
#include <blaze/math/DenseRow.h>


#define Blaze_NO_DEBUG

template<class TScalar>
class BlazeMatrixWrapper;

template<class TScalar>
class BlazeVectorWrapper;



template<class TScalar>
BlazeMatrixWrapper<TScalar> operator-(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);
template<class TScalar>
BlazeMatrixWrapper<TScalar> operator+(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);
template<class TScalar>
BlazeMatrixWrapper<TScalar> operator*(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);

template<class TScalar>
BlazeMatrixWrapper<TScalar> operator/(const BlazeMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
template<class TScalar>
BlazeMatrixWrapper<TScalar> operator*(const BlazeMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
template<class TScalar>
BlazeMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const BlazeMatrixWrapper<TScalar> & rightMatrix);

template<class TScalar>
BlazeVectorWrapper<TScalar> operator*(BlazeMatrixWrapper<TScalar> const &leftMatrix, BlazeVectorWrapper<TScalar> const &rightVector);

template<class TScalar>
TScalar vnl_trace(BlazeMatrixWrapper<TScalar> const& M);
template<class TScalar>
BlazeMatrixWrapper<TScalar> vnl_inverse(BlazeMatrixWrapper<TScalar> const& M);
template<class TScalar>
BlazeMatrixWrapper<TScalar> inverse_sympd(BlazeMatrixWrapper<TScalar> const& M);

template <class TScalar>
std::ostream & operator<<(std::ostream& os, BlazeMatrixWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical matrix class (Blaze).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The BlazeMatrixWrapper class is a wrapping for mathematical matrix class, templated by type of element
 *              and using the Blaze library.
 */
template<class TScalar>
class BlazeMatrixWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend BlazeMatrixWrapper<TScalar> operator-<>(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);
	friend BlazeMatrixWrapper<TScalar> operator+<>(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);
	friend BlazeMatrixWrapper<TScalar> operator*<>(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right);

	friend BlazeMatrixWrapper<TScalar> operator/<>(const BlazeMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
	friend BlazeMatrixWrapper<TScalar> operator*<>(const BlazeMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
	friend BlazeMatrixWrapper<TScalar> operator*<>(TScalar const& leftScalar, const BlazeMatrixWrapper<TScalar> & rightMatrix);

	friend BlazeVectorWrapper<TScalar> operator*<>(BlazeMatrixWrapper<TScalar> const &leftMatrix, BlazeVectorWrapper<TScalar> const &rightVector);

	friend TScalar vnl_trace<>(BlazeMatrixWrapper<TScalar> const& M);
	friend BlazeMatrixWrapper<TScalar> vnl_inverse<>(BlazeMatrixWrapper<TScalar> const& M);
	friend BlazeMatrixWrapper<TScalar> inverse_sympd<>(BlazeMatrixWrapper<TScalar> const& M);

	friend std::ostream & operator<<<>(std::ostream& os, BlazeMatrixWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename blaze::DynamicVector<TScalar, blaze::columnVector> BlazeVectorType;
	/// Matrix type.
	typedef typename blaze::DynamicMatrix<TScalar, blaze::rowMajor> BlazeMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline BlazeMatrixWrapper() : m_Matrix() {}

	explicit inline BlazeMatrixWrapper(unsigned r, unsigned c, TScalar const& v0) : m_Matrix(r, c, v0) {}

	explicit inline BlazeMatrixWrapper(unsigned r, unsigned c) : m_Matrix(r, c) {}

	inline BlazeMatrixWrapper(const BlazeMatrixType & m) : m_Matrix(m) {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return number of rows
	inline unsigned rows()    const { return m_Matrix.rows(); }

	//: Return number of columns
	// A synonym for cols()
	inline unsigned columns()  const { return m_Matrix.columns(); }

	//: Return number of columns
	// A synonym for columns()
	inline unsigned cols()    const { return m_Matrix.columns(); }

	//: Return number of elements
	// This equals rows() * cols()
	inline unsigned size()    const { return (m_Matrix.rows()*m_Matrix.columns()); }

	//: Access the contiguous block storing the elements in the matrix row-wise. O(1).
	// 1d array, row-major order.
	inline TScalar const* data_block() const { exit(42); return 0; }

	//: Access the contiguous block storing the elements in the matrix row-wise. O(1).
	// 1d array, row-major order.
	inline TScalar      * data_block() { exit(42); return 0; }

	//: Resize to r rows by c columns. Old data lost.
	// Returns true if size changed.
	inline bool set_size(unsigned r, unsigned c) { m_Matrix.resize(r, c); return true; }

	//: Set all elements of matrix to specified value.
	// Complexity $O(r.c)$
	inline BlazeMatrixWrapper<TScalar>& fill(TScalar const& scalar) { m_Matrix = scalar; return *this; }

	//: Get n rows beginning at rowstart
	inline BlazeMatrixWrapper<TScalar> get_n_rows(unsigned rowstart, unsigned n) const {
		return BlazeMatrixWrapper<TScalar>(blaze::submatrix(m_Matrix, rowstart, 0, n, m_Matrix.columns()));
	}

	//: Get n columns beginning at colstart
	inline BlazeMatrixWrapper<TScalar> get_n_columns(unsigned colstart, unsigned n) const {
		return BlazeMatrixWrapper<TScalar>(blaze::submatrix(m_Matrix, 0, colstart, m_Matrix.rows(), n));
	}

	//: Set columns to those in M, starting at starting_column, then return *this.
	inline BlazeMatrixWrapper<TScalar> & set_columns(unsigned starting_column, BlazeMatrixWrapper<TScalar> const& M) {
		blaze::DenseSubmatrix<BlazeMatrixType>  submatrix = blaze::submatrix(m_Matrix, 0, starting_column, m_Matrix.rows(), M.columns());
		m_Matrix = M.m_Matrix;
		return *this;
	}

	//: Get a vector equal to the given row
	inline BlazeVectorWrapper<TScalar> get_row(unsigned r) const { return BlazeVectorWrapper<TScalar>(blaze::trans(blaze::row(m_Matrix, r))); }

	//: Get a vector equal to the given column
	inline BlazeVectorWrapper<TScalar> get_column(unsigned c) const { return BlazeVectorWrapper<TScalar>(blaze::column(m_Matrix, c)); }

	//: Set the i-th row
	inline void set_row(unsigned i, BlazeVectorWrapper<TScalar> const& v) { blaze::row(m_Matrix, i) = blaze::trans(v.toBlaze()); }

	//: Set the elements of the i'th row to value, then return *this.
	inline BlazeMatrixWrapper<TScalar> & set_row(unsigned i, TScalar value ) { blaze::row(m_Matrix, i) = value; return *this; }

	//: Set j-th column to v
	inline void set_column(unsigned i, BlazeVectorWrapper<TScalar> const& v) { blaze::column(m_Matrix, i) = v.toBlaze(); }

	//: Set this matrix to an identity matrix
	//  Abort if the matrix is not square
	inline void set_identity() {
		m_Matrix = 0.0;
		for( size_t i=0; i<m_Matrix.rows(); i++ )
			m_Matrix(i,i) = 1.0;
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return Frobenius norm of matrix (sqrt of sum of squares of its elements)
	inline TScalar frobenius_norm() const {
		TScalar result = 0.0;
//		for( size_t i=0UL; i<m_Matrix.rows(); ++i )
//		   for( blaze::CompressedMatrix<int,rowMajor>::Iterator it=m_Matrix.begin(i); it!=m_Matrix.end(i); ++it )
//			  result += it->value()*it->value();
		for( size_t j=0; j<m_Matrix.columns(); j++ )
			for( size_t i=0; i<m_Matrix.rows(); i++ )
			  result += m_Matrix(i,j)*m_Matrix(i,j);
		return result;
	}

	//: Return transpose
	inline BlazeMatrixWrapper<TScalar> transpose() const { return BlazeMatrixWrapper<TScalar>(blaze::trans(m_Matrix)); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar const & operator()(long r, long c) const { return m_Matrix(r, c); }
	inline TScalar       & operator()(long r, long c)       { return m_Matrix(r, c); }

	inline BlazeMatrixWrapper<TScalar>& operator/=(TScalar rhsScalar)
							{ m_Matrix /= rhsScalar; return *this; }
	inline BlazeMatrixWrapper<TScalar>& operator*=(TScalar rhsScalar)
							{ m_Matrix *= rhsScalar; return *this; }
	inline BlazeMatrixWrapper<TScalar>& operator+=(TScalar rhsScalar)
							{ m_Matrix += rhsScalar; return *this; }

	inline BlazeMatrixWrapper<TScalar>& operator*=(BlazeMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix *= rhsMatrix.m_Matrix; return *this; }
	inline BlazeMatrixWrapper<TScalar>& operator+=(BlazeMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix += rhsMatrix.m_Matrix; return *this; }
	inline BlazeMatrixWrapper<TScalar>& operator-=(BlazeMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix -= rhsMatrix.m_Matrix; return *this; }

	inline BlazeMatrixWrapper<TScalar> operator-() const { return BlazeMatrixWrapper<TScalar>(- m_Matrix); }



private :

	BlazeMatrixType m_Matrix;

};
//:  public vnl_vector<TScalar> {}; /* class BlazeVectorWrapper */


#include "BlazeMatrixWrapper_friend.h"


#endif /* _BlazeMatrixWrapper_h */

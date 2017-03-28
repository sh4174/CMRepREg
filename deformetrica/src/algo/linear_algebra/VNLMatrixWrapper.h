/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _VNLMatrixWrapper_h
#define _VNLMatrixWrapper_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_cross.h"



template<class TScalar>
class VNLMatrixWrapper;

template<class TScalar>
class VNLVectorWrapper;



template<class TScalar>
VNLMatrixWrapper<TScalar> operator-(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);
template<class TScalar>
VNLMatrixWrapper<TScalar> operator+(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);
template<class TScalar>
VNLMatrixWrapper<TScalar> operator*(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);

template<class TScalar>
VNLMatrixWrapper<TScalar> operator/(const VNLMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
template<class TScalar>
VNLMatrixWrapper<TScalar> operator*(const VNLMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
template<class TScalar>
VNLMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const VNLMatrixWrapper<TScalar> & rightMatrix);

template<class TScalar>
VNLVectorWrapper<TScalar> operator*(VNLMatrixWrapper<TScalar> const &leftMatrix, VNLVectorWrapper<TScalar> const &rightVector);

template<class TScalar>
TScalar vnl_trace(VNLMatrixWrapper<TScalar> const& M);
template<class TScalar>
VNLMatrixWrapper<TScalar> vnl_inverse(VNLMatrixWrapper<TScalar> const& M);
template<class TScalar>
VNLMatrixWrapper<TScalar> inverse_sympd(VNLMatrixWrapper<TScalar> const& M);

template <class TScalar>
std::ostream & operator<<(std::ostream& os, VNLMatrixWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical matrix class (VNL).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The VNLMatrixWrapper class is a wrapping for mathematical matrix class, templated by type of element
 *              and using the VNL library.
 */
template<class TScalar>
class VNLMatrixWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend VNLMatrixWrapper<TScalar> operator-<>(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);
	friend VNLMatrixWrapper<TScalar> operator+<>(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);
	friend VNLMatrixWrapper<TScalar> operator*<>(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right);

	friend VNLMatrixWrapper<TScalar> operator/<>(const VNLMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar);
	friend VNLMatrixWrapper<TScalar> operator*<>(const VNLMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar);
	friend VNLMatrixWrapper<TScalar> operator*<>(TScalar const& leftScalar, const VNLMatrixWrapper<TScalar> & rightMatrix);

	friend VNLVectorWrapper<TScalar> operator*<>(VNLMatrixWrapper<TScalar> const &leftMatrix, VNLVectorWrapper<TScalar> const &rightVector);

	friend TScalar vnl_trace<>(VNLMatrixWrapper<TScalar> const& M);
	friend VNLMatrixWrapper<TScalar> vnl_inverse<>(VNLMatrixWrapper<TScalar> const& M);
	friend VNLMatrixWrapper<TScalar> inverse_sympd<>(VNLMatrixWrapper<TScalar> const& M);

	friend std::ostream & operator<<<>(std::ostream& os, VNLMatrixWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline VNLMatrixWrapper() : m_Matrix() {}

	explicit inline VNLMatrixWrapper(unsigned r, unsigned c, TScalar const& v0) : m_Matrix(r, c, v0) {}

	explicit inline VNLMatrixWrapper(unsigned r, unsigned c) : m_Matrix(r, c) {}

	inline VNLMatrixWrapper(const vnl_matrix<TScalar> & m) : m_Matrix(m) {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Return number of rows
	inline unsigned rows()    const { return m_Matrix.rows(); }

	/// Return number of columns
	// A synonym for cols()
	inline unsigned columns()  const { return m_Matrix.columns(); }

	/// Return number of columns
	// A synonym for columns()
	inline unsigned cols()    const { return m_Matrix.cols(); }

	/// Return number of elements
	// This equals rows() * cols()
	inline unsigned size()    const { return m_Matrix.size(); }

	/// Access the contiguous block storing the elements in the matrix row-wise. O(1).
	// 1d array, row-major order.
	TScalar const* data_block() const { return m_Matrix.data_block(); }

	/// Access the contiguous block storing the elements in the matrix row-wise. O(1).
	// 1d array, row-major order.
	TScalar      * data_block() { return m_Matrix.data_block(); }

	/// Resize to r rows by c columns. Old data lost.
	// Returns true if size changed.
	bool set_size(unsigned r, unsigned c) { return m_Matrix.set_size(r, c); }

	/// Set all elements of matrix to specified value.
	// Complexity $O(r.c)$
	VNLMatrixWrapper<TScalar>& fill(TScalar const& scalar) { VNLMatrixWrapper<TScalar>(m_Matrix.fill(scalar)); return *this; }

	/// Get n rows beginning at rowstart
	VNLMatrixWrapper<TScalar> get_n_rows(unsigned rowstart, unsigned n) const {
		return VNLMatrixWrapper<TScalar>(m_Matrix.get_n_rows(rowstart, n));
	}

	/// Get n columns beginning at colstart
	VNLMatrixWrapper<TScalar> get_n_columns(unsigned colstart, unsigned n) const {
		return VNLMatrixWrapper<TScalar>(m_Matrix.get_n_columns(colstart, n));
	}

	/// Set columns to those in M, starting at starting_column, then return *this.
	VNLMatrixWrapper<TScalar> & set_columns(unsigned starting_column, VNLMatrixWrapper<TScalar> const& M) {
		m_Matrix.set_columns(starting_column, M.m_Matrix);
		return *this;
	}

	/// Get a vector equal to the given row
	inline VNLVectorWrapper<TScalar> get_row(unsigned r) const { return VNLVectorWrapper<TScalar>(m_Matrix.get_row(r)); }

	/// Get a vector equal to the given column
	inline VNLVectorWrapper<TScalar> get_column(unsigned c) const { return VNLVectorWrapper<TScalar>(m_Matrix.get_column(c)); }

	/// Set the i-th row
	inline void set_row(unsigned i, VNLVectorWrapper<TScalar> const& v) { m_Matrix.set_row(i, v.toVNL()); }

	/// Set the elements of the i'th row to value, then return *this.
	inline VNLMatrixWrapper<TScalar> & set_row(unsigned i, TScalar value ) { m_Matrix.set_row(i, value); return *this; }

	/// Set j-th column to v
	inline void set_column(unsigned i, VNLVectorWrapper<TScalar> const& v) { m_Matrix.set_column(i, v.toVNL()); }

	/// Set this matrix to an identity matrix
	//  Abort if the matrix is not square
	inline void set_identity() { m_Matrix.set_identity(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Return Frobenius norm of matrix (sqrt of sum of squares of its elements)
	inline TScalar frobenius_norm() const { return m_Matrix.frobenius_norm(); }

	/// Return transpose
	inline VNLMatrixWrapper<TScalar> transpose() const { return VNLMatrixWrapper<TScalar>(m_Matrix.transpose()); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar const & operator()(unsigned r, unsigned c) const { return m_Matrix(r, c); }
	inline TScalar       & operator()(unsigned r, unsigned c)       { return m_Matrix(r, c); }

	inline VNLMatrixWrapper<TScalar>& operator/=(TScalar rhsScalar)
							{ m_Matrix /= rhsScalar; return *this; }
	inline VNLMatrixWrapper<TScalar>& operator*=(TScalar rhsScalar)
							{ m_Matrix *= rhsScalar; return *this; }
	inline VNLMatrixWrapper<TScalar>& operator+=(TScalar rhsScalar)
							{ m_Matrix += rhsScalar; return *this; }

	inline VNLMatrixWrapper<TScalar>& operator*=(VNLMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix *= rhsMatrix.m_Matrix; return *this; }
	inline VNLMatrixWrapper<TScalar>& operator+=(VNLMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix += rhsMatrix.m_Matrix; return *this; }
	inline VNLMatrixWrapper<TScalar>& operator-=(VNLMatrixWrapper<TScalar> const& rhsMatrix)
							{ m_Matrix -= rhsMatrix.m_Matrix; return *this; }

	inline VNLMatrixWrapper<TScalar> operator-() const { return VNLMatrixWrapper<TScalar>(- m_Matrix); }



private :

	vnl_matrix<TScalar> m_Matrix;

};
///  public vnl_vector<TScalar> {}; /* class VNLVectorWrapper */


#include "VNLMatrixWrapper_friend.h"


#endif /* _VNLMatrixWrapper_h */

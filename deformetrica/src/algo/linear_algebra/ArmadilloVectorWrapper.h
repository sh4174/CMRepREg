/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _ArmadilloVectorWrapper_h
#define _ArmadilloVectorWrapper_h

#include <armadillo>

#include <cstdlib>



template<class TScalar>
class ArmadilloMatrixWrapper;

template<class TScalar>
class ArmadilloVectorWrapper;



template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator*(ArmadilloMatrixWrapper<TScalar> const &leftMatrix, ArmadilloVectorWrapper<TScalar> const &rightVector);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator-(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator+(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right);

template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator+(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator-(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator/(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const ArmadilloVectorWrapper<TScalar> & rightVector);

template<class TScalar>
std::ostream & operator<<(std::ostream& os, ArmadilloVectorWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical vector class (Armadillo).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The ArmadilloVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the Armadillo library.
 */
template<class TScalar>
class ArmadilloVectorWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend ArmadilloVectorWrapper<TScalar> operator*<>(ArmadilloMatrixWrapper<TScalar> const &leftMatrix, ArmadilloVectorWrapper<TScalar> const &rightVector);

	friend ArmadilloVectorWrapper<TScalar> operator-<>(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right);
	friend ArmadilloVectorWrapper<TScalar> operator+<>(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right);

	friend ArmadilloVectorWrapper<TScalar> operator+<>(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend ArmadilloVectorWrapper<TScalar> operator-<>(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend ArmadilloVectorWrapper<TScalar> operator/<>(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

	friend ArmadilloVectorWrapper<TScalar> operator*<>(const TScalar & leftScalar, const ArmadilloVectorWrapper<TScalar> & rightVector);

	/// Overload of the stream operator.
	friend std::ostream & operator<<<>(std::ostream& os, ArmadilloVectorWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename arma::Col<TScalar> ArmadilloVectorType;
	/// Matrix type.
	typedef typename arma::Mat<TScalar> ArmadilloMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Default constructor.
	inline ArmadilloVectorWrapper() : m_Vector() {}

 	/// Creates vector of \e len elements.
 	explicit inline ArmadilloVectorWrapper(unsigned len) : m_Vector(len) { }

 	/// Creates vector of \e len elements, all set to \e v0.
 	explicit inline ArmadilloVectorWrapper(unsigned len, TScalar const& v0) : m_Vector(len) { m_Vector.fill(v0); }

	/// Special constructor (do not use it in Deformetrica!).
	inline ArmadilloVectorWrapper(const ArmadilloVectorType & v) : m_Vector(v) {}

	/// Access directly to the Armadillo vector (do not use it in Deformetrica!).
	inline ArmadilloVectorType const & toArmadillo() const { return m_Vector; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the length, number of elements, dimension of this vector.
	inline unsigned size() const { return m_Vector.n_rows; }

	/// Resizes to \e n elements.
	inline bool set_size(unsigned n) { m_Vector.set_size(n); return true; }

	/// Gets value at element \e i.
	inline TScalar get(unsigned int i) const {
		return m_Vector.at(i);
	}

	/// Sets all elements of matrix to specified value.
	inline ArmadilloVectorWrapper<TScalar>& fill(TScalar const& scalar) { m_Vector.fill(scalar); return *this; }

	/// Converts a (column) vector to a matrix of size \e rows x \e cols.
	inline ArmadilloMatrixWrapper<TScalar> convert_to_matrix_row_wise(unsigned rows, unsigned cols) const {
		ArmadilloMatrixType result(m_Vector);
		result.reshape(cols, rows);
		return ArmadilloMatrixWrapper<TScalar>(arma::trans(result));
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns magnitude (length) of vector.
	inline TScalar magnitude() const { return arma::norm(m_Vector, 2); }

	/// Returns sum of squares of elements.
	inline TScalar squared_magnitude() const { TScalar result = arma::norm(m_Vector, 2); return result*result; }

	/// Returns sum of elements.
	inline TScalar sum() const { return arma::sum(m_Vector); }

	/// Smallest value.
	inline TScalar min_value() const { return m_Vector.min(); }

	/// Largest value.
	inline TScalar max_value() const { return m_Vector.max(); }

	/// Access the contiguous block storing the elements in the vector.
	inline TScalar const* memptr() const { return m_Vector.memptr(); }

	/// Access the contiguous block storing the elements in the vector.
	inline TScalar      * memptr() { return m_Vector.memptr(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns reference to the element at index \e i.
	inline TScalar       & operator()(unsigned int i) { return m_Vector(i); }
	/// Returns reference to the element at index \e i.
	inline TScalar const & operator()(unsigned int i) const { return m_Vector(i); }

	/// Returns reference to the element at index \e i.
	/// \warning Please use operator() in Deformetrica instead of operator[].
	inline TScalar       & operator[](unsigned int i) { return m_Vector(i); }
	/// Returns reference to the element at index \e i.
	/// \warning Please use operator() in Deformetrica instead of operator[].
	inline TScalar const & operator[](unsigned int i) const { return m_Vector(i); }

	/// Scalar multiplication of lhs vector in situ by \e rhsScalar.
	inline ArmadilloVectorWrapper<TScalar>& operator*=(TScalar rhsScalar) { m_Vector *= rhsScalar; return *this; }
	/// Scalar division of lhs vector in situ by \e rhsScalar.
	inline ArmadilloVectorWrapper<TScalar>& operator/=(TScalar rhsScalar) { m_Vector /= rhsScalar; return *this; }

	/// Adds \e rhs to lhs vector in situ.
	inline ArmadilloVectorWrapper<TScalar>& operator+=(ArmadilloVectorWrapper<TScalar> const& rhs) { m_Vector += rhs.m_Vector; return *this; }

	/// Unary minus operator.
	inline ArmadilloVectorWrapper<TScalar> operator-() const { return ArmadilloVectorWrapper<TScalar>(- m_Vector); }
	/// Unary plus operator.
	inline ArmadilloVectorWrapper<TScalar> operator+() const { return ArmadilloVectorWrapper<TScalar>(+ m_Vector); }

	/// Scalar multiplication of lhs vector by \e scalar.
	inline ArmadilloVectorWrapper<TScalar> operator*(TScalar scalar) const { return ArmadilloVectorWrapper<TScalar>(m_Vector*scalar); }
	/// Scalar division of lhs vector by \e scalar.
	inline ArmadilloVectorWrapper<TScalar> operator/(TScalar scalar) const { return ArmadilloVectorWrapper<TScalar>(m_Vector/scalar); }



private :

	ArmadilloVectorType m_Vector;

};



template<class TScalar>
inline TScalar dot_product  (ArmadilloVectorWrapper<TScalar> const & v1, ArmadilloVectorWrapper<TScalar> const & v2) {
	return arma::dot(v1.toArmadillo(), v2.toArmadillo());
}

template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> cross_3d  (ArmadilloVectorWrapper<TScalar> const & v1, ArmadilloVectorWrapper<TScalar> const & v2) {
	return ArmadilloVectorWrapper<TScalar>(arma::cross(v1.toArmadillo(), v2.toArmadillo()));
}


#include "ArmadilloVectorWrapper_friend.h"



#endif /* _ArmadilloVectorWrapper_h */

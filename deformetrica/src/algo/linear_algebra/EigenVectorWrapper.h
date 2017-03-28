/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _EigenVectorWrapper_h
#define _EigenVectorWrapper_h

#include <Eigen/Eigen>
#include <Eigen/Core>

#include <cstdlib>

#define EIGEN_NO_DEBUG


template<class TScalar>
class EigenMatrixWrapper;

template<class TScalar>
class EigenVectorWrapper;



template<class TScalar>
EigenVectorWrapper<TScalar> operator*(EigenMatrixWrapper<TScalar> const &leftMatrix, EigenVectorWrapper<TScalar> const &rightVector);
template<class TScalar>
EigenVectorWrapper<TScalar> operator-(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right);
template<class TScalar>
EigenVectorWrapper<TScalar> operator+(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right);

template<class TScalar>
EigenVectorWrapper<TScalar> operator+(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
EigenVectorWrapper<TScalar> operator-(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
EigenVectorWrapper<TScalar> operator/(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

template<class TScalar>
EigenVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const EigenVectorWrapper<TScalar> & rightVector);

template<class TScalar>
std::ostream & operator<<(std::ostream& os, EigenVectorWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical vector class (Eigen).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The EigenVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the Eigen library.
 */
template<class TScalar>
class EigenVectorWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend EigenVectorWrapper<TScalar> operator*<>(EigenMatrixWrapper<TScalar> const &leftMatrix, EigenVectorWrapper<TScalar> const &rightVector);

	friend EigenVectorWrapper<TScalar> operator-<>(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right);
	friend EigenVectorWrapper<TScalar> operator+<>(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right);

	friend EigenVectorWrapper<TScalar> operator+<>(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend EigenVectorWrapper<TScalar> operator-<>(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend EigenVectorWrapper<TScalar> operator/<>(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

	friend EigenVectorWrapper<TScalar> operator*<>(const TScalar & leftScalar, const EigenVectorWrapper<TScalar> & rightVector);

	friend std::ostream & operator<<<>(std::ostream& os, EigenVectorWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline EigenVectorWrapper() : m_Vector() {}

	explicit inline EigenVectorWrapper(unsigned len, TScalar const& v0) : m_Vector(len) { m_Vector.fill(v0); }

	explicit inline EigenVectorWrapper(unsigned len) : m_Vector(len) {}

	inline EigenVectorWrapper(const Eigen::Matrix<TScalar, Eigen::Dynamic, 1> & v) : m_Vector(v) {}

	inline Eigen::Matrix<TScalar, Eigen::Dynamic, 1> const & toEigen() const { return m_Vector; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return the length, number of elements, dimension of this vector.
	inline unsigned size() const { return m_Vector.innerSize(); }

	//: Resize to n elements.
	// This is a destructive resize, in that the old data is lost if size() != \a n before the call.
	// If size() is already \a n, this is a null operation.
	inline bool set_size(unsigned n) { m_Vector.resize(n); return true; }

	//: Get value at element i
	inline TScalar get(unsigned int i) const {
		return m_Vector(i);
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return magnitude (length) of vector
	inline TScalar magnitude() const { return m_Vector.norm(); }

	//: Return sum of squares of elements
	inline TScalar squared_magnitude() const { return m_Vector.squaredNorm(); }

	//: Smallest value
	inline TScalar min_value() const { return m_Vector.minCoeff(); }

	//: Largest value
	inline TScalar max_value() const { return m_Vector.maxCoeff(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar       & operator()(unsigned int i) { return m_Vector(i); }
	inline TScalar const & operator()(unsigned int i) const { return m_Vector(i); }

	inline TScalar       & operator[](unsigned int i) { return m_Vector[i]; }
	inline TScalar const & operator[](unsigned int i) const { return m_Vector[i]; }

	inline EigenVectorWrapper<TScalar>& operator*=(TScalar rhsScalar) { m_Vector *= rhsScalar; return *this; }
	inline EigenVectorWrapper<TScalar>& operator/=(TScalar rhsScalar) { m_Vector /= rhsScalar; return *this; }

	inline EigenVectorWrapper<TScalar>& operator+=(EigenVectorWrapper<TScalar> const& rhs) {
		m_Vector += rhs.m_Vector;
		return *this;
	}

	inline EigenVectorWrapper<TScalar> operator*(TScalar scalar) const { return EigenVectorWrapper<TScalar>(m_Vector*scalar); }



private :

	Eigen::Matrix<TScalar, Eigen::Dynamic, 1> m_Vector;

};
//:  public vnl_vector<TScalar> {}; /* class EigenVectorWrapper */



template<class TScalar>
inline TScalar dot_product  (EigenVectorWrapper<TScalar> const & v1, EigenVectorWrapper<TScalar> const & v2) {
	return v1.toEigen().dot(v2.toEigen());
}

template<class TScalar>
inline EigenVectorWrapper<TScalar> vnl_cross_3d  (EigenVectorWrapper<TScalar> const & v1, EigenVectorWrapper<TScalar> const & v2) {
	EigenVectorWrapper<TScalar> result(3);
	result(0) = v1(1) * v2(2) - v1(2) * v2(1);
	result(1) = v1(2) * v2(0) - v1(0) * v2(2);
	result(2) = v1(0) * v2(1) - v1(1) * v2(0);
	return result;
}


#include "EigenVectorWrapper_friend.h"



#endif /* _EigenVectorWrapper_h */

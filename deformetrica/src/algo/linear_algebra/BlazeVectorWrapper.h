/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _BlazeVectorWrapper_h
#define _BlazeVectorWrapper_h

#include <blaze/Math.h>

#include <cstdlib>

#define Blaze_NO_DEBUG


template<class TScalar>
class BlazeMatrixWrapper;

template<class TScalar>
class BlazeVectorWrapper;



template<class TScalar>
BlazeVectorWrapper<TScalar> operator*(BlazeMatrixWrapper<TScalar> const &leftMatrix, BlazeVectorWrapper<TScalar> const &rightVector);
template<class TScalar>
BlazeVectorWrapper<TScalar> operator-(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right);
template<class TScalar>
BlazeVectorWrapper<TScalar> operator+(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right);

template<class TScalar>
BlazeVectorWrapper<TScalar> operator+(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
BlazeVectorWrapper<TScalar> operator-(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
BlazeVectorWrapper<TScalar> operator/(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

template<class TScalar>
BlazeVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const BlazeVectorWrapper<TScalar> & rightVector);

template<class TScalar>
std::ostream & operator<<(std::ostream& os, BlazeVectorWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical vector class (Blaze).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The BlazeVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the Blaze library.
 */
template<class TScalar>
class BlazeVectorWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend BlazeVectorWrapper<TScalar> operator*<>(BlazeMatrixWrapper<TScalar> const &leftMatrix, BlazeVectorWrapper<TScalar> const &rightVector);

	friend BlazeVectorWrapper<TScalar> operator-<>(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right);
	friend BlazeVectorWrapper<TScalar> operator+<>(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right);

	friend BlazeVectorWrapper<TScalar> operator+<>(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend BlazeVectorWrapper<TScalar> operator-<>(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend BlazeVectorWrapper<TScalar> operator/<>(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

	friend BlazeVectorWrapper<TScalar> operator*<>(const TScalar & leftScalar, const BlazeVectorWrapper<TScalar> & rightVector);

	friend std::ostream & operator<<<>(std::ostream& os, BlazeVectorWrapper<TScalar> const& rhs);



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

	inline BlazeVectorWrapper() : m_Vector() {}

	explicit inline BlazeVectorWrapper(unsigned len, TScalar const& v0) : m_Vector(len, v0) {}

	explicit inline BlazeVectorWrapper(unsigned len) : m_Vector(len) {}

	inline BlazeVectorWrapper(const BlazeVectorType & v) : m_Vector(v) {}

	inline BlazeVectorType const & toBlaze() const { return m_Vector; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return the length, number of elements, dimension of this vector.
	inline unsigned size() const { return m_Vector.size(); }

	//: Resize to n elements.
	// This is a destructive resize, in that the old data is lost if size() != \a n before the call.
	// If size() is already \a n, this is a null operation.
	inline bool set_size(unsigned n) { m_Vector.resize(n); return true; }

	//: Get value at element i
	inline TScalar get(unsigned int i) const {
		return m_Vector[i];
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return magnitude (length) of vector
	inline TScalar magnitude() const {
		// rows = columns ?
		TScalar result = 0.0;
		for( size_t i=0; i<m_Vector.size(); i++ )
			result += m_Vector[i]*m_Vector[i];
		return sqrt(result);
	}

	//: Return sum of squares of elements
	inline TScalar squared_magnitude() const {
		TScalar result = 0.0;
			for( size_t i=0; i<m_Vector.size(); i++ )
			  result += m_Vector[i]*m_Vector[i];
		return result;
	}

	//: Smallest value
	inline TScalar min_value() const { return blaze::min(m_Vector); }

	//: Largest value
	inline TScalar max_value() const { return blaze::max(m_Vector); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar       & operator()(unsigned int i) { return m_Vector[i]; }
	inline TScalar const & operator()(unsigned int i) const { return m_Vector[i]; }

	inline TScalar       & operator[](unsigned int i) { return m_Vector[i]; }
	inline TScalar const & operator[](unsigned int i) const { return m_Vector[i]; }

	inline BlazeVectorWrapper<TScalar>& operator*=(TScalar rhsScalar) { m_Vector *= rhsScalar; return *this; }
	inline BlazeVectorWrapper<TScalar>& operator/=(TScalar rhsScalar) { m_Vector /= rhsScalar; return *this; }

	inline BlazeVectorWrapper<TScalar>& operator+=(BlazeVectorWrapper<TScalar> const& rhs) {
		m_Vector += rhs.m_Vector;
		return *this;
	}

	inline BlazeVectorWrapper<TScalar> operator*(TScalar scalar) const { return BlazeVectorWrapper<TScalar>(m_Vector*scalar); }



private :

	BlazeVectorType m_Vector;

};
//:  public vnl_vector<TScalar> {}; /* class BlazeVectorWrapper */



template<class TScalar>
inline TScalar dot_product  (BlazeVectorWrapper<TScalar> const & v1, BlazeVectorWrapper<TScalar> const & v2) {
	return (v1.toBlaze(), v2.toBlaze());
}

template<class TScalar>
inline BlazeVectorWrapper<TScalar> vnl_cross_3d  (BlazeVectorWrapper<TScalar> const & v1, BlazeVectorWrapper<TScalar> const & v2) {
/*	BlazeVectorWrapper<TScalar> result(3);
	result(0) = v1(1) * v2(2) - v1(2) * v2(1);
	result(1) = v1(2) * v2(0) - v1(0) * v2(2);
	result(2) = v1(0) * v2(1) - v1(1) * v2(0);
	return result;
*/
	return BlazeVectorWrapper<TScalar>(v1.toBlaze() % v2.toBlaze());
}


#include "BlazeVectorWrapper_friend.h"



#endif /* _BlazeVectorWrapper_h */

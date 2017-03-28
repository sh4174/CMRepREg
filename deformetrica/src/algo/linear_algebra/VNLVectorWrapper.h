/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _VNLVectorWrapper_h
#define _VNLVectorWrapper_h

#include "vnl/vnl_vector.h"
#include "vnl/vnl_cross.h"

#include <cstdlib>



template<class TScalar>
class VNLMatrixWrapper;

template<class TScalar>
class VNLVectorWrapper;



template<class TScalar>
VNLVectorWrapper<TScalar> operator*(VNLMatrixWrapper<TScalar> const &leftMatrix, VNLVectorWrapper<TScalar> const &rightVector);
template<class TScalar>
VNLVectorWrapper<TScalar> operator-(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right);
template<class TScalar>
VNLVectorWrapper<TScalar> operator+(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right);

template<class TScalar>
VNLVectorWrapper<TScalar> operator+(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
VNLVectorWrapper<TScalar> operator-(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
VNLVectorWrapper<TScalar> operator/(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
template<class TScalar>
VNLVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const VNLVectorWrapper<TScalar> & rightVector);

template<class TScalar>
std::ostream & operator<<(std::ostream& os, VNLVectorWrapper<TScalar> const& rhs);



/**
 *  \brief      Mathematical vector class (VNL).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The VNLVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the VNL library.
 */
template<class TScalar>
class VNLVectorWrapper {

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Friend methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	friend VNLVectorWrapper<TScalar> operator*<>(VNLMatrixWrapper<TScalar> const &leftMatrix, VNLVectorWrapper<TScalar> const &rightVector);

	friend VNLVectorWrapper<TScalar> operator-<>(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right);
	friend VNLVectorWrapper<TScalar> operator+<>(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right);

	friend VNLVectorWrapper<TScalar> operator+<>(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend VNLVectorWrapper<TScalar> operator-<>(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);
	friend VNLVectorWrapper<TScalar> operator/<>(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar);

	friend VNLVectorWrapper<TScalar> operator*<>(const TScalar & leftScalar, const VNLVectorWrapper<TScalar> & rightVector);

	friend std::ostream & operator<<<>(std::ostream& os, VNLVectorWrapper<TScalar> const& rhs);



public :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline VNLVectorWrapper() : m_Vector() {}

	explicit inline VNLVectorWrapper(unsigned len, TScalar const& v0) : m_Vector(len, v0) {}

	explicit inline VNLVectorWrapper(unsigned len) : m_Vector(len) {}

	/// Special constructor (do not use it in Deformetrica!).
	inline VNLVectorWrapper(const vnl_vector<TScalar> & v) : m_Vector(v) {}

	/// Access directly to the VNL vector (do not use it in Deformetrica!).
	inline vnl_vector<TScalar> const & toVNL() const { return m_Vector; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Methods :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return the length, number of elements, dimension of this vector.
	inline unsigned size() const { return m_Vector.size(); }

	//: Resize to n elements.
	// This is a destructive resize, in that the old data is lost if size() != \a n before the call.
	// If size() is already \a n, this is a null operation.
	inline bool set_size(unsigned n) { return m_Vector.set_size(n); }

	//: Get value at element i
	inline TScalar get(unsigned int i) const {
		return m_Vector.get(i);
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operations :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//: Return magnitude (length) of vector
	inline TScalar magnitude() const { return m_Vector.magnitude(); }

	//: Return sum of squares of elements
	inline TScalar squared_magnitude() const { return m_Vector.squared_magnitude(); }

	//: Smallest value
	TScalar min_value() const { return m_Vector.min_value(); }

	//: Largest value
	TScalar max_value() const { return m_Vector.max_value(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Operator overloading :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar       & operator()(unsigned int i) { return m_Vector(i); }
	inline TScalar const & operator()(unsigned int i) const { return m_Vector(i); }

	inline TScalar       & operator[](unsigned int i) { return m_Vector[i]; }
	inline TScalar const & operator[](unsigned int i) const { return m_Vector[i]; }

	inline VNLVectorWrapper<TScalar>& operator*=(TScalar rhsScalar) { m_Vector *= rhsScalar; return *this; }
	inline VNLVectorWrapper<TScalar>& operator/=(TScalar rhsScalar) { m_Vector /= rhsScalar; return *this; }

	inline VNLVectorWrapper<TScalar>& operator+=(VNLVectorWrapper<TScalar> const& rhs) {
		m_Vector += rhs.m_Vector;
		return *this;
	}

	inline VNLVectorWrapper<TScalar> operator*(TScalar scalar) const { return VNLVectorWrapper<TScalar>(m_Vector*scalar); }



private :

	vnl_vector<TScalar> m_Vector;

};
//:  public vnl_vector<TScalar> {}; /* class VNLVectorWrapper */



template<class TScalar>
TScalar dot_product  (VNLVectorWrapper<TScalar> const & v1, VNLVectorWrapper<TScalar> const & v2) {
	TScalar result = 0.0;
	for (size_t i=0; i<v1.size(); i++)
		result += v1[i]*v2[i];
	return result;
}

template<class TScalar>
inline VNLVectorWrapper<TScalar> vnl_cross_3d  (VNLVectorWrapper<TScalar> const & v1, VNLVectorWrapper<TScalar> const & v2) {
	VNLVectorWrapper<TScalar> result(3);
	result(0) = v1(1) * v2(2) - v1(2) * v2(1);
	result(1) = v1(2) * v2(0) - v1(0) * v2(2);
	result(2) = v1(0) * v2(1) - v1(1) * v2(0);
	return result;
}


#include "VNLVectorWrapper_friend.h"



#endif /* _VNLVectorWrapper_h */

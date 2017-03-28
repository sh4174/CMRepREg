/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _BlazeMatrixWrapper_friend_h
#define _BlazeMatrixWrapper_friend_h

#ifndef _BlazeMatrixWrapper_h
#error Do not include BlazeMatrixWrapper_friend.h : include BlazeMatrixWrapper.h instead
#endif

#include "BlazeVectorWrapper.h"
#include <iostream>



//
// Left Matrix / Right Matrix :
//
template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator-(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right) {
	return BlazeMatrixWrapper<TScalar>(left.m_Matrix - right.m_Matrix);
}

template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator+(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right) {
	return BlazeMatrixWrapper<TScalar>(left.m_Matrix + right.m_Matrix);
}

template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator*(const BlazeMatrixWrapper<TScalar> & left, const BlazeMatrixWrapper<TScalar> & right) {
	return BlazeMatrixWrapper<TScalar>(left.m_Matrix * right.m_Matrix);
}



//
// (Left-Right) Matrix / Scalar :
//

template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator/(const BlazeMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar) {
	return BlazeMatrixWrapper<TScalar>(leftMatrix.m_Matrix / rightScalar);
}


template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator*(const BlazeMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar) {
	return BlazeMatrixWrapper<TScalar>(leftMatrix.m_Matrix * rightScalar);
}

template<class TScalar>
inline BlazeMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const BlazeMatrixWrapper<TScalar> & rightMatrix) {
	return BlazeMatrixWrapper<TScalar>(leftScalar * rightMatrix.m_Matrix);
}



//
// Left Matrix / Right Vector :
//
template<class TScalar>
BlazeVectorWrapper<TScalar> operator*(BlazeMatrixWrapper<TScalar> const &leftMatrix, BlazeVectorWrapper<TScalar> const &rightVector) {
	return BlazeVectorWrapper<TScalar>(leftMatrix.m_Matrix * rightVector.m_Vector);
}



//
// Other operations :
//
template <class TScalar>
TScalar vnl_trace(BlazeMatrixWrapper<TScalar> const& M) {
//	if(M.cols()!=M.rows())
	TScalar result = 0.0;
	for(size_t i=0; i<M.rows(); i++)
		result += M(i,i);
	return result;
}

template <class TScalar>
BlazeMatrixWrapper<TScalar> vnl_inverse(BlazeMatrixWrapper<TScalar> const& M) {
	exit(42);
//	return BlazeMatrixWrapper<TScalar>(M.m_Matrix.inverse());
}

template <class TScalar>
BlazeMatrixWrapper<TScalar> inverse_sympd(BlazeMatrixWrapper<TScalar> const& M) {
	exit(42);
//	return BlazeMatrixWrapper<TScalar>(M.m_Matrix.inverse());
}


template <class TScalar>
std::ostream & operator<<(std::ostream& os, BlazeMatrixWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Matrix);
}




#endif /* _BlazeMatrixWrapper_friend_h */

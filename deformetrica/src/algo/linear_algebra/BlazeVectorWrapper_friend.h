/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _BlazeVectorWrapper_friend_h
#define _BlazeVectorWrapper_friend_h

#ifndef _BlazeVectorWrapper_h
#error Do not include BlazeVectorWrapper_friend.h : include BlazeVectorWrapper.h instead
#endif

#include "BlazeVectorWrapper.h"
#include <iostream>


//
// Left Vector / Right Vector :
//
template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator-(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right) {
	return BlazeVectorWrapper<TScalar>(left.m_Vector - right.m_Vector);
}

template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator+(const BlazeVectorWrapper<TScalar> & left, const BlazeVectorWrapper<TScalar> & right) {
	return BlazeVectorWrapper<TScalar>(left.m_Vector + right.m_Vector);
}



//
// Left Vector / Right Scalar :
//
template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator+(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return (BlazeVectorWrapper<TScalar>(leftVector.m_Vector) + rightScalar);//*BlazeVectorWrapper<TScalar>(leftVector.size(), 1.0));
}

template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator-(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return (BlazeVectorWrapper<TScalar>(leftVector.m_Vector) - rightScalar);//*BlazeVectorWrapper<TScalar>(leftVector.size(), 1.0));
}

template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator/(const BlazeVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return BlazeVectorWrapper<TScalar>(leftVector.m_Vector / rightScalar);
}



//
// Left Scalar / Right Vector :
//
template<class TScalar>
inline BlazeVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const BlazeVectorWrapper<TScalar> & rightVector) {
	return BlazeVectorWrapper<TScalar>(leftScalar * rightVector.m_Vector);
}

//
// Other operations :
//
template <class TScalar>
std::ostream & operator<<(std::ostream& os, BlazeVectorWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Vector);
}



#endif /* _BlazeVectorWrapper_friend_h */

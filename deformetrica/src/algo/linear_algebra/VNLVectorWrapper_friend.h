/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _VNLVectorWrapper_friend_h
#define _VNLVectorWrapper_friend_h

#ifndef _VNLVectorWrapper_h
#error Do not include VNLVectorWrapper_friend.h : include VNLVectorWrapper.h instead
#endif

#include "VNLVectorWrapper.h"
#include <iostream>


//
// Left Vector / Right Vector :
//
template<class TScalar>
inline VNLVectorWrapper<TScalar> operator-(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right) {
	return VNLVectorWrapper<TScalar>(left.m_Vector - right.m_Vector);
}

template<class TScalar>
inline VNLVectorWrapper<TScalar> operator+(const VNLVectorWrapper<TScalar> & left, const VNLVectorWrapper<TScalar> & right) {
	return VNLVectorWrapper<TScalar>(left.m_Vector + right.m_Vector);
}



//
// Left Vector / Right Scalar :
//
template<class TScalar>
inline VNLVectorWrapper<TScalar> operator+(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return VNLVectorWrapper<TScalar>(leftVector.m_Vector + rightScalar);
}

template<class TScalar>
inline VNLVectorWrapper<TScalar> operator-(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return VNLVectorWrapper<TScalar>(leftVector.m_Vector - rightScalar);
}

template<class TScalar>
inline VNLVectorWrapper<TScalar> operator/(const VNLVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return VNLVectorWrapper<TScalar>(leftVector.m_Vector / rightScalar);
}



//
// Left Scalar / Right Vector :
//
template<class TScalar>
inline VNLVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const VNLVectorWrapper<TScalar> & rightVector) {
	return VNLVectorWrapper<TScalar>(leftScalar * rightVector.m_Vector);
}

//
// Other operations :
//
template <class TScalar>
std::ostream & operator<<(std::ostream& os, VNLVectorWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Vector);
}



#endif /* _VNLVectorWrapper_friend_h */

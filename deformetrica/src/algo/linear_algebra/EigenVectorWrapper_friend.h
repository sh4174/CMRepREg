/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _EigenVectorWrapper_friend_h
#define _EigenVectorWrapper_friend_h

#ifndef _EigenVectorWrapper_h
#error Do not include EigenVectorWrapper_friend.h : include EigenVectorWrapper.h instead
#endif

#include "EigenVectorWrapper.h"
#include <iostream>


//
// Left Vector / Right Vector :
//
template<class TScalar>
inline EigenVectorWrapper<TScalar> operator-(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right) {
	return EigenVectorWrapper<TScalar>(left.m_Vector - right.m_Vector);
}

template<class TScalar>
inline EigenVectorWrapper<TScalar> operator+(const EigenVectorWrapper<TScalar> & left, const EigenVectorWrapper<TScalar> & right) {
	return EigenVectorWrapper<TScalar>(left.m_Vector + right.m_Vector);
}



//
// Left Vector / Right Scalar :
//
template<class TScalar>
inline EigenVectorWrapper<TScalar> operator+(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return (EigenVectorWrapper<TScalar>(leftVector.m_Vector) + rightScalar*EigenVectorWrapper<TScalar>(leftVector.size(), 1.0));
}

template<class TScalar>
inline EigenVectorWrapper<TScalar> operator-(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return (EigenVectorWrapper<TScalar>(leftVector.m_Vector) - rightScalar*EigenVectorWrapper<TScalar>(leftVector.size(), 1.0));
}

template<class TScalar>
inline EigenVectorWrapper<TScalar> operator/(const EigenVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return EigenVectorWrapper<TScalar>(leftVector.m_Vector / rightScalar);
}



//
// Left Scalar / Right Vector :
//
template<class TScalar>
inline EigenVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const EigenVectorWrapper<TScalar> & rightVector) {
	return EigenVectorWrapper<TScalar>(leftScalar * rightVector.m_Vector);
}

//
// Other operations :
//
template <class TScalar>
std::ostream & operator<<(std::ostream& os, EigenVectorWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Vector);
}



#endif /* _EigenVectorWrapper_friend_h */

/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _ArmadilloVectorWrapper_friend_h
#define _ArmadilloVectorWrapper_friend_h

#ifndef _ArmadilloVectorWrapper_h
#error Do not include ArmadilloVectorWrapper_friend.h : include ArmadilloVectorWrapper.h instead
#endif

#include "ArmadilloVectorWrapper.h"
#include <iostream>


//
// Left Vector / Right Vector :
//
template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator-(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right) {
	return ArmadilloVectorWrapper<TScalar>(left.m_Vector - right.m_Vector);
}

template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator+(const ArmadilloVectorWrapper<TScalar> & left, const ArmadilloVectorWrapper<TScalar> & right) {
	return ArmadilloVectorWrapper<TScalar>(left.m_Vector + right.m_Vector);
}



//
// Left Vector / Right Scalar :
//
template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator+(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return ArmadilloVectorWrapper<TScalar>(leftVector.m_Vector + rightScalar);
}

template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator-(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return ArmadilloVectorWrapper<TScalar>(leftVector.m_Vector - rightScalar);
}

template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator/(const ArmadilloVectorWrapper<TScalar> & leftVector, const TScalar & rightScalar) {
	return ArmadilloVectorWrapper<TScalar>(leftVector.m_Vector / rightScalar);
}



//
// Left Scalar / Right Vector :
//
template<class TScalar>
inline ArmadilloVectorWrapper<TScalar> operator*(const TScalar & leftScalar, const ArmadilloVectorWrapper<TScalar> & rightVector) {
	return ArmadilloVectorWrapper<TScalar>(leftScalar * rightVector.m_Vector);
}

//
// Other operations :
//
template <class TScalar>
std::ostream & operator<<(std::ostream& os, ArmadilloVectorWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Vector.t());
}



#endif /* _ArmadilloVectorWrapper_friend_h */

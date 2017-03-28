/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _EigenMatrixWrapper_friend_h
#define _EigenMatrixWrapper_friend_h

#ifndef _EigenMatrixWrapper_h
#error Do not include EigenMatrixWrapper_friend.h : include EigenMatrixWrapper.h instead
#endif

#include "EigenVectorWrapper.h"
#include <iostream>



//
// Left Matrix / Right Matrix :
//
template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator-(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right) {
	return EigenMatrixWrapper<TScalar>(left.m_Matrix - right.m_Matrix);
}

template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator+(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right) {
	return EigenMatrixWrapper<TScalar>(left.m_Matrix + right.m_Matrix);
}

template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator*(const EigenMatrixWrapper<TScalar> & left, const EigenMatrixWrapper<TScalar> & right) {
	return EigenMatrixWrapper<TScalar>(left.m_Matrix * right.m_Matrix);
}



//
// (Left-Right) Matrix / Scalar :
//

template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator/(const EigenMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar) {
	return EigenMatrixWrapper<TScalar>(leftMatrix.m_Matrix / rightScalar);
}


template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator*(const EigenMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar) {
	return EigenMatrixWrapper<TScalar>(leftMatrix.m_Matrix * rightScalar);
}

template<class TScalar>
inline EigenMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const EigenMatrixWrapper<TScalar> & rightMatrix) {
	return EigenMatrixWrapper<TScalar>(leftScalar * rightMatrix.m_Matrix);
}



//
// Left Matrix / Right Vector :
//
template<class TScalar>
EigenVectorWrapper<TScalar> operator*(EigenMatrixWrapper<TScalar> const &leftMatrix, EigenVectorWrapper<TScalar> const &rightVector) {
	return EigenVectorWrapper<TScalar>(leftMatrix.m_Matrix * rightVector.m_Vector);
}



//
// Other operations :
//
template <class TScalar>
TScalar vnl_trace(EigenMatrixWrapper<TScalar> const& M) {
	return M.m_Matrix.trace();
}

template <class TScalar>
EigenMatrixWrapper<TScalar> vnl_inverse(EigenMatrixWrapper<TScalar> const& M) {
	return EigenMatrixWrapper<TScalar>(M.m_Matrix.inverse());
}

template <class TScalar>
EigenMatrixWrapper<TScalar> inverse_sympd(EigenMatrixWrapper<TScalar> const& M) {
	return EigenMatrixWrapper<TScalar>(M.m_Matrix.inverse());
}


template <class TScalar>
std::ostream & operator<<(std::ostream& os, EigenMatrixWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Matrix);
}




#endif /* _EigenMatrixWrapper_friend_h */

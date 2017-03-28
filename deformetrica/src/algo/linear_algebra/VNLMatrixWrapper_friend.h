/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _VNLMatrixWrapper_friend_h
#define _VNLMatrixWrapper_friend_h

#ifndef _VNLMatrixWrapper_h
#error Do not include VNLMatrixWrapper_friend.h : include VNLMatrixWrapper.h instead
#endif

#include "VNLVectorWrapper.h"
#include <iostream>

#include "vnl/vnl_trace.h"
#include "vnl/algo/vnl_matrix_inverse.h"



//
// Left Matrix / Right Matrix :
//
template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator-(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right) {
	return VNLMatrixWrapper<TScalar>(left.m_Matrix - right.m_Matrix);
}

template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator+(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right) {
	return VNLMatrixWrapper<TScalar>(left.m_Matrix + right.m_Matrix);
}

template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator*(const VNLMatrixWrapper<TScalar> & left, const VNLMatrixWrapper<TScalar> & right) {
	return VNLMatrixWrapper<TScalar>(left.m_Matrix * right.m_Matrix);
}



//
// (Left-Right) Matrix / Scalar :
//

template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator/(const VNLMatrixWrapper<TScalar> & leftMatrix, int const& rightScalar) {
	return VNLMatrixWrapper<TScalar>(leftMatrix.m_Matrix / rightScalar);
}


template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator*(const VNLMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar) {
	return VNLMatrixWrapper<TScalar>(leftMatrix.m_Matrix * rightScalar);
}

template<class TScalar>
inline VNLMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const VNLMatrixWrapper<TScalar> & rightMatrix) {
	return VNLMatrixWrapper<TScalar>(leftScalar * rightMatrix.m_Matrix);
}



//
// Left Matrix / Right Vector :
//
template<class TScalar>
VNLVectorWrapper<TScalar> operator*(VNLMatrixWrapper<TScalar> const &leftMatrix, VNLVectorWrapper<TScalar> const &rightVector) {
	return VNLVectorWrapper<TScalar>(leftMatrix.m_Matrix * rightVector.m_Vector);
}



//
// Other operations :
//
template <class TScalar>
TScalar vnl_trace(VNLMatrixWrapper<TScalar> const& M) {
	return vnl_trace(M.m_Matrix);
}


template <class TScalar>
VNLMatrixWrapper<TScalar> vnl_inverse(VNLMatrixWrapper<TScalar> const& M) {
	return VNLMatrixWrapper<TScalar>(vnl_matrix_inverse<TScalar>(M.m_Matrix));
}

template <class TScalar>
VNLMatrixWrapper<TScalar> inverse_sympd(VNLMatrixWrapper<TScalar> const& M) {
	return VNLMatrixWrapper<TScalar>(vnl_matrix_inverse<TScalar>(M.m_Matrix));
}

template <class TScalar>
std::ostream & operator<<(std::ostream& os, VNLMatrixWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Matrix);
}




#endif /* _VNLMatrixWrapper_friend_h */

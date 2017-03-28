/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#ifndef _ArmadilloMatrixWrapper_friend_h
#define _ArmadilloMatrixWrapper_friend_h

#ifndef _ArmadilloMatrixWrapper_h
#error Do not include ArmadilloMatrixWrapper_friend.h : include ArmadilloMatrixWrapper.h instead
#endif

#include "ArmadilloVectorWrapper.h"
#include <iostream>




//
// Left Matrix / Right Matrix :
//
template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator-(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right) {
	return ArmadilloMatrixWrapper<TScalar>(left.m_Matrix - right.m_Matrix);
}

template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator+(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right) {
	return ArmadilloMatrixWrapper<TScalar>(left.m_Matrix + right.m_Matrix);
}

template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator*(const ArmadilloMatrixWrapper<TScalar> & left, const ArmadilloMatrixWrapper<TScalar> & right) {
	return ArmadilloMatrixWrapper<TScalar>(left.m_Matrix * right.m_Matrix);
}



//
// (Left-Right) Matrix / Scalar :
//
template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator*(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, TScalar const& rightScalar) {
	return ArmadilloMatrixWrapper<TScalar>(leftMatrix.m_Matrix * rightScalar);
}

template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator*(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix) {
	return ArmadilloMatrixWrapper<TScalar>(leftScalar * rightMatrix.m_Matrix);
}


template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator/(const ArmadilloMatrixWrapper<TScalar> & leftMatrix, unsigned int const& rightScalar) {
	return ArmadilloMatrixWrapper<TScalar>(leftMatrix.m_Matrix / rightScalar);
}

template<class TScalar>
inline ArmadilloMatrixWrapper<TScalar> operator/(TScalar const& leftScalar, const ArmadilloMatrixWrapper<TScalar> & rightMatrix) {
	return ArmadilloMatrixWrapper<TScalar>(leftScalar / rightMatrix.m_Matrix);
}


//
// Left Matrix / Right Vector :
//
template<class TScalar>
ArmadilloVectorWrapper<TScalar> operator*(ArmadilloMatrixWrapper<TScalar> const &leftMatrix, ArmadilloVectorWrapper<TScalar> const &rightVector) {
	return ArmadilloVectorWrapper<TScalar>(leftMatrix.m_Matrix * rightVector.m_Vector);
}



//
// Other operations :
//

template <class TScalar>
ArmadilloMatrixWrapper<TScalar> diagonal_matrix(unsigned N, TScalar const & value) {
	typename ArmadilloMatrixWrapper<TScalar>::ArmadilloMatrixType result;
	return ArmadilloMatrixWrapper<TScalar>( result.eye(N, N)*value );
}

template <class TScalar>
TScalar trace(ArmadilloMatrixWrapper<TScalar> const& M) {
	return arma::trace(M.m_Matrix);
}

template <class TScalar>
ArmadilloMatrixWrapper<TScalar> inverse(ArmadilloMatrixWrapper<TScalar> const& M) {
	return ArmadilloMatrixWrapper<TScalar>(arma::inv(M.m_Matrix));
}

template <class TScalar>
ArmadilloMatrixWrapper<TScalar> inverse_sympd(ArmadilloMatrixWrapper<TScalar> const& M) {
	return ArmadilloMatrixWrapper<TScalar>(arma::inv_sympd(M.m_Matrix));
}

template <class TScalar>
ArmadilloVectorWrapper<TScalar> eigenvalues_sym(ArmadilloMatrixWrapper<TScalar> const& M) {
	//	arma::Col<TScalar> eigen_values;
	typename ArmadilloVectorWrapper<TScalar>::ArmadilloVectorType eigen_values;
	arma::eig_sym(eigen_values, M.m_Matrix);
	return ArmadilloVectorWrapper<TScalar>(eigen_values);
}


template <class TScalar>
std::ostream & operator<<(std::ostream& os, ArmadilloMatrixWrapper<TScalar> const& rhs) {
	return (os << rhs.m_Matrix);
}




#endif /* _ArmadilloMatrixWrapper_friend_h */

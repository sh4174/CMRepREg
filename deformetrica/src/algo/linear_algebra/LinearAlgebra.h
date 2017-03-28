/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LinearAlgebra_h
#define _LinearAlgebra_h

/* C++11
#include "VectorInterface.h"
#include "MatrixInterface.h"
*/

#undef DEFORMETRICA_WITH_LINALG_VNL
#undef DEFORMETRICA_WITH_LINALG_EIGEN
#undef DEFORMETRICA_WITH_LINALG_BLAZE

//#define DEFORMETRICA_WITH_LINALG_VNL
//#define DEFORMETRICA_WITH_LINALG_EIGEN
//#define DEFORMETRICA_WITH_LINALG_BLAZE
//#define EIGEN_NO_DEBUG


#ifdef DEFORMETRICA_WITH_LINALG_VNL
	#include "VNLVectorWrapper.h"
	#include "VNLMatrixWrapper.h"
#elif defined(DEFORMETRICA_WITH_LINALG_EIGEN)
	#include "EigenVectorWrapper.h"
	#include "EigenMatrixWrapper.h"
#elif defined(DEFORMETRICA_WITH_LINALG_BLAZE)
	#include "BlazeVectorWrapper.h"
	#include "BlazeMatrixWrapper.h"
#else
	#include "ArmadilloVectorWrapper.h"
	#include "ArmadilloMatrixWrapper.h"
#endif


/**
 *  \brief      Linear algebra class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The LinearAlgebra class contains only typedef for Matrix and Vector selecting the best linear algebra library suited for Deformetrica.
 */
template<class TScalar>
class LinearAlgebra {

public :

#ifdef DEFORMETRICA_WITH_LINALG_VNL
	/// (VNL) Vector type.
	typedef VNLVectorWrapper<TScalar> Vector;
	/// (VNL) Matrix type.
	typedef VNLMatrixWrapper<TScalar> Matrix;
#elif defined(DEFORMETRICA_WITH_LINALG_EIGEN)
	/// (Eigen) Vector type.
	typedef EigenVectorWrapper<TScalar> Vector;
	/// (Eigen) Matrix type.
	typedef EigenMatrixWrapper<TScalar> Matrix;
#elif defined(DEFORMETRICA_WITH_LINALG_BLAZE)
	/// (Blaze) Vector type.
	typedef BlazeVectorWrapper<TScalar> Vector;
	/// (Blaze) Matrix type.
	typedef BlazeMatrixWrapper<TScalar> Matrix;
#else
	/// (Armadillo) Vector type.
	typedef ArmadilloVectorWrapper<TScalar> Vector;
	/// (Armadillo) Matrix type.
	typedef ArmadilloMatrixWrapper<TScalar> Matrix;
#endif

};

/* C++11
#ifdef DEFORMETRICA_WITH_LINALG_VNL
	template <class TScalar >
	using Vector = VectorInterface<TScalar, VNLVectorWrapper>;
	//typedef VectorInterface<double, VNLVectorWrapper> VectorTest;
#endif*/


#endif /* _LinearAlgebra_h */

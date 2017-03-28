/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _KernelType_h
#define _KernelType_h


///	Possible type of kernels.
typedef enum
{
	null, 		/*!< Default value. */
	Exact,		/*!< Kernel with exact computation (see ExactKernel). */
#ifdef USE_CUDA
	CUDAExact,	/*!< Kernel with exact computation on GPU (see ExactKernel). */
#endif
	FGT,		/*!< Kernel with Fast Gauss Transform computation (see FGTKernel). */
	P3M			/*!< Kernel with linearly spaced grid computation (see P3MKernel). */
} KernelEnumType;


#endif /* _KernelType_h */

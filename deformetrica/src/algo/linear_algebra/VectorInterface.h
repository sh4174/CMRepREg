/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _VectorInterface_h
#define _VectorInterface_h

template <class TScalar, template<class> class Wrapper>
class VectorInterface : public Wrapper<TScalar> {

};


#endif /* _VectorInterface_h */

#ifndef _CauchyFunction_h
#define _CauchyFunction_h

#include "RadialFunction.h"



template < typename TYPE >
class CauchyFunction : public RadialFunction<TYPE>
{
    TYPE Sigma, ooSigma2, ooSigma4;

    public:

    CauchyFunction() { }

    CauchyFunction(TYPE sigma)
    {
        Sigma = sigma;
        ooSigma2 = 1.0/(Sigma*Sigma);
        ooSigma4 = 1.0/(ooSigma2*ooSigma2);
    }

    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        return 1.0/(1.0+r2*ooSigma2);
    }

    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        TYPE u = 1.0+r2*ooSigma2;
        return - ooSigma2 / (u*u);
    }

    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        TYPE u = 1.0+r2*ooSigma2;
        return 2.0 * ooSigma4 / (u*u*u);
    }

    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        TYPE u = 1.0/(1.0+r2*ooSigma2);
        *d1 = - ooSigma2 * u * u;
        *d2 = - 2.0 * ooSigma2 * *d1 * u;
    }

}; /* class CauchyFunction */



#endif /* _CauchyFunction_h */

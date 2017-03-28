#ifndef _GpuConv2D_h
#define _GpuConv2D_h


///////////////////////////////////////
///// Conv2D ////////////////////////////
///////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuEvalConv2D(KER Ker, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);



template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuEvalConv2D(TYPE sigma, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);



////////////////////////////////////////
///// GRAD1 Conv2D ///////////////////////
////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGrad1Conv2D(KER Ker, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);


template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGrad1Conv2D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);



///////////////////////////////////////
////////// GRAD Conv2D //////////////////
///////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradConv2D(KER Ker,
        TYPE* alpha_h, TYPE* x_h, TYPE* beta_h, TYPE* gamma_h,
         int nx);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradConv2D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* beta_h, TYPE* gamma_h, int nx);



////////////////////////////////////////////
////////// GRAD DIFF Conv2D //////////////////
////////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradDiffConv2D(KER Ker,
        TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h,
         int nx);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradDiffConv2D(TYPE sigma, TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h, int nx);



////////////////////////////////////////////
////////// DIFF Conv2D ///////////////////////
////////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuDiffConv2D(KER Ker,
        TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h,
        int nx);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuDiffConv2D(TYPE sigma, TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h, int nx);



#endif /* _GpuConv2D_h */

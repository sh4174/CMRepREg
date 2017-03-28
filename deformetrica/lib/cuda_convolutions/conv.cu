#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#include "GpuConv1D.cu"
#include "GpuConv2D.cu"


#define DIMPOINT 2 // dimension of the ambiant space : for curve ==2 or 3; for surface ==3
#define DIMVECT 2 // dimension of the object : for curve ==2; for surface ==3


// Compilation instruction : nvcc -arch=sm_?? conv.cu -o conv


/////////////////////////////////////
//               MAIN
/////////////////////////////////////


int main()
{

	/*---------*/
	/*  init   */
	/*---------*/

	float* x_h;
	float* y_h;
	float* beta_h;

	// Dimension of the problem 
	int nx = 100; // nbr of "row" (== nbre of point in the first object) may vary from 100 to 100000
	int ny = 1000000; // nbr of "column" (== nbre of point in the second object) may vary from 2000 to 100000 

	// arbitrary data
	x_h = (float *) malloc (nx * DIMPOINT * sizeof(float));
	for (int i=0;i<nx*DIMPOINT;i++){x_h[i] = (float)i/((float)nx*3.0f);}

	y_h = (float *) malloc (ny * DIMPOINT * sizeof(float));
	for (int i=0;i<ny*DIMPOINT;i++){y_h[i] = logf( ((float)i+1.0f)/((float)ny*3.0f));}
	
	beta_h = (float *) malloc (ny * DIMVECT * sizeof(float));
	for (int i=0;i<ny*DIMVECT;i++){beta_h[i] = cosf((float)i*(float)i/((float)ny*3.0f));}


	/*----------*/
	/*  kernel  */
	/*----------*/

	printf("\n");

	///////////
	// Conv2 //
	///////////

	clock_t tic2 = clock();
	// pointer to output
	float* gamma_h2;
	gamma_h2 = (float *) malloc (nx * DIMVECT * sizeof(float));	
	GaussGpuEvalConv2D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h2, nx, ny);


	clock_t toc2 = clock();
	printf("Conv2d with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc2 - tic2) / CLOCKS_PER_SEC);

	///////////////
	// Conv2 bis //
	///////////////

	clock_t tic2_bis = clock();
	// pointer to output
	float* gamma_h2_bis;
	gamma_h2_bis = (float *) malloc (nx * DIMVECT * sizeof(float));	
	GaussGpuEvalConv2D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h2_bis, nx, ny);
	clock_t toc2_bis = clock();
	printf("Conv2d bis with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc2_bis - tic2_bis) / CLOCKS_PER_SEC);

	///////////////
	// Conv2 ter //
	///////////////

	clock_t tic2_ter = clock();
	// pointer to output
	float* gamma_h2_ter;
	gamma_h2_ter = (float *) malloc (nx * DIMVECT * sizeof(float));	
	GaussGpuEvalConv2D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h2_ter, nx, ny);
	clock_t toc2_ter = clock();
	printf("Conv2d ter with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc2_ter - tic2_ter) / CLOCKS_PER_SEC);

	///////////
	// Conv1 //
	///////////

	clock_t tic1 = clock();
	// pointer to output
	float* gamma_h1;
	gamma_h1 = (float *) malloc (nx * DIMVECT * sizeof(float));
	GaussGpuEvalConv1D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h1, nx, ny);
	clock_t toc1 = clock();
	printf("Conv1d with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc1 - tic1) / CLOCKS_PER_SEC);

	///////////////
	// Conv1 bis //
	///////////////

	clock_t tic1_bis = clock();
	// pointer to output
	float* gamma_h1_bis;
	gamma_h1_bis = (float *) malloc (nx * DIMVECT * sizeof(float));	
	GaussGpuEvalConv1D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h1_bis, nx, ny);
	clock_t toc1_bis = clock();
	printf("Conv1d bis with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc1_bis - tic1_bis) / CLOCKS_PER_SEC);

	///////////////
	// Conv1 ter //
	///////////////

	clock_t tic1_ter = clock();
	// pointer to output
	float* gamma_h1_ter;
	gamma_h1_ter = (float *) malloc (nx * DIMVECT * sizeof(float));	
	GaussGpuEvalConv1D<float,DIMPOINT,DIMVECT>(0.5,x_h, y_h, beta_h, gamma_h1_ter, nx, ny);
	clock_t toc1_ter = clock();
	printf("Conv1d ter with nx=%d and ny=%d took %f seconds\n",nx,ny ,(double)(toc1_ter - tic1_ter) / CLOCKS_PER_SEC);


	/*--------*/
	/* Output */
	/*--------*/

	float err =0;
	for (int i=0;i<DIMVECT*nx;i++){
	       float erR = fabs((gamma_h1_bis[i] - gamma_h2_bis[i]) / gamma_h1_bis[i]) ;
	       if (erR > err) {
	          err = erR;
	       }		  
	}

	printf("Max relative error : %f \n\n",err);
	return 0;

}




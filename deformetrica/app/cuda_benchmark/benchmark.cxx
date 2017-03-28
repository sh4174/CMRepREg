#include "DeformetricaConfig.h"

#include <iostream>

#include "LinearAlgebra.h"

#include "KernelFactory.h"
#include "KernelType.h"

#include "SimpleTimer.h"

#include "readMatrixDLM.txx"

using namespace std;

int main(int argc, char** argv)
{
#ifdef USE_CUDA

	#ifdef USE_DOUBLE_PRECISION
		typedef double TScalar;
		std::cout << "(Computations are in double precision)" << std::endl << std::endl;
	#else
		typedef float TScalar;
		std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
	#endif
	typedef LinearAlgebra<TScalar>::Vector VectorType;
	typedef LinearAlgebra<TScalar>::Matrix MatrixType;

	const unsigned int Dimension = 3;
	const unsigned int WeightDimension = 1*Dimension;

	int N = 5000;
	int M = 5000;
	TScalar sigma = 0.25;

	if (argc != 4 )
	{
		cerr << "Usage: " << argv[0] << " sigma N M " << endl;
		cerr << "Parameters are set to default values." << endl;
	}
	else
	{
		sigma = atof(argv[1]);
		N = atoi(argv[2]);
		M = atoi(argv[3]);
	}
	cout << "Sigma = " << sigma << ", N = " << N << ", M = " << M << endl;

	MatrixType X(M, Dimension, 0.0);
	MatrixType alpha(M, Dimension, 0.0);
	MatrixType sources(N, Dimension, 0.0);
	MatrixType weights(N, WeightDimension, 0.0);
	MatrixType xi(N, WeightDimension, 0.0);

	for(int j = 0; j < Dimension; j++) {
		for(int i = 0; i < N; i++)
		{
			sources(i, j) = cos((i+j)*(i+j));
		}
	}
	for(int j = 0; j < Dimension; j++) {
		for(int i = 0; i < M; i++)
		{
			alpha(i, j) = cos((i+j)*(i+j));
		}
	}
	for(int j = 0; j < WeightDimension; j++) {
		for(int i = 0; i < N; i++)
		{
			weights(i, j) = sin((i+j)*(i+j));
		}
	}
	for(int j = 0; j < Dimension; j++) {
		for(int i = 0; i < M; i++)
		{
			X(i, j) = sqrt(1+i+j)/sqrt(M*N);
			xi(i, j) = cos((i+j)*(i+j));
		}
	}

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef KernelFactoryType::KernelBaseType KernelType;

	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObjCPU = kFactory->CreateKernelObject(sources, weights, sigma, Exact);

	KernelType* kernelObjCUDA = kFactory->CreateKernelObject(CUDAExact);//sources, weights, sigma);
	kernelObjCUDA->SetKernelWidth(sigma);

	//
	// Initialization of the kernel :
	//
	SimpleTimer initGPU;
	kernelObjCUDA->SetSources(MatrixType(3, 3, 1.0));
	kernelObjCUDA->SetWeights(MatrixType(3, 3, 1.0));
	kernelObjCUDA->Convolve(MatrixType(3, 3, 1.0));
	TScalar timeInitGPU = initGPU.GetElapsedCPUTimeInSecondsOnly();
	cout << endl << "GPU Initialization : " << timeInitGPU << "s" << endl;
	// End initialization :

	kernelObjCUDA->SetSources(sources);
	kernelObjCUDA->SetWeights(weights);

	cout << "-------------------------" << endl;
	cout << "Convolution" << endl;
	cout << "-------------------------" << endl;

	SimpleTimer timerGPU;
	MatrixType result_cuda = kernelObjCUDA->Convolve(X);
	TScalar timeGPU = timerGPU.GetElapsedCPUTimeInSecondsOnly();
	cout << "GPU : " << timeGPU << "s" << endl;

	SimpleTimer timerCPU;
	MatrixType result_cpu = kernelObjCPU->Convolve(X);
	TScalar timeCPU = timerCPU.GetElapsedCPUTimeInSecondsOnly();
	cout << "CPU : " << timeCPU << "s" << endl;


	TScalar diff = 0.0;
	TScalar error = 0.0;
	TScalar norm = 0.0;
	TScalar temp = 0.0;
	for(int j = 0; j < WeightDimension; j++) {
		for(int i = 0; i < result_cpu.rows(); i++)
		{
			diff = result_cpu(i, j) - result_cuda(i, j);
			error += sqrt(diff*diff);
			temp = result_cpu(i, j);
			norm += sqrt(temp*temp);
		}
	}
//	cout << "CPU = " << result_cpu << endl;
//	cout << "GPU = " << result_cuda << endl;
/*	for(int j = 0; j < WeightDimension; j++) {
		for(int i = 0; i < 10; i++)
		{
			cout << "CPU = " << result_cpu(i, j) << endl;
			cout << "GPU = " << result_cuda(i, j) << endl;
		}
	}*/
	error /= result_cpu.rows();
	norm /= result_cpu.rows();

	cout << "Erreur moyenne (relative): " << error << " (" << error/norm << ")" << endl;
	cout << "Gain : " << timeCPU/timeGPU << " x" << endl << endl;

	cout << "-------------------------" << endl;
	cout << "Convolution with gradient" << endl;
	cout << "-------------------------" << endl;

	SimpleTimer timerGPU_grad;
	MatrixType result_cuda_grad = kernelObjCUDA->ConvolveGradient(X, alpha);
	TScalar timeGPU_grad = timerGPU_grad.GetElapsedCPUTimeInSecondsOnly();
	cout << "GPU : " << timeGPU_grad << "s" << endl;

	SimpleTimer timerCPU_grad;
	MatrixType result_cpu_grad = kernelObjCPU->ConvolveGradient(X, alpha);
	TScalar timeCPU_grad = timerCPU_grad.GetElapsedCPUTimeInSecondsOnly();
	cout << "CPU : " << timeCPU_grad << "s" << endl;

	diff = 0.0;
	error = 0.0;
	norm = 0.0;
	temp = 0.0;
	for(int j = 0; j < result_cpu_grad.cols(); j++)
	{
		for(int i = 0; i < result_cpu_grad.rows(); i++)
		{
			diff = result_cpu_grad(i, j) - result_cuda_grad(i, j);
			error += sqrt(diff*diff);
			temp = result_cpu_grad(i, j);
			norm += sqrt(temp*temp);
		}
	}
	error /= result_cpu_grad.rows();
	norm /= result_cpu_grad.rows();

	cout << "Erreur moyenne (relative): " << error << " (" << error/norm << ")" << endl;
	cout << "Gain : " << timeCPU_grad/timeGPU_grad << " x" << endl << endl;




	cout << "-------------------------" << endl;
	cout << "Convolution with hessian" << endl;
	cout << "-------------------------" << endl;

	SimpleTimer timerGPU_hess;
	MatrixType result_cuda_hess = kernelObjCUDA->ConvolveSpecialHessian(xi);
	TScalar timeGPU_hess = timerGPU_hess.GetElapsedCPUTimeInSecondsOnly();
	cout << "GPU : " << timeGPU_hess << "s" << endl;

	SimpleTimer timerCPU_hess;
	MatrixType result_cpu_hess = kernelObjCPU->ConvolveSpecialHessian(xi);
	TScalar timeCPU_hess = timerCPU_hess.GetElapsedCPUTimeInSecondsOnly();
	cout << "CPU : " << timeCPU_hess << "s" << endl;

	diff = 0.0;
	error = 0.0;
	norm = 0.0;
	temp = 0.0;
	for(int j = 0; j < result_cpu_hess.cols(); j++)
	{
		for(int i = 0; i < result_cpu_hess.rows(); i++)
		{
			diff = result_cpu_hess(i, j) - result_cuda_hess(i, j);
			error += sqrt(diff*diff);
			temp = result_cpu_hess(i, j);
			norm += sqrt(temp*temp);
		}
	}
	error /= result_cpu_hess.rows();
	norm /= result_cpu_hess.rows();

	cout << "Erreur moyenne (relative): " << error << " (" << error/norm << ")" << endl;
	cout << "Gain : " << timeCPU_hess/timeGPU_hess << " x" << endl << endl;

#endif

	return 0;
}

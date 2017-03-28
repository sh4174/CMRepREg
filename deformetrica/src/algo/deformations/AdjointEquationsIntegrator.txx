#ifndef _AdjointEquationsIntegrator_txx
#define _AdjointEquationsIntegrator_txx

#include "SimpleTimer.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
AdjointEquationsIntegrator<TScalar, Dimension>
::AdjointEquationsIntegrator() :
m_IsLandmarkPoints(false), m_IsImagePoints(false), m_HasJumps(false), m_KernelObj1(NULL), m_KernelObj2(NULL), m_KernelObj3(NULL), m_KernelObj4(NULL)
{
}



template <class TScalar, unsigned int Dimension>
AdjointEquationsIntegrator<TScalar, Dimension>
::~AdjointEquationsIntegrator()
{
	
}

/*
template <class TScalar, unsigned int Dimension>
void 
AdjointEquationsIntegrator<TScalar, Dimension>
::SetInitialConditionsLandmarkPoints(MatrixType& M)
{
	m_InitialConditionsLandmarkPoints.resize(1);
	m_InitialConditionsLandmarkPoints[0] = M;
}

template <class TScalar, unsigned int Dimension>
void 
AdjointEquationsIntegrator<TScalar, Dimension>
::SetInitialConditionsImagePoints(MatrixType& M)
{
	m_InitialConditionsImagePoints.resize(1);
	m_InitialConditionsImagePoints[0] = M;
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::Update()
{
	m_XiPosT.resize(m_NumberOfTimePoints);
	m_XiMomT.resize(m_NumberOfTimePoints);
	m_ThetaT.resize(m_NumberOfTimePoints);
	m_EtaT.resize(m_NumberOfTimePoints);
	
	if (m_IsLandmarkPoints)
		this->IntegrateAdjointOfLandmarkPointsEquations();
	
	if (m_IsImagePoints)
	{
		if (m_ComputeTrueInverseFlow)
			this->IntegrateAdjointOfImagePointsBackward();
		else
			this->IntegrateAdjointOfImagePointsForward();
		
	}
	
	this->IntegrateAdjointOfDiffeoParametersEquations();
}


template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::IntegrateAdjointOfLandmarkPointsEquations()
{
	// Propagate theta backward
	m_ThetaT.resize(m_NumberOfTimePoints); // <--- Initialize each matrix at 0 ?
		
	m_ThetaT[m_NumberOfTimePoints-1] = m_ListInitialConditionsLandmarkPoints[0];

	// Initial condition is zero for regression
	if (m_HasJumps)
	{
		m_ThetaT[m_NumberOfTimePoints-1].fill(0.0);
	}
	
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);
	int subjIndex = m_JumpTimes.size()-1;

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	for (long t = m_NumberOfTimePoints-1; t > 0; t--)
	{

		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		if (m_HasJumps)
		{
			if (t == m_JumpTimes[subjIndex])
			{
				m_ThetaT[t] += m_ListInitialConditionsLandmarkPoints[subjIndex];
			}
		}
		
		MatrixType dTheta = kernelObj->ConvolveGradient(m_LandmarkPointsT[t], m_ThetaT[t]);

		m_ThetaT[t-1] = m_ThetaT[t] + dTheta * dt; // the plus is correct! dTheta should be negative.

		// Heun's method
		if (m_UseImprovedEuler)
		{

			kernelObj->SetSources(m_PosT[t-1]);
			kernelObj->SetWeights(m_MomT[t-1]);

			if (m_HasJumps)
			{
				if (t == m_JumpTimes[subjIndex])
				{
					//m_ThetaT[t-1] += m_ListInitialConditionsLandmarkPoints[subjIndex];
					subjIndex--;
					if (subjIndex < 0)
						subjIndex = 0;
				}
			}
			
			MatrixType dTheta2 = kernelObj->ConvolveGradient(m_LandmarkPointsT[t-1], m_ThetaT[t-1]);
			m_ThetaT[t-1] = m_ThetaT[t] + (dTheta + dTheta2) * (dt * 0.5);
		}

	}

	delete kernelObj;
}


// This computes \dot{\eta}(t) = -\partial_1 G(t) \eta(t)
template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::IntegrateAdjointOfImagePointsBackward()
{
	
	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// 
	// // The path of the image points over time
	// YT.resize(m_NumberOfTimePoints);
	// for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
	// {
	// 	//		YT[t] = GridFunctionsType::UpsampleImagePoints(tempLII->GetImage(), tempLII->GetDownSampledWorkingImage(), m_InverseMapsT[t]);
	// 	YT[t] = tempLII->UpSampleImageMap(m_InverseMapsT[t]);
	// }

	MatrixType& YT0 = m_ImagePointsT[0];

	// long numVoxels = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);
	int subjIndex = m_JumpTimes.size()-1;

	// The aux variable Eta over time
	m_EtaT.resize(m_NumberOfTimePoints);
	// The source term for integration. There is a + sign here contrary to in TransportAlongGeodesicForward
	m_EtaT[m_NumberOfTimePoints-1] = m_ListInitialConditionsImagePoints[0];

	// Initial condition is zero for regression
	if (m_HasJumps)
	{
		m_EtaT[m_NumberOfTimePoints-1].fill(0.0);
	}

	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	// Integrate backwards in time
	for (long t = (m_NumberOfTimePoints - 1); t>=1; t--)
	{
		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		// The velocity is always v_t(y0)
		// MatrixType VtY0 = kernelObj->Convolve(m_DownSampledY1);
		MatrixType VtY0 = kernelObj->Convolve(YT0);

		// This gives us a Dim*Dim matrix at every pixel of the original image
		// std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(m_DownSampledY1);
		std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(YT0);


		// We need to have eta(t) in the form of an image, so we can compute the jacobian
		std::vector<ImageTypePointer> etaTImage;
		etaTImage.resize(Dimension);
		for (unsigned int d = 0; d < Dimension; d++)
		{
			// Create an image from each Dimension
			ImageTypePointer img = GridFunctionsType::VectorToImage(m_FullResolutionImage, m_EtaT[t].get_column(d));
			etaTImage[d] = img;
		}

		// Compute jacobian of eta(t)
		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
		std::vector<ImageTypePointer> gradImages;
		gradImages.resize(Dimension*Dimension);
		unsigned int indx = 0;
		// Compute the gradient in all directions
		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
		{
			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
			{
				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
				derivf->SetInput(etaTImage[dim1]);
				derivf->SetDirection(dim2);
				derivf->SetOrder(1);
				derivf->SetUseImageSpacingOn();
				derivf->Update();

				gradImages[indx] = derivf->GetOutput();
				indx++;
			}
		}

		// Now we have all the information we need to compute dEta, but we need to loop over all the pixels
		// and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
		// IteratorType it(m_DownSampledWorkingImage, m_DownSampledWorkingImage->GetLargestPossibleRegion());
		// int matrixIndex = 0;
		// int numVoxelsDS = m_DownSampledWorkingImage->GetLargestPossibleRegion().GetNumberOfPixels();
		// MatrixType dEta(numVoxelsDS, Dimension, 0);
		IteratorType it(m_FullResolutionImage, m_FullResolutionImage->GetLargestPossibleRegion());
		int matrixIndex = 0;
		int numVoxels = m_FullResolutionImage->GetLargestPossibleRegion().GetNumberOfPixels();
		MatrixType dEta(numVoxels, Dimension, 0);

		// Loop over the grid to construct jacobian matrices and compute dY
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			typedef typename ImageType::IndexType ImageIndexType;
			ImageIndexType imageIndex = it.GetIndex();
			MatrixType jacobian(Dimension, Dimension);

			unsigned int tempMatrixIndex = 0;
			// Build the (DimensionxDimension) jacobian matrix
			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
			{
				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
				{
					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
					tempMatrixIndex++;
				}
			}

			// The trace of the first term
			TScalar traceValue = trace(firstTerm[matrixIndex]);

			// Get the (Dimension) vector corresponding to this points v_t(y0)
			VectorType VtY0k = VtY0.get_row(matrixIndex);

			// The final value for this pixel of dEta
			// JAMES
			// dEta.set_row(matrixIndex, traceValue*EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);
			// STANLEY
			dEta.set_row(matrixIndex, -traceValue*m_EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);


			matrixIndex++;
		}

		
		if (m_HasJumps)
		{
			if (t == m_JumpTimes[subjIndex])
			{
				dEta += m_ListInitialConditionsImagePoints[subjIndex];
				subjIndex--;
				if (subjIndex < 0)
					subjIndex = 0;
			}
		}

		// JAMES
		// EtaT[t-1] = EtaT[t] + dEta * dt;
		// STANLEY (backward integration: the differential changes sign)
		// also need to upsample image map
		// MatrixType dEta_HD = this->UpSampleImageMap(dEta);
		// EtaT[t-1] = EtaT[t] - dEta_HD * dt;
		m_EtaT[t-1] = m_EtaT[t] - dEta * dt;


	}

	// We have computed Eta(t), the backwards propagator actually needs -(d_yp Y)^t Eta(t)
	for (long t = 0; t<m_NumberOfTimePoints; t++)
	{
		// We need to have Eta(t) in the form of an image, so we can compute the jacobian
		std::vector<ImageTypePointer> YTImage;
		YTImage.resize(Dimension);
		for (unsigned int d = 0; d < Dimension; d++)
		{
			// Create an image from each Dimension
			ImageTypePointer img = GridFunctionsType::VectorToImage(m_FullResolutionImage, m_ImagePointsT[t].get_column(d));
			YTImage[d] = img;
		}

		// Compute jacobian of eta(t)
		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
		std::vector<ImageTypePointer> gradImages;
		gradImages.resize(Dimension*Dimension);
		unsigned int indx = 0;
		// Compute the gradient in all directions
		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
		{
			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
			{
				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
				derivf->SetInput(YTImage[dim1]);
				derivf->SetDirection(dim2);
				derivf->SetOrder(1);
				derivf->SetUseImageSpacingOn();
				derivf->Update();

				gradImages[indx] = derivf->GetOutput();
				indx++;
			}
		}

		// Now we need to loop over all the pixels and construct the 3x3 jacobian matrix
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
		// JAMES
		// IteratorType it(m_DownSampledWorkingImage, m_DownSampledWorkingImage->GetLargestPossibleRegion());
		// STANLEY: YT is high res image
		IteratorType it(m_FullResolutionImage, m_FullResolutionImage->GetLargestPossibleRegion());
		int matrixIndex = 0;

		// Loop over the grid to construct jacobian matrices and compute dY
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			typedef typename ImageType::IndexType ImageIndexType;
			ImageIndexType imageIndex = it.GetIndex();
			MatrixType jacobian(Dimension, Dimension);

			unsigned int tempMatrixIndex = 0;
			// Build the (DimensionxDimension) jacobian matrix
			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
			{
				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
				{
					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
					tempMatrixIndex++;
				}
			}

			// The final value for this pixel of dEta
			//JAMES
			// EtaT[t].set_row(matrixIndex, jacobian.transpose()*EtaT[t].get_row(matrixIndex));
			//STANLEY
			m_EtaT[t].set_row(matrixIndex, -jacobian.transpose()*m_EtaT[t].get_row(matrixIndex));
			matrixIndex++;
		}
	}


	// Deporte dans IntegrateAdjointDiffeoParametersEquations
	// // STANLEY
	// // the backward propagator needs YT = Y0 for all t
	// for (int t = 1; t < m_NumberOfTimePoints; t++)
	// {
	// 	YT[t] = YT[0];
	// }

    delete kernelObj;

}



template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::IntegrateAdjointOfImagePointsForward()
{

	typedef GridFunctions<TScalar, Dimension> GridFunctionsType;

	// long numVoxels = image->GetLargestPossibleRegion().GetNumberOfPixels();

	// Propagate eta forward
	m_EtaT.resize(m_NumberOfTimePoints);
	// STANLEY: PREVIOUS
	// EtaT[0] = Eta0;
	// STANLEY: NOW (the sign is taken into account here)


	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	for (long t = 0; t < (m_NumberOfTimePoints - 1); t++)
	{
		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		/*
		// Grad mom splatted at CP locations, with convolutions evaluated at y(t-1)
		// This computes alpha_p \nabla_1 K(y_k, c_p)
		std::vector<MatrixType> gradMom = kernelObj->ConvolveGradient(YT[t]);

		// This computes \eta_k^t alpha_p /nabla_1 K(y_k, c_p)
		MatrixType dEta(numVoxels, Dimension, 0);
		for (unsigned int j = 0; j < numVoxels; j++)
		{
		dEta.set_row(j, gradMom[j].transpose() * EtaT[t].get_row(j));
		}
		*/

		// Grad mom splatted at CP locations, with convolutions evaluated at y(t-1)
		// This computes \eta_k^t alpha_p /nabla_1 K(y_k, c_p)
		MatrixType dEta = kernelObj->ConvolveGradient(m_ImagePointsT[t], m_EtaT[t]);
		// assert(dEta.rows() == numVoxels);
		// assert(dEta.columns() == Dimension);

		m_EtaT[t + 1] = m_EtaT[t] - dEta * dt;

		// Heun's method
		if (m_UseImprovedEuler)
		{
			kernelObj->SetSources(m_PosT[t + 1]);
			kernelObj->SetWeights(m_MomT[t + 1]);

			/*
			gradMom = kernelObj->ConvolveGradient(YT[t + 1]);

			MatrixType dEta2(numVoxels, Dimension, 0);
			for (unsigned int j = 0; j < numVoxels; j++)
			{
			dEta2.set_row(j, gradMom[j].transpose() * EtaT[t + 1].get_row(j));
			}
			*/
			MatrixType dEta2 = kernelObj->ConvolveGradient(m_ImagePointsT[t + 1], m_EtaT[t + 1]);
			// assert(dEta2.rows() == numVoxels);
			// assert(dEta2.columns() == Dimension);

			m_EtaT[t + 1] = m_EtaT[t] - (dEta + dEta2) * (dt * 0.5);
		}

	}

	// std::cout << "EtaT[end] = " << std::endl;
	// std::cout << EtaT[numTimePoints-1].get_n_rows(0,20) << std::endl;


	// typedef GridFunctions<TScalar, Dimension> GridFunctionsType;
	// ImageTypePointer testT = GridFunctionsType::VectorToImage(m_Image, YT[19].get_column(0));
	// typedef itk::ImageFileWriter<ImageType> WriterType;
	//   typename WriterType::Pointer writer = WriterType::New();
	//   writer->SetInput(testT);
	//   writer->SetFileName("YT.vtk");
	//   writer->Update();

	delete kernelObj;
	

}





	
template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::IntegrateAdjointOfDiffeoParametersEquations()

{	
	long numCP = m_PosT[0].rows();

	// Initialization:  Xi(1) = 0 or Xi[T-1] = 0
	MatrixType zeroM(numCP, Dimension, 0);

	// m_XiPosT.resize(m_numTimePoints);
	// m_XiMomT.resize(m_numTimePoints);
	for (long t = 0; t < m_NumberOfTimePoints; t++)
	{
		m_XiPosT[t] = zeroM;
		m_XiMomT[t] = zeroM;
	}

	// Set up kernel objects
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	m_KernelObj1 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj1->SetKernelWidth(m_KernelWidth);
	m_KernelObj2 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj2->SetKernelWidth(m_KernelWidth);
	m_KernelObj3 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj3->SetKernelWidth(m_KernelWidth);
//	m_KernelObj4 = kFactory->CreateKernelObject(m_KernelType);
//	m_KernelObj4->SetKernelWidth(m_KernelWidth);

	// Backward propagation
	//TScalar dt = 1.0 / (numTimePoints - 1.0);
	TScalar dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);

	for (long t = m_NumberOfTimePoints-1; t >= 1; t--)
	{
		MatrixType dPos;
		MatrixType dMom;
		this->ComputeUpdateAt(t, dPos, dMom);

		//cout<<"dPos "<<dPos<<"\n";
		//cout<<"dMom "<<dMom<<"\n";

		//cout<<"dt "<<dt<<"\n";
		
		m_XiPosT[t-1] = m_XiPosT[t] + dPos*dt;
		m_XiMomT[t-1] = m_XiMomT[t] + dMom*dt;

		//cout<<"m_XiMomT[t] "<<m_XiMomT[t]<<"\n";

		//cout<<"m_XiPosT[t-1] "<<m_XiPosT[t-1]<<"\n";
		//cout<<"m_XiMomT[t-1] "<<m_XiMomT[t-1]<<"\n";
		//exit(0);

		// Heun's method
		if (m_UseImprovedEuler)
		{
			MatrixType dPos2;
			MatrixType dMom2;
			this->ComputeUpdateAt(t-1, dPos2, dMom2);

			m_XiPosT[t-1] = m_XiPosT[t] + (dPos + dPos2) * (dt * 0.5);
			m_XiMomT[t-1] = m_XiMomT[t] + (dMom + dMom2) * (dt * 0.5);
		}

		// std::cout << "numTimePoints = " << numTimePoints << std::endl;
		// std::cout << " t = " << t << std::endl;
		// std::cout << m_XiPosT[t-1] << std::endl;
		// std::cout << m_XiMomT[t-1] << std::endl;
		// std::cout << std::endl;

	} // for t

	//cout<<"m_XiMomT[0] "<<m_XiMomT[0]<<"\n";

	delete m_KernelObj1;
	delete m_KernelObj2;
	delete m_KernelObj3;

}








////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
AdjointEquationsIntegrator<TScalar, Dimension>
::ComputeUpdateAt(unsigned int s, MatrixType& dPos, MatrixType &dMom)
{

	long numCP = m_PosT[0].rows();

	KernelType* momXiPosKernelObj = m_KernelObj1;
	KernelType* etaKernelObj = m_KernelObj2;
	KernelType* tmpKernelObj = m_KernelObj3;


	// Concatenate landmark and image points, as well as their adjoint variables. Save time in convolution.
	int nbOfLandmarkPoints = m_IsLandmarkPoints?m_LandmarkPointsT[0].rows():0;
	int nbOfImagePoints = m_IsImagePoints?m_ImagePointsT[0].rows():0;
	int nbTotalPoints = nbOfLandmarkPoints + nbOfImagePoints;

	MatrixType ConcatenatedPoints(nbTotalPoints, Dimension);
	MatrixType ConcatenatedVectors(nbTotalPoints, Dimension);
	
	for (int r = 0; r < nbOfLandmarkPoints; r++)
	{
		ConcatenatedPoints.set_row(r, m_LandmarkPointsT[s].get_row(r));
		ConcatenatedVectors.set_row(r, m_ThetaT[s].get_row(r));
	}
	for (int r = 0; r < nbOfImagePoints; r++)
	{
		// Be careful: if ComputeTrueInverseFlow, vectors m_EtaT[s] are always attached to the fixed points m_ImagePointsT[0]!
		ConcatenatedPoints.set_row(nbOfLandmarkPoints + r, m_ComputeTrueInverseFlow?m_ImagePointsT[0].get_row(r):m_ImagePointsT[s].get_row(r));
		ConcatenatedVectors.set_row(nbOfLandmarkPoints + r, m_EtaT[s].get_row(r));
	}
	
	// Term 1
	// SimpleTimer timer;
	etaKernelObj->SetSources(ConcatenatedPoints);
	etaKernelObj->SetWeights(ConcatenatedVectors);

	/*
	// SimpleTimer subtimer;
	MatrixList gradEta = etaKernelObj->ConvolveGradient(PositionsT[s]);

	// std::cout << "\tsub-term 1.1: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// subtimer.Start();


	MatrixType dXi1(numCP, Dimension, 0);
	for (unsigned int i = 0; i < numCP; i++)
	{
	dXi1.set_row(i, gradEta[i].transpose() * MomentasT[s].get_row(i));
	}

	// std::cout << "\tsub-term 1.2: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	*/

	MatrixType dXi1 = etaKernelObj->ConvolveGradient(m_PosT[s], m_MomT[s]);
	// assert(dXi1.rows() == numCP);
	// assert(dXi1.columns() == Dimension);

	// std::cout << "term 1: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();
	// Term 2
	MatrixType dXi2 = etaKernelObj->Convolve(m_PosT[s]);

	// std::cout << "term 2: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();
	// Term 3
	MatrixType AXiPos(numCP, Dimension*2, 0);
	AXiPos.set_columns(0, m_MomT[s]);
	AXiPos.set_columns(Dimension, m_XiPosT[s]);

	momXiPosKernelObj->SetSources(m_PosT[s]);
	momXiPosKernelObj->SetWeights(AXiPos);

	// subtimer.Start();
	MatrixType kAXiPos = momXiPosKernelObj->Convolve(m_PosT[s]);

	MatrixList gradAXiPos =	momXiPosKernelObj->ConvolveGradient(m_PosT[s]);

	// std::cout << "\tsub term 3.1: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// subtimer.Start();

	MatrixType dXi3(numCP, Dimension, 0);
	for (unsigned int i = 0; i < numCP; i++)
	{
		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		MatrixType gradXiPos_i = gradAXiPos[i].get_n_rows(Dimension, Dimension);
		dXi3.set_row(i,gradMom_i.transpose() * m_XiPosT[s].get_row(i) + gradXiPos_i.transpose() * m_MomT[s].get_row(i));
	}

	// std::cout << "\tsub term 3.2: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// 
	// std::cout << "term 3: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();
	// Term 5
	MatrixType dXi5 = kAXiPos.get_n_columns(Dimension, Dimension);

	// std::cout << "term 5: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();
	// Term 6
	MatrixType dXi6(numCP, Dimension, 0);

	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		MatrixType W(numCP, Dimension, 0.0);
		for (unsigned int i = 0; i < numCP; i++)
			W.set_row(i, m_XiMomT[s](i,dim) * m_MomT[s].get_row(i));

		tmpKernelObj->SetSources(m_PosT[s]);
		tmpKernelObj->SetWeights(W);
		MatrixType tmpgrad = tmpKernelObj->ConvolveGradient(m_PosT[s], dim);

		dXi6 += tmpgrad;
	}

	for (unsigned int i = 0; i < numCP; i++)
	{
		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		dXi6.set_row(i, dXi6.get_row(i) - gradMom_i * m_XiMomT[s].get_row(i));
	}

	// std::cout << "term 6: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();

	// Term 4
	/*
	MatrixType dXi4(numCP, Dimension, 0);

	// dXi4 first term
	// timer.Start();
	tmpKernelObj->SetSources(PositionsT[s]);
	tmpKernelObj->SetWeights(MomentasT[s]);
	std::vector< std::vector<MatrixType> > HessA = tmpKernelObj->ConvolveHessian(PositionsT[s]);
	for (unsigned int i = 0; i < numCP; i++)
	{
	VectorType XiMomi = m_XiMomT[s].get_row(i);
	MatrixType Auxi(Dimension, Dimension, 0.0);
	for (unsigned int dim = 0; dim < Dimension; dim++)
	Auxi -= HessA[i][dim] * MomentasT[s](i,dim);

	dXi4.set_row(i, Auxi * XiMomi);
	}
	// dXi4 second term
	std::vector< MatrixType > H(numCP);
	MatrixType Zeros(Dimension, Dimension, 0.0);
	for (int i = 0; i < numCP; i++)
	H[i] = Zeros;

	for (unsigned r = 0; r < Dimension; r++)
	{
	MatrixType Wd(numCP,Dimension,0.0);
	for (int i = 0; i < numCP; i++)
	Wd.set_row(i, MomentasT[s].get_row(i) * m_XiMomT[s](i,r) );

	tmpKernelObj->SetWeights(Wd);
	MatrixType Hrr = tmpKernelObj->ConvolveHessian(PositionsT[s],r,r);
	for (int i = 0; i < numCP; i ++)
	for (unsigned int q = 0; q < Dimension; q++)
	H[i](r,q) += Hrr(i,q);

	for (unsigned q = r+1; q < Dimension; q++)
	{
	MatrixType W(numCP, 2*Dimension, 0.0);
	for (int i = 0; i < numCP; i++)
	for (unsigned int d = 0; d < Dimension; d++)
	{
	W(i,d) = Wd(i,d);
	W(i, d + Dimension) = MomentasT[s](i,d) * m_XiMomT[s](i,q);
	}

	tmpKernelObj->SetWeights(W);
	MatrixType Hrq = tmpKernelObj->ConvolveHessian(PositionsT[s],r,q);
	for (int i = 0; i < numCP; i ++)
	for (unsigned int d = 0; d < Dimension; d++)
	{
	H[i](q,d) += Hrq(i,d);
	H[i](r,d) += Hrq(i,d + Dimension);
	}
	}		
	}

	for (int i = 0; i < numCP; i++)
	dXi4.set_row(i, dXi4.get_row(i) + H[i] * MomentasT[s].get_row(i));


	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();
	*/

	// Alexandre :
	// timer.Start();
	tmpKernelObj->SetSources(m_PosT[s]);
	tmpKernelObj->SetWeights(m_MomT[s]);
	MatrixType dXi4 = tmpKernelObj->ConvolveSpecialHessian(m_XiMomT[s]);
	assert(dXi4.rows() == numCP);
	assert(dXi4.columns() == Dimension);
	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;



	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
	// timer.Start();

	dPos = dXi1 + dXi3 + dXi4;
	dMom = dXi2 + dXi5 + dXi6;


	// std::cout << "dXi1 = "<< std::endl;
	// std::cout << dXi1 << std::endl;
	// std::cout << std::endl;
	// std::cout << "dXi2 = " << std::endl;
	// std::cout << dXi2 << std::endl;
	// std::cout << std::endl;
	// std::cout << "dXi3 = " << std::endl;
	// std::cout << dXi3 << std::endl;
	// std::cout << std::endl;
	// std::cout << "dXi4 = " << std::endl;
	// std::cout << dXi4 << std::endl;
	// std::cout << std::endl;
	// std::cout << "dXi5 = " << std::endl;
	// std::cout << dXi5 << std::endl;
	// std::cout << std::endl;
	// std::cout << "dXi6 = " << std::endl;
	// std::cout << dXi6 << std::endl;
	// std::cout << std::endl;


}





// template <class TScalar, unsigned int Dimension>
// void
// AdjointEquationsIntegrator<TScalar, Dimension>
// ::ComputeUpdateAt(unsigned int s, MatrixType& dPos, MatrixType &dMom)
// {
// 
// 	MatrixList PositionsT = m_Def->GetTrajectoryPositions();
// 	MatrixList MomentasT = m_Def->GetTrajectoryMomentas();
// 
// 	long numCP = PositionsT[0].rows();
// 
// 	KernelType* momXiPosKernelObj = m_KernelObj1;
// 	KernelType* etaKernelObj = m_KernelObj2;
// 	KernelType* tmpKernelObj = m_KernelObj3;
// 
// 	// Term 1
// 	// SimpleTimer timer;
// 	etaKernelObj->SetSources(m_PointsT[s]);
// 	etaKernelObj->SetWeights(m_VectorsT[s]);
// 
// 	/*
// 	// SimpleTimer subtimer;
// 	MatrixList gradEta = etaKernelObj->ConvolveGradient(PositionsT[s]);
// 
// 	// std::cout << "\tsub-term 1.1: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// subtimer.Start();
// 
// 
// 	MatrixType dXi1(numCP, Dimension, 0);
// 	for (unsigned int i = 0; i < numCP; i++)
// 	{
// 	dXi1.set_row(i, gradEta[i].transpose() * MomentasT[s].get_row(i));
// 	}
// 
// 	// std::cout << "\tsub-term 1.2: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	*/
// 
// 	MatrixType dXi1 = etaKernelObj->ConvolveGradient(PositionsT[s], MomentasT[s]);
// 	assert(dXi1.rows() == numCP);
// 	assert(dXi1.columns() == Dimension);
// 
// 	// std::cout << "term 1: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 	// Term 2
// 	MatrixType dXi2 = etaKernelObj->Convolve(PositionsT[s]);
// 
// 	// std::cout << "term 2: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 	// Term 3
// 	MatrixType AXiPos(numCP, Dimension*2, 0);
// 	AXiPos.set_columns(0, MomentasT[s]);
// 	AXiPos.set_columns(Dimension, m_XiPosT[s]);
// 
// 	momXiPosKernelObj->SetSources(PositionsT[s]);
// 	momXiPosKernelObj->SetWeights(AXiPos);
// 
// 	// subtimer.Start();
// 	MatrixType kAXiPos = momXiPosKernelObj->Convolve(PositionsT[s]);
// 
// 	MatrixList gradAXiPos =	momXiPosKernelObj->ConvolveGradient(PositionsT[s]);
// 
// 	// std::cout << "\tsub term 3.1: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// subtimer.Start();
// 
// 	MatrixType dXi3(numCP, Dimension, 0);
// 	for (unsigned int i = 0; i < numCP; i++)
// 	{
// 		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
// 		MatrixType gradXiPos_i = gradAXiPos[i].get_n_rows(Dimension, Dimension);
// 		dXi3.set_row(i,gradMom_i.transpose() * m_XiPosT[s].get_row(i) + gradXiPos_i.transpose() * MomentasT[s].get_row(i));
// 	}
// 
// 	// std::cout << "\tsub term 3.2: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// 
// 	// std::cout << "term 3: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 	// Term 5
// 	MatrixType dXi5 = kAXiPos.get_n_columns(Dimension, Dimension);
// 
// 	// std::cout << "term 5: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 	// Term 6
// 	MatrixType dXi6(numCP, Dimension, 0);
// 
// 	// MARCEL
// 	// grad(p) = sum_q grad_q xi_mom(q) * mom(p)
// 	//   for (unsigned int p = 0; p < Dimension; p++)
// 	//   {
// 	//     MatrixType W_pq(numCP, Dimension, 0);
// 	//     for (unsigned int q = 0; q < Dimension; q++)
// 	//     {
// 	//       for (unsigned int i = 0; i < numCP; i++)
// 	//         W_pq(i, q) = MomentasT[s](i, p) * m_XiMomT[s](i, q);
// 	//     }
// 	//
// 	//     tmpKernelObj->SetSources(PositionsT[s]);
// 	//     tmpKernelObj->SetWeights(W_pq);
// 	//
// 	//     // Sum of grad_q xi_mom(q) * mom(p) so use diagonal elements
// 	// /*
// 	//     MatrixList gradTmp = tmpKernelObj->ConvolveGradient(PositionsT[s]);
// 	//     for (unsigned int i = 0; i < numCP; i++)
// 	//       for (unsigned int q = 0; q < Dimension; q++)
// 	//         dXi6(i, p) += gradTmp[i](q, q);
// 	// */
// 	//
// 	//     for (unsigned int q = 0; q < Dimension; q++)
// 	//     {
// 	// 			subtimer.Start();
// 	//       VectorType grad_qq =
// 	//         tmpKernelObj->ConvolveGradient(PositionsT[s], q, q);
// 	// 			std::cout << "\tsub term: " << subtimer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	//       dXi6.set_column(p, dXi6.get_column(p) + grad_qq);
// 	//     }
// 	//   }
// 
// 	//STANLEY
// 	for (unsigned int dim = 0; dim < Dimension; dim++)
// 	{
// 		MatrixType W(numCP, Dimension, 0.0);
// 		for (unsigned int i = 0; i < numCP; i++)
// 			W.set_row(i, m_XiMomT[s](i,dim) * MomentasT[s].get_row(i));
// 
// 		tmpKernelObj->SetSources(PositionsT[s]);
// 		tmpKernelObj->SetWeights(W);
// 		MatrixType tmpgrad = tmpKernelObj->ConvolveGradient(PositionsT[s], dim);
// 
// 		dXi6 += tmpgrad;
// 	}
// 
// 	for (unsigned int i = 0; i < numCP; i++)
// 	{
// 		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
// 		dXi6.set_row(i, dXi6.get_row(i) - gradMom_i * m_XiMomT[s].get_row(i));
// 	}
// 
// 	// std::cout << "term 6: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 
// 	// Term 4
// 	/*
// 	MatrixType dXi4(numCP, Dimension, 0);
// 
// 	// MARCEL
// 	//   for (unsigned int i = 0; i < numCP; i++)
// 	//   {
// 	//     VectorType ai = MomentasT[s].get_row(i);
// 	//     VectorType ei = m_XiMomT[s].get_row(i);
// 	// 
// 	//     MatrixType xi(1, Dimension, 0.0);
// 	//     for (unsigned int c = 0; c < Dimension; c++)
// 	//       xi(0, c) = PositionsT[s](i, c);
// 	// 
// 	//     MatrixType Hq(Dimension, Dimension, 0.0);
// 	//     for (unsigned int q = 0; q < Dimension; q++)
// 	//     {
// 	//       for (unsigned int r = 0; r < Dimension; r++)
// 	//         for (unsigned int c = 0; c < Dimension; c++)
// 	//         {
// 	//           VectorType h_rc = momXiPosKernelObj->ConvolveHessian(
// 	//             xi, q, r, c);
// 	//           Hq(c, r) = h_rc[0];
// 	//         }
// 	// 			//MARCEL
// 	// 			// dXi4.set_row(i, (Hq * ei) * (ai[q] * -1.0));
// 	//       	// STANLEY
// 	// 			dXi4.set_row(i, dXi4.get_row(i) - (Hq * ei) * ai[q]);
// 	//     }
// 	//   }
// 	// 
// 	// // DEBUG
// 	// // MatrixType dXi4_1(numCP, Dimension, 0);
// 	// // dXi4_1 = dXi4;
// 	// // std::cout << "dXi4 intermediate = " << std::endl << dXi4 << std::endl;
// 	// 
// 	//   MatrixType XqAc(numCP, Dimension*Dimension, 0);
// 	//   for (unsigned int i = 0; i < numCP; i++)
// 	//   {
// 	//     for (unsigned int q = 0; q < Dimension; q++)
// 	//       for (unsigned int c = 0; c < Dimension; c++)
// 	//         XqAc(i, q*Dimension+c) = m_XiMomT[s](i, q) * MomentasT[s](i, c);
// 	//   }
// 	// 
// 	//   tmpKernelObj->SetSources(PositionsT[s]);
// 	//   tmpKernelObj->SetWeights(XqAc);
// 	// 
// 	//   // dXi4 += M * ai
// 	//   // M(r, c) = \sum_j \del_{q,r} K xi_q a_c
// 	//   for (unsigned int i = 0; i < numCP; i++)
// 	//   {
// 	//     MatrixType xi(1, Dimension, 0.0);
// 	//     for (unsigned int c = 0; c < Dimension; c++)
// 	//       xi(0, c) = PositionsT[s](i, c);
// 	// 
// 	//     MatrixType M_i(Dimension, Dimension, 0.0);
// 	//     for (unsigned int q = 0; q < Dimension; q++)
// 	//     {
// 	//       for (unsigned int r = 0; r < Dimension; r++)
// 	//         for (unsigned int c = 0; c < Dimension; c++)
// 	//         {
// 	// 				VectorType h_qr = tmpKernelObj->ConvolveHessian(
// 	//             xi, q*Dimension+c, r, q);
// 	//           M_i(r, c) += h_qr[0];
// 	//         }
// 	//     }
// 	// 
// 	//     VectorType ai = MomentasT[s].get_row(i);
// 	//     dXi4.set_row(i, dXi4.get_row(i) + M_i*ai);
// 	//   }
// 
// 
// 
// 	// STANLEY
// 	// dXi4 first term
// 	// timer.Start();
// 	tmpKernelObj->SetSources(PositionsT[s]);
// 	tmpKernelObj->SetWeights(MomentasT[s]);
// 	std::vector< std::vector<MatrixType> > HessA = tmpKernelObj->ConvolveHessian(PositionsT[s]);
// 	for (unsigned int i = 0; i < numCP; i++)
// 	{
// 	VectorType XiMomi = m_XiMomT[s].get_row(i);
// 	MatrixType Auxi(Dimension, Dimension, 0.0);
// 	for (unsigned int dim = 0; dim < Dimension; dim++)
// 	Auxi -= HessA[i][dim] * MomentasT[s](i,dim);
// 
// 	dXi4.set_row(i, Auxi * XiMomi);
// 	}
// 	// dXi4 second term
// 	std::vector< MatrixType > H(numCP);
// 	MatrixType Zeros(Dimension, Dimension, 0.0);
// 	for (int i = 0; i < numCP; i++)
// 	H[i] = Zeros;
// 
// 	for (unsigned r = 0; r < Dimension; r++)
// 	{
// 	MatrixType Wd(numCP,Dimension,0.0);
// 	for (int i = 0; i < numCP; i++)
// 	Wd.set_row(i, MomentasT[s].get_row(i) * m_XiMomT[s](i,r) );
// 
// 	tmpKernelObj->SetWeights(Wd);
// 	MatrixType Hrr = tmpKernelObj->ConvolveHessian(PositionsT[s],r,r);
// 	for (int i = 0; i < numCP; i ++)
// 	for (unsigned int q = 0; q < Dimension; q++)
// 	H[i](r,q) += Hrr(i,q);
// 
// 	for (unsigned q = r+1; q < Dimension; q++)
// 	{
// 	MatrixType W(numCP, 2*Dimension, 0.0);
// 	for (int i = 0; i < numCP; i++)
// 	for (unsigned int d = 0; d < Dimension; d++)
// 	{
// 	W(i,d) = Wd(i,d);
// 	W(i, d + Dimension) = MomentasT[s](i,d) * m_XiMomT[s](i,q);
// 	}
// 
// 	tmpKernelObj->SetWeights(W);
// 	MatrixType Hrq = tmpKernelObj->ConvolveHessian(PositionsT[s],r,q);
// 	for (int i = 0; i < numCP; i ++)
// 	for (unsigned int d = 0; d < Dimension; d++)
// 	{
// 	H[i](q,d) += Hrq(i,d);
// 	H[i](r,d) += Hrq(i,d + Dimension);
// 	}
// 	}		
// 	}
// 
// 	for (int i = 0; i < numCP; i++)
// 	dXi4.set_row(i, dXi4.get_row(i) + H[i] * MomentasT[s].get_row(i));
// 
// 
// 	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 	*/
// 
// 	// Alexandre :
// 	// timer.Start();
// 	tmpKernelObj->SetSources(PositionsT[s]);
// 	tmpKernelObj->SetWeights(MomentasT[s]);
// 	MatrixType dXi4 = tmpKernelObj->ConvolveSpecialHessian(m_XiMomT[s]);
// 	assert(dXi4.rows() == numCP);
// 	assert(dXi4.columns() == Dimension);
// 	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 
// 
// 
// 	// std::cout << "term 4: " << timer.GetElapsedCPUTimeInSecondsOnly() << std::endl;
// 	// timer.Start();
// 
// 	dPos = dXi1 + dXi3 + dXi4;
// 	dMom = dXi2 + dXi5 + dXi6;
// 
// 
// 	// std::cout << "dXi1 = "<< std::endl;
// 	// std::cout << dXi1 << std::endl;
// 	// std::cout << std::endl;
// 	// std::cout << "dXi2 = " << std::endl;
// 	// std::cout << dXi2 << std::endl;
// 	// std::cout << std::endl;
// 	// std::cout << "dXi3 = " << std::endl;
// 	// std::cout << dXi3 << std::endl;
// 	// std::cout << std::endl;
// 	// std::cout << "dXi4 = " << std::endl;
// 	// std::cout << dXi4 << std::endl;
// 	// std::cout << std::endl;
// 	// std::cout << "dXi5 = " << std::endl;
// 	// std::cout << dXi5 << std::endl;
// 	// std::cout << std::endl;
// 	// std::cout << "dXi6 = " << std::endl;
// 	// std::cout << dXi6 << std::endl;
// 	// std::cout << std::endl;
// 
// 
// }



#endif /* _AdjointEquationsIntegrator_txx */

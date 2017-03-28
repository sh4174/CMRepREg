/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _estimateGeodesicRegression_txx
#define _estimateGeodesicRegression_txx

#include "itkVersion.h"
#include "vtkVersion.h"

#include "LinearAlgebra.h"

#include "KernelFactory.h"

#include "Regression.h"

#include "RegressionEstimator.h"

#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "DeformableObject.h"

#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

#include "DeformableObjectReader.h"

#include "writeMatrixDLM.txx"

#include "SimpleTimer.h"

#if ITK_VERSION_MAJOR >= 4
    #include "itkFFTWGlobalConfiguration.h"
#endif

#include "itksys/SystemTools.hxx"


#include <cstring>
#include <iostream>
#include <sstream>

/** 
 estimateGeodesicRegression()
 Manages the estimation of sparse geodesic regression.  Reads and organizes input and handles output of results.
 
 @param[in] paramDiffeos An object containing diffeo parameters
 @param[in] numObservations The number of observations per object
 @param[in] numObjects The number of unique objects (eg 2 for multiple observations of left/right hemisphere)
 @param[in] paramObjectsList A vector of objects containing parameter values for each object
 @param[in] templatefnList A vector of vectors of string paths to data defining the inital template for each observation and object
 @param[in] observationfnList A vector of vectors of string paths to data defining the observations for each observation and object
 @param[in] observationTimesList A vector of vectors of times for each observation and object
 */
template<unsigned int Dimension>
void estimateGeodesicRegression(SparseDiffeoParameters* paramDiffeos, int numObservations, int numObjects, 
				std::vector<DeformableObjectParameters::Pointer>& paramObjectsList, const std::vector<char*>& templatefnList, 
				const std::vector<std::vector<char*> >& observationfnList, const std::vector<std::vector<double> >& observationTimesList)
{
	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

#ifdef USE_DOUBLE_PRECISION
	typedef double TScalar;
	std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
	typedef float TScalar;
	std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

	std::cout << "Sparse diffeomorphic atlas estimation\n===" << std::endl;
	std::cout << "ITK version " << itk::Version::GetITKVersion() << std::endl;
	std::cout << "VTK version " << vtkVersion::GetVTKVersion() << std::endl;
	std::cout << "Number of objects: " << numObjects << std::endl;
	std::cout << "Number of observations: " << numObservations << std::endl;
	std::cout << "\n===\n" << std::endl;
	std::cout << "Deformation Parameters:" << std::endl;
	paramDiffeos->PrintSelf(std::cout);
	std::cout << "\n===\n" << std::endl;
	for (int i = 0; i < numObjects; i++)
	{
		std::cout << "Object " << i << ": " << std::endl;
		std::cout << "template file: " << templatefnList[i] << std::endl;
		paramObjectsList[i]->PrintSelf(std::cout);
		std::cout << "\n===\n" << std::endl;
	}

#if ITK_VERSION_MAJOR >= 4
	itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;

	typedef Diffeos<TScalar, Dimension> DiffeosType;

	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	typedef std::vector<DeformableObjectType*> DeformableObjectList;
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

	typedef Regression<TScalar, Dimension> RegressionType;
	
	typedef RegressionEstimator< TScalar, Dimension> RegressionEstimatorType;

	// Set up the kernel factory - can't use any Update() before setting this up!
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
	kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

	// Create the deformation object
	DiffeosType* def = new DiffeosType();
	def->SetT0(paramDiffeos->GetT0());
	def->SetTN(paramDiffeos->GetTN());
	def->SetKernelWidth(paramDiffeos->GetKernelWidth());
	def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
	def->SetPaddingFactor( paramDiffeos->GetP3MPaddingFactor() ); // to define the bounding box
	def->UseImprovedEuler();

	if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0)
	{
		def->SetKernelType(P3M);
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "fgt") == 0)
	{
		def->SetKernelType(FGT);
	}
#ifdef USE_CUDA
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0)
	{
		def->SetKernelType(CUDAExact);
	}
#endif
	else
	{
		if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
			std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
		def->SetKernelType(Exact);
	}

	// Create template and target objects
	DeformableObjectList templateObjectList(numObjects);
	std::vector<DeformableObjectList> targetObjectList(numObservations);
	
	// For CMRep Regression Only
	std::vector<DeformableObjectList> targetCMRepObjectList(numObservations);

	for (int s = 0; s < numObservations; s++)
	{
		DeformableObjectList Aux(numObjects);
		targetObjectList[s] = Aux;

		DeformableObjectList Aux2( numObjects );
		targetCMRepObjectList[ s ] = Aux2;
	}

	// Scan the list of objects
	typedef DeformableObjectReader<TScalar,Dimension> ObjectReaderType;
	
	std::vector<unsigned int> timeIndices;
	timeIndices.resize(numObservations);

	for (int i = 0; i < numObjects; i++)
	{
		ObjectReaderType* reader = new ObjectReaderType();
		reader->SetObjectParameters(paramObjectsList[i]);

		reader->SetFileName(templatefnList[i]);
		reader->Update();

		templateObjectList[i] = reader->GetOutput();

		for (int s = 0; s < numObservations; s++)
		{
			cout << paramObjectsList[ i ]->GetDeformableObjectType() << endl;
			ObjectReaderType* reader = new ObjectReaderType();
			
			if( paramObjectsList[ i ]->GetDeformableObjectType() == "CMRep" )
			{
				ifstream datFile( observationfnList[ i ][ s ] );
				string vtkBndPath;
				getline( datFile, vtkBndPath );

				string vtkMedPath;
				getline( datFile, vtkMedPath );

				paramObjectsList[ i ]->SetDeformableObjectType( "NonOrientedSurfaceMesh" );
				reader->SetObjectParameters(paramObjectsList[i]);
				reader->SetFileName(const_cast<char*> ( vtkBndPath.c_str() ));
				reader->Update();
				targetObjectList[s][i] = reader->GetOutput();


				ObjectReaderType* CMRepReader = new ObjectReaderType();
				CMRepReader->SetObjectParameters(paramObjectsList[i]);
				CMRepReader->SetFileName( const_cast<char*> ( vtkMedPath.c_str() ) );
				CMRepReader->Update();
								
				targetCMRepObjectList[ s ][ i ] = CMRepReader->GetOutput();
				paramObjectsList[ i ]->SetDeformableObjectType( "CMRep" );
			}
			else
			{
				reader->SetObjectParameters(paramObjectsList[i]);
				reader->SetFileName(observationfnList[i][s]);
				reader->Update();
				targetObjectList[s][i] = reader->GetOutput();
			}

			// Get time point information
			TScalar timept = observationTimesList[i][s];

			// Figure out which time index the time point corresponds to
			unsigned int T = paramDiffeos->GetNumberOfTimePoints();
			TScalar tn = paramDiffeos->GetTN();
			TScalar t0 = paramDiffeos->GetT0();
			int timeIndex = int((timept-t0)*((T-1)/(tn-t0))+0.5f);

			std::cout<<"Time = "<<timept<<"   Time index = "<<timeIndex<<"\n";
			timeIndices[s] = timeIndex;
		}

	}
	
	// Creating multi-objects for template and targets
	typename std::vector< DeformableMultiObjectType* > target(numObservations);
	for (unsigned int s = 0; s < numObservations; s++)
	{
		target[s] = new DeformableMultiObjectType();
		target[s]->SetObjectList(targetObjectList[s]);
		if ( paramObjectsList[ 0 ]->GetDeformableObjectType() == "CMRep" )
		{
			target[ s ]->SetCMRepObjectList( targetCMRepObjectList[s] );
		}

		target[s]->Update();
	}

	DeformableMultiObjectType* templateObjects = new DeformableMultiObjectType();
	templateObjects->SetObjectList(templateObjectList);
	templateObjects->Update();

	RegressionType* regression = new RegressionType();
	RegressionEstimatorType* regressionBuilder = new RegressionEstimatorType();

	regressionBuilder->SetTimeIndices(timeIndices);
	regression->SetSparsityPrior(paramDiffeos->GetSparsityPrior());
	regressionBuilder->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
	
	if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientdescent") == 0)
	{
		regressionBuilder->SetGradientDescent();
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "ista") == 0)
	{
		regressionBuilder->SetISTA();
	}
	else
	{
		if ( (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "f_ista") != 0) &&
			(itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fista") != 0) )
		std::cerr << "Unknown optimization method type: defaulting to FISTA" << std::endl;
		
		regressionBuilder->SetFISTA();
	}

	regression->SetTemplate(templateObjects);

	std::string CP_fn = paramDiffeos->GetInitialCPPosition_fn();
	regression->SetControlPoints(CP_fn);	
	regression->SetCPSpacing( paramDiffeos->GetInitialCPSpacing() );
	regression->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
	regression->SetDiffeos(def);
	
	regression->Update();

	regressionBuilder->SetTargetList(target);
	regressionBuilder->SetRegression(regression);
	paramDiffeos->FreezeCP()?regressionBuilder->SetFreezeCP():regressionBuilder->UnsetFreezeCP();
	paramDiffeos->FreezeTemplate()?regressionBuilder->SetFreezeTemplate():regressionBuilder->UnsetFreezeTemplate();
	regressionBuilder->SetMaxIterations(paramDiffeos->GetMaxIterations());
	regressionBuilder->SetSaveEveryNIters(paramDiffeos->GetSaveEveryNIters());
	regressionBuilder->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
	regressionBuilder->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
	regressionBuilder->SetInitialStepMultiplier(paramDiffeos->GetInitialStepMultiplier());

	VectorType DataSigmaSquared(numObjects, 0.0);
	for (unsigned int i = 0; i < numObjects; i++)
	{
		TScalar DSS = paramObjectsList[i]->GetDataSigma();
		DataSigmaSquared(i) = DSS*DSS;
	}
	regression->SetDataSigmaSquared(DataSigmaSquared);

	std::string RegressionName = paramDiffeos->GetAtlasName_fn();
	if (RegressionName.empty())
		RegressionName = "Regression";

	std::vector< std::string > objectsName(numObjects);
	std::vector< std::string > objectsNameExtension(numObjects);

	for (int i = 0; i < numObjects; i++)
	{
		std::string sourcefn;
		sourcefn.assign(templatefnList[i]);
		int length = sourcefn.size();
		int index = length - 1;
		while( (sourcefn[index] != '.') && (index > 0) )
			index--;

		if (index == 0)
			throw std::runtime_error("template file name has no extension");

		int index2 = index - 1;
		while ( (sourcefn[index2] != '/' && (index2 > 0)) )
			index2--;
	
		if (index2 > 0)
			index2 += 1;
				
		objectsName[i] = sourcefn.substr(index2,index-index2);
		objectsNameExtension[i] = sourcefn.substr(index,length-index);
	}
	
	regressionBuilder->SetTemplateObjectsName(objectsName, objectsNameExtension);
	regressionBuilder->SetRegressionName(RegressionName);

	// create timer
	SimpleTimer timer;

	regressionBuilder->Update();

	timer.Stop();

	std::cout << "Regression estimation took "
			<< timer.GetElapsedHours() << " hours, "
			<< timer.GetElapsedMinutes() << " minutes, "
			<< timer.GetElapsedSeconds() << " seconds"
			<< std::endl;


	std::cout << "Write output files..." << std::endl;
		
	regressionBuilder->WriteOutput();
	std::cout << "... done" << std::endl;

	delete regression;
	delete regressionBuilder;

}

#endif

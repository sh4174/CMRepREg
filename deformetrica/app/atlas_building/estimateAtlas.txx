#ifndef _estimateAtlas_txx
#define _estimateAtlas_txx

#include "itkVersion.h"
#include "vtkVersion.h"

#include "LinearAlgebra.h"

#include "KernelFactory.h"

#include "Atlas.h"
#include "DeterministicAtlas.h"
#include "BayesianAtlas.h"


#include "AbstractAtlasEstimator.h"
#include "DeterministicAtlasEstimator.h"
#include "BayesianAtlasEstimator.h"
#include "MultiClassBayesianAtlasEstimator.h"

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

template <unsigned int Dimension>
void estimateAtlas(
		SparseDiffeoParameters* paramDiffeos,
		int numSubjects,
		int numObjects,
		std::vector<DeformableObjectParameters::Pointer>& paramObjectsList,
		const std::vector<char*>& templatefnList,
		const std::vector< std::vector<char*> >& subjectfnList)
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
	std::cout << "Number of subjects: " << numSubjects << std::endl;
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

	typedef Atlas<TScalar, Dimension> AtlasType;
	typedef DeterministicAtlas<TScalar, Dimension> DeterministicAtlasType;
	typedef BayesianAtlas<TScalar, Dimension> BayesianAtlasType;

	typedef AbstractAtlasEstimator< TScalar, Dimension> AtlasEstimatorType;
	typedef DeterministicAtlasEstimator<TScalar, Dimension> DeterministicAtlasEstimatorType;
	typedef BayesianAtlasEstimator<TScalar, Dimension> BayesianAtlasEstimatorType;
    typedef MultiClassBayesianAtlasEstimator<TScalar, Dimension> MultiClassBayesianAtlasEstimatorType;
    
	

	// Set up the kernel factory - can't use any Update() before setting this up!
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
	kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());
	

	// create template and target objects
	DeformableObjectList templateObjectList(numObjects);
	std::vector<DeformableObjectList> targetObjectList(numSubjects);
	for (int s = 0; s < numSubjects; s++)
	{
		DeformableObjectList Aux(numObjects);
		targetObjectList[s] = Aux;
	}

	// // PIETRO
	// VectorType DataSigmaSquared_HyperParameter(numObjects, 0.0);
	// VectorType DataSigmaSquared_Prior(numObjects, 0.0);
	// VectorType DataSigmaSquared(numObjects, 0.0);

	// scan the list of objects
	typedef DeformableObjectReader<TScalar,Dimension> ObjectReaderType;

	for (int i = 0; i < numObjects; i++)
	{
		for (int s = 0; s < numSubjects; s++)
		{
			ObjectReaderType* reader = new ObjectReaderType();
			reader->SetObjectParameters(paramObjectsList[i]);

			reader->SetFileName(subjectfnList[i][s]);
			reader->Update();

			targetObjectList[s][i] = reader->GetOutput();
		}
				
		ObjectReaderType* reader = new ObjectReaderType();
		reader->SetObjectParameters(paramObjectsList[i]);

		reader->SetFileName(templatefnList[i]);
		reader->Update();

		templateObjectList[i] = reader->GetOutput();

		// // PIETRO
		// if (paramDiffeos->BayesianFramework())
		// {
		// 	DataSigmaSquared_HyperParameter[i] = (TScalar) (reader->GetWeightData());
		// 	DataSigmaSquared_Prior[i] = (TScalar) (reader->GetPriorData());
		//
		// 	if (DataSigmaSquared_HyperParameter[i]<=0)
		// 	{
		// 		std::cout << "WARNING! The data weight is equal to 0 or negative!" << std::endl;
		// 		throw std::runtime_error("If you want to use the Bayesian Framework it is not a good idea to set Weight Data to zero...");
		// 	}
		//
		// 	if (DataSigmaSquared_Prior[i]<=0)
		// 	{
		// 		std::cout << "WARNING! The data prior is equal to 0 or negative!" << std::endl;
		// 		throw std::runtime_error("If you want to use the Bayesian Framework it is not a good idea to set Prior Data to zero...");
		// 	}
		// }
		// else
		// {
		// 	DataSigmaSquared[i] = (TScalar) (reader->GetDataSigmaSquared()); // so it is easier to check that the data sigma squared is well computed in the bayesian framework
		//
		// 	if (DataSigmaSquared[i]<=0)
		// 	{
		// 		std::cout << "WARNING! The data sigma squared is equal to 0 or negative!" << std::endl;
		// 		throw std::runtime_error("If you want to use the Deterministic Framework it is not a good idea to set Data Sigma Squared to zero...");
		// 	}
		// }

	} // end scanning objects
	
	// creating multi-objects for template and targets
	typename std::vector< DeformableMultiObjectType* > target(numSubjects);
	for (unsigned int s = 0; s < numSubjects; s++)
	{
		target[s] = new DeformableMultiObjectType();
		target[s]->SetObjectList(targetObjectList[s]);
		target[s]->Update();
	}

	DeformableMultiObjectType* templateObjects = new DeformableMultiObjectType();
	templateObjects->SetObjectList(templateObjectList);
	templateObjects->Update();
	
	// to save memory space
	for (unsigned int i = 0; i < numObjects; i++)
	{
		for (unsigned int s = 0; s < numSubjects; s++)
			delete targetObjectList[s][i];
		
		delete templateObjectList[i];
	}
	
	// warnings
    bool isOfBayesianType = ((itksys::SystemTools::Strucmp(paramDiffeos->GetAtlasType().c_str(),"BAYESIAN") == 0)||
                             (itksys::SystemTools::Strucmp(paramDiffeos->GetAtlasType().c_str(),"MULTICLASSBAYESIAN") == 0));
	if ( (numSubjects == 1) && (!paramDiffeos->FreezeTemplate()) )
	{
		std::cout << "Warning: for what is obviously a matching problem, source object will be froozen" << std::endl;
		paramDiffeos->SetFreezeTemplate();
	}
	else if ( (paramDiffeos->FreezeTemplate()) && isOfBayesianType)
	{
		std::cout << "Warning: template freeze is not implemented for Bayesian framework: template will be updated!" << std::endl;
		paramDiffeos->UnsetFreezeTemplate();
	}
	
	// create the deformation object
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
                std::cout << "Working on Cuda Exact" << std::endl;
		def->SetKernelType(CUDAExact);
	}
#endif
	else
	{
		if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
			std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
		def->SetKernelType(Exact);
	}

	
	// create the atlas, whose parameters will be optimized, as well as the optimization scheme
	AtlasType* atlas;
	AtlasEstimatorType* atlasbuilder;
	
	if (isOfBayesianType)
	{
		BayesianAtlasType* BayesianAtlas = new BayesianAtlasType();
        
		BayesianAtlas->SetCovarianceMomenta_HyperParameter(
			numSubjects * paramDiffeos->GetCovarianceMomenta_Normalized_Hyperparameter() );
		BayesianAtlas->SetCovarianceMomenta_Prior_Inverse( paramDiffeos->GetCovarianceMomenta_Prior_Inverse_fn() );
		
		BayesianAtlas->SetNoiseDimension( templateObjects->GetDimensionOfDiscretizedObjects() );

		VectorType DataSigmaSquared_Hyperparameter(numObjects, 0.0);
		VectorType DataSigmaSquared_Prior(numObjects, 0.0);
		std::vector<int> dim = templateObjects->GetDimensionOfDiscretizedObjects();
		for (unsigned int i = 0; i < numObjects; i++)
		{
			DataSigmaSquared_Hyperparameter(i) = paramObjectsList[i]->GetDataSigma_Normalized_Hyperparameter() * dim[i] * numSubjects;
			DataSigmaSquared_Prior(i) = paramObjectsList[i]->GetDataSigma_Prior();
			
			// std::cout << "Object " << i << std::endl;
			// std::cout << "Number points grid object: " << dim[i] << std::endl;
			// DataSigmaSquared_HyperParameter_Normalized[i] = DataSigmaSquared_HyperParameter[i]*dim[i]*numSubjects;
			// std::cout << "Data Sigma Weight is: " << DataSigmaSquared_HyperParameter[i] << " and it is normalized to: " << DataSigmaSquared_HyperParameter_Normalized[i] << std::endl;
			// std::cout << "Data Sigma Prior is: " << DataSigmaSquared_Prior[i] << " and it is not normalized " << std::endl;
			// DataSigmaSquared_Prior_Normalized[i]=DataSigmaSquared_Prior[i];
		}
		BayesianAtlas->SetDataSigmaSquared_HyperParameter(DataSigmaSquared_Hyperparameter);
		BayesianAtlas->SetDataSigmaSquared_Prior(DataSigmaSquared_Prior);
        
        if (itksys::SystemTools::Strucmp(paramDiffeos->GetAtlasType().c_str(),"BAYESIAN") == 0)
        {
            BayesianAtlasEstimatorType* BayesianAtlasBuilder = new BayesianAtlasEstimatorType();
            
            BayesianAtlasBuilder->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
            atlasbuilder = BayesianAtlasBuilder;
        }
        else
        {
            MultiClassBayesianAtlasEstimatorType* MultiClassBayesianAtlasBuilder = new MultiClassBayesianAtlasEstimatorType();
            MultiClassBayesianAtlasBuilder->SetMaximumNumberOfClasses(paramDiffeos->GetMaximumNumberOfClasses());
            atlasbuilder = MultiClassBayesianAtlasBuilder;
        }

		atlas = BayesianAtlas;
		
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetAtlasType().c_str(),"DETERMINISTIC") == 0)
	{
		DeterministicAtlasType* DeterministicAtlas = new DeterministicAtlasType();
		DeterministicAtlasEstimatorType* DeterministicAtlasBuilder = new DeterministicAtlasEstimatorType();
		
		DeterministicAtlas->SetSparsityPrior(paramDiffeos->GetSparsityPrior());
		DeterministicAtlasBuilder->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
		
		if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientdescent") == 0)
		{
			DeterministicAtlasBuilder->SetGradientDescent();
		}
		else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "ista") == 0)
		{
			DeterministicAtlasBuilder->SetISTA();
		}
		else
		{
			if ( (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "f_ista") != 0) &&
				(itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fista") != 0) )
			std::cerr << "Unknown optimization method type: defaulting to FISTA" << std::endl;
			
			DeterministicAtlasBuilder->SetFISTA();
		}
		
		std::string CMI_fn = paramDiffeos->GetCovarianceMomentaInverse_fn();
		if (CMI_fn.empty())
			DeterministicAtlas->SetRKHSNormForRegularization();
		else
			DeterministicAtlas->SetCovarianceMomentaInverse( CMI_fn );
				
		VectorType DataSigmaSquared(numObjects, 0.0);
		for (unsigned int i = 0; i < numObjects; i++)
		{
			TScalar DSS = paramObjectsList[i]->GetDataSigma();
			DataSigmaSquared(i) = DSS*DSS;
		}
		DeterministicAtlas->SetDataSigmaSquared(DataSigmaSquared);

		
		atlas = DeterministicAtlas;
		atlasbuilder = DeterministicAtlasBuilder;
	}
    else
    {
        throw std::runtime_error("Unknown atlas type: currently supported types are: Deterministic, Bayesian or MulticlassBayesian");
    }
	
	
	atlas->SetTemplate(templateObjects);
	std::string CP_fn = paramDiffeos->GetInitialCPPosition_fn();
	atlas->SetControlPoints(CP_fn);	
	
	
	atlas->SetCPSpacing( paramDiffeos->GetInitialCPSpacing() );
	atlas->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
	atlas->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads());
	atlas->SetDiffeos(def);
	
	atlas->Update();

	// Purely informative
	MatrixType DataDomain = atlas->GetBoundingBox(); // include template domain and control point positions
	std::cout << "Working domain only Template: origin =  [" << DataDomain.get_column(0) << "] length = [" << DataDomain.get_column(1) - DataDomain.get_column(0) << "]" << std::endl;
	for (int s = 0; s < numSubjects; s++)
	{
		MatrixType BB = target[s]->GetBoundingBox();
		//std::cout << "BB " << s << " = " << BB << std::endl;
		for (int d = 0; d < Dimension; d++)
		{
			DataDomain(d,0) = ( DataDomain(d,0)<BB(d,0)?DataDomain(d,0):BB(d,0) );
			DataDomain(d,1) = ( DataDomain(d,1)>BB(d,1)?DataDomain(d,1):BB(d,1) );
		}
	}
		
	std::cout << "Working domain: origin =  [" << DataDomain.get_column(0) << "] length = [" << DataDomain.get_column(1) - DataDomain.get_column(0) << "]" << std::endl; 

	atlasbuilder->SetTargetList(target);
	atlasbuilder->SetAtlas(atlas);
	// atlasbuilder->SetDiffeos(def);
	// atlasbuilder->SetInitialCPSpacing(InitialCPSpacing);
	paramDiffeos->FreezeCP()?atlasbuilder->SetFreezeCP():atlasbuilder->UnsetFreezeCP();
	paramDiffeos->FreezeTemplate()?atlasbuilder->SetFreezeTemplate():atlasbuilder->UnsetFreezeTemplate();
	// atlasbuilder->SetInitialCPPosition(paramDiffeos->GetInitialCPPosition_fn()); // null string means no InitialCPPosition given (will be initialized with a regular lattice)
	// atlasbuilder->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
	// atlasbuilder->SetSparsityPrior(paramDiffeos->GetSparsityPrior());
	// if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientdescent") == 0)
	// {
	// 	atlasbuilder->SetGradientDescent();
	// }
	// else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "ista") == 0)
	// {
	// 	atlasbuilder->SetISTA();
	// }
	// else
	// {
	// 	if ( (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "f_ista") != 0) &&
	// 			(itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fista") != 0) )
	// 		std::cerr << "Unknown optimization method type: defaulting to FISTA" << std::endl;
	// 		atlasbuilder->SetFISTA();
	// }
//	paramDiffeos->AdaptiveDataSigma()?atlasbuilder->SetAdaptiveDataSigma():atlasbuilder->UnsetAdaptiveDataSigma();
	// atlasbuilder->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
	atlasbuilder->SetMaxIterations(paramDiffeos->GetMaxIterations());
	atlasbuilder->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
	atlasbuilder->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
	atlasbuilder->SetInitialStepMultiplier(paramDiffeos->GetInitialStepMultiplier());
	// atlasbuilder->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads());
	// atlasbuilder->UseGradientDescent();

	// create output names for saving
	// All saving name will start with AtlasName, so this is the place to add the path to a particular folder!
	///TODO : add the possibility to choose a folder

	// PIETRO
	std::string AtlasName = paramDiffeos->GetAtlasName_fn();
	if (AtlasName.empty())
		AtlasName = "Atlas";

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
	
	atlasbuilder->SetTemplateObjectsName(objectsName, objectsNameExtension);
	atlasbuilder->SetAtlasName(AtlasName);


	// create timer
	SimpleTimer timer;

	atlasbuilder->Update();

	timer.Stop();

	std::cout << "Atlas estimation took "
			<< timer.GetElapsedHours() << " hours, "
			<< timer.GetElapsedMinutes() << " minutes, "
			<< timer.GetElapsedSeconds() << " seconds"
			<< std::endl;


	std::cout << "Write output files..." << std::endl;
		
	atlasbuilder->WriteOutput();
	atlasbuilder->WriteAtlasToSubjectDeformations();
	
	std::cout << "... done" << std::endl;


	if (isOfBayesianType)
	{
		std::cout << "Write xml files that could be used later on to match the atlas to new data" << std::endl;
		VectorType DSS = atlasbuilder->GetAtlas()->GetDataSigmaSquared();
		for (int i = 0; i < numObjects; i++)
		{
			std::ostringstream oss;
			oss << AtlasName << "_" << paramObjectsList[i]->GetXMLFileName() << std::ends;
			paramObjectsList[i]->SetDataSigma( sqrt(DSS(i)) );
			bool b = writeDeformableObjectParametersXML(oss.str().c_str(), paramObjectsList[i]);
			if (b)
				std::cout << "XML file " << oss.str() << " created" << std::endl;
			else
				std::cout << "Warning: failed writing XMLfile " << oss.str() << std::endl;
				
				
		}
	}


/*
	// std::vector< MatrixType> Mom = atlasbuilder->GetMomenta();
	//
	// std::ostringstream ossMom;
	// ossMom << AtlasName << "_InitialMomentas.txt" << std::ends;
	// writeMultipleMatrixDLM<TScalar>(ossMom.str().c_str(),Mom);
	//
	// std::vector<std::string> outname(numObjects);
	// std::vector<std::string> outext(numObjects);
	//
	// for (unsigned int s = 0; s < numSubjects; s++)
	// {
	// 	for (int i = 0; i < numObjects; i++)
	// 	{
	// 		std::string sourcefn;
	// 		sourcefn.assign(templatefnList[i]);
	// 		int length = sourcefn.size();
	// 		int index = length - 1;
	// 		while( (sourcefn[index] != '.') && (index>=0) )
	// 			index--;
	//
	// 		std::ostringstream oss;
	// 		oss << sourcefn.substr(0,index) << "_to_subject_" << s;
	//
	// 		outname[i] = oss.str();
	// 		outext[i] = sourcefn.substr(index,length-index);
	//
	// 	}
	// 	atlasbuilder->GetAtlas()->WriteAtlasDeformation(Mom[s], outname, outext);
	// }

	//std::cout << "Write Flow BEGIN" << std::endl;
	// atlasbuilder->WriteFlow(outname, outext);
	//std::cout << "Write Flow END" << std::endl;

	// //std::cout << "Write Template BEGIN" << std::endl;
	// TemplateType* templ = atlasbuilder->GetTemplate();
	// templ->WriteMultiObject(outfn);
	// //std::cout << "Write template END" << std::endl;
*/
	delete atlas;
	delete atlasbuilder;

}

#endif

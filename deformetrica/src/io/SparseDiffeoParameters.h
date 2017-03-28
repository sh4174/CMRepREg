/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoParameters_h
#define _SparseDiffeoParameters_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>

#include <limits>
#include <cstddef>


class SparseDiffeoParameters : public itk::Object
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef SparseDiffeoParameters Self;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	itkNewMacro(Self);

	// Make sure all values are OK
	virtual bool CheckValues();

	void PrintSelf(std::ostream& os);

	// Compulsory parameters
	itkGetMacro(KernelWidth, double);
	itkSetMacro(KernelWidth, double);

	// Optional parameters...
	// ... for the diffeos
	itkGetMacro(KernelType, std::string);
	itkSetMacro(KernelType, std::string);

	itkGetMacro(NumberOfTimePoints, unsigned int);
	itkSetMacro(NumberOfTimePoints, unsigned int);

	itkGetMacro(T0, double);itkSetMacro(T0, double);
	itkGetMacro(TN, double);itkSetMacro(TN, double);

	itkGetMacro(InitialCPSpacing, double);
	itkSetMacro(InitialCPSpacing, double);

	itkGetMacro(P3MWorkingSpacingRatio, double);
	itkSetMacro(P3MWorkingSpacingRatio, double);

	itkGetMacro(P3MPaddingFactor, double);
	itkSetMacro(P3MPaddingFactor, double);


	// ... for the objective function
	itkGetMacro(SparsityPrior, double);
	itkSetMacro(SparsityPrior, double);


//	// ... for the gradient descent
//	inline void SetUseFISTA() { m_UseFISTA = true; }
//	inline void UnsetUseFISTA() { m_UseFISTA = false; }
//	inline bool UseFISTA() { return m_UseFISTA; }

	// ... for the optimization method
	itkGetMacro(OptimizationMethodType, std::string);
	itkSetMacro(OptimizationMethodType, std::string);

//	inline void SetAdaptiveDataSigma() { m_AdaptiveDataSigma = true; }
//	inline void UnsetAdaptiveDataSigma() { m_AdaptiveDataSigma = false; }
//	inline bool AdaptiveDataSigma() { return m_AdaptiveDataSigma; }

	//  ********************** CHANGE *******************
	inline void SetCovarianceMomenta_Normalized_Hyperparameter(double d) { m_CovarianceMomenta_Normalized_Hyperparameter = d; }
	inline double GetCovarianceMomenta_Normalized_Hyperparameter() { return m_CovarianceMomenta_Normalized_Hyperparameter; }

//	inline void SetBayesianFramework() { m_BayesianFramework = true; }
//	inline void UnsetBayesianFramework() { m_BayesianFramework = false; }
//	inline bool BayesianFramework() { return m_BayesianFramework; }
    itkGetMacro(AtlasType, std::string);
    itkSetMacro(AtlasType, std::string);

    itkGetMacro(MaximumNumberOfClasses, int);
    itkSetMacro(MaximumNumberOfClasses, int);

	itkGetMacro(CovarianceMomentaInverse_fn, std::string);
	itkSetMacro(CovarianceMomentaInverse_fn, std::string);

	itkGetMacro(CovarianceMomenta_Prior_Inverse_fn, std::string);
	itkSetMacro(CovarianceMomenta_Prior_Inverse_fn, std::string);

	itkGetMacro(AtlasName_fn, std::string);
	itkSetMacro(AtlasName_fn, std::string);

	//  ********************** CHANGE *******************

	itkGetMacro(InitialCPPosition_fn, std::string);
	itkSetMacro(InitialCPPosition_fn, std::string);

	itkGetMacro(InitialMomenta_fn, std::string);
	itkSetMacro(InitialMomenta_fn, std::string);

	inline void SetCPsAtShapePoints(bool yesNo) {m_CPsAtShapePoints = yesNo; }
	inline bool GetCpsAtShapePoints() {return m_CPsAtShapePoints; }

	inline void SetFreezeCP() { m_FreezeCP = true; }
	inline void UnsetFreezeCP() { m_FreezeCP = false; }
	inline bool FreezeCP() { return m_FreezeCP; }


	inline void SetFreezeTemplate() { m_FreezeTemplate = true; }
	inline void UnsetFreezeTemplate() { m_FreezeTemplate = false; }
	inline bool FreezeTemplate() { return m_FreezeTemplate; }

	itkGetMacro(SaveEveryNIters, unsigned int);
	void SetSaveEveryNIters(unsigned int n)
	{
		if (n<1)
		{
			std::cout<<"Invalid parameter <save-every-n-iters>, not saving intermediate results.\n";
			// A very large value means no temporary output is saved
			m_SaveEveryNIters = std::numeric_limits<unsigned int>::max();
		}
		else
		{
			m_SaveEveryNIters = n;
		}
	}

 	itkGetMacro(SmoothingKernelWidthRatio, double);
	itkSetMacro(SmoothingKernelWidthRatio, double);

	itkGetMacro(MaxIterations, unsigned int);
	itkSetMacro(MaxIterations, unsigned int);

	itkGetMacro(MaxLineSearchIterations, unsigned int);
	itkSetMacro(MaxLineSearchIterations, unsigned int);

	itkGetMacro(StepExpand, double);
	itkSetMacro(StepExpand, double);

	itkGetMacro(StepShrink, double);
	itkSetMacro(StepShrink, double);

	itkGetMacro(AdaptiveTolerance, double);
	itkSetMacro(AdaptiveTolerance, double);

	itkGetMacro(InitialStepMultiplier, double);
	itkSetMacro(InitialStepMultiplier, double);

	itkGetMacro(NumberOfThreads, unsigned int);
	itkSetMacro(NumberOfThreads, unsigned int);

	inline void SetComputeTrueInverseFlow(){ m_ComputeTrueInverseFlow = true; }
	inline void UnsetComputeTrueInverseFlow(){ m_ComputeTrueInverseFlow = false; }
	inline bool ComputeTrueInverseFlow() { return m_ComputeTrueInverseFlow; }


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SparseDiffeoParameters();

	~SparseDiffeoParameters();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	double m_KernelWidth;

	std::string m_KernelType;
	double m_T0;
	double m_TN;
	unsigned int m_NumberOfTimePoints;
	double m_InitialCPSpacing;
	double m_P3MWorkingSpacingRatio;
	double m_P3MPaddingFactor;

	bool m_ComputeTrueInverseFlow;

	double m_SparsityPrior;

	std::string m_OptimizationMethodType;
//	bool m_AdaptiveDataSigma;

	//  ********************** CHANGE *******************
	double m_CovarianceMomenta_Normalized_Hyperparameter;
//	bool m_BayesianFramework;
    std::string m_AtlasType;
    int m_MaximumNumberOfClasses;
	std::string m_CovarianceMomenta_Prior_Inverse_fn; // Inverse of the Prior of the Cov Momenta in a Bayesian Framework
	std::string m_CovarianceMomentaInverse_fn; // Matrix to use in place of the kernel in a Deterministic Framework
	std::string m_AtlasName_fn;
	//  ********************** CHANGE *******************

	bool m_CPsAtShapePoints;					// Is it used in Deformetrica ?
	bool m_FreezeCP;
	bool m_FreezeTemplate;
	std::string m_InitialCPPosition_fn;
	std::string m_InitialMomenta_fn;

	unsigned int m_SaveEveryNIters;

	double m_SmoothingKernelWidthRatio;

	unsigned int m_MaxIterations;
	unsigned int m_MaxLineSearchIterations;

	double m_StepExpand;
	double m_StepShrink;
	double m_AdaptiveTolerance;
	double m_InitialStepMultiplier;

	unsigned int m_NumberOfThreads;


}; /* class SparseDiffeoParameters */

#endif /* _SparseDiffeoParameters_h */

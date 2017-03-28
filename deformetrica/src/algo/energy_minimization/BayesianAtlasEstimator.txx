/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BayesianAtlasEstimator_txx
#define _BayesianAtlasEstimator_txx

#include "AbstractAtlasEstimator.h"

#include "KernelFactory.h"

#include "Diffeos.h"

#include "writeMatrixDLM.txx"

#include <iostream>
#include <iomanip>
#include <sstream>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
BayesianAtlasEstimator<TScalar, Dimension>
::BayesianAtlasEstimator() : Superclass()
 {
	 // // convert m_Atlas of type Atlas.h into BayesianAtlas.h
	 // m_Atlas = dynamic_cast<AtlasType*>(Superclass::m_Atlas);

 }



template <class TScalar, unsigned int Dimension>
BayesianAtlasEstimator<TScalar, Dimension>
::~BayesianAtlasEstimator()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
BayesianAtlasEstimator<TScalar, Dimension>
::Update()
 {
	
 	if ( (m_InitialMomentas.size() == Superclass::m_NumberOfSubjects) && (m_InitialMomentas[0].rows() == m_NumberOfCPs) && (m_InitialMomentas[0].columns() == Dimension) )
	{
		std::cout << "Using predefined set of momenta" << std::endl;
	}
	else
	{
		if (m_InitialMomentas.size() > 0)
			std::cout << "Warning: initial momenta file has incompatible number of subjects and/or vectors. Initial momenta reset to zero" << std::endl;

		m_InitialMomentas.resize(Superclass::m_NumberOfSubjects);
		for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
		{
			m_InitialMomentas[s].set_size(m_NumberOfCPs, Dimension);
			m_InitialMomentas[s].fill(0.0);
		}
	}
	
	m_ValuesHistory.resize(Superclass::m_MaxIterations + 1);
	
	MatrixType ControlPoints = m_Atlas->GetControlPoints();
	MatrixList TemplateData = m_Atlas->GetTemplateData();
	VectorType DataSigmaSquared = m_Atlas->GetDataSigmaSquared();
	MatrixType CovarianceMomenta_Inverse = m_Atlas->GetCovarianceMomentaInverse();
		
	MatrixType& X = ControlPoints;
	MatrixList& A = m_InitialMomentas;
	MatrixList& T = TemplateData;
	VectorType& DSS = DataSigmaSquared;
	MatrixType& CovMomInv = CovarianceMomenta_Inverse;
	
	this->GradientDescent(X, A, T, DSS, CovMomInv);
	
	m_Atlas->SetControlPoints(ControlPoints);
	m_Atlas->SetTemplateData(TemplateData);
	m_Atlas->SetDataSigmaSquared(DataSigmaSquared);
	m_Atlas->SetCovarianceMomentaInverse(CovarianceMomenta_Inverse);
	
		
 }
 


template <class TScalar, unsigned int Dimension>
void
BayesianAtlasEstimator<TScalar, Dimension>
::WriteOutput(std::string AtlasName)
 {
	 
	Superclass::WriteOutput(AtlasName);
			
	std::ostringstream oss;

	MatrixList Mom = this->GetMomenta();
	oss << AtlasName << "_InitialMomentas.txt" << std::ends;
	writeMultipleMatrixDLM<TScalar>(oss.str().c_str(), Mom);

 }



template <class TScalar, unsigned int Dimension>
void
BayesianAtlasEstimator<TScalar, Dimension>
::WriteAtlasToSubjectDeformations()
{
	
	MatrixList Mom = this->GetMomenta();
	
	for (unsigned int s = 0; s < Mom.size(); s++)
	{
		std::vector< std::string > tempfn(Superclass::m_NumberOfObjects);
		for (unsigned int i = 0; i < Superclass::m_NumberOfObjects; i++)
		{
			std::ostringstream oss;
			oss << Superclass::m_TemplateObjectsName[i] << "_to_subject_" << s;
			tempfn[i] = oss.str();
		}
		
		m_Atlas->WriteAtlasDeformation(Mom[s], Superclass::m_AtlasName, tempfn, Superclass::m_TemplateObjectsNameExtension);
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
template <class TScalar, unsigned int Dimension>
void
BayesianAtlasEstimator<TScalar, Dimension>
::GradientDescent(MatrixType& X, MatrixList& A, MatrixList& T, VectorType& DSS, MatrixType& CovMomInv)
 {
	 
	MatrixType Xnew;
	MatrixList Anew(Superclass::m_NumberOfSubjects);
	MatrixList Tnew(Superclass::m_NumberOfObjects);

	// m_GradPos.set_size(m_NumberOfCPs, Dimension);
	// m_GradMom.resize(m_NumberOfSubjects);
	// m_GradTemplate_L2.resize(m_NumberOfObjects);
	// for (int i = 0; i < m_NumberOfObjects; i++)
	// 	m_GradTemplate_L2[i].set_size(T[i].rows(), T[i].columns());

	TScalar stepXA, stepT;
	
	m_Atlas->SetControlPoints(X);
	m_Atlas->SetTemplateData(T);
	m_Atlas->SetDataSigmaSquared(DSS);
	m_Atlas->SetCovarianceMomentaInverse(CovMomInv);
	
	// PIETRO, it does not work since it is protected
	//m_Atlas->UpdateDeformationAndKernelDataDomain(Superclass::m_TargetList);

	// Get the value of the likelihood for the initial momentas and set DSS and CovMomInv, which minimizes the likelihood for these values of the initial momenta
	m_ValuesHistory[0] = m_Atlas->ComputeLikelihoodWithOptimalParameters(A, Superclass::m_TargetList);
	DSS = m_ValuesHistory[0]->GetDataSigmaSquared();
	CovMomInv = m_ValuesHistory[0]->GetCovMomInverse();
	m_ValuesHistory[0]->PrintIter(0);
	
	TScalar lsqRef = m_ValuesHistory[0]->GetLikelihood();
	unsigned int iterRef = 0;
	

	for (unsigned int iter = 0; iter < Superclass::m_MaxIterations; iter++)
	{

		m_Atlas->SetControlPoints(X);
		m_Atlas->SetTemplateData(T);
		m_Atlas->SetDataSigmaSquared(DSS);
		m_Atlas->SetCovarianceMomentaInverse(CovMomInv);


		if (!(iter % 1))
		{
			std::ostringstream oss;
			oss << "Iter" << std::setfill('0') << std::setw(4) << iter + 1 << "_" << Superclass::m_AtlasName;
			this->WriteOutput(oss.str());
		}

	
		m_Atlas->ComputeLikelihoodGradient(A, Superclass::m_TargetList, m_GradMom, m_GradPos, m_GradTempL);
		// DSS = m_ValuesHistory[0]->GetOptimalDataSigmaSquared();
		// CovMomInv = m_ValuesHistory[0]->GetOptimalCovMomInverse();
		// m_ValuesHistory[0]->PrintIter(0);
				
		// this->ThreadedComputeGradient(X, A, T);

		// DEBUG
		//std::cout << "GradPos = " << std::endl << m_GradPos << std::endl;
		//for (int s = 0; s < m_NumberOfSubjects; s++)
		//	std::cout << "GradMom[" <<  s << "] = " << std::endl << m_GradMom[s] << std::endl;
		//std::cout << std::endl;


		if (iter == 0)
		{
			TScalar maxGrad = 1e-20;
			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			{
				TScalar g = 0;
				g = m_GradPos.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
				for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
				{
					g = m_GradMom[s].get_row(i).magnitude();
					if (g > maxGrad)
						maxGrad = g;
				}
			}

			//std::cout << "maxGrad = " << maxGrad << std::endl;

			TScalar initStepXA = Superclass::m_InitialStepMultiplier * fabs(m_ValuesHistory[0]->GetLikelihood()) / maxGrad / maxGrad;
			TScalar initStepT = initStepXA;

			stepXA = initStepXA;
			stepT = initStepT;
		}

		// std::vector<int> activeCPMask(m_NumberOfCPs, 1);
		bool foundMin = false;

		for (unsigned int li = 0; li < Superclass::m_MaxLineSearchIterations; li++)
		{
			std::cout << "stepsizeXA = " << stepXA << "\tstepsizeT = " << stepT << std::endl;

			this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);
			//std::cout << "Xnew =\n " << Xnew << std::endl;
			//std::cout << "X =\n " << X << std::endl;
			//std::cout << "Anew 0=\n " << Anew[0] << std::endl;
			//std::cout << "A 0=\n " << A[0] << std::endl;
			//std::cout << "Tnew 0=\n " << Tnew[0] << std::endl;
			//std::cout << "T 0=\n " << T[0] << std::endl;
			m_Atlas->SetControlPoints(Xnew);
			m_Atlas->SetTemplateData(Tnew);
			m_ValuesHistory[iter+1] = m_Atlas->ComputeLikelihoodWithOptimalParameters(Anew, Superclass::m_TargetList);


			TScalar Q = lsqRef - m_ValuesHistory[iter+1]->GetLikelihood();
			//std::cout << "Q = " << Q << std::endl;

			if (Q >= 0) //(J1 < Jcurr)
			{
				foundMin = true;
				break;
			}
			else
			{
				// case test 1:
				stepXA *= Superclass::m_AdaptiveShrink;
				
				this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);
				m_Atlas->SetControlPoints(Xnew);
				m_Atlas->SetTemplateData(Tnew);
				FunctionalValuesType* valuesTest1 = m_Atlas->ComputeLikelihoodWithOptimalParameters(Anew, Superclass::m_TargetList);

				TScalar Q1 = lsqRef - valuesTest1->GetLikelihood();

				//std::cout << "test case 1 = " << Q1 << std::endl;

				// case test 2:
				stepXA /= Superclass::m_AdaptiveShrink;
				stepT *= Superclass::m_AdaptiveShrink;
				
				this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);
				m_Atlas->SetControlPoints(Xnew);
				m_Atlas->SetTemplateData(Tnew);
				FunctionalValuesType* valuesTest2 = m_Atlas->ComputeLikelihoodWithOptimalParameters(Anew, Superclass::m_TargetList);

				TScalar Q2 = lsqRef - valuesTest2->GetLikelihood();
				//std::cout << "test case 2 = " << Q2 << std::endl;

				if ( (Q1 >= 0) || (Q2 >= 0) ) //( (Jtest1 < Jcurr) || (Jtest2 < Jcurr) )
				{
					if ( Q1 >= Q2)
					{
						stepXA *= Superclass::m_AdaptiveShrink;
						stepT /= Superclass::m_AdaptiveShrink;
						this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);
						m_ValuesHistory[iter+1] = valuesTest1->Clone();
					}
					else
					{
						m_ValuesHistory[iter+1] = valuesTest2->Clone();
					}

					foundMin = true;
					break;
				}
				else
				{
					stepXA *= Superclass::m_AdaptiveShrink;
				}

				delete valuesTest1;
				delete valuesTest2;

			}
		} // for li


		if (foundMin)
		{
			X = Xnew;
			A = Anew;
			T = Tnew;
			DSS = m_ValuesHistory[iter+1]->GetDataSigmaSquared();
			CovMomInv = m_ValuesHistory[iter+1]->GetCovMomInverse();

			stepXA *= Superclass::m_AdaptiveExpand;
			stepT *= Superclass::m_AdaptiveExpand;
			lsqRef = m_ValuesHistory[iter+1]->GetLikelihood();
		}

		if (!foundMin)
		{
			// Loop terminated without finding smaller functional
			std::cout << " number of loops exceeded " << std::endl;
			break;
		}

		if ( (m_ValuesHistory[iter]->GetLikelihood() - m_ValuesHistory[iter+1]->GetLikelihood())
				< Superclass::m_AdaptiveTolerance * (m_ValuesHistory[iterRef]->GetLikelihood() - m_ValuesHistory[iter+1]->GetLikelihood()) )
			break;

		// Jcurr = J1;

		m_ValuesHistory[iter+1]->PrintIter(iter+1);

	} // for iter
 }



template <class TScalar, unsigned int Dimension>
void
BayesianAtlasEstimator<TScalar, Dimension>
::GradientDescentStep(
		MatrixType& Xtest, const MatrixType& X,
		MatrixList& Atest, const MatrixList& A,
		MatrixList& Ttest, const MatrixList& T,
		TScalar stepXA, TScalar stepT)
{

	if (Superclass::m_UpdateCP)
		Xtest = X - m_GradPos*stepXA;
	else
		Xtest = X;

	for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
		Atest[s] = A[s] - m_GradMom[s] * stepXA;

	for (int i = 0; i < Superclass::m_NumberOfObjects; i++)
		Ttest[i] = T[i] - m_GradTempL[i] * stepT;
		
	return;

}



#endif /* _BayesianAtlasEstimator_txx */

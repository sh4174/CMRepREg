/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeterministicAtlasEstimator_txx
#define _DeterministicAtlasEstimator_txx

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
DeterministicAtlasEstimator<TScalar, Dimension>
::DeterministicAtlasEstimator() : Superclass()
 {
	 m_OptimizationMethod = null;
 }



template <class TScalar, unsigned int Dimension>
DeterministicAtlasEstimator<TScalar, Dimension>
::~DeterministicAtlasEstimator()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeterministicAtlasEstimator<TScalar, Dimension>
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
		
	MatrixType& X = ControlPoints;
	MatrixList& A = m_InitialMomentas;
	MatrixList& T = TemplateData;
	
	
	// run optimization method
	if ( (m_OptimizationMethod == GradientDescent) && (m_Atlas->GetSparsityPrior() > 1e-10) )
	{
		std::cout << "Warning: gradient descent does not work with L^1 penalty term. Optimization method switched to ISTA" << std::endl;
		m_OptimizationMethod = ISTA;
	}
		
	if (m_OptimizationMethod == F_ISTA)
		this->FISTA(X, A, T);
	else
		this->GradientDescentAndISTA(X, A, T);
	
	if (Superclass::m_UpdateCP)
		m_Atlas->SetControlPoints(ControlPoints);
	if (Superclass::m_UpdateTemplate)
		m_Atlas->SetTemplateData(TemplateData);
 }
 


template <class TScalar, unsigned int Dimension>
void
DeterministicAtlasEstimator<TScalar, Dimension>
::WriteOutput(std::string AtlasName)
{
	Superclass::WriteOutput(AtlasName);
			
	MatrixList Mom = this->GetMomenta();

	std::ostringstream ossMom;
	ossMom << AtlasName << "_InitialMomentas.txt" << std::ends;

	writeMultipleMatrixDLM<TScalar>(ossMom.str().c_str(),Mom);
}



template <class TScalar, unsigned int Dimension>
void
DeterministicAtlasEstimator<TScalar, Dimension>
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
DeterministicAtlasEstimator<TScalar, Dimension>
::GradientDescentAndISTA(MatrixType& X, MatrixList& A, MatrixList& T)
 {
	 
	MatrixType Xnew;
	MatrixList Anew(Superclass::m_NumberOfSubjects);
	MatrixList Tnew(Superclass::m_NumberOfObjects);
	
	TScalar stepXA, stepT;
	
	if (Superclass::m_UpdateCP)
		m_Atlas->SetControlPoints(X);
	if (Superclass::m_UpdateTemplate)
		m_Atlas->SetTemplateData(T);
	
	// Get the value of the functional for the initial momentas
	m_ValuesHistory[0] = m_Atlas->ComputeFunctional(A, Superclass::m_TargetList);
	m_ValuesHistory[0]->PrintIter(0);
	
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();
	unsigned int iterRef = 0;
	

	for (unsigned int iter = 0; iter < Superclass::m_MaxIterations; iter++)
	{

		if (!(iter % 3))
		{
			std::ostringstream oss;
			oss << "Iter" << std::setfill('0') << std::setw(4) << iter + 1 << "_" << Superclass::m_AtlasName;
			this->WriteOutput(oss.str());
		}

		if (Superclass::m_UpdateCP)
			m_Atlas->SetControlPoints(X);
		if (Superclass::m_UpdateTemplate)
			m_Atlas->SetTemplateData(T);
	
		if (Superclass::m_UpdateTemplate)
			m_Atlas->ComputeFunctionalGradient(A, Superclass::m_TargetList, m_GradMom, m_GradPos, m_GradTempL_L2, m_GradTempL_Sob);
		else
			m_Atlas->ComputeFunctionalGradient(A, Superclass::m_TargetList, m_GradMom, m_GradPos, m_GradTempL_L2, false);
				
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

			TScalar initStepXA = Superclass::m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			TScalar initStepT = initStepXA;

			stepXA = initStepXA;
			stepT = initStepT;
		}

		std::vector<int> activeCPMask(m_NumberOfCPs, 1);
		bool foundMin = false;

		for (unsigned int li = 0; li < Superclass::m_MaxLineSearchIterations; li++)
		{
		 	std::cout << "stepsizeXA = " << stepXA;
			if (Superclass::m_UpdateTemplate)
				std::cout << "\tstepsizeT = " << stepT;
			std::cout << std::endl;

			this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
			
			if (Superclass::m_UpdateCP)
				m_Atlas->SetControlPoints(Xnew);
			if (Superclass::m_UpdateTemplate)
				m_Atlas->SetTemplateData(Tnew);
			
			
			m_ValuesHistory[iter+1] = m_Atlas->ComputeFunctional(Anew, Superclass::m_TargetList);


			TScalar Q = lsqRef - m_ValuesHistory[iter+1]->GetTotalL2Value();
			if (m_OptimizationMethod == ISTA)
					Q += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);


			if (Q >= 0) //(J1 < Jcurr)
			{
				foundMin = true;
				break;
			}
			else
			{
				// case test 1:
				stepXA *= Superclass::m_AdaptiveShrink;
				
				this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
				
				if (Superclass::m_UpdateCP)
					m_Atlas->SetControlPoints(Xnew);
				if (Superclass::m_UpdateTemplate)
					m_Atlas->SetTemplateData(Tnew);

				FunctionalValuesType* valuesTest1 = m_Atlas->ComputeFunctional(Anew, Superclass::m_TargetList);

				TScalar Q1 = lsqRef - valuesTest1->GetTotalL2Value();
				if (m_OptimizationMethod == ISTA)
						Q1 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

				// case test 2:
				TScalar Q2 = -1.0;
				stepXA /= Superclass::m_AdaptiveShrink;
				FunctionalValuesType* valuesTest2;
				if (Superclass::m_UpdateTemplate)
				{
					stepT *= Superclass::m_AdaptiveShrink;
				
					this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
					m_Atlas->SetControlPoints(Xnew);
					m_Atlas->SetTemplateData(Tnew);
					valuesTest2 = m_Atlas->ComputeFunctional(Anew, Superclass::m_TargetList);

					Q2 = lsqRef - valuesTest2->GetTotalL2Value();
					if (m_OptimizationMethod == ISTA)
							Q2 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

					//std::cout << "test case 2 = " << Q2 << std::endl;
				}

				
				if ( (Q1 >= 0) || (Q2 >= 0) ) //( (Jtest1 < Jcurr) || (Jtest2 < Jcurr) )
				{
					if ( Q1 >= Q2 )
					{
						stepXA *= Superclass::m_AdaptiveShrink;
						stepT /= Superclass::m_AdaptiveShrink;
						this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
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
				if (Superclass::m_UpdateTemplate)
					 delete valuesTest2;

			}
		} // for li


		if (foundMin)
		{
			X = Xnew;
			A = Anew;
			T = Tnew;

			stepXA *= Superclass::m_AdaptiveExpand;
			stepT *= Superclass::m_AdaptiveExpand;
			lsqRef = m_ValuesHistory[iter+1]->GetTotalL2Value();
		}

		if (!foundMin)
		{
			// Loop terminated without finding smaller functional
			std::cout << " number of loops exceeded " << std::endl;
			break;
		}

		if ( (m_ValuesHistory[iter]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue())
				< Superclass::m_AdaptiveTolerance * (m_ValuesHistory[iterRef]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue()) )
			break;

		// Jcurr = J1;

		m_ValuesHistory[iter+1]->PrintIter(iter+1);

	} // for iter
 }



template <class TScalar, unsigned int Dimension>
void
DeterministicAtlasEstimator<TScalar, Dimension>
::FISTA(MatrixType& X, MatrixList& A, MatrixList& T)
 {
	 
	MatrixType Xtest;
 	MatrixList Atest(Superclass::m_NumberOfSubjects);
 	MatrixList Ttest(Superclass::m_NumberOfObjects);
	
 	// m_GradPos.set_size(m_NumberOfCPs, Dimension);
 	// m_GradMom.resize(m_NumberOfSubjects);
 	// m_GradTemplate_L2.resize(m_NumberOfObjects);
 	// for (int i = 0; i < m_NumberOfObjects; i++)
 	// 	m_GradTemplate_L2[i].set_size(T[i].rows(), T[i].columns());

 	TScalar stepXA, stepT;
	
	if (Superclass::m_UpdateCP)
		m_Atlas->SetControlPoints(X);
	if (Superclass::m_UpdateTemplate)
		m_Atlas->SetTemplateData(T);
	
 	// Get the value of the functional for the initial momentas	
 	m_ValuesHistory[0] = m_Atlas->ComputeFunctional(A, Superclass::m_TargetList);	
 	m_ValuesHistory[0]->PrintIter(0);
	
 	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();
	
	MatrixType Xprev = X;
	MatrixList Aprev = A;
	MatrixList Tprev = T;

	std::vector<int> activeCPMask(m_NumberOfCPs, 1);
	std::vector<int> activeCPMask_test(m_NumberOfCPs, 1);

	TScalar Qarray[Superclass::m_UpdateTemplate?4:2];
	unsigned int freezeDirectionCounter = 0;
	unsigned int iterRef = 0;

	TScalar tau = 1.0;
 	for (unsigned int iter = 0; iter < Superclass::m_MaxIterations; iter++)
 	{

 		if (!(iter % 3))
 		{
			std::ostringstream oss;
			oss << "Iter" << std::setfill('0') << std::setw(4) << iter + 1 << "_" << Superclass::m_AtlasName;
			this->WriteOutput(oss.str());
 		}

		if (Superclass::m_UpdateCP)
			m_Atlas->SetControlPoints(X);
		if (Superclass::m_UpdateTemplate)
 			m_Atlas->SetTemplateData(T);
	
		if (Superclass::m_UpdateTemplate)
 			m_Atlas->ComputeFunctionalGradient(A, Superclass::m_TargetList, m_GradMom, m_GradPos, m_GradTempL_L2, m_GradTempL_Sob);
		else
			m_Atlas->ComputeFunctionalGradient(A, Superclass::m_TargetList, m_GradMom, m_GradPos, m_GradTempL_L2, 0);
		
				
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

 			TScalar initStepXA = Superclass::m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
 			TScalar initStepT = initStepXA;

 			stepXA = initStepXA;
 			stepT = initStepT;
 		}

		bool minimtest = false;
		unsigned int lineIter = 0;
		for (; lineIter < Superclass::m_MaxLineSearchIterations; lineIter++)
		{
		 	std::cout << "stepsizeXA = " << stepXA;
			if (Superclass::m_UpdateTemplate)
				std::cout << "\tstepsizeT = " << stepT;
			std::cout << std::endl;

		 	this->GradientDescentStep(Xtest, X, Atest, A, Ttest, T, activeCPMask_test, stepXA, stepT);

			if (Superclass::m_UpdateCP)
				m_Atlas->SetControlPoints(Xtest);
			if (Superclass::m_UpdateTemplate)
				m_Atlas->SetTemplateData(Ttest);

			m_ValuesHistory[iter+1] = m_Atlas->ComputeFunctional(Atest, Superclass::m_TargetList);

			TScalar Q_XAT = lsqRef -  m_ValuesHistory[iter+1]->GetTotalL2Value() + this->_QdiffTerm(
					Xtest, X, Atest, A, Ttest, T, stepXA, stepT);

		 	if (Q_XAT >= 0) //(J1 < Jcurr)
		 	{
		 		minimtest = true;
		 		break;
		 	}
		 	else
		 	{
				TScalar Q_X1, Q_X2, Q_T1, Q_T2;
				FunctionalValuesType *values_XAprime1, *values_XAprime2, *values_Tprime1, *values_Tprime2;

		 		// case test 0:
		 		stepXA *= Superclass::m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime1(m_NumberOfCPs, 1);

				MatrixType X_prime1;
				MatrixList A_prime1(Superclass::m_NumberOfSubjects);

		 		this->GradientDescentStep(X_prime1, X, A_prime1, A, Ttest, T, activeCPMask_prime1, stepXA, stepT);
				
				if (Superclass::m_UpdateCP)
		 			m_Atlas->SetControlPoints(X_prime1);
				if (Superclass::m_UpdateTemplate)
		 			m_Atlas->SetTemplateData(Ttest);
		 		
				values_XAprime1 = m_Atlas->ComputeFunctional(A_prime1, Superclass::m_TargetList);

				Q_X1 = lsqRef - values_XAprime1->GetTotalL2Value() + this->_QdiffTerm(
						X_prime1, X, A_prime1, A, Ttest, T, stepXA, stepT);

		 				//std::cout << "test case 0 = " << Q_X1 << std::endl;

		 				// case test 1:
				stepXA *= Superclass::m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime2(m_NumberOfCPs, 1);

				MatrixType X_prime2;
				MatrixList A_prime2(Superclass::m_NumberOfSubjects);

				this->GradientDescentStep(X_prime2, X, A_prime2, A, Ttest, T, activeCPMask_prime2, stepXA, stepT);
				
				if (Superclass::m_UpdateCP)
 					m_Atlas->SetControlPoints(X_prime2);
				if (Superclass::m_UpdateTemplate)
					m_Atlas->SetTemplateData(Ttest);
 				
				values_XAprime2 = m_Atlas->ComputeFunctional(A_prime2, Superclass::m_TargetList);

				Q_X2 = lsqRef - values_XAprime2->GetTotalL2Value() + this->_QdiffTerm(
						X_prime2, X, A_prime2, A, Ttest, T, stepXA, stepT);

		 				//std::cout << "test case 1 = " << Q2 << std::endl;

				// test case 2
				stepXA /= (Superclass::m_AdaptiveShrink * Superclass::m_AdaptiveShrink);
				MatrixList T_prime1(Superclass::m_NumberOfObjects);
				MatrixList T_prime2(Superclass::m_NumberOfObjects);
				if (Superclass::m_UpdateTemplate)
				{
					stepT *= Superclass::m_AdaptiveShrink;

					this->GradientDescentStep(Xtest, X, Atest, A, T_prime1, T, activeCPMask_test, stepXA, stepT);
					
					m_Atlas->SetControlPoints(Xtest);
					m_Atlas->SetTemplateData(T_prime1);
			 			
					values_Tprime1 = m_Atlas->ComputeFunctional(Atest, Superclass::m_TargetList);

					Q_T1 = lsqRef - values_Tprime1->GetTotalL2Value() + this->_QdiffTerm(
							Xtest, X, Atest, A, T_prime1, T, stepXA, stepT);


					// test case 3
					stepT *= Superclass::m_AdaptiveShrink;
					this->GradientDescentStep(Xtest, X, Atest, A, T_prime2, T, activeCPMask_test, stepXA, stepT);
					
	 				m_Atlas->SetControlPoints(Xtest);
	 				m_Atlas->SetTemplateData(T_prime2);
	 				values_Tprime2 = m_Atlas->ComputeFunctional(Atest, Superclass::m_TargetList);

					Q_T2 = lsqRef - values_Tprime2->GetTotalL2Value() + this->_QdiffTerm(
							Xtest, X, Atest, A, T_prime2, T, stepXA, stepT);
				}


				// Select best case
				Qarray[0] = Q_X1;
				Qarray[1] = Q_X2;
				if (Superclass::m_UpdateTemplate)
				{
					Qarray[2] = Q_T1;
					Qarray[3] = Q_T2;
				}

				int ind_max = 0;
				TScalar Q_max = Qarray[0];
				for (int j = 1; j < ( (Superclass::m_UpdateTemplate)?4:2 ); j++)
				{
					if (Qarray[j] > Q_max)
					{
						ind_max = j;
						Q_max = Qarray[j];
					}
				}


		 		if ( Q_max >= 0 )
				{
					minimtest = true;
					if (ind_max == 0)
					{
						// std::cout << " Q_X1" << std::endl;
						Xtest = X_prime1;
						Atest = A_prime1;
						activeCPMask_test = activeCPMask_prime1;
						m_ValuesHistory[iter+1] = values_XAprime1->Clone();

						stepXA *= Superclass::m_AdaptiveShrink;
						stepT /= (Superclass::m_AdaptiveShrink * Superclass::m_AdaptiveShrink);
					}
					else if (ind_max == 1)
					{
						// std::cout << " Q_X2" << std::endl;
						Xtest = X_prime2;
						Atest = A_prime2;
						activeCPMask_test = activeCPMask_prime2;
						m_ValuesHistory[iter+1] = values_XAprime2->Clone();

						stepXA *= (Superclass::m_AdaptiveShrink * Superclass::m_AdaptiveShrink);
						stepT /= (Superclass::m_AdaptiveShrink * Superclass::m_AdaptiveShrink);
					}
					else if (ind_max == 2)
					{
						// std::cout << " Q_T1" << std::endl;
						Ttest = T_prime1;
						m_ValuesHistory[iter+1] = values_Tprime1->Clone();

						stepT /= Superclass::m_AdaptiveShrink;
					}
					else
					{
						// std::cout << " Q_T2" << std::endl;
						Ttest = T_prime2;
						m_ValuesHistory[iter+1] = values_Tprime2->Clone();
					}
				} // if (Q_max > 0)
				else // update failed, continue line search from X - stepX*gradX, A - stepA*gradA, T - stepT*gradT
				{
					stepXA *= Superclass::m_AdaptiveShrink;
					stepT /= Superclass::m_AdaptiveShrink;
				}

				delete values_XAprime1;
				delete values_XAprime2;
				if (Superclass::m_UpdateTemplate)
				{
					delete values_Tprime1;
					delete values_Tprime2;					
				}

			 }  // if (Q_XAT >= 0)

			if (minimtest)
				break;

		} // end line search


		if (minimtest)
		{
			TScalar tau_next = 1.0;

			if ( (lineIter == 0)
					&& ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze)) )
			{
				// stepsize increased
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau / Superclass::m_AdaptiveExpand)) / 2.0;
				if (freezeDirectionCounter > NbIterFreeze)
					freezeDirectionCounter = 0;
			}
			else // stepsize decreased
			{
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau)) / 2.0;
				freezeDirectionCounter++;
			}

			TScalar tau_scale = (tau-1.0) / tau_next;
			if (Superclass::m_UpdateCP)
				X = Xtest + (Xtest - Xprev) * tau_scale;
			for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
				 A[s] = Atest[s] + (Atest[s] - Aprev[s]) * tau_scale;
			if (Superclass::m_UpdateTemplate)
			{
				for (unsigned int i = 0; i < Superclass::m_NumberOfObjects; i++)
					T[i] = Ttest[i] + (Ttest[i] - Tprev[i]) * tau_scale;				
			}
			activeCPMask = activeCPMask_test;

			Xprev = Xtest;
			Aprev = Atest;
			Tprev = Ttest;

			tau = tau_next;

			// Check if we're in increasing trend
			if ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze))
			{
				stepXA *= Superclass::m_AdaptiveExpand;
				stepT *= Superclass::m_AdaptiveExpand;
			}

			// Update the reference L2 part of the functional
			if (Superclass::m_UpdateCP)
				m_Atlas->SetControlPoints(X);
			if (Superclass::m_UpdateTemplate)
				m_Atlas->SetTemplateData(T);
			
			FunctionalValuesType* values = m_Atlas->ComputeFunctional(A, Superclass::m_TargetList);
			if (values->IsOutOfBox()) //m_Def->OutOfBox() || m_Template->OutOfBox())
			{
				std::cerr << "Out of box: needs to restart FISTA from this point." << std::endl;
				X = Xtest; A = Atest; 
				if (Superclass::m_UpdateTemplate)
					T = Ttest;
				
				lsqRef = m_ValuesHistory[iter+1]->GetTotalL2Value();
			}
			else
				lsqRef = values->GetTotalL2Value();

			// these values are not necessarily decreasing, whereas Functional(Xprev) is.
			m_ValuesHistory[iter+1]->PrintIter(iter+1);

			int NbActiveCP = 0;
			for (int i = 0; i < m_NumberOfCPs; i++)
				NbActiveCP += (activeCPMask[i] != 0);
			std::cout << NbActiveCP << " active control points" << std::endl;
		}

		if (!minimtest) // Inner loop terminated without finding smaller functional
			break;

		TScalar deltaF_curr = m_ValuesHistory[iter]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue();
		TScalar deltaF_ref = m_ValuesHistory[iterRef]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue();
		// if ( (deltaF_curr > 0) && (deltaF_curr < m_AdaptiveTolerance*deltaF_ref) )
		if (fabs(deltaF_curr) < Superclass::m_AdaptiveTolerance*fabs(deltaF_ref))
		{
			std::cout << "Tolerance BREAK" << std::endl;
			std::cout << "FINAL VALUES: ";
			m_ValuesHistory[iter+1]->PrintIter(iter+1);
			break;
		}

	} // end iter

	X = this->_maskMatrix(X, activeCPMask);
	for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
		A[s] = this->_maskMatrix(A[s], activeCPMask);

 }




template <class TScalar, unsigned int Dimension>
typename DeterministicAtlasEstimator<TScalar, Dimension>::MatrixType
DeterministicAtlasEstimator<TScalar, Dimension>
::_maskMatrix(const MatrixType& M, const std::vector<int>& mask)
 {
	unsigned int numElements = 0;
	for (unsigned int i = 0; i < mask.size(); i++)
		if (mask[i] != 0)
			numElements++;

	if (numElements == M.rows())
		return M;

	MatrixType Z(numElements, M.columns());

	unsigned int r = 0;
	for (unsigned int i = 0; i < M.rows(); i++)
	{
		if (mask[i] != 0)
			Z.set_row(r++, M.get_row(i));
	}

	return Z;
 }



template <class TScalar, unsigned int Dimension>
void
DeterministicAtlasEstimator<TScalar, Dimension>
::GradientDescentStep(
		MatrixType& Xtest, const MatrixType& X,
		MatrixList& Atest, const MatrixList& A,
		MatrixList& Ttest, const MatrixList& T,
		std::vector< int >& activeCPMask,
		TScalar stepXA, TScalar stepT)
{

	if (Superclass::m_UpdateCP)
		Xtest = X - m_GradPos*stepXA;
	else
		Xtest = X;

	if (m_Atlas->GetSparsityPrior() < 1e-10)
	{
		activeCPMask.assign(m_NumberOfCPs, 1);
		for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
			Atest[s] = A[s] - m_GradMom[s] * stepXA;
	}
	else
	{
		activeCPMask.assign(m_NumberOfCPs, 0); // gives for each CP the number of non-zero momenta that it carries
		for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
			Atest[s] = this->SoftThresholdUpdate(A[s], m_GradMom[s], stepXA, activeCPMask);
	}

	if (Superclass::m_UpdateTemplate)
	{
		for (int i = 0; i < Superclass::m_NumberOfObjects; i++)
			Ttest[i] = T[i] - m_GradTempL_Sob[i]*stepT;		
	}
	
	return;

}



template <class TScalar, unsigned int Dimension>
typename DeterministicAtlasEstimator<TScalar, Dimension>::MatrixType
DeterministicAtlasEstimator<TScalar, Dimension>
::SoftThresholdUpdate(const MatrixType& X, MatrixType& gradX, TScalar step, std::vector<int>& mask)
 {
	MatrixType Xnew = X - step*gradX;
	TScalar thres = step * m_Atlas->GetSparsityPrior();

	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
	{
		TScalar ai = Xnew.get_row(i).magnitude();
		TScalar ai_thres = ai;
		// soft-thresholding
		if (ai > thres)
			ai_thres -= thres;
		else if (ai < -thres)
			ai_thres += thres;
		else
			ai_thres = 0;

		// Hard-thresholding
		// if (fabs(ai) < thres)
		// 	ai_thres = 0;

		mask[i] += (ai_thres != 0);

		if (ai != 0.0)
			Xnew.set_row(i, Xnew.get_row(i) * (ai_thres / ai));
	}

	return Xnew;
 }



template <class TScalar, unsigned int Dimension>
TScalar
DeterministicAtlasEstimator<TScalar, Dimension>
::_QdiffTerm(
		const MatrixType& Xtest, const MatrixType& X,
		const MatrixList& Atest, const MatrixList& A,
		const MatrixList& Ttest, const MatrixList& T,
		TScalar stepXA, TScalar stepT)
 {
			
 	// for variables updated with a standard gradient step
 	TScalar valT = 0.0;
	if (Superclass::m_UpdateTemplate)
	{
	 	for (int i = 0; i < Superclass::m_NumberOfObjects; i++)
	 	{
	 		for (unsigned int d = 0; d < Dimension; d++)
	 		{
	 			valT += dot_product(m_GradTempL_L2[i].get_column(d),m_GradTempL_Sob[i].get_column(d)); 
	 		}
	 	}
	 	valT *= ( -0.5 * stepT );		
	}

	TScalar valX = 0.0;
 	if (Superclass::m_UpdateCP)
	{
		valX = m_GradPos.frobenius_norm();
	 	valX *= ( -0.5 * stepXA * valX );		
	}

 	// for variables updated with a soft-thresholded gradient step
 	MatrixList diffA(Superclass::m_NumberOfSubjects);
 	for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
 		diffA[s] = Atest[s] - A[s];

 	// Compute first term of the form (xtest - x)^t * Grad
 	TScalar termVal1 = 0.0;
 	for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
 	{
 		MatrixType& diffA_s = diffA[s];
 		const MatrixType& gradA_s = m_GradMom[s];
 		for (unsigned int i = 0; i < m_NumberOfCPs; i++)
 			termVal1 += dot_product(diffA_s.get_row(i), gradA_s.get_row(i));
 	}

 	// Compute second term of the form ||xtest - x||
 	TScalar termVal2 = 0.0;
 	for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
 	{
 		MatrixType& diffA_s = diffA[s];
 		for (unsigned int i = 0; i < m_NumberOfCPs; i++)
 			termVal2 += dot_product(diffA_s.get_row(i), diffA_s.get_row(i));
 	}

 	// sum all this
 	TScalar val = termVal1 + termVal2 / (2.0 * stepXA) + valT + valX;

 	return val;
}



#endif /* _DeterministicAtlasEstimator_txx */

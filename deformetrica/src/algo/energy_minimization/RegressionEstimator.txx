/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _RegressionEstimator_txx
#define _RegressionEstimator_txx

#include "RegressionEstimator.h"
 
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
RegressionEstimator<TScalar, Dimension>
::RegressionEstimator()
{
	this->m_SaveEveryNIters = 65535;	// by default don't save during optimization
	this->SetMaxIterations(100);
	m_MaxLineSearchIterations = 20;

	m_AdaptiveExpand = 1.2;
	m_AdaptiveShrink = 0.5;
	m_AdaptiveTolerance = 1e-4;

	m_InitialStepMultiplier = 1.0;

	m_UpdateCP = true;		
	m_UpdateTemplate = true;
	m_OptimizationMethod = null;
}



template <class TScalar, unsigned int Dimension>
RegressionEstimator<TScalar, Dimension>
::~RegressionEstimator()
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
RegressionEstimator<TScalar, Dimension>
::WriteOutput()
{
	this->WriteOutput(this->m_RegressionName);
}


template <class TScalar, unsigned int Dimension>
void
RegressionEstimator<TScalar, Dimension>
::WriteOutput(std::string regressionName)
{

	MatrixType mom = this->GetMomenta();
	std::vector< std::string > tempfn(this->m_NumberOfObjects);
	for (unsigned int i = 0; i < this->m_NumberOfObjects; i++)
	{
		std::ostringstream oss;
		oss << this->m_TemplateObjectsName[i] << "_trajectory_";
		tempfn[i] = oss.str();
	}
	
	m_Regression->WriteTemplateFlow(mom, regressionName, tempfn, this->m_TemplateObjectsNameExtension);

	m_Regression->WriteRegressionParameters(mom, regressionName);
}


template <class TScalar, unsigned int Dimension>
void
RegressionEstimator<TScalar, Dimension>
::Update()
{
	if ( (m_InitialMomenta.rows() == m_NumberOfCPs) && (m_InitialMomenta.columns() == Dimension) )
	{
		std::cout << "Using predefined set of momenta" << std::endl;
	}
	else
	{
		if (m_InitialMomenta.rows() > 0)
			std::cout << "Warning: initial momenta file has incompatible number of vectors. Initial momenta reset to zero" << std::endl;

		m_InitialMomenta.set_size(m_NumberOfCPs, Dimension);
		m_InitialMomenta.fill(0.0);
		
	}
	
	m_ValuesHistory.resize(this->m_MaxIterations + 1);
	
	MatrixType ControlPoints = m_Regression->GetControlPoints();
	MatrixList TemplateData = m_Regression->GetTemplateData();
		
	MatrixType& X = ControlPoints;
	MatrixType& A = m_InitialMomenta;
	MatrixList& T = TemplateData;
	
	// Run optimization method
	if ( (m_OptimizationMethod == GradientDescent) && (m_Regression->GetSparsityPrior() > 1e-10) )
	{
		std::cout << "Warning: gradient descent does not work with L^1 penalty term. Optimization method switched to ISTA" << std::endl;
		m_OptimizationMethod = ISTA;
	}
		
	if (m_OptimizationMethod == F_ISTA)
		this->FISTA(X, A, T);
	else
		this->GradientDescentAndISTA(X, A, T);
	
	if (this->m_UpdateCP)
		m_Regression->SetControlPoints(ControlPoints);
	if (this->m_UpdateTemplate)
	{
		m_Regression->SetTemplateData(TemplateData);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
template <class TScalar, unsigned int Dimension>
void
RegressionEstimator<TScalar, Dimension>
::GradientDescentAndISTA(MatrixType& X, MatrixType& A, MatrixList& T)
{
	MatrixType Xnew;
	MatrixType Anew;
	MatrixList Tnew(this->m_NumberOfObjects);
	
	TScalar stepXA, stepT;
	
	if (this->m_UpdateCP)
		m_Regression->SetControlPoints(X);
	if (this->m_UpdateTemplate)
		m_Regression->SetTemplateData(T);
	
	// Get the value of the functional for the initial momentas
	m_ValuesHistory[0] = m_Regression->ComputeFunctional(A, this->m_TargetList, this->m_TimeIndices);
	m_ValuesHistory[0]->PrintIter(0);
	
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();
	unsigned int iterRef = 0;

	for (unsigned int iter = 0; iter < this->m_MaxIterations; iter++)
	{
		if ((!((iter+1) % m_SaveEveryNIters)) && (iter != 0))
 		{
			std::ostringstream oss;
			oss << "Iter" << std::setfill('0') << std::setw(4) << iter << "_" << this->m_RegressionName;
			this->WriteOutput(oss.str());
 		}

		if (this->m_UpdateCP)
			m_Regression->SetControlPoints(X);
		if (this->m_UpdateTemplate)
			m_Regression->SetTemplateData(T);
	
		if (this->m_UpdateTemplate)
			m_Regression->ComputeFunctionalGradient(A, this->m_TargetList, this->m_TimeIndices, m_GradMom, m_GradPos, m_GradTempL_L2, m_GradTempL_Sob);
		else
			m_Regression->ComputeFunctionalGradient(A, this->m_TargetList, this->m_TimeIndices, m_GradMom, m_GradPos, m_GradTempL_L2, false);
		
		if (iter == 0)
		{
			TScalar maxGrad = 1e-20;
			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			{
				TScalar g = 0;
				g = m_GradPos.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
				g = m_GradMom.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
			}

			TScalar initStepXA = this->m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			TScalar initStepT = initStepXA;

			stepXA = initStepXA;
			stepT = initStepT;
		}

		std::vector<int> activeCPMask(m_NumberOfCPs, 1);
		bool foundMin = false;

		for (unsigned int li = 0; li < this->m_MaxLineSearchIterations; li++)
		{
		 	std::cout << "stepsizeXA = " << stepXA;
			if (this->m_UpdateTemplate)
				std::cout << "\tstepsizeT = " << stepT;
			std::cout << std::endl;

			this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
			
			if (this->m_UpdateCP)
				m_Regression->SetControlPoints(Xnew);
			if (this->m_UpdateTemplate)
				m_Regression->SetTemplateData(Tnew);
			
			
			m_ValuesHistory[iter+1] = m_Regression->ComputeFunctional(Anew, this->m_TargetList, this->m_TimeIndices);


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
				stepXA *= this->m_AdaptiveShrink;
				
				this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
				
				if (this->m_UpdateCP)
					m_Regression->SetControlPoints(Xnew);
				if (this->m_UpdateTemplate)
					m_Regression->SetTemplateData(Tnew);

				FunctionalValuesType* valuesTest1 = m_Regression->ComputeFunctional(Anew, this->m_TargetList, this->m_TimeIndices);

				TScalar Q1 = lsqRef - valuesTest1->GetTotalL2Value();
				if (m_OptimizationMethod == ISTA)
					Q1 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

				// case test 2:
				TScalar Q2 = -1.0;
				stepXA /= this->m_AdaptiveShrink;
				FunctionalValuesType* valuesTest2;
				if (this->m_UpdateTemplate)
				{
					stepT *= this->m_AdaptiveShrink;
				
					this->GradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
					m_Regression->SetControlPoints(Xnew);
					m_Regression->SetTemplateData(Tnew);
					valuesTest2 = m_Regression->ComputeFunctional(Anew, this->m_TargetList, this->m_TimeIndices);

					Q2 = lsqRef - valuesTest2->GetTotalL2Value();
					if (m_OptimizationMethod == ISTA)
						Q2 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

					//std::cout << "test case 2 = " << Q2 << std::endl;
				}

				
				if ( (Q1 >= 0) || (Q2 >= 0) ) //( (Jtest1 < Jcurr) || (Jtest2 < Jcurr) )
				{
					if ( Q1 >= Q2 )
					{
						stepXA *= this->m_AdaptiveShrink;
						stepT /= this->m_AdaptiveShrink;
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
					stepXA *= this->m_AdaptiveShrink;
				}

				delete valuesTest1;
				if (this->m_UpdateTemplate)
					 delete valuesTest2;

			}
		} // for li


		if (foundMin)
		{
			X = Xnew;
			A = Anew;
			T = Tnew;

			stepXA *= this->m_AdaptiveExpand;
			stepT *= this->m_AdaptiveExpand;
			lsqRef = m_ValuesHistory[iter+1]->GetTotalL2Value();
		}

		if (!foundMin)
		{
			// Loop terminated without finding smaller functional
			std::cout << " number of loops exceeded " << std::endl;
			break;
		}

		if ( (m_ValuesHistory[iter]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue())
				< this->m_AdaptiveTolerance * (m_ValuesHistory[iterRef]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue()) )
			break;

		// Jcurr = J1;

		m_ValuesHistory[iter+1]->PrintIter(iter+1);

	} // for iter
}



template <class TScalar, unsigned int Dimension>
void
RegressionEstimator<TScalar, Dimension>
::FISTA(MatrixType& X, MatrixType& A, MatrixList& T)
{
	MatrixType Xtest;
 	MatrixType Atest;
 	MatrixList Ttest(this->m_NumberOfObjects);
	
 	TScalar stepXA, stepT;
	
	if (this->m_UpdateCP)
		m_Regression->SetControlPoints(X);
	if (this->m_UpdateTemplate)
		m_Regression->SetTemplateData(T);

 	// Get the value of the functional for the initial momentas	
 	m_ValuesHistory[0] = m_Regression->ComputeFunctional(A, this->m_TargetList, this->m_TimeIndices);	
	m_ValuesHistory[0]->PrintIter(0);

 	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();
	
	MatrixType Xprev = X;
	MatrixType Aprev = A;
	MatrixList Tprev = T;

	std::vector<int> activeCPMask(m_NumberOfCPs, 1);
	std::vector<int> activeCPMask_test(m_NumberOfCPs, 1);

	TScalar Qarray[this->m_UpdateTemplate?4:2];
	unsigned int freezeDirectionCounter = 0;
	unsigned int iterRef = 0;

	TScalar tau = 1.0;
 	for (unsigned int iter = 0; iter < this->m_MaxIterations; iter++)
 	{
		if ((!((iter+1) % m_SaveEveryNIters)) && (iter != 0))
 		{
			std::ostringstream oss;
			oss << "Iter" << std::setfill('0') << std::setw(4) << iter << "_" << this->m_RegressionName;
			this->WriteOutput(oss.str());
 		}

		if (this->m_UpdateCP)
			m_Regression->SetControlPoints(X);
		if (this->m_UpdateTemplate)
 			m_Regression->SetTemplateData(T);
	
		if (this->m_UpdateTemplate)
 			m_Regression->ComputeFunctionalGradient(A, this->m_TargetList, this->m_TimeIndices, m_GradMom, m_GradPos, m_GradTempL_L2, m_GradTempL_Sob);
		else
			m_Regression->ComputeFunctionalGradient(A, this->m_TargetList, this->m_TimeIndices, m_GradMom, m_GradPos, m_GradTempL_L2, 0);

 		if (iter == 0)
 		{
 			TScalar maxGrad = 1e-20;
 			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
 			{
 				TScalar g = 0;
 				g = m_GradPos.get_row(i).magnitude();
 				if (g > maxGrad)
 					maxGrad = g;
 				g = m_GradMom.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
 				
 			}

 			TScalar initStepXA = this->m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;		
			TScalar initStepT = initStepXA;

 			stepXA = initStepXA;
 			stepT = initStepT;
 		}

		bool minimtest = false;
		unsigned int lineIter = 0;
		for (; lineIter < this->m_MaxLineSearchIterations; lineIter++)
		{
		 	std::cout << "stepsizeXA = " << stepXA;
			if (this->m_UpdateTemplate)
				std::cout << "\tstepsizeT = " << stepT;
			std::cout << std::endl;

			
		 	this->GradientDescentStep(Xtest, X, Atest, A, Ttest, T, activeCPMask_test, stepXA, stepT);

			if (this->m_UpdateCP)
				m_Regression->SetControlPoints(Xtest);
			if (this->m_UpdateTemplate)
				m_Regression->SetTemplateData(Ttest);

			m_ValuesHistory[iter+1] = m_Regression->ComputeFunctional(Atest, this->m_TargetList, this->m_TimeIndices);

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
		 		stepXA *= this->m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime1(m_NumberOfCPs, 1);

				MatrixType X_prime1;
				MatrixType A_prime1;

		 		this->GradientDescentStep(X_prime1, X, A_prime1, A, Ttest, T, activeCPMask_prime1, stepXA, stepT);
				
				if (this->m_UpdateCP)
		 			m_Regression->SetControlPoints(X_prime1);
				if (this->m_UpdateTemplate)
		 			m_Regression->SetTemplateData(Ttest);
		 		
				values_XAprime1 = m_Regression->ComputeFunctional(A_prime1, this->m_TargetList, this->m_TimeIndices);

				Q_X1 = lsqRef - values_XAprime1->GetTotalL2Value() + this->_QdiffTerm(
						X_prime1, X, A_prime1, A, Ttest, T, stepXA, stepT);

		 				//std::cout << "test case 0 = " << Q_X1 << std::endl;

 				// case test 1:
				stepXA *= this->m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime2(m_NumberOfCPs, 1);

				MatrixType X_prime2;
				MatrixType A_prime2;

				this->GradientDescentStep(X_prime2, X, A_prime2, A, Ttest, T, activeCPMask_prime2, stepXA, stepT);
				
				if (this->m_UpdateCP)
 					m_Regression->SetControlPoints(X_prime2);
				if (this->m_UpdateTemplate)
					m_Regression->SetTemplateData(Ttest);
 				
				values_XAprime2 = m_Regression->ComputeFunctional(A_prime2, this->m_TargetList, this->m_TimeIndices);

				Q_X2 = lsqRef - values_XAprime2->GetTotalL2Value() + this->_QdiffTerm(
						X_prime2, X, A_prime2, A, Ttest, T, stepXA, stepT);

		 				//std::cout << "test case 1 = " << Q2 << std::endl;

				// test case 2
				stepXA /= (this->m_AdaptiveShrink * this->m_AdaptiveShrink);
				MatrixList T_prime1(this->m_NumberOfObjects);
				MatrixList T_prime2(this->m_NumberOfObjects);
				if (this->m_UpdateTemplate)
				{
					stepT *= this->m_AdaptiveShrink;

					this->GradientDescentStep(Xtest, X, Atest, A, T_prime1, T, activeCPMask_test, stepXA, stepT);
					
					m_Regression->SetControlPoints(Xtest);
					m_Regression->SetTemplateData(T_prime1);
			 			
					values_Tprime1 = m_Regression->ComputeFunctional(Atest, this->m_TargetList, this->m_TimeIndices);

					Q_T1 = lsqRef - values_Tprime1->GetTotalL2Value() + this->_QdiffTerm(
							Xtest, X, Atest, A, T_prime1, T, stepXA, stepT);


					// test case 3
					stepT *= this->m_AdaptiveShrink;
					this->GradientDescentStep(Xtest, X, Atest, A, T_prime2, T, activeCPMask_test, stepXA, stepT);
					
	 				m_Regression->SetControlPoints(Xtest);
	 				m_Regression->SetTemplateData(T_prime2);
	 				values_Tprime2 = m_Regression->ComputeFunctional(Atest, this->m_TargetList, this->m_TimeIndices);

					Q_T2 = lsqRef - values_Tprime2->GetTotalL2Value() + this->_QdiffTerm(
							Xtest, X, Atest, A, T_prime2, T, stepXA, stepT);
				}


				// Select best case
				Qarray[0] = Q_X1;
				Qarray[1] = Q_X2;
				if (this->m_UpdateTemplate)
				{
					Qarray[2] = Q_T1;
					Qarray[3] = Q_T2;
				}

				int ind_max = 0;
				TScalar Q_max = Qarray[0];
				for (int j = 1; j < ( (this->m_UpdateTemplate)?4:2 ); j++)
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

						stepXA *= this->m_AdaptiveShrink;
						stepT /= (this->m_AdaptiveShrink * this->m_AdaptiveShrink);
					}
					else if (ind_max == 1)
					{
						// std::cout << " Q_X2" << std::endl;
						Xtest = X_prime2;
						Atest = A_prime2;
						activeCPMask_test = activeCPMask_prime2;
						m_ValuesHistory[iter+1] = values_XAprime2->Clone();

						stepXA *= (this->m_AdaptiveShrink * this->m_AdaptiveShrink);
						stepT /= (this->m_AdaptiveShrink * this->m_AdaptiveShrink);
					}
					else if (ind_max == 2)
					{
						// std::cout << " Q_T1" << std::endl;
						Ttest = T_prime1;
						m_ValuesHistory[iter+1] = values_Tprime1->Clone();

						stepT /= this->m_AdaptiveShrink;
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
					stepXA *= this->m_AdaptiveShrink;
					stepT /= this->m_AdaptiveShrink;
				}

				delete values_XAprime1;
				delete values_XAprime2;
				if (this->m_UpdateTemplate)
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

			if ( (lineIter == 0) && ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze)) )
			{
				// stepsize increased
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau / this->m_AdaptiveExpand)) / 2.0;
				if (freezeDirectionCounter > NbIterFreeze)
					freezeDirectionCounter = 0;
			}
			else // stepsize decreased
			{
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau)) / 2.0;
				freezeDirectionCounter++;
			}

			TScalar tau_scale = (tau-1.0) / tau_next;
			if (this->m_UpdateCP)
				X = Xtest + (Xtest - Xprev) * tau_scale;
			A = Atest + (Atest - Aprev) * tau_scale;

			if (this->m_UpdateTemplate)
			{
				for (unsigned int i = 0; i < this->m_NumberOfObjects; i++)
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
				stepXA *= this->m_AdaptiveExpand;
				stepT *= this->m_AdaptiveExpand;
			}

			// Update the reference L2 part of the functional
			if (this->m_UpdateCP)
				m_Regression->SetControlPoints(X);
			if (this->m_UpdateTemplate)
				m_Regression->SetTemplateData(T);
			
			FunctionalValuesType* values = m_Regression->ComputeFunctional(A, this->m_TargetList, this->m_TimeIndices);

			if (values->IsOutOfBox()) //m_Def->OutOfBox() || m_Template->OutOfBox())
			{
				std::cerr << "Out of box: needs to restart FISTA from this point." << std::endl;
				X = Xtest; A = Atest; 
				if (this->m_UpdateTemplate)
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
		if (fabs(deltaF_curr) < this->m_AdaptiveTolerance*fabs(deltaF_ref))
		{
			std::cout << "Tolerance BREAK" << std::endl;
			std::cout << "FINAL VALUES: ";
			m_ValuesHistory[iter+1]->PrintIter(iter+1);
			break;
		}

	} // end iter

	X = this->_maskMatrix(X, activeCPMask);
	A = this->_maskMatrix(A, activeCPMask);
}




template <class TScalar, unsigned int Dimension>
typename RegressionEstimator<TScalar, Dimension>::MatrixType
RegressionEstimator<TScalar, Dimension>
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
RegressionEstimator<TScalar, Dimension>
::GradientDescentStep(
		MatrixType& Xtest, const MatrixType& X,
		MatrixType& Atest, const MatrixType& A,
		MatrixList& Ttest, const MatrixList& T,
		std::vector< int >& activeCPMask,
		TScalar stepXA, TScalar stepT)
{
	if (this->m_UpdateCP)
		Xtest = X - m_GradPos*stepXA;
	else
		Xtest = X;

	if (m_Regression->GetSparsityPrior() < 1e-10)
	{
		activeCPMask.assign(m_NumberOfCPs, 1);
		Atest = A - m_GradMom * stepXA;
	}
	else
	{
		activeCPMask.assign(m_NumberOfCPs, 0); // gives for each CP the number of non-zero momenta that it carries
		Atest = this->SoftThresholdUpdate(A, m_GradMom, stepXA, activeCPMask);
	}

	if (this->m_UpdateTemplate)
	{
		for (int i = 0; i < this->m_NumberOfObjects; i++)
			Ttest[i] = T[i] - m_GradTempL_Sob[i]*stepT;		
	}
	
	return;

}



template <class TScalar, unsigned int Dimension>
typename RegressionEstimator<TScalar, Dimension>::MatrixType
RegressionEstimator<TScalar, Dimension>
::SoftThresholdUpdate(const MatrixType& X, MatrixType& gradX, TScalar step, std::vector<int>& mask)
{
	MatrixType Xnew = X - step*gradX;
	TScalar thres = step * m_Regression->GetSparsityPrior();

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
RegressionEstimator<TScalar, Dimension>
::_QdiffTerm(const MatrixType& Xtest, const MatrixType& X,
	     const MatrixType& Atest, const MatrixType& A,
	     const MatrixList& Ttest, const MatrixList& T,
	     TScalar stepXA, TScalar stepT)
{
	// For variables updated with a standard gradient step
 	TScalar valT = 0.0;
	if (this->m_UpdateTemplate)
	{
	 	for (int i = 0; i < this->m_NumberOfObjects; i++)
	 	{
	 		for (unsigned int d = 0; d < Dimension; d++)
	 		{
	 			valT += dot_product(m_GradTempL_L2[i].get_column(d),m_GradTempL_Sob[i].get_column(d)); 
	 		}
	 	}
	 	valT *= ( -0.5 * stepT );		
	}

	TScalar valX = 0.0;
 	if (this->m_UpdateCP)
	{
		valX = m_GradPos.frobenius_norm();
	 	valX *= ( -0.5 * stepXA * valX );		
	}

 	// for variables updated with a soft-thresholded gradient step
 	MatrixType diffA;
 	diffA = Atest - A;

 	// Compute first term of the form (xtest - x)^t * Grad
 	TScalar termVal1 = 0.0;
 	MatrixType& diffA_s = diffA;
 	const MatrixType& gradA_s = m_GradMom;
 	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
 		termVal1 += dot_product(diffA_s.get_row(i), gradA_s.get_row(i));
 	
 	// Compute second term of the form ||xtest - x||
 	TScalar termVal2 = 0.0;
 	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
 		termVal2 += dot_product(diffA_s.get_row(i), diffA_s.get_row(i));
 	
 	// sum all this
 	TScalar val = termVal1 + termVal2 / (2.0 * stepXA) + valT + valX;

 	return val;
}



#endif

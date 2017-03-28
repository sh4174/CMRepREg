/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeterministicAtlas_txx
#define _DeterministicAtlas_txx

#include "DeterministicAtlas.h"
#include "DeformableMultiObject.h"

#include "KernelFactory.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeterministicAtlas<TScalar, Dimension>
::DeterministicAtlas() : Superclass()
{
	this->SetDeterministicAtlasType();
	m_UseRKHSNormForRegularization = true;
	m_SparsityPrior = 0.0;
}



template <class TScalar, unsigned int Dimension>
DeterministicAtlas<TScalar, Dimension>
::~DeterministicAtlas()
{}



template <class TScalar, unsigned int Dimension>
DeterministicAtlas<TScalar, Dimension>
::DeterministicAtlas(const DeterministicAtlas& other) : Superclass(other)
{
	m_UseRKHSNormForRegularization = other.m_UseRKHSNormForRegularization;
	m_SparsityPrior = other.m_SparsityPrior;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeterministicAtlas<TScalar, Dimension>
::Update()
{
	Superclass::Update();

	if ( !m_UseRKHSNormForRegularization )
	{
		if ( Superclass::m_CovMomInverse.rows() == 0 )
		{
			std::cout << "Warning: no covariance momenta matrix set for computing regularization term. We will use the RKHS norm for computing the regularization term" << std::endl;
			m_UseRKHSNormForRegularization = true;
		}
		if ( Superclass::m_CovMomInverse.rows() != Dimension*Superclass::m_ControlPoints.rows() )
		{
		    std::cout << " Warning: the size of the covariance momenta matrix that you set does not match the number of control points. Matrix set to identity" << std::endl;
			Superclass::m_CovMomInverse.set_size(Superclass::m_ControlPoints.rows()*Dimension, Superclass::m_ControlPoints.rows()*Dimension);
		 	Superclass::m_CovMomInverse.set_identity();
		}
	}
}



template <class TScalar, unsigned int Dimension>
typename DeterministicAtlas<TScalar, Dimension>::FunctionalValuesType*
DeterministicAtlas<TScalar, Dimension>
::ComputeFunctional(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target)
{
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of Momentas and target multi-object mismatch in DeterministicAtlas::ComputeFunctional");

	this->UpdateDeformationAndKernelDataDomain(target);

	int numSubjects = Momentas.size();

	FunctionalValuesType* out = new FunctionalValuesType();

	// compute \forall k\in Objects \forall s\in Subjects \norm(\phi^{Momentas_s}(template_k) - target_{k,s})^2 / (2.0 * DataSigmaSquared_k^2)
	std::vector< std::vector< TScalar > > Residuals(numSubjects);

	bool oob = this->ComputeResiduals(Momentas, target, Residuals);
	if (oob)
	{
		out->SetOutOfBox();
		return out;
	}

	//cout<<"Printing residuals...\n";
	//for(int i=0; i<Residuals[0].size(); ++i)
  	//	std::cout << Residuals[0][i] << '\n';


	// DataTerm = sum of residuals
	TScalar DataTerm = 0.0;
	for (unsigned int s = 0; s < numSubjects; s++)
	{
		for (unsigned int i = 0; i < Superclass::m_NumberOfObjects; i++)
			 DataTerm += Residuals[s][i] / (2.0 * Superclass::m_DataSigmaSquared(i)) ;
	}

	//cout<<"Data sigma squared "<<Superclass::m_DataSigmaSquared(0)<<"\n";

	// Regularity term
	TScalar RegTerm = 0.0;
	if ( m_UseRKHSNormForRegularization ) // use the RKHS norm
	{
		typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
		typedef typename KernelFactoryType::KernelBaseType KernelType;
		KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	    KernelType* momKernelObj  = kfac->CreateKernelObject(Superclass::m_Def->GetKernelType());
		momKernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());
		momKernelObj->SetSources(Superclass::m_ControlPoints);

		for (unsigned int s = 0; s < numSubjects; s++)
		{
			momKernelObj->SetWeights(Momentas[s]);

			MatrixType kMom = momKernelObj->Convolve(Superclass::m_ControlPoints);

			for (unsigned int i = 0; i < Superclass::m_ControlPoints.rows(); i++)
				RegTerm += dot_product(kMom.get_row(i), Momentas[s].get_row(i));
		}

		delete momKernelObj;
	}
	else // covariance matrix given
	{
		for (unsigned int s = 0; s < numSubjects; s++)
		{
			VectorType Moms = this->Vectorize(Momentas[s]);
			RegTerm += dot_product( Moms, Superclass::m_CovMomInverse * Moms );
		}

	}
	RegTerm *= 0.5;

	// Sparsity Prior
	TScalar SparsityTerm = 0.0;
	for (unsigned int s = 0; s < numSubjects; s++)
	{
		for (unsigned int i = 0; i < Superclass::m_ControlPoints.rows(); i++)
			SparsityTerm += Momentas[s].get_row(i).magnitude();
	}
	SparsityTerm *= m_SparsityPrior;

	//cout<<"Data term = "<<DataTerm<<"\n";
	//cout<<"Reg term = "<<RegTerm<<"\n";
	//cout<<"Residuals = "<<Residuals<<"\n";

	out->SetDataTerm(DataTerm);
	out->SetRegularityTerm(RegTerm);
	out->SetSparsityTerm(SparsityTerm);
	out->SetResiduals(Residuals);

	return out;

}



template <class TScalar, unsigned int Dimension>
void
DeterministicAtlas<TScalar, Dimension>
::ComputeFunctionalGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> target,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL_L2, MatrixList& GradTempL_Sob)
{
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of Momentas and target multi-object mismatch in DeterministicAtlas::ComputeFunctionalGradient");

	this->UpdateDeformationAndKernelDataDomain(target);
	
	this->ComputeDataTermGradient(Momentas, target, GradPos, GradMom, GradTempL_L2, GradTempL_Sob);

	this->AddGradientRegularityTerm(GradPos, GradMom, Momentas);
}


template <class TScalar, unsigned int Dimension>
void
DeterministicAtlas<TScalar, Dimension>
::ComputeFunctionalGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> target,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL, bool do_SobolevGradient)
{
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of Momentas and target multi-object mismatch in DeterministicAtlas::ComputeFunctionalGradient");
	
	this->UpdateDeformationAndKernelDataDomain(target);

	this->ComputeDataTermGradient(Momentas, target, GradPos, GradMom, GradTempL, do_SobolevGradient);

	this->AddGradientRegularityTerm(GradPos, GradMom, Momentas);

}


template <class TScalar, unsigned int Dimension>
void
DeterministicAtlas<TScalar, Dimension>
::AddGradientRegularityTerm(MatrixType GradPos, MatrixList GradMom, const MatrixList& Momentas)
{

	int numSubjects = Momentas.size();

	// gradient of the regularity term
	if ( m_UseRKHSNormForRegularization ) // use the RKHS norm
	{
		typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
		typedef typename KernelFactoryType::KernelBaseType KernelType;
		KernelFactoryType* kfac = KernelFactoryType::Instantiate();

		KernelType* momKernelObj  = kfac->CreateKernelObject(Superclass::m_Def->GetKernelType());
		momKernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());
		momKernelObj->SetSources(Superclass::m_ControlPoints);

		for (unsigned int s = 0; s < numSubjects; s++)
		{
			momKernelObj->SetWeights(Momentas[s]);
			GradMom[s] += momKernelObj->Convolve(Superclass::m_ControlPoints);

			GradPos += momKernelObj->ConvolveGradient(Superclass::m_ControlPoints, Momentas[s]);
		}

		delete momKernelObj;
	}
	else // covariance matrix given
	{
		for (unsigned int s = 0; s < numSubjects; s++)
		{
			VectorType aux = Superclass::m_CovMomInverse * this->Vectorize(Momentas[s]);
			GradMom[s] += this->VectorToMatrix( aux );
		}
	}

}



#endif /* _DeterministicAtlas_txx */

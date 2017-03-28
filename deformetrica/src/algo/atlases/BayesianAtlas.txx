/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BayesianAtlas_txx
#define _BayesianAtlas_txx

#include "BayesianAtlas.h"
#include "DeformableMultiObject.h"

#include "KernelFactory.h"
#include "ExactKernel.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
BayesianAtlas<TScalar, Dimension>
::BayesianAtlas() : Superclass()
{
	this->SetBayesianAtlasType();

	m_CovMom_HyperParameter = 1.0;
	
	m_CovMom_Prior_Inverse.set_size(0,0);
	
	m_DataSigmaSquared_HyperParameter.set_size(Superclass::m_NumberOfObjects);
	m_DataSigmaSquared_HyperParameter.fill(1.0);
	
	m_DataSigmaSquared_Prior.set_size(Superclass::m_NumberOfObjects);
	m_DataSigmaSquared_Prior.fill(1.0);
}


template <class TScalar, unsigned int Dimension>
BayesianAtlas<TScalar, Dimension>
::~BayesianAtlas()
{}



template <class TScalar, unsigned int Dimension>
BayesianAtlas<TScalar, Dimension>
::BayesianAtlas(const BayesianAtlas& other) : Superclass(other)
{
	m_CovMom_HyperParameter = other.m_CovMom_HyperParameter;
 	m_CovMom_Prior_Inverse = other.m_CovMom_Prior_Inverse;
 	m_DataSigmaSquared_HyperParameter = other.m_DataSigmaSquared_HyperParameter;
 	m_DataSigmaSquared_Prior = other.m_DataSigmaSquared_Prior;
 	m_NoiseDimension = other.m_NoiseDimension;
}

template <class TScalar, unsigned int Dimension>
void
BayesianAtlas<TScalar, Dimension>
::SetCovarianceMomenta_Prior_Inverse(std::string fn)
{
	if (strlen(fn.c_str()))
	{
		m_CovMom_Prior_Inverse = readMatrixDLM<TScalar>(fn.c_str());
		std::cout << "Using a prior for covariance momenta of size " << m_CovMom_Prior_Inverse.rows() << " x " << m_CovMom_Prior_Inverse.columns() << std::endl;
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
BayesianAtlas<TScalar, Dimension>
::Update()
 {
	 Superclass::Update();

	 if (m_CovMom_Prior_Inverse.rows() == 0)
	 {
		 std::cout << "no prior for momenta covariance matrix given: set to kernel matrix with deformation kernel width = " << Superclass::m_Def->GetKernelWidth() << std::endl;
		 
		 this->InitializeCovMomPrior();
	 }

	Superclass::m_CovMomInverse.set_size(Superclass::m_ControlPoints.rows()*Dimension, Superclass::m_ControlPoints.rows()*Dimension);
	
	// We need this "initialization", since we call Atlas::ComputeResiduals before computing the optimal data sigma squared
//	Superclass::m_DataSigmaSquared.set_size(Superclass::m_NumberOfObjects);
//	Superclass::m_DataSigmaSquared = 1.0; // All sigma squared equal to 1
	 
}


template <class TScalar, unsigned int Dimension>
void
BayesianAtlas<TScalar, Dimension>
::WriteAtlasParameters(const std::string& AtlasName)
{

	Superclass::WriteAtlasParameters(AtlasName);

	// write data sigma
	VectorType DSS = this->GetDataSigmaSquared();
	MatrixType DS(1,DSS.size());
	for (int i = 0; i < DSS.size(); i++)
		DS(0,i) = sqrt(DSS(i));
	
	std::ostringstream ossDSS;
	ossDSS << AtlasName << "_DataSigma.txt" << std::ends;
	writeMatrixDLM<TScalar>(ossDSS.str().c_str(), DS);
	
	// write inverse of optimal covariance momenta matrix
	MatrixType CovMomInv = this->GetCovarianceMomentaInverse();
	std::ostringstream ossCMI;
	ossCMI << AtlasName << "_CovarianceMomentaInverse.txt" << std::ends;
	writeMatrixDLM<TScalar>(ossCMI.str().c_str(), CovMomInv);

}



template <class TScalar, unsigned int Dimension>
void
BayesianAtlas<TScalar, Dimension>
::InitializeCovMomPrior()
{
	typedef ExactKernel<TScalar, Dimension> KernelType;
	    KernelType* ker = new KernelType();
	ker->SetKernelWidth( Superclass::m_Def->GetKernelWidth() );

	// be careful: this should be consistent with the way this->Vectorize works!
	int NumberOfCPs = Superclass::m_ControlPoints.rows();
	m_CovMom_Prior_Inverse.set_size(NumberOfCPs*Dimension, NumberOfCPs*Dimension);
	m_CovMom_Prior_Inverse.fill(0.0);

	for (unsigned int i = 0; i < NumberOfCPs; i++)
	{
		for (unsigned int j = 0; j < NumberOfCPs; j++)
		{
			VectorType CPi = Superclass::m_ControlPoints.get_row(i);
			VectorType CPj = Superclass::m_ControlPoints.get_row(j);

			TScalar k = ker->EvaluateKernel(CPi, CPj);
			MatrixType Aux = diagonal_matrix(Dimension, k);

			m_CovMom_Prior_Inverse.update(Aux, Dimension*i, Dimension*j);
			m_CovMom_Prior_Inverse.update(Aux, Dimension*j, Dimension*i);
		}
	}
	delete ker;
	
	
	// m_CovMom_Prior_Inverse.set_size(Superclass::m_ControlPoints.rows()*Dimension, Superclass::m_ControlPoints.rows()*Dimension);
	// m_CovMom_Prior_Inverse.fill(0.0);
	// m_CovMom_Prior_Inverse.fill_diagonal(1.0);
	
	// std::cout << log( determinant(m_CovMom_Prior_Inverse) ) << std::endl;
	
}


template <class TScalar, unsigned int Dimension>
typename BayesianAtlas<TScalar, Dimension>::FunctionalValuesType*
BayesianAtlas<TScalar, Dimension>
::ComputeLikelihoodWithOptimalParameters(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target)
{
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLikelihood");

	this->UpdateDeformationAndKernelDataDomain(target);
		
	int numSubjects = Momentas.size();

	FunctionalValuesType* out = new FunctionalValuesType();
	
	// compute \forall k\in Objects \forall s\in Subjects \norm(\phi^{Momentas_s}(template_k) - target_{k,s})^2
    std::vector< std::vector< TScalar > > Residuals(numSubjects);
	bool oob = this->ComputeResiduals(Momentas, target, Residuals);

	
	if (oob)
	{
		out->SetOutOfBox();
		return out;
	}

	VectorType OptimalDataSigmaSquared;
	MatrixType OptimalCovMomInverse;
    std::vector< std::vector< TScalar> > DataTerms;
    VectorType DataPriors, CovMomTerms;
    TScalar CovMomPrior;

//	TScalar llh_ds = this->ComputeOptimalDataSigmaSquared(Residuals, OptimalDataSigmaSquared, DataTerms, DataPriors);
//	TScalar llh_cm = this->ComputeOptimalCovMomInverse(Momentas, OptimalCovMomInverse, CovMomTerms, CovMomPrior);
    this->ComputeOptimalDataSigmaSquared(Residuals, OptimalDataSigmaSquared, DataTerms, DataPriors);
    this->ComputeOptimalCovMomInverse(Momentas, OptimalCovMomInverse, CovMomTerms, CovMomPrior);
	
	out->SetDataSigmaSquared(OptimalDataSigmaSquared);
	out->SetCovMomInverse(OptimalCovMomInverse);
    out->SetResiduals(Residuals);
	out->SetDataTerms(DataTerms);
    out->SetDataPriors(DataPriors);
    out->SetCovMomTerms(CovMomTerms);
    out->SetCovMomPrior(CovMomPrior);
    out->SetLikelihoodWithPriors();
	
    // comment TScalar llh_ds = ... and TScalar llh_cm = ... when commenting this line!
//    std::cout << "sanity check = is " << llh_ds + llh_cm << " equal to " << out->GetLikelihood() << std::endl;
    
	return out;
}


template <class TScalar, unsigned int Dimension>
typename BayesianAtlas<TScalar, Dimension>::FunctionalValuesType*
BayesianAtlas<TScalar, Dimension>
::ComputeLikelihoodWithoutPriors(const MatrixList& Momentas, const std::vector<DeformableMultiObjectType* > target)
{
    if (Momentas.size() != target.size())
        throw std::runtime_error("Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLikelihood");
    
    FunctionalValuesType* out = new FunctionalValuesType();

    this->UpdateDeformationAndKernelDataDomain(target);

    // compute \forall k\in Objects \forall s\in Subjects \norm(\phi^{Momentas_s}(template_k) - target_{k,s})^2
    std::vector< std::vector< TScalar > > Residuals;
    bool oob = this->ComputeResiduals(Momentas, target, Residuals);
    
    if (oob)
    {
        out->SetOutOfBox();
        return out;
    }
    
//    // Compute the log-determinant of the Optimal momenta covariance matrix
//    TScalar Log_det_Cov_Mom = 0.0;
//    if (NumberOfSubjects > 0)
//    {
//        vnl_symmetric_eigensystem<TScalar> CovMomInverse_eig(Superclass::m_CovMomInverse);
//        vnl_diag_matrix<TScalar> Lambda = CovMomInverse_eig.D;
//        VectorType temp = Lambda.diagonal();
//        TScalar LE =  temp.size();
//        
//        for (unsigned int jj=0; jj<LE; jj++)
//        {
//            if (temp[jj]<=0)
//                throw std::runtime_error("Eigenvalue of the initial momenta covariance matrix smaller or equal to zero!");
//            
//            Log_det_Cov_Mom -= log(temp[jj]);
//        }
//    }
//    
//    VectorType CovMomTerms(NumberOfSubjects);
//    for (unsigned int s = 0; s < NumberOfSubjects; s++)
//    {
//        VectorType Moms = this->Vectorize( Momentas[s] );
//        CovMomTerms(s) = dot_product( Moms, Superclass::m_CovMomInverse * Moms) + Log_det_Cov_Mom;
//    }
//    CovMomTerms *= 0.5;
//
//    std::vector< TScalar >DataTerms(NumberOfSubjects);
//    for (unsigned int s = 0; s < NumberOfSubjects; s++)
//        DataTerms[s].resize(NumberOfObjects);
//    
//    for (int i = 0; i < NumberOfObjects; i++)
//    {
//        TScalar DSS_i = Superclass::m_DataSigmaSquared(i);
//        for (unsigned int s = 0; s < NumberOfSubjects; s++)
//        {
//            DataTerms[s][i] = 0.5*( Residuals[s][i]/DSS_i + m_NoiseDimension[i]*log(DSS_i) );
//        }
//    }

    
    VectorType CovMomTerms;
    std::vector< std::vector< TScalar > > DataTerms;
    
    this->ComputeCovMomTerms(Momentas, CovMomTerms);
    this->ComputeDataTerms(Residuals, DataTerms);
    
    out->SetDataSigmaSquared(Superclass::m_DataSigmaSquared);
    out->SetCovMomInverse(Superclass::m_CovMomInverse);
    out->SetResiduals(Residuals);
    out->SetDataTerms(DataTerms);
    out->SetCovMomTerms(CovMomTerms);
    out->UnsetLikelihoodWithPriors();
    
    return out;

}



template <class TScalar, unsigned int Dimension>
TScalar
BayesianAtlas<TScalar, Dimension>
::ComputeCovMomTerms(const MatrixList& Momentas, VectorType& CovMomTerms)
{

    int NumberOfSubjects = Momentas.size();
    
    // Compute the log-determinant of the Optimal momenta covariance matrix
    TScalar Log_det_Cov_Mom = 0.0;
    if (NumberOfSubjects > 0)
    {
/////////////////////////////////// Before LinAlg ////////////////////////////////////////////
//        vnl_symmetric_eigensystem<TScalar> CovMomInverse_eig(Superclass::m_CovMomInverse);
//        vnl_diag_matrix<TScalar> Lambda = CovMomInverse_eig.D;
//        VectorType temp = Lambda.diagonal();
//        TScalar LE =  temp.size();
//
//        for (unsigned int jj=0; jj<LE; jj++)
//        {
//            if (temp[jj]<=0)
//                throw std::runtime_error("Eigenvalue of the initial momenta covariance matrix smaller or equal to zero!");
//
//            Log_det_Cov_Mom -= log(temp[jj]);
//        }
//////////////////////////////////////////////////////////////////////////////////////////////
        VectorType temp = eigenvalues_sym(Superclass::m_CovMomInverse);
        TScalar LE =  temp.size();

        for (unsigned int jj=0; jj<LE; jj++)
        {
            if (temp(jj)<=0)
                throw std::runtime_error("Eigenvalue of the initial momenta covariance matrix smaller or equal to zero!");

            Log_det_Cov_Mom -= log(temp(jj));
        }

    }
    
    CovMomTerms.set_size(NumberOfSubjects);
    for (unsigned int s = 0; s < NumberOfSubjects; s++)
    {
        VectorType Moms = this->Vectorize( Momentas[s] );
        CovMomTerms(s) = dot_product( Moms, Superclass::m_CovMomInverse * Moms) + Log_det_Cov_Mom;
    }
    CovMomTerms *= 0.5;
    
    return CovMomTerms.sum();
    
}


template <class TScalar, unsigned int Dimension>
TScalar
BayesianAtlas<TScalar, Dimension>
::ComputeDataTerms(const std::vector< std::vector< TScalar> >& Residuals, std::vector< std::vector<TScalar > >& DataTerms)
{
    
    int NumberOfSubjects = Residuals.size();
    int NumberOfObjects = Superclass::m_NumberOfObjects;
    
    DataTerms.resize(NumberOfSubjects);
    for (unsigned int s = 0; s < NumberOfSubjects; s++)
        DataTerms[s].resize(NumberOfObjects);
    
    TScalar out = 0.0;
    for (int i = 0; i < NumberOfObjects; i++)
    {
        TScalar DSS_i = Superclass::m_DataSigmaSquared(i);
        for (unsigned int s = 0; s < NumberOfSubjects; s++)
        {
            DataTerms[s][i] = 0.5*( Residuals[s][i]/DSS_i + m_NoiseDimension[i]*log(DSS_i) );
            out+= DataTerms[s][i];
        }
    }
    
    return out;
}




template <class TScalar, unsigned int Dimension>
TScalar
BayesianAtlas<TScalar, Dimension>
::ComputeOptimalCovMomInverse(const MatrixList& Momentas, MatrixType& OptimalCovMomInverse, VectorType& CovMomTerms, TScalar& CovMomPrior)
{
	int NumberOfSubjects = Momentas.size();
	int NumberOfCPs = Superclass::m_ControlPoints.rows();
    
    // Compute \sum_i \alpha_i\alpha_i^T
	MatrixType AlphaPart(NumberOfCPs*Dimension, NumberOfCPs*Dimension, 0.0);
	for (unsigned int s = 0; s < NumberOfSubjects; s++)
	{
        if ( Momentas[s].rows() != NumberOfCPs )
            throw std::runtime_error("Number of Control points and size of momenta mismatch in BayesianAtlas::ComputeOptimalCovMomInverse");

		MatrixType alpha_temp(NumberOfCPs*Dimension , 1, 0.0);
		alpha_temp.set_column(0, this->Vectorize(Momentas[s]));

		AlphaPart += ( alpha_temp * alpha_temp.transpose() );
		
	}

/////////////////////////////////// Before LinAlg ////////////////////////////////////////////
//    // Compute the inverse of the optimal momenta covariance matrix
//	MatrixType diag(NumberOfCPs*Dimension, NumberOfCPs*Dimension, 0.0);
//	diag.fill_diagonal(m_CovMom_HyperParameter);
//
//	MatrixType aux = vnl_matrix_inverse<TScalar>(diag + m_CovMom_Prior_Inverse*AlphaPart);
//	aux *= (m_CovMom_HyperParameter + NumberOfSubjects);
//
//	OptimalCovMomInverse =  aux * m_CovMom_Prior_Inverse;

//	// Compute the log-determinant of the Optimal momenta covariance matrix
//    vnl_symmetric_eigensystem<TScalar> OptimalCovMomInverse_eig(OptimalCovMomInverse);
//	vnl_diag_matrix<TScalar> Lambda = OptimalCovMomInverse_eig.D;
//	VectorType temp = Lambda.diagonal();
//	TScalar LE =  temp.size();
//	TScalar Log_det_Cov_Mom = 0.0;

//	for (unsigned int jj=0; jj<LE; jj++)
//	{
//		if (temp[jj]<=0)
//			throw std::runtime_error("Eigenvalue of the initial momenta covariance matrix smaller or equal to zero!");
//
//		Log_det_Cov_Mom -= log(temp[jj]);
//	}
/////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute the inverse of the optimal momenta covariance matrix
	MatrixType aux = inverse( diagonal_matrix(NumberOfCPs*Dimension, m_CovMom_HyperParameter) + m_CovMom_Prior_Inverse*AlphaPart );
	aux *= (m_CovMom_HyperParameter + NumberOfSubjects);
	
	OptimalCovMomInverse =  aux * m_CovMom_Prior_Inverse;

	// Compute the log-determinant of the Optimal momenta covariance matrix
    VectorType temp = eigenvalues_sym(OptimalCovMomInverse);
	TScalar LE =  temp.size();
	TScalar Log_det_Cov_Mom = 0.0;

	for (unsigned int jj=0; jj<LE; jj++)
	{
		if (temp[jj]<=0)
			throw std::runtime_error("Eigenvalue of the initial momenta covariance matrix smaller or equal to zero!");

		Log_det_Cov_Mom -= log(temp[jj]);
	}
    
	//std::cout << "log( vnl_determinant(OptimalCovMomInverse) ) = " << log( vnl_determinant(OptimalCovMomInverse) ) << std::endl;
	//std::cout << "Log_det_Cov_Mom = " << Log_det_Cov_Mom << std::endl;

//	TScalar PartialLikelihood_CovMom = m_CovMom_HyperParameter * trace(aux) + (NumberOfSubjects + m_CovMom_HyperParameter) * Log_det_Cov_Mom;
    CovMomPrior = 0.5 * m_CovMom_HyperParameter * (trace(aux) + Log_det_Cov_Mom);
    
    CovMomTerms.set_size(NumberOfSubjects);
	for (unsigned int s = 0; s < NumberOfSubjects; s++)
	{
		VectorType Moms = this->Vectorize( Momentas[s] );
        CovMomTerms(s) = dot_product( Moms, OptimalCovMomInverse * Moms) + Log_det_Cov_Mom;
	}
    CovMomTerms *= 0.5;
    
    TScalar PartialLikelihood_CovMom = CovMomTerms.sum() + CovMomPrior;
	return PartialLikelihood_CovMom;

}

template <class TScalar, unsigned int Dimension>
TScalar
BayesianAtlas<TScalar, Dimension>
::ComputeOptimalDataSigmaSquared(const std::vector< std::vector< TScalar> >& Residuals, VectorType& OptimalDataSigmaSquared, std::vector< std::vector<TScalar > >& DataTerms, VectorType& DataPriors)
{
	int NumberOfObjects = Superclass::m_Template->GetNumberOfObjects();
	int NumberOfSubjects = Residuals.size();

    // Compute optimal data sigma squared
	OptimalDataSigmaSquared.set_size(NumberOfObjects);
	OptimalDataSigmaSquared.fill(0.0);

	for (int s = 0; s < NumberOfSubjects; s++)
	{
		if (Residuals[s].size() != NumberOfObjects)
			throw std::runtime_error("Number of objects in Residual Norms and Template mismatch in BayesianAtlas::ComptueOptimalDataSigmaSquared");
				
		for (int i = 0; i < NumberOfObjects; i++)
			OptimalDataSigmaSquared(i) += Residuals[s][i];
	}

    
	for (int i = 0; i < NumberOfObjects; i++)
		OptimalDataSigmaSquared(i) = (OptimalDataSigmaSquared[i] + m_DataSigmaSquared_HyperParameter[i] * m_DataSigmaSquared_Prior[i]) / (m_DataSigmaSquared_HyperParameter[i] + m_NoiseDimension[i]*NumberOfSubjects);
	

    // Compute DataTerms and DataPriors
    TScalar PartialLikelihood_DataSigmaSquared = 0.0;
    DataPriors.set_size(NumberOfObjects);
    DataTerms.resize(NumberOfSubjects);
    for (unsigned int s = 0; s < NumberOfSubjects; s++)
        DataTerms[s].resize(NumberOfObjects);

    for (int i = 0; i < NumberOfObjects; i++)
	{
        TScalar OptimalDSS_i = OptimalDataSigmaSquared(i);
        for (unsigned int s = 0; s < NumberOfSubjects; s++)
        {
            DataTerms[s][i] = 0.5*( Residuals[s][i]/OptimalDSS_i + m_NoiseDimension[i]*log(OptimalDSS_i) );
        }
        DataPriors(i) = 0.5 * m_DataSigmaSquared_HyperParameter[i] * (
                            log(OptimalDSS_i) + m_DataSigmaSquared_Prior[i]/OptimalDSS_i );
        
        PartialLikelihood_DataSigmaSquared += (m_DataSigmaSquared_HyperParameter[i] + m_NoiseDimension[i]*NumberOfSubjects) * (log(OptimalDataSigmaSquared(i)) + 1 );
	}
	
	PartialLikelihood_DataSigmaSquared *= 0.5;

//    // sanity check
//    TScalar total = 0.0;
//    for (int i = 0; i < NumberOfObjects; i++)
//    {
//        for (unsigned int s = 0; s < NumberOfSubjects; s++)
//            total += DataTerms[s][i];
//
//        total+= DataPriors(i);
//    }
//    std::cout << "PartialLikelihood_DataSigmaSquared: is " << PartialLikelihood_DataSigmaSquared << " equal to " << total << " ?" << std::endl;

	return PartialLikelihood_DataSigmaSquared;
}



template <class TScalar, unsigned int Dimension>
void
BayesianAtlas<TScalar, Dimension>
::ComputeLikelihoodGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> target,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL)
{
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLikelihood");
	
	this->UpdateDeformationAndKernelDataDomain(target);
	
	int numSubjects = Momentas.size();

	this->ComputeDataTermGradient(Momentas, target, GradPos, GradMom, GradTempL, true);

	for (unsigned int s = 0; s < numSubjects; s++)
	{
		// Be careful: it should have been updated Superclass::m_CovMomInverse
		// STANLEY
		// VectorType aux = 2.0 * ( Superclass::m_CovMomInverse * this->Vectorize( Momentas[s] ) );
		// PIETRO
		VectorType aux = Superclass::m_CovMomInverse * this->Vectorize( Momentas[s] ) ;
// Before LinAlg:		GradMom[s] += this->VectorToMatrix( aux );
		GradMom[s] += aux.convert_to_matrix_row_wise(aux.size()/Dimension, Dimension);
	}
}


#endif /* _BayesianAtlas_txx */

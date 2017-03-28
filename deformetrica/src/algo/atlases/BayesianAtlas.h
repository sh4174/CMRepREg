/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BayesianAtlas_h
#define _BayesianAtlas_h

#include "Atlas.h"
#include "BayesianAtlasFunctionalValues.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      Atlas object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    A child class of the class Atlas for Bayesian atlas estimation.
 */
template <class TScalar, unsigned int Dimension>
class BayesianAtlas : public Atlas<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Atlas Type
	typedef Atlas<TScalar, Dimension> Superclass;
	
	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;
	/// List of matrices type.
	typedef typename Superclass::MatrixList MatrixList;

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	/// List of deformable objects type.
	typedef std::vector<DeformableObjectType*> DeformableObjectList;

	/// Multi-Deformable object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	
	/// Values of the different terms in the log-likelihood
	typedef BayesianAtlasFunctionalValues<TScalar> FunctionalValuesType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	BayesianAtlas();
	/// Copy constructor.
	BayesianAtlas(const BayesianAtlas& other);

	/// Makes a copy of the object.
	BayesianAtlas* Clone() { return new BayesianAtlas(*this); }

	virtual ~BayesianAtlas();
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets HyperParameters of the covariance matrices to \e d.
	inline void SetCovarianceMomenta_HyperParameter(TScalar d) { m_CovMom_HyperParameter = d; }
	
	/// TODO .
	inline void SetCovarianceMomenta_Prior_Inverse(MatrixType M) { m_CovMom_Prior_Inverse = M; }
	
	/// TODO
	inline void SetCovarianceMomenta_Prior_Inverse(std::string fn);
	
	/// TODO .
	inline void SetDataSigmaSquared_HyperParameter(VectorType d) { m_DataSigmaSquared_HyperParameter = d; }
	
	/// TODO .
	inline void SetDataSigmaSquared_Prior(VectorType V) { m_DataSigmaSquared_Prior = V; }

	/// TODO .
	inline void SetNoiseDimension(std::vector<int> V) { m_NoiseDimension = V; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	virtual void Update();
	
	/// TODO .
	virtual void WriteAtlasParameters(const std::string& AtlasName);
	
    /**
     *  \brief  Compute the optimal atlas parameters (data sigma squared and covariance matrix of momentas), which maximizes the likelihood
     *
     *  \details    For a given set of \e Momentas and \e observations, compute the momenta covariance matrix and the data sigma squared that minimize the likelihood and return the values of the different terms of the log-likelihood
     *
     *  \param[in]    Momentas    The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
     *  \param[in]    observations The target objects\f$S\f$.
     */
	FunctionalValuesType* ComputeLikelihoodWithOptimalParameters(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > observations);
    
    /// Compute the log-likelihood for deformation parameters \e Momentas and target objects /e observations except the prior terms
    FunctionalValuesType* ComputeLikelihoodWithoutPriors(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > observations);

    
    /**
     *  \brief  Compute the terms of the log-likelihood, which involves the momenta covariance matrix exept the priors
     *
     *  \details Compute the terms of the likelihood, which involves the momenta covariance matrix except the priors, using the current value of Superclass::m_CovMomInverse
     *
     *  \param[in]    Momentas    The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
    *  \param[out]   CovMomTerms  Part of the log-likelihood with the variance terms: \f$ 0.5*\alpha_i^T OptimalCovMomInverse \alpha_i + 0.5*log(det(OptimalCovMomInverse))\f$ for each subject index \f$i\f$
     * Return \f$\sum_{i=1}^{NumberOfSubjects} CovMomTerms(i) \f$
     */
    TScalar ComputeCovMomTerms(const MatrixList& Momentas, VectorType& CovMomTerms);
    
    

    /**
     *  \brief  Compute the terms of the log-likelihood, which involves the variance of noise except the priors
     *
     *  \details Compute the terms of the log-likelihood, which involves the variance of noise except the priors using the current value of Superclass::m_DataSigmaSquared
     *
     *  \param[in]    Residuals    The values of the residuals \f$\left\Vert\phi^{\alpha_i}(X_{0,k}) - S_k\right\Vert_{W}^2 \f$ for each subject index \f$i/f$ and each object index \f$k\f$
     *  \param[out]   DataTerms  Part of the log-likelihood with the variance terms: \f$ 0.5*Residuals_{i,k}/OptimalDataSigmaSquared(k) + 0.5*dim_k log(OptimalDataSigmaSquared(k))\f$ for each object index \f$k\f$ and each subject index \f$i\f$
     * Return \sum_{k=1}^{NumberOfObjects} \left( \sum_{i=1}^{NumberOfSubjects} DataTerms(i,k)\right)
     */
    TScalar ComputeDataTerms(const std::vector< std::vector< TScalar> >& Residuals, std::vector< std::vector<TScalar > >& DataTerms);

    
    
    /**
     *  \brief  Compute the optimal covariance matrix of momentas, which maximizes the likelihood
     *
     *  \details Compute the optimal covariance matrix and the \e CovMomTerms and \e CovMomPrior and returns the sum of the terms in the log-likelihood corresponding the momenta latent variables
     *
     *  \param[in]    Momentas    The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
     *  \param[out]   OptimalCovMomInverse  The optimal value of the momenta covariance matrix
     *  \param[out]   CovMomTerms  Part of the log-likelihood with the variance terms: \f$ 0.5*\alpha_i^T OptimalCovMomInverse \alpha_i + 0.5*log(det(OptimalCovMomInverse))\f$ for each subject index \f$i\f$
     *  \param[out]   CovMomPrior  Part of the log-likelihood with the priors of the variance terms: \f$ 0.5*m_CovMom_Hyperparameter* (Trace(m_CovMom_Prior_Inverse^T OptimalCovMomInverse) + log(det(OptimalCovMomInverse)) )\f$
     * Return \f$CovMomPriors + \sum_{i=1}^{NumberOfSubjects} CovMomTerms(i) \f$
     */
    TScalar ComputeOptimalCovMomInverse(const MatrixList& Momentas, MatrixType& OptimalCovMomInverse, VectorType& CovMomTerms, TScalar& CovMomPrior);
    
    /**
     *  \brief  Compute the optimal data sigma squared for each object, which maximizes the likelihood
     *
     *  \details Compute the optimal data sigma squared and the \e DataTerms and \e DataPrior and returns the sum of the terms in the log-likelihood corresponding the object noise latent variables
     *
     *  \param[in]    Residuals    The values of the residuals \f$\left\Vert\phi^{\alpha_i}(X_{0,k}) - S_k\right\Vert_{W}^2 \f$ for each subject index \f$i/f$ and each object index \f$k\f$
     *  \param[out]   OptimalDataSigmaSquared  The optimal value of the data sigma squared for each object
     *  \param[out]   DataTerms  Part of the log-likelihood with the variance terms: \f$ 0.5*Residuals_{i,k}/OptimalDataSigmaSquared(k) + 0.5*dim_k log(OptimalDataSigmaSquared(k))\f$ for each object index \f$k\f$ and each subject index \f$i\f$
     *  \param[out]   DataPrior  Part of the log-likelihood with the priors of the variance terms: \f$ 0.5*m_DataSigmaSquared_HyperParameter* (m_DataSigmaSquared_Prior(k) / OptimalDataSigmaSquared(k)) + log(OptimalDataSigmaSquared(k)) )\f$ for each object index \f$k\f$
     * Return \sum_{k=1}^{NumberOfObjects} \left(DataPrior(k) + \sum_{i=1}^{NumberOfSubjects} DataTerms(i,k)\right)
     */
    TScalar ComputeOptimalDataSigmaSquared(const std::vector< std::vector< TScalar> >& Residuals, VectorType& OptimalDataSigmaSquared, std::vector< std::vector<TScalar > >& DataTerms, VectorType& DataPriors);

	
	/// For given values of the \e Momentas and \e observations, compute the gradient of the log-likelihood with respect to the momenta in \e GradMom, the positions of the control points in \e GradPos and template objects in \e GradTempL
	void ComputeLikelihoodGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> observations,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL);
	
	// void Write(std::string AtlasName, std::vector<std::string>& TemplateObjectsNames);

protected:

	/// TODO .
	void InitializeCovMomPrior();	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Protected attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// TODO .
	TScalar m_CovMom_HyperParameter;
	
	/// TODO .
	MatrixType m_CovMom_Prior_Inverse;
	
	/// TODO .
	// STANLEY
	//TScalar m_DataSigmaSquared_HyperParameter;
	// PIETRO
	VectorType m_DataSigmaSquared_HyperParameter;
	
	/// TODO .
	VectorType m_DataSigmaSquared_Prior;
	
	/// TODO .
	std::vector<int> m_NoiseDimension;


}; /* class BayesianAtlas */


#ifndef MU_MANUAL_INSTANTIATION
#include "BayesianAtlas.txx"
#endif


#endif /* _BayesianAtlas_h */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BayesianAtlasFunctionalValues_h
#define _BayesianAtlasFunctionalValues_h

#include "LinearAlgebra.h"
#include <vector>

/**
 *  \brief 		Values of the functional during bayesian atlas estimation.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details   The BayesianAtlasFunctionalValues class stores the values of the different
 *				parts of the functional at a given iteration. The functional is :
 // *				\f[
 // *				\functional =
 // *				\sum_{i=1}^{N} \left\lbrace \sum_{k=1}^{\nbobj}\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2 +
 // *				\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i +
 // *				\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i} \right\rbrace
 // *				\f]
 // *				where :
 // *				\li \f$\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2\f$ denotes the fidelity-to-data term of the k-th object of the i-th subject (BayesianAtlasFunctionalValues::m_DataTerm) ;
 // *				\li \f$\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i \f$ denotes the regularity term of the i-th deformation, which deforms the template objects to the objects of the i-th subject (BayesianAtlasFunctionalValues::m_Regularity) ;
 // *				\li \f$\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i}\f$ denotes the \f$\Lone\f$ prior (BayesianAtlasFunctionalValues::m_Sparsity) . \n
 // *				This class goes with the SparseDiffeoAtlasEstimator class.
 */
template <class TScalar>
class BayesianAtlasFunctionalValues
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type (std).
	typedef std::vector<TScalar> SDTVectorType;
	
	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;

	

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Constructor with initialization of the number of deformable objects
	/// (resp. the number of subjects) to \e nobj (resp. \e to nsubj).
	BayesianAtlasFunctionalValues()
	{
		m_OutOfBox = false;
        m_LikelihoodWithPriors = true;
	}

	/// Copy constructor.
	BayesianAtlasFunctionalValues(const BayesianAtlasFunctionalValues& other)
	{
		
		m_OutOfBox = other.m_OutOfBox;
        m_LikelihoodWithPriors = other.m_LikelihoodWithPriors;
		m_DataSigmaSquared = other.m_DataSigmaSquared;
		m_CovMomInverse = other.m_CovMomInverse;
		m_Residuals = other.m_Residuals;
        m_CovMomTerms = other.m_CovMomTerms;
        m_CovMomPrior = other.m_CovMomPrior;
        m_DataTerms = other.m_DataTerms;
        m_DataPriors = other.m_DataPriors;
        
		
	}

	/// Makes a copy of the object.
	BayesianAtlasFunctionalValues* Clone() { return new BayesianAtlasFunctionalValues(*this); }

	~BayesianAtlasFunctionalValues() {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Get the values of the data sigma squared of each object used in the computation of the different terms of the log-likelihood
	inline VectorType GetDataSigmaSquared() const { return m_DataSigmaSquared; }
	/// Set the values of the data sigma squared of each object used in the computation of the different terms of the log-likelihood
	inline void SetDataSigmaSquared(const VectorType t) { m_DataSigmaSquared = t; }
	
	/// Get the momenta covariance matrix used in the computation of the different terms of the log-likelihood
	inline MatrixType GetCovMomInverse() const { return m_CovMomInverse; }
	/// Set the momenta covariance matrix used in the computation of the different terms of the log-likelihood
	inline void SetCovMomInverse(const MatrixType t) { m_CovMomInverse = t; }
    
    /// Set the values of the residuals \f$\left\Vert\phi^{\alpha_i}(X_{0,k}) - S_k\right\Vert_{W}^2 \f$ for each subject index \f$i/f$ and each object index \f$k\f$
    inline void SetResiduals(const std::vector<SDTVectorType> t) { m_Residuals = t; }
    /// Get the values of the residuals \f$\left\Vert\phi^{\alpha_i}(X_{0,k}) - S_k\right\Vert_{W}^2 \f$ for each subject index \f$i/f$ and each object index \f$k\f$
    inline std::vector< SDTVectorType > GetResiduals() const { return m_Residuals; }

    /// Set the value of the part of the log-likelihood with the variance terms: \f$ 0.5*Residuals_{i,k}/DataSigmaSquared(k) + 0.5*dim_k log(DataSigmaSquared(k))\f$ for each object index \f$k\f$ and each subject index \f$i\f$
    inline void SetDataTerms(const std::vector<SDTVectorType> t) { m_DataTerms = t; }
    
    /// Set the value of the part of the log-likelihood with the priors of the variance terms: \f$ 0.5*DataSigmaSquared_HyperParameter* (DataSigmaSquared_Prior(k) / DataSigmaSquared(k)) + log(DataSigmaSquared(k)) )\f$ for each object index \f$k\f$
    inline void SetDataPriors(const VectorType t) { m_DataPriors = t; }
    
    /// Set the value of the part of the log-likelihood with the variance terms: \f$ 0.5*\alpha_i^T CovMomInverse \alpha_i + 0.5*log(det(CovMomInverse))\f$ for each subject index \f$i\f$
    inline void SetCovMomTerms(const VectorType t) { m_CovMomTerms = t; }
    
    /// Set the value of the part of the log-likelihood with the priors of the variance terms: \f$ 0.5*CovMom_Hyperparameter* (Trace(CovMom_Prior_Inverse^T CovMomInverse) + log(det(CovMomInverse)) )\f$
    inline void SetCovMomPrior(const TScalar t) { m_CovMomPrior = t; }
    
    /// Return true if any point of the trajectory is out of box, false otherwise.
    inline bool IsOutOfBox() { return m_OutOfBox; }
    /// Set the boolean OutOfBox to true.
    inline void SetOutOfBox() { m_OutOfBox = true; }
    
    // Set flag to say that the likelihood should be computed with priors
    inline void SetLikelihoodWithPriors() { m_LikelihoodWithPriors = true; }
    // Set flag to say that the likelihood should be computed omitting prior terms
    inline void UnsetLikelihoodWithPriors() { m_LikelihoodWithPriors = false; }





	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /// Return the value of the likelihood.
    TScalar GetLikelihood()
    {
        TScalar likelihood = this->GetLikelihoodWithoutPriorTerms();
        
        likelihood += ( (m_LikelihoodWithPriors)?(this->GetPriorTerms()):0.0 );

        return likelihood;
    }

    /// Compute the value of the likelihood without prior terms
    TScalar GetLikelihoodWithoutPriorTerms()
    {
        int NumberOfSubjects = m_DataTerms.size();
        
        TScalar likelihood = 0.0;
        for (unsigned int s = 0; s < NumberOfSubjects; s++)
            likelihood += this->GetLikelihoodWithoutPriorTermsForSubject(s);
        
        
        return likelihood;
    }
    
    /// Compute the contribution of the s-th subject to the likelihood (prior terms are not summed over the subjects and therefore do not count in this computation)
    TScalar GetLikelihoodWithoutPriorTermsForSubject(unsigned int s)
    {
        TScalar likelihood = m_CovMomTerms(s);
        for (unsigned int i = 0;  i < m_DataTerms[s].size(); i++)
            likelihood += m_DataTerms[s][i];
        
        return likelihood;
    }

    std::vector< TScalar > GetLikelihoodWithoutPriorTermsForAllSubjects()
    {
        int NumberOfSubjects = m_DataTerms.size();
        std::vector< TScalar > out(NumberOfSubjects);
        for (unsigned int s = 0; s < NumberOfSubjects; s++)
        {
            out[s] = GetLikelihoodWithoutPriorTermsForSubject(s);
        }
        return out;
    }
    
    
    TScalar GetPriorTerms()
    {
        if (!m_LikelihoodWithPriors)
        {
            throw std::runtime_error("could not get prior terms in BayesianAtlasFunctionalValues");
        }
        return (m_CovMomPrior + m_DataPriors.sum());
    }
   
    
    /// Sum the residuals over the subjects, for each object.
    VectorType GetResidualsPerObjects()
    {
        int NumberOfSubjects = m_Residuals.size();
        int NumberOfObjects = m_Residuals[0].size();
        VectorType ResidualsPerObjects(NumberOfObjects, 0.0);
        for (unsigned int s = 0; s < NumberOfSubjects; s++)
        {
            for (unsigned int i = 0; i < NumberOfObjects; i++)
                ResidualsPerObjects[i] += m_Residuals[s][i];
        }
        return ResidualsPerObjects;
    }

    

	/// Displays the different parts of the functional with an indication of iteration \e iter.
	void PrintIter(const int iter)
	{
		std::cout << "Iter " << iter << "  >> - log_Likelihood = " << this->GetLikelihood() << std::endl;
		std::cout << "\t DataSigmaSquared = " << m_DataSigmaSquared << std::endl;
		std::cout << "\t Residual Norms Per Objects = " << this->GetResidualsPerObjects() << std::endl;
//		std::cout << "\t log(Determinant Covariance Matrix) = " << m_Log_det_Cov_Mom << std::endl;
        std::cout << "\t CovMomTerms = " << m_CovMomTerms.sum() << std::endl;
	}



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// See AbstractDeformations::m_OutOfBox for details.
	bool m_OutOfBox;

    /// flag to know whether including prior terms in the likelihood or not
    bool m_LikelihoodWithPriors;
    
    /// Values of the data sigma squared used for computing the different terms of the likelihood
    VectorType m_DataSigmaSquared;
    
    /// Inverse of the momenta covariance matrix used for computing the different terms of the likelihood
    MatrixType m_CovMomInverse;
    
	/// Residuals per subject per object
	std::vector<SDTVectorType> m_Residuals;
    
    /// Terms of the likelihood that involve the covariance of the momentas
    VectorType  m_CovMomTerms;
    
    /// Term of the likelihood that corresponds to the prior of the momenta covariance matrix
    TScalar m_CovMomPrior;
    
    /// Terms of the likelihood that involve the variance of the object's noise
    std::vector<SDTVectorType> m_DataTerms;
    
    /// Term of the likelihood that corresponds to the prior of the noise variance (for each object)
    VectorType m_DataPriors;
  

}; /* class BayesianAtlasFunctionalValues */


#endif /* _BayesianAtlasFunctionalValues_h */

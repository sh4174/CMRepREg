/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeterministicAtlasFunctionalValues_h
#define _DeterministicAtlasFunctionalValues_h

#include "LinearAlgebra.h"
#include <vector>

/**
 *  \brief      Values of the functional during deterministic atlas estimation.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The DeterministicAtlasFunctionalValues class stores the values of the different
 *              parts of the functional at a given iteration. The functional is :
 // *				\f[
 // *				\functional =
 // *				\sum_{i=1}^{N} \left\lbrace \sum_{k=1}^{\nbobj}\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2 +
 // *				\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i +
 // *				\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i} \right\rbrace
 // *				\f]
 // *				where :
 // *				\li \f$\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2\f$ denotes the fidelity-to-data term of the k-th object of the i-th subject (DeterministicAtlasFunctionalValues::m_DataTerm) ;
 // *				\li \f$\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i \f$ denotes the regularity term of the i-th deformation, which deforms the template objects to the objects of the i-th subject (DeterministicAtlasFunctionalValues::m_Regularity) ;
 // *				\li \f$\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i}\f$ denotes the \f$\Lone\f$ prior (DeterministicAtlasFunctionalValues::m_Sparsity) . \n
 // *				This class goes with the SparseDiffeoAtlasEstimator class.
 */
template <class TScalar>
class DeterministicAtlasFunctionalValues
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type (std).
	typedef std::vector<TScalar> STDVectorType;
	
	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Constructor with initialization of the number of deformable objects
	/// (resp. the number of subjects) to \e nobj (resp. \e to nsubj).
	DeterministicAtlasFunctionalValues()
	{
		m_OutOfBox = false;
	}

	/// Copy constructor.
	DeterministicAtlasFunctionalValues(const DeterministicAtlasFunctionalValues& other)
	{
		
		m_OutOfBox = other.m_OutOfBox;
		m_Residuals = other.m_Residuals;
		m_ResidualsPerObjects = other.m_ResidualsPerObjects;
		m_DataTerm = other.m_DataTerm;
		m_RegTerm = other.m_RegTerm;
		m_SparsityTerm = other.m_SparsityTerm;
		
	}

	/// Makes a copy of the object.
	DeterministicAtlasFunctionalValues* Clone() { return new DeterministicAtlasFunctionalValues(*this); }

	~DeterministicAtlasFunctionalValues() {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the total value of the functional.
	inline TScalar GetTotalValue() { return (m_DataTerm + m_RegTerm + m_SparsityTerm); }
	/// Returns the value of the L2 part of the functional.
	inline TScalar GetTotalL2Value() { return (m_DataTerm + m_RegTerm); }
	
	/// Returns the fidelity-to-data term.
	inline TScalar GetDataTerm() const { return m_DataTerm; }
	/// Sets the fidelity-to-data term to \e t.
	inline void SetDataTerm(TScalar t) { m_DataTerm = t; }
	
	/// Sets the regularity term to \e t.
	inline void SetRegularityTerm(TScalar t) { m_RegTerm = t; }

	/// Sets the sparsity term to \e t.
	inline void SetSparsityTerm(TScalar t) { m_SparsityTerm = t; }
	
	/// TODO .
	inline void SetResiduals(const std::vector<STDVectorType> t)
	{
		 m_Residuals = t;
		 int numSubjects = t.size();
		 int numObjects = t[0].size();
		 m_ResidualsPerObjects.set_size(numObjects);
		 m_ResidualsPerObjects.fill(0.0);
		 for (unsigned int s = 0; s < numSubjects; s++)
		 {
			 for (unsigned int i = 0; i < numObjects; i++)
				 m_ResidualsPerObjects[i] += m_Residuals[s][i] / numSubjects;
		 }		 
	}
	
	/// Returns true if any point of the trajectory is out of box, false otherwise.
	inline bool IsOutOfBox() { return m_OutOfBox; }
	/// Sets the boolean OutOfBox to true.
	inline void SetOutOfBox() { m_OutOfBox = true; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Displays the different parts of the functional with an indication of iteration \e iter.
	void PrintIter(const int iter)
	{
		std::cout << "Iter " << iter << "  >> Functional = " << this->GetTotalValue() <<
					"\t (Data Term = " << m_DataTerm << 
					"\t Regularity Term  = " << m_RegTerm <<
					"\t Sparsity Prior = " << m_SparsityTerm << " )" << std::endl;
		//std::cout << "\tResidual Norms per Objects (divided by the number of subjects) = " << m_ResidualsPerObjects << std::endl;
	}



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// If an OutOfBox exception occurs during computation of the functional values.
	bool m_OutOfBox;

	/// Vector of residual norms per objects and per subjects.
	std::vector<STDVectorType> m_Residuals;
	
	/// Vector of residual norms per objects (sum over the subjects of \e m_Residuals).
	VectorType m_ResidualsPerObjects;

	/// Value of the fidelity-to-data term.
	TScalar m_DataTerm;
	
	/// Value of the regularity term.
	TScalar m_RegTerm;
	
	/// Value of the sparsity term.
	TScalar m_SparsityTerm;



}; /* class DeterministicAtlasFunctionalValues */


#endif /* _DeterministicAtlasFunctionalValues_h */

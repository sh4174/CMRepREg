/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _RegressionFunctionalValues_h
#define _RegressionFunctionalValues_h

#include "LinearAlgebra.h"
#include <vector>

/**
 *  \brief      Values of the functional during regression estimation.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The RegressionFunctionalValues class stores the values of the different
 *              parts of the functional at a given iteration. The functional is :
 *
 */
template <class TScalar>
class RegressionFunctionalValues
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
	RegressionFunctionalValues()
	{
		m_OutOfBox = false;
	}

	/// Copy constructor.
	RegressionFunctionalValues(const RegressionFunctionalValues& other)
	{
		
		m_OutOfBox = other.m_OutOfBox;
		m_Residuals = other.m_Residuals;
		m_ResidualsPerObjects = other.m_ResidualsPerObjects;
		m_DataTerm = other.m_DataTerm;
		m_RegTerm = other.m_RegTerm;
		m_SparsityTerm = other.m_SparsityTerm;
		
	}

	/// Makes a copy of the object.
	RegressionFunctionalValues* Clone() { return new RegressionFunctionalValues(*this); }

	~RegressionFunctionalValues() {}



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
		//std::cout << "\tResidual Norms per Objects (divided by the number of observations) = " << m_ResidualsPerObjects << std::endl;
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



}; /* class RegressionFunctionalValues */


#endif /* _RegressionFunctionalValues_h */

/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeterministicAtlas_h
#define _DeterministicAtlas_h

#include "Atlas.h"
#include "DeterministicAtlasFunctionalValues.h"

#include <vector>

#include <cstring>
#include <iostream>
#include <sstream>


/**
 *  \brief      Atlas object class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    A child class of the class Atlas for deterministic atlas estimation.
 */
template <class TScalar, unsigned int Dimension>
class DeterministicAtlas : public Atlas<TScalar, Dimension>
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
	typedef DeterministicAtlasFunctionalValues<TScalar> FunctionalValuesType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeterministicAtlas();
	/// Copy constructor.
	DeterministicAtlas(const DeterministicAtlas& other);

	/// Makes a copy of the object.
	DeterministicAtlas* Clone() { return new DeterministicAtlas(*this); }

	virtual ~DeterministicAtlas();
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// use the norm of the RKHS for computing the regularity term. Otherwise, use the pre-defined matrix \e m_CovMom_Inverse
	inline void SetRKHSNormForRegularization() { m_UseRKHSNormForRegularization = true; }
	/// use the  pre-defined matrix \e m_CovMom_Inverse for computing the regularity term
	inline void UnSetRKHSNormForRegularization() { m_UseRKHSNormForRegularization = false; }

	/// Returns the sparsity prior.
	inline TScalar GetSparsityPrior() const { return m_SparsityPrior; }
	/// Sets the sparsity prior to \e d.
	inline void SetSparsityPrior(const TScalar d) { m_SparsityPrior = d; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void Update();
	
	/// TODO .
	FunctionalValuesType* ComputeFunctional(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > observations);
	
	/// TODO .
	/// \warning    Computes the gradient of the L2 part of the functional only!
	void ComputeFunctionalGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> observations,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL_L2, MatrixList& GradTempL_Sob);

	void ComputeFunctionalGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> observations,
			MatrixList& GradMom, MatrixType& GradPos, MatrixList& GradTempL, bool do_SobolevGradient);
			
	void AddGradientRegularityTerm(MatrixType GradPos, MatrixList GradMom, const MatrixList& Momentas);
	

protected:

	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Protected attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// flag to select the method to compute the regularity term
	bool m_UseRKHSNormForRegularization;

	/// Coefficient added to the \f$L^1\f$ penalty.
	TScalar m_SparsityPrior;
	

}; /* class DeterministicAtlas */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeterministicAtlas.txx"
#endif


#endif /* _DeterministicAtlas_h */

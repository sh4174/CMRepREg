/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BayesianAtlasEstimator_h
#define _BayesianAtlasEstimator_h

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "LinearAlgebra.h"
#include <vector>

#include "readMatrixDLM.txx"

#include "Diffeos.h"


#include "DeformableMultiObject.h"
#include "BayesianAtlas.h"

#include "BayesianAtlasFunctionalValues.h"

/**
 *  \brief      Bayesian Atlas construction class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The BayesianAtlasEstimator class enables to build atlases from collections of objects configurations.\n \n
 *              Given a population of M objects and N subjects, it estimates a common template complex T made of M template objects,
 *              the noise variance of each object and the covariance matrix of the momentum vectors.
 */
template <class TScalar, unsigned int Dimension>
class BayesianAtlasEstimator : public AbstractAtlasEstimator< TScalar, Dimension >
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// AbstractAtlasEstimator
	typedef AbstractAtlasEstimator< TScalar, Dimension > Superclass;

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	/// List of matrices type.
	typedef std::vector<MatrixType> MatrixList;

	/// Deformable multi-object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	/// List of deformable multi-objects type.
	typedef std::vector<DeformableMultiObjectType*> DeformableMultiObjectList;
	/// Abstract Atlas type.
	typedef Atlas<TScalar, Dimension> AtlasType;
	/// Bayesian Atlas type.
	typedef BayesianAtlas<TScalar, Dimension> BayesianAtlasType;
	
	/// Functional values type.
	typedef BayesianAtlasFunctionalValues<TScalar> FunctionalValuesType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	BayesianAtlasEstimator();

	virtual ~BayesianAtlasEstimator();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the atlas.
	virtual AtlasType* GetAtlas() const { return m_Atlas; }
	/// Sets atlas to \e obj.
	virtual void SetAtlas(AtlasType* obj)
	{
		if (! obj->IsBayesian())
			throw std::runtime_error("Atlas given in Bayesian Atlas Estimator is not Bayesian");
		
		m_Atlas = dynamic_cast<BayesianAtlasType*>(obj);
		m_NumberOfCPs = m_Atlas->GetControlPoints().rows();
	}

	/// Returns the list of the momentas.
	inline MatrixList GetMomenta() const { return m_InitialMomentas; }
	/// Initializes the momentas by reading the file \e fn.
	inline void SetInitialMomentas(std::string fn)
	{
		if (strlen(fn.c_str())) // null file name means no file to be read
			m_InitialMomentas = readMultipleMatrixDLM<TScalar>(fn.c_str());
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Performs the optimization method .
	virtual void Update();
	
	virtual void WriteOutput(std::string AtlasName);
	
	virtual void WriteAtlasToSubjectDeformations();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 *	\brief		Performs line search gradient descent.
	 *
	 *	\details	This method implements gradient descent with a line search strategy.
	 *  		\li	gradient descent step with current stepsize is accepted if energy
	 *  			is decreased (\f$f(x^{n+1}) - f(x^n) < 0\f$) ;
	 *  \param[in, out]	X	Control points to be optimized or not (depending on
	 *						BayesianAtlasEstimator::m_freezeCP).
	 *  \param[in, out]	A	Momentas to be optimized.
	 *  \param[in, out]	T	Template to be optimized.
	 *  \param[in, out]	DSS		vector of noise variance of each object to be optimized
	 *  \param[in, out]	CovMom	covariance matrix of the momentum vectors to be optimized
	
	 */
	void GradientDescent(MatrixType& X, MatrixList& A, MatrixList& T, VectorType& DSS, MatrixType& CovMom);


	/**
	 *  \brief       Computes a gradient descent step for test values of the variables and step sizes.
	 *
	 *  \details     If there is no sparsity prior and if control points can move, this method computes :
	 *               \f[ x^{n+1} \leftarrow x^n - \tau \left(\nabla_{x}\functional\right) ,\f]
	 *               \f[ \alpha^{n+1} \leftarrow \alpha^n - \tau \left(\nabla_{\alpha}\functional\right), \f]
	 *               where \f$\tau\f$ is the current step-size of the gradient descent. \n
	 *               If the sparsity prior is positive, then the method SoftThresholdUpdate() is called
	 *               to update the momentas using a soft-thresholding function.
	 *               If control points can't move, then \f$x^{n+1} \leftarrow x^n\f$.
	 *
	 *  \param[out]   Xtest    Updated value of the control points \f$\left(x_i^{n+1}\right)_i\f$.
	 *  \param[in]    X        Current value of the control points \f$\left(x_i^n\right)_i\f$.
	 *  \param[out]   Atest    Updated value of the momentas \f$\left(\alpha_i^{n+1}\right)_i\f$.
	 *  \param[in]    A        Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 *  \param[out]   Ttest    Updated value of the template \f$\left(x_0\right)^{n+1}\f$.
	 *  \param[in]    T        Current value of the template \f$\left(x_0\right)^n\f$.
	 *  \param[in]    stepXA   Step-size for momentas and control point positions update.
	 *  \param[in]    stepT    Step-size for template update.
	 */
	void GradientDescentStep(
			MatrixType& Xtest, const MatrixType& X,
			MatrixList& Atest, const MatrixList& A,
			MatrixList& Ttest, const MatrixList& T,
			TScalar stepXA, TScalar stepT);


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// The atlas.
	BayesianAtlasType* m_Atlas;

	/// Gradient of the functional w.r.t. the CP position.
	MatrixType m_GradPos;

	/// Gradient of the functional w.r.t. momentas.
	MatrixList m_GradMom;
	
	/// Gradient w.r.t. template objects.
	MatrixList m_GradTempL;

	/// List of matrices of the coordinates of the momenta (List size: NumberOfSubjects, Matrices size: NumberOfCPs x Dimension).
	MatrixList m_InitialMomentas;
    
    /// Number of control points in the atlas
    int m_NumberOfCPs;

	/// Vector containing a history of the the values of the functional during optimization method.
	std::vector< FunctionalValuesType* > m_ValuesHistory;


protected:


}; /* class BayesianAtlasEstimator */


#ifndef MU_MANUAL_INSTANTIATION
#include "BayesianAtlasEstimator.txx"
#endif


#endif /* _BayesianAtlasEstimator_h */

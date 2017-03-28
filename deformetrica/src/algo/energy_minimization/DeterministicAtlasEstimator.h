/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeterministicAtlasEstimator_h
#define _DeterministicAtlasEstimator_h

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "LinearAlgebra.h"
#include <vector>

#include "readMatrixDLM.txx"

#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "DeterministicAtlas.h"

#include "DeterministicAtlasFunctionalValues.h"


/**
 *  \brief      Deterministic Atlas construction class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The DeterministicAtlasEstimator class enables to build atlases from collections of objects configurations.\n \n
 *              Given a population of M objects and N subjects, it estimates a common template complex T made of M template objects,
 *              the noise variance of each object and the covariance matrix of the momentum vectors.
 */
template <class TScalar, unsigned int Dimension>
class DeterministicAtlasEstimator : public AbstractAtlasEstimator< TScalar, Dimension >
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
	/// Abstract Atlas type
	typedef Atlas<TScalar, Dimension> AtlasType;
	/// Atlas type type.
	typedef DeterministicAtlas<TScalar, Dimension> DeterministicAtlasType;

	/// Functional values type.
	typedef DeterministicAtlasFunctionalValues<TScalar> FunctionalValuesType;

	/// Type of optimization method.
	typedef enum
	{
		null,			/*!< Null value. */
		GradientDescent,/*!< Usual line search gradient descent (see GradientDescentAndISTA()). */
		ISTA,			/*!< Line search gradient descent compatible with \f$L^1\f$ penalty terms (see GradientDescentAndISTA()). */
		F_ISTA,			/*!< Faster gradient scheme built on ISTA and Nesterov's scheme (see FISTA()). */
	} OptimizationMethodType;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeterministicAtlasEstimator();

	virtual ~DeterministicAtlasEstimator();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the atlas.
	virtual AtlasType* GetAtlas() const { return m_Atlas; }
	/// Sets the atlas to \e obj.
	virtual void SetAtlas(AtlasType* obj)
	{
		if (! obj->IsDeterministic())
			throw std::runtime_error("Atlas given in Deterministic Atlas Estimator is not deterministic");
		
		m_Atlas = dynamic_cast<DeterministicAtlasType*>(obj);
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

	/// Selects the Gradient Descent optimization method.
	inline void SetGradientDescent() { m_OptimizationMethod = GradientDescent; }
	/// Selects the ISTA algorithm.
	inline void SetISTA() { m_OptimizationMethod = ISTA; }
	/// Selects the FISTA algorithm.
	inline void SetFISTA() { m_OptimizationMethod = F_ISTA; }
	

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
	 *	\brief		Perform line search gradient descent.
	 *
	 *	\details	This method implements gradient descent with a line search strategy.
	 *  		\li	gradient descent step with current stepsize is accepted if energy
	 *  			is decreased (\f$f(x^{n+1}) - f(x^n) < 0\f$) ;
	 *	\param[in, out]	X	Control points to be optimized or not (depending on
	 *						DeterministicAtlasEstimator::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 *	\param[in, out]	T	Template to be optimized.	
	 */
	void GradientDescentAndISTA(MatrixType& X, MatrixList& A, MatrixList& T);
	
	/**
	 *	\brief		Perform FISTA.
	 *
	 *	\details	This optimization method combines the ISTA method (condition for accepting
	 *				gradient step) with an improved gradient descent scheme introduced by Nesterov :\n
	 *  		\li	If m_sparsityPrior = 0.0, the Nesterov scheme is applied as explained in :\n
	 *				A method of solving a convex programming problem with convergence rate O (1/k2),
	 *				Y. Nesterov, Soviet Mathematics Doklady 27 (2), 372-376 (1983) ;
	 *  		\li	If m_SparsityPrior > 0, the method implements the FISTA method as explained in :\n
	 *  			A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems,
	 *  			A. Beck & M. Teboulle, SIAM J. Imaging Sciences 2(1), 183--202 (2009) .
	 *
	 *	\param[in, out]	X	Control points positions (optimized or not, depending on
	 *						SparseDiffeoAtlasEstimator::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 *	\param[in, out]	T	Template to be optimized.
	 */
	void FISTA(MatrixType& X, MatrixList& A, MatrixList& T);


	/**
	 * \brief		Computes a gradient descent step for test values of the variables and step sizes
	 *
	 * \details		If there is no sparsity prior and if control points can move, this method computes :
	 *				\f[ x^{n+1} \leftarrow x^n - \tau \left(\nabla_{x}\functional\right) ,\f]
	 *				\f[ \alpha^{n+1} \leftarrow \alpha^n - \tau \left(\nabla_{\alpha}\functional\right), \f]
	 *				where \f$\tau\f$ is the current step-size of the gradient descent. \n
	 *				If the sparsity prior is positive, then the method SoftThresholdUpdate() is called
	 *				to update the momentas using a soft-thresholding function. If control points can't move, then \f$x^{n+1} \leftarrow x^n\f$.
	 *
	 * \param[out]	Xtest			Updated value of the control points \f$\left(x_i^{n+1}\right)_i\f$.
	 * \param[in]	X				Current value of the control points \f$\left(x_i^n\right)_i\f$.
	 * \param[out]	Atest			Updated value of the momentas \f$\left(\alpha_i^{n+1}\right)_i\f$.
	 * \param[in]	A				Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[out]	Ttest			Updated value of the template \f$\left(x_0\right)^{n+1}\f$.
	 * \param[in]	T				Current value of the template \f$\left(x_0\right)^n\f$.
	 * \param[in]	activeCPMask	Mask of active control points (those carrying non-zero momenta).
	 * \param[in]	stepXA			Step-size for momentas and control point positions update.
	 * \param[in]	stepT			Step-size for template update.
	 */
	void GradientDescentStep(
			MatrixType& Xtest, const MatrixType& X,
			MatrixList& Atest, const MatrixList& A,
			MatrixList& Ttest, const MatrixList& T,
			std::vector< int >& activeCPMask,
			TScalar stepXA, TScalar stepT);
			
			
	/// Masks matrix rows. Rows whose index are masked (i.e. mask[i] = 1) are removed from matrix \e M.
	MatrixType _maskMatrix(const MatrixType& M, const std::vector<int>& mask);

	/**
	 * \brief		Soft-Threshold momentas updates.
	 *
	 * \details		This method soft-threshold the gradient descent step of the momentas. Has no effect if m_sparsityPrior = 0.
	 * 				The momentas are updated according to :
	 * 				\f[
	 * 				\alpha^{n+1}_{0,p} \leftarrow S_{\tau\sparsityprior}\left(\norm{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}\right)
	 * 				\frac{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}{\norm{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}},
	 * 				\f]
	 * 				where \f$\tau\f$ is the current step-size of the gradient descent and \f$S\f$ is
	 * 				the thresholding function \f$S_{\lambda}(x)=\max(0,x-\lambda)+\min(0,x+\lambda)\f$.
	 *
	 * \param[in]	X		Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[in]	gradX	Gradient of the functional at current iteration.
	 * \param[in]	step	Descent step size for the momentas.
	 * \param[out]	mask	Mask of active control points.
	 * \return		The soft-thresholded momentas.
	 */
	MatrixType SoftThresholdUpdate(const MatrixType& X, MatrixType& gradX, TScalar step, std::vector<int>& mask);

	/**
	 * \brief	Computes the condition for accepting gradient descent step for the F/ISTA methods.
	 *
	 * \details	F/ISTA methods require the current stepsize to satisfy
	 * 			\f$f(x^{n+1}) - f(x^n) < \verb#_QdiffTerm#(x^{n+1},x^{n})\f$ :
	 * 		\li	If m_sparsityPrior = 0.0, \f$\verb#_QdiffTerm# = -0.5 * stepsize * ||\nabla_{x^n} f||^2\f$ ;
	 * 		\li	If m_sparsityPrior > 0.0, _QdiffTerm() takes into account the soft-thresholded values of \f$x^n\f$.
	 *
	 * See details in \n
	 * A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems,
	 * A. Beck & M. Teboulle, SIAM J. Imaging Sciences 2(1), 183--202 (2009)
	 *
	 * \param[out]	Xtest	Updated value of the control points \f$\left(x_i^{n+1}\right)_i\f$.
	 * \param[in]	X		Current value of the control points \f$\left(x_i^n\right)_i\f$.
	 * \param[out]	Atest	Updated value of the momentas \f$\left(\alpha_i^{n+1}\right)_i\f$.
	 * \param[in]	A		Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[out]	Ttest	Updated value of the template \f$\left(x_0\right)^{n+1}\f$.
	 * \param[in]	T		Current value of the template \f$\left(x_0\right)^n\f$.
	 * \param[in]	stepXA	Current stepsize for the control points and the momentas.
	 * \param[in]	stepT	Current stepsize for the template.
	 * \return		The value of QdiffTerm.
	 */
	TScalar _QdiffTerm(
			const MatrixType& Xtest, const MatrixType& X,
			const MatrixList& Atest, const MatrixList& A,
			const MatrixList& Ttest, const MatrixList& T,
			TScalar stepXA, TScalar stepT);


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// The atlas
	DeterministicAtlasType* m_Atlas;

	/// Gradient of the functional w.r.t. the CP position.
	MatrixType m_GradPos;

	/// Gradient of the functional w.r.t. momentas
	MatrixList m_GradMom;
	
	/// Gradient w.r.t. template objects using the L2 metric
	MatrixList m_GradTempL_L2;
	
	/// Gradient w.r.t. template objects using the Sobolev metric (convolution of the L2 metric with a smoothing kernel)
	MatrixList m_GradTempL_Sob;

	/// List of matrices of the coordinates of the momenta (List size: NumberOfSubjects, Matrices size: NumberOfCPs x Dimension).
	MatrixList m_InitialMomentas;
	
    /// Number of control points in the atlas
    int m_NumberOfCPs;
    
	/// Type of optimization method
	OptimizationMethodType m_OptimizationMethod;
	
	/// FISTA parameter. Can't increase stepsize during NbIterFreeze iterations.
	static const unsigned int NbIterFreeze = 3;

	/// Vector containing a history of the the values of the functional during optimization method.
	std::vector< FunctionalValuesType* > m_ValuesHistory;
	

protected:


}; /* class DeterministicAtlasEstimator */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeterministicAtlasEstimator.txx"
#endif


#endif /* _DeterministicAtlasEstimator_h */

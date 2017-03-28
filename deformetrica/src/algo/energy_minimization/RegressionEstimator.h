/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _RegressionEstimator_h
#define _RegressionEstimator_h

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "LinearAlgebra.h"
#include <vector>

#include "readMatrixDLM.txx"

#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "Regression.h"

#include "RegressionFunctionalValues.h"

/**
 *  \brief      Regression estimation
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The Regression class enables to compute the optimal parameters of for a given set
 *              of time indexed shape data (i.e. the observations).
 */
template <class TScalar, unsigned int Dimension>
class RegressionEstimator
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////


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
	/// Regression type type.
	typedef Regression<TScalar, Dimension> RegressionType;

	/// Functional values type.
	typedef RegressionFunctionalValues<TScalar> FunctionalValuesType;

	/// Type of optimization method.
	typedef enum
	{
		null,			/*!< Null value. */
		GradientDescent,	/*!< Usual line search gradient descent (see GradientDescentAndISTA()). */
		ISTA,			/*!< Line search gradient descent compatible with \f$L^1\f$ penalty terms (see GradientDescentAndISTA()). */
		F_ISTA,			/*!< Faster gradient scheme built on ISTA and Nesterov's scheme (see FISTA()). */
	} OptimizationMethodType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	RegressionEstimator();

	~RegressionEstimator();

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline std::vector< unsigned int > GetTimeIndices() const { return m_TimeIndices; }
	inline void SetTimeIndices(const std::vector<unsigned int> timeIndices) { m_TimeIndices = timeIndices; }

	/// Sets the list of targets to \e obj.
	inline void SetTargetList(DeformableMultiObjectList& obj)
	{
		m_TargetList = obj;
		m_NumberOfTimepoints = obj.size();
		m_NumberOfObjects = obj[0]->GetNumberOfObjects();
	}

	/// Sets the regression to \e obj.
	void SetRegression(RegressionType* obj)
	{
	 	m_Regression = obj;
	 	m_NumberOfCPs = m_Regression->GetControlPoints().rows(); // this number of CP is not supposed to change during optimization
	}
	/// Returns the output regression.
	RegressionType* GetRegression() { return m_Regression; }

	/// Returns the list of the momentas.
	inline MatrixType GetMomenta() const { return m_InitialMomenta; }
	/// Initializes the momentas by reading the file \e fn.
	inline void SetInitialMomentas(std::string fn)
	{
		if (strlen(fn.c_str())) // null file name means no file to be read
			m_InitialMomenta = readMatrixDLM<TScalar>(fn.c_str());
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
	
	/// Disables optimization in control point positions.
	inline void SetFreezeCP() { m_UpdateCP = false; }
	/// Enables optimization in control point positions.
	inline void UnsetFreezeCP() { m_UpdateCP = true; }

	/// Disables optimization of the baseline shape
	inline void SetFreezeTemplate() { m_UpdateTemplate = false; }
	/// Enables optimization of the baseline shape
	inline void UnsetFreezeTemplate() { m_UpdateTemplate = true; }

	inline void SetSaveEveryNIters(unsigned int n) { m_SaveEveryNIters = n; }

	/// Sets the maximum of descent iterations to \e n.
	inline void SetMaxIterations(unsigned int n) { m_MaxIterations = n; }

	/// Sets the maximum number of iterations to \e n.
	inline void SetMaxLineSearchIterations(unsigned int n) { m_MaxLineSearchIterations = n; }

	/// Sets the expand parameter to \e d.
	inline void SetAdaptiveExpand(TScalar d) { m_AdaptiveExpand = d; }

	/// Sets the shrink parameter to \e d.
	inline void SetAdaptiveShrink(TScalar d) { m_AdaptiveShrink = d; }

	/// Sets the tolerance parameter to \e d.
	inline void SetAdaptiveTolerance(TScalar d) { m_AdaptiveTolerance = d; }

	/// Sets the initial step to \e d.
	inline void SetInitialStepMultiplier(TScalar d) { m_InitialStepMultiplier = d; }
	
	/// Sets atlas name for saving.
	inline void SetRegressionName(std::string RegressionName) { m_RegressionName = RegressionName; }
	
	/// Sets template objects name for saving.
	inline void SetTemplateObjectsName(std::vector< std::string >& ObjectsName, std::vector< std::string >& ObjectsNameExt)
	{
		 m_TemplateObjectsName = ObjectsName;
		 m_TemplateObjectsNameExtension = ObjectsNameExt; 
	}
	

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Performs the optimization method.
	void Update();
	
	/// Writes output 
	void WriteOutput(std::string regressionName);
	void WriteOutput();
	

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
	void GradientDescentAndISTA(MatrixType& X, MatrixType& A, MatrixList& T);
	
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
	void FISTA(MatrixType& X, MatrixType& A, MatrixList& T);


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
	void GradientDescentStep(MatrixType& Xtest, const MatrixType& X,
				 MatrixType& Atest, const MatrixType& A,
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
	TScalar _QdiffTerm(const MatrixType& Xtest, const MatrixType& X,
			   const MatrixType& Atest, const MatrixType& A,
			   const MatrixList& Ttest, const MatrixList& T,
			   TScalar stepXA, TScalar stepT);

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Pointer on the regression (baseline + control points)
	RegressionType* m_Regression;

	/// \f$\nbsubj\f$ collections of deformable target objects.
	DeformableMultiObjectList m_TargetList;

	/// Time indices of the targets
	std::vector< unsigned int > m_TimeIndices;

	/// Regression name.
	std::string m_RegressionName;
	
	/// Template Objects Name.
	std::vector< std::string > m_TemplateObjectsName;
	
	/// Extension of Objects Name.
	std::vector < std::string > m_TemplateObjectsNameExtension;

	/// Gradient of the functional w.r.t. the CP position.
	MatrixType m_GradPos;

	/// Gradient of the functional w.r.t. momentas
	MatrixType m_GradMom;
	
	/// Gradient w.r.t. template objects using the L2 metric
	MatrixList m_GradTempL_L2;
	
	/// Gradient w.r.t. template objects using the Sobolev metric (convolution of the L2 metric with a smoothing kernel)
	MatrixList m_GradTempL_Sob;

	/// Matrix of the coordinates of the momenta (Matrix size: NumberOfCPs x Dimension).
	MatrixType m_InitialMomenta;
	
	/// Number of control points in the Regression.
	int m_NumberOfCPs;

	/// Number of deformable objects of the target and the source.
	int m_NumberOfObjects;

	/// Number of timepoints.
	int m_NumberOfTimepoints;	

	/// Determines if control points are optimized or not.
	bool m_UpdateCP;
	
	/// Determines if baseline shapes are optimized
	bool m_UpdateTemplate;
    
	/// Type of optimization method
	OptimizationMethodType m_OptimizationMethod;
	
	/// FISTA parameter. Can't increase stepsize during NbIterFreeze iterations.
	static const unsigned int NbIterFreeze = 3;

	/// Vector containing a history of the the values of the functional during optimization method.
	std::vector< FunctionalValuesType* > m_ValuesHistory;

	unsigned int m_SaveEveryNIters;

	/// Maximum number of descent iterations.
	unsigned int m_MaxIterations;

	/// Maximum number of iterations in the line search procedure.
	unsigned int m_MaxLineSearchIterations;

	/// Shrink parameter for the line search.
	TScalar m_AdaptiveShrink;
	/// Expand parameter for the line search.
	TScalar m_AdaptiveExpand;
	/// The algorithm stops when \f$F(end-1)-F(end) < \verb#m_AdaptiveTolerance# * \left( F(0)-F(end) \right)\f$.
	TScalar m_AdaptiveTolerance;

	/// Initial step multiplier.
	TScalar m_InitialStepMultiplier;


}; /* class RegressionEstimator */


#ifndef MU_MANUAL_INSTANTIATION
#include "RegressionEstimator.txx"
#endif


#endif /* _RegressionEstimator_h */

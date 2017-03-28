/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _MultiClassBayesianAtlasEstimator_h
#define _MultiClassBayesianAtlasEstimator_h

#include "LinearAlgebra.h"
#include <vector>

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "readMatrixDLM.txx"

#include "Diffeos.h"

#include "Atlas.h"
#include "DeformableMultiObject.h"
#include "BayesianAtlas.h"

#include "BayesianAtlasFunctionalValues.h"

/**
 *  \brief      Multi Class Bayesian Atlas construction class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The MultiClassBayesianAtlasEstimator class enables to build atlases from collections of objects configurations.\n \n
 *              Given a population of M objects and N subjects, it clusters data and estimates for each cluster a common template complex T
 *              made of M template objects,
 *              the noise variance of each object and the covariance matrix of the momentum vectors also for each cleuas.
 */
template <class TScalar, unsigned int Dimension>
class MultiClassBayesianAtlasEstimator : public AbstractAtlasEstimator< TScalar, Dimension >
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
    /// List of vectors type.
    typedef std::vector<VectorType> VectorList;

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
    
    MultiClassBayesianAtlasEstimator();
    
    virtual ~MultiClassBayesianAtlasEstimator();
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Return first atlas.
    virtual AtlasType* GetAtlas() const
    {
        std::cout << "Warning: only atlas in first cluster is returned in MultiClassBayesianAtlasEstimator" << std::endl;
        
        AtlasType* out = m_AtlasList[0];
        return out;
    }
    /// Return Atlas List
    std::vector< AtlasType* > GetAtlasList() const { return m_AtlasList; }
    /// Return i-th atlas from the list
    AtlasType* GetAtlas(unsigned int i) const
    {
        if (i < m_MaximumNumberOfClasses)
            return m_AtlasList[i];
        else
            throw std::runtime_error("Class index out of bound");
    }
    
    /// Set atlas to \e obj.
    virtual void SetAtlas(AtlasType* exampleAtlas)
    {
        if (!exampleAtlas->IsBayesian())
            throw std::runtime_error("Atlas given in Bayesian Atlas Estimator is not Bayesian");
        
        m_ExampleAtlas = dynamic_cast<BayesianAtlasType*>(exampleAtlas);
        Superclass::m_NumberOfCPs = exampleAtlas->GetControlPoints().rows();
    }
    
    /// Returns the list of the momentas.
    inline std::vector<MatrixList> GetMomentaList() const { return m_InitialMomentas; }
//    /// Initializes the momentas by reading the file \e fn.
//    inline void SetInitialMomentas(std::string fn)
//    {
//        
//        std::cout << "Warning: Setting initial momenta is not permitted in MulticlassBayesianAtlasEstimator\n Initialization discarded" << std::endl;
//        return;
////        if (strlen(fn.c_str())) // null file name means no file to be read
////            m_InitialMomentas = readMultipleMatrixDLM<TScalar>(fn.c_str());
//    }

    /// Set maximum number of classes
    inline void SetMaximumNumberOfClasses(int i) { m_MaximumNumberOfClasses = i; }

    /// Get maximum number of classes
    inline unsigned int GetMaximumNumberOfClasses() const { return m_MaximumNumberOfClasses; }

    /// Get probability of each class
    inline VectorType GetClassProbability() const { return m_ClassProba; }

    /// Compute probabilities of each class given the class index of each subject in \e m_ClassProba
    void ComputeClassProbabilities();


    
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
     *  \param[in, out]	X	List of control points to be optimized or not (depending on
     *						MultiClassBayesianAtlasEstimator::m_freezeCP).
     *  \param[in, out]	A	List of momentas to be optimized.
     *  \param[in, out]	T	List of template to be optimized.
     *  \param[in, out]	DSS	List of vector of noise variance of each object to be optimized
     *  \param[in, out]	CovMom	List of covariance matrix of the momentum vectors to be optimized
     *  \param[in, out] SubjectClass class label for each subject
     
     */
    void GradientDescent(MatrixList& X, std::vector<MatrixList>& A, std::vector<MatrixList>& T, VectorList& DSS, MatrixList& CovMom, VectorType& SubjectClass);
    
    
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


    /// Compute Class Probabilities from vector assigning classes to subjects
    VectorType ComputeClassProbabilities(VectorType SubjectClass);


    /// Split a vector of objects (\e vector) in two vectors: the one containing the objects that are assigned to the cluster \e c and the other the remaining objects
    template<class ObjectType> bool Split(unsigned int c, const VectorType& SubjectClass, const std::vector<ObjectType>& vector, std::vector<ObjectType>& vectorInCluster, std::vector<ObjectType>& vectorOther);

    /// Converse method of Split
    template<class ObjectType> void Combine(unsigned int c, const VectorType& SubjectClass, std::vector<ObjectType>& vector, const std::vector<ObjectType>& vectorInCluster, const std::vector<ObjectType>& vectorOther);

    /// Convert vectors of type VNL to std::vector< TScalar >
    std::vector< TScalar > VNLtoSTDVector( const VectorType& in);



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// The example atlas.
    BayesianAtlasType* m_ExampleAtlas;

    /// The atlas List
    std::vector<BayesianAtlasType*> m_AtlasList;

    /// Maximum Number of Classes
    int m_MaximumNumberOfClasses;

    /// Probability of each class
    VectorType m_ClassProba;

    /// Class label of each observation
    VectorType m_SubjectClass;

    /// Gradient of the functional w.r.t. the CP position.
    MatrixType m_GradPos;
    
    /// Gradient of the functional w.r.t. momentas.
    MatrixList m_GradMom;
    
    /// Gradient w.r.t. template objects.
    MatrixList m_GradTempL;
    
    /// for each cluster, this is a list of matrices of the momenta coordinates (List size: NumberOfClusters times NumberOfSubjects, Matrices size: NumberOfCPs x Dimension).
    std::vector<MatrixList> m_InitialMomentas;
    
    /// Vector containing a history of the the values of the functional during optimization method.
    VectorType m_ValuesHistory;
//    std::vector< FunctionalValuesType* > m_ValuesHistory;


    
}; /* class MultiClassBayesianAtlasEstimator */


#ifndef MU_MANUAL_INSTANTIATION
#include "MultiClassBayesianAtlasEstimator.txx"
#endif


#endif /* _MultiClassBayesianAtlasEstimator_h */

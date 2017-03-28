/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _MultiClassBayesianEstimator_txx
#define _MultiClassBayesianEstimator_txx

#include "AbstractAtlasEstimator.h"

#include "KernelFactory.h"

#include "Diffeos.h"

#include "writeMatrixDLM.txx"

#include <iostream>
#include <iomanip>
#include <sstream>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::MultiClassBayesianAtlasEstimator() : Superclass()
{
    // // convert m_Atlas of type Atlas.h into BayesianAtlas.h
    // m_Atlas = dynamic_cast<AtlasType*>(Superclass::m_Atlas);
    m_MaximumNumberOfClasses = 1;
    
}



template <class TScalar, unsigned int Dimension>
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::~MultiClassBayesianAtlasEstimator()
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::Update()
{
    
    // Duplicate the example atlas, as initialization of each class
    m_AtlasList.resize(m_MaximumNumberOfClasses);
    for (unsigned int i = 0; i < m_MaximumNumberOfClasses; i++)
    {
        m_AtlasList[i] = m_ExampleAtlas->Clone();
    }
    
    // class labels are randomly assigned to each subject
    m_SubjectClass.set_size( Superclass::m_NumberOfSubjects );
    for (unsigned int i = 0; i < Superclass::m_NumberOfSubjects; i++) {
        m_SubjectClass[i] = rand() % m_MaximumNumberOfClasses;
    }
//    m_ClassProba.set_size( m_MaximumNumberOfClasses );
//    this->ComputeClassProbabilities( m_SubjectClass );
    
//    // Set initial momenta
//    if ( (m_InitialMomentas.size() == Superclass::m_NumberOfSubjects) && (m_InitialMomentas[0].rows() == Superclass::m_NumberOfCPs) && (m_InitialMomentas[0].columns() == Dimension) )
//    {
//        std::cout << "Using predefined set of momenta" << std::endl;
//    }
//    else
//    {
//        if (m_InitialMomentas.size() > 0)
//            std::cout << "Warning: initial momenta file has incompatible number of subjects and/or vectors. Initial momenta reset to zero" << std::endl;
//        
//        m_InitialMomentas.resize(Superclass::m_NumberOfSubjects);
//        for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//        {
//            m_InitialMomentas[s].set_size(Superclass::m_NumberOfCPs, Dimension);
//            m_InitialMomentas[s].fill(0.0);
//        }
//    }

    m_InitialMomentas.resize( m_MaximumNumberOfClasses );
    for (unsigned int i = 0; i < m_MaximumNumberOfClasses; i++)
    {
        m_InitialMomentas[i].resize( Superclass::m_NumberOfSubjects );
        for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
        {
            m_InitialMomentas[i][s].set_size(m_AtlasList[i]->GetControlPoints().rows(), Dimension);
            m_InitialMomentas[i][s].fill(0.0);
        }
    }
    
    
    m_ValuesHistory.set_size(Superclass::m_MaxIterations + 1);
    
    MatrixList ControlPoints( m_MaximumNumberOfClasses );
    std::vector<MatrixList> TemplateData(m_MaximumNumberOfClasses);
    VectorList DataSigmaSquared(m_MaximumNumberOfClasses);
    MatrixList CovarianceMomenta_Inverse(m_MaximumNumberOfClasses);

    for (unsigned int i = 0; i < m_MaximumNumberOfClasses; i++) {
        ControlPoints[i] = m_AtlasList[i]->GetControlPoints();
        TemplateData[i] = m_AtlasList[i]->GetTemplateData();
        DataSigmaSquared[i] = m_AtlasList[i]->GetDataSigmaSquared();
        CovarianceMomenta_Inverse[i] = m_AtlasList[i]->GetCovarianceMomentaInverse();
    }
    
    
    // Parameters of the model to be estimated
    MatrixList& X = ControlPoints;
    std::vector<MatrixList>& A = m_InitialMomentas;
    std::vector<MatrixList>& T = TemplateData;
    VectorList& DSS = DataSigmaSquared;
    std::vector<MatrixType>& CovMomInv = CovarianceMomenta_Inverse;
    VectorType& SubjectClass = m_SubjectClass;
    
    this->GradientDescent(X, A, T, DSS, CovMomInv, SubjectClass);
    
    
    for (unsigned int i = 0; i < m_MaximumNumberOfClasses; i++)
    {
        m_AtlasList[i]->SetControlPoints( ControlPoints[i] );
        m_AtlasList[i]->SetTemplateData( TemplateData[i] );
        m_AtlasList[i]->SetDataSigmaSquared( DataSigmaSquared[i] );
        m_AtlasList[i]->SetCovarianceMomentaInverse( CovarianceMomenta_Inverse[i] );
    }
    m_ClassProba = this->ComputeClassProbabilities( m_SubjectClass );
    
}




template <class TScalar, unsigned int Dimension>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::WriteOutput(std::string AtlasName)
{
    
    if ( ( AtlasName.size() == 0 ) || ( Superclass::m_TemplateObjectsName.size() == 0 ) || ( Superclass::m_TemplateObjectsNameExtension.size()) == 0 )
        throw std::runtime_error("No file names given");
    
    if ( ( Superclass::m_TemplateObjectsName.size() != Superclass::m_NumberOfObjects ) || ( Superclass::m_TemplateObjectsNameExtension.size() != Superclass::m_NumberOfObjects ) )
        throw std::runtime_error("Number of objects and objects'name mismatch in MultiClassBayesianAtlasEstimator::WriteOutput");
    
    std::vector<MatrixList> Mom = this->GetMomentaList();
    m_ClassProba = this->ComputeClassProbabilities( m_SubjectClass );

    for (unsigned int c = 0; c < m_MaximumNumberOfClasses; c++)
    {
        if (m_ClassProba(c) > 1e-20)
        {
            std::ostringstream oss;
            oss << AtlasName << "_cluster_" << c;
            
            if (Superclass::m_UpdateTemplate)
                m_AtlasList[c]->WriteTemplateData(oss.str(), Superclass::m_TemplateObjectsName, Superclass::m_TemplateObjectsNameExtension);
            
            m_AtlasList[c]->WriteAtlasParameters(oss.str());
            
            oss << "_InitialMomentas.txt" << std::ends;
            writeMultipleMatrixDLM<TScalar>(oss.str().c_str(), Mom[c]);
        }
    }

}



template <class TScalar, unsigned int Dimension>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::WriteAtlasToSubjectDeformations()
{

    std::cout << "Warning: WriteAtlasToSubjectDeformations not yet implemented in MultiClassBayesianAtlasEstimator" << std::endl;
//    std:vector<MatrixList> Mom = this->GetMomenta();
//    
//    for (unsigned int s = 0; s < Mom.size(); s++)
//    {
//        std::vector< std::string > tempfn(Superclass::m_NumberOfObjects);
//        for (unsigned int i = 0; i < Superclass::m_NumberOfObjects; i++)
//        {
//            std::ostringstream oss;
//            oss << Superclass::m_TemplateObjectsName[i] << "_to_subject_" << s;
//            tempfn[i] = oss.str();
//        }
//        
//        m_Atlas->WriteAtlasDeformation(Mom[s], Superclass::m_AtlasName, tempfn, Superclass::m_TemplateObjectsNameExtension);
//    }
}


template <class TScalar, unsigned int Dimension>
typename MultiClassBayesianAtlasEstimator<TScalar, Dimension>::VectorType
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::ComputeClassProbabilities(VectorType SubjectClass)
{
    VectorType out(m_MaximumNumberOfClasses, 0.0);
    for (unsigned int i = 0; i < SubjectClass.size(); i++)
    {
        out( SubjectClass(i) ) += 1.0;
    }
    out /= SubjectClass.size();
    
    return out;
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

//template <class TScalar, unsigned int Dimension>
//void
//MultiClassBayesianAtlasEstimator<TScalar, Dimension>
//::GradientDescent(MatrixList& X, std::vector<MatrixList>& A, std::vector<MatrixList>& T, VectorList& DSS, MatrixList& CovMomInv, VectorType& SubjectClass)
//{
//    
//    MatrixType Xnew;
//    MatrixList Anew( Superclass::m_NumberOfSubjects );
//    MatrixList Tnew( Superclass::m_NumberOfObjects );
//    
//    std::vector<TScalar> stepXA( m_MaximumNumberOfClasses );
//    std::vector<TScalar> stepT( m_MaximumNumberOfClasses );
//    
////    std::vector< std::vector< FunctionalValuesType* > > LikelihoodWithoutPriorsHistory(Superclass::m_MaxIterations +1 );
//    
//    for (unsigned int c = 0; c < m_MaximumNumberOfClasses; c++)
//    {
//        m_AtlasList[c]->SetControlPoints(X[c]);
//        m_AtlasList[c]->SetTemplateData(T[c]);
//        m_AtlasList[c]->SetDataSigmaSquared(DSS[c]);
//        m_AtlasList[c]->SetCovarianceMomentaInverse(CovMomInv[c]);
//    }
//    
//    // Update Class Proba according to the initial values in SubjectClass
//    VectorType ClassProba = this->ComputeClassProbabilities(SubjectClass);
//
//    
//    // reference values
//    TScalar totalLikelihood = 0.0;
//    std::vector<TScalar> lsqRef_T(m_MaximumNumberOfClasses, 0.0);
//    std::vector< std::vector<TScalar> > lsqRef_A(m_MaximumNumberOfClasses);
//    for (int c = 0; c < m_MaximumNumberOfClasses; c++)
//    {
//        lsqRef_A[c].resize(Superclass::m_NumberOfSubjects, 0.0);
//    }
//    
//    // Compute the optimal atlas parameters for each cluster (data sigma squared and momenta covariance matrix)
//    for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//    {
//        std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
//        MatrixList MomentaInCluster, MomentaOther;
//        
//        bool isEmpty = this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
//        this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
//        
//        if (isEmpty)
//        {
//            std::cout << "Cluster " << clusterIndex << " is empty" << std::endl;
//        }
//        
//        // Compute optimal values of data sigma and Convariance matrices for each cluster (using observations in this cluster only)
//        FunctionalValuesType* valuesInCluster = m_AtlasList[clusterIndex]->ComputeLikelihoodWithOptimalParameters(MomentaInCluster, targetsInCluster);
//
//        // Save these optimal values
//        DSS[clusterIndex] = valuesInCluster->GetDataSigmaSquared();
//        CovMomInv[clusterIndex] = valuesInCluster->GetCovMomInverse();
//        m_AtlasList[clusterIndex]->SetDataSigmaSquared(DSS[clusterIndex]);
//        m_AtlasList[clusterIndex]->SetCovarianceMomentaInverse(CovMomInv[clusterIndex]);
//        
//        // Compute the reference values of the functional (lsqRef_T, lsqRef_A and totallikelihood)
//        lsqRef_T[clusterIndex] = valuesInCluster->GetLikelihoodWithoutPriorTerms();
//        
//        int numSubjectsInCluster = MomentaInCluster.size();
//        totalLikelihood += valuesInCluster->GetLikelihood() + numSubjectsInCluster * log(ClassProba[clusterIndex]);
//
//        FunctionalValuesType* valuesOther = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(MomentaOther, targetsOther);
//
//        std::vector<TScalar> LikelihoodInCluster = valuesInCluster->GetLikelihoodWithoutPriorTermsForAllSubjects();
//        std::vector<TScalar> LikelihoodOther = valuesOther->GetLikelihoodWithoutPriorTermsForAllSubjects();
//        std::vector<TScalar> Likelihood;
//        this->Combine(clusterIndex, SubjectClass, Likelihood, LikelihoodInCluster, LikelihoodOther);
//        for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//        {
//            lsqRef_A[clusterIndex][s] = Likelihood[s];
//        }
//
//        std::cout << "Cluster " << clusterIndex << " ";
//        valuesInCluster->PrintIter(0);
//    }
//    m_ValuesHistory[0] = totalLikelihood;
//    
//    std::cout << "Iter 0 - total likelihood: " << totalLikelihood << std::endl << std::endl;
//    
//    
//    // ATTENTION: CHANGER LE TYPE DE m_GradPos , m_GradTempL, m_GradMom
//    
//    
//    
//    // gradient descent begins
//    for (unsigned int iter = 0; iter < Superclass::m_MaxIterations; iter++)
//    {
//        
//        bool foundMin = false;
//        
//        if (!(iter % 1))
//        {
//            std::ostringstream oss;
//            oss << "Iter" << std::setfill('0') << std::setw(4) << iter + 1 << "_" << Superclass::m_AtlasName;
//            this->WriteOutput(oss.str());
//        }
//
//        std::vector< FunctionalValuesType* > LikelihoodWithoutPriors_Ref(m_MaximumNumberOfClasses);
//        
//        
//        // compute the gradient w.r.t template shape and momenta
//        for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//        {
//
//            // Compute likelihood gradient separating terms from observations assigned to the current cluster from others assigned to another cluster
//            MatrixList MomentaInCluster, MomentaOther;
//            MatrixList gradMomInCluster, gradMomOther, Aux2;
//            MatrixType Aux;
//            std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
//            this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
//            this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
//            
//            m_AtlasList[clusterIndex]->ComputeLikelihoodGradient(MomentaInCluster, targetsInCluster, gradMomInCluster, m_GradPos[clusterIndex], m_GradTempL[clusterIndex]);
//            m_AtlasList[clusterIndex]->ComputeLikelihoodGradient(MomentaOther, targetsOther, gradMomOther, Aux, Aux2);
//            
//            this->Combine(clusterIndex, SubjectClass, m_GradMom[clusterIndex], gradMomInCluster, gradMomOther);
//        }
//        
//        // determine step size at first iteration
//        if (iter == 0)
//        {
//            TScalar maxGrad = 1e-20;
//            TScalar g = 0;
//            for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//            {
//                for (unsigned int i = 0; i < m_AtlasList[clusterIndex]->GetControlPoints().rows(); i++)
//                {
//                    g = m_GradPos[clusterIndex].get_row(i).magnitude();
//                    if (g > maxGrad)
//                        maxGrad = g;
//                    for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//                    {
//                        g = m_GradMom[clusterIndex][s].get_row(i).magnitude();
//                        if (g > maxGrad)
//                            maxGrad = g;
//                    }
//                }
//
//            }
//            
//            TScalar initStepXA = Superclass::m_InitialStepMultiplier * fabs(totalLikelihood) / maxGrad / maxGrad; // It is not clear whether totalLikelihood is here the good normalization for gradMom...
//            TScalar initStepT = initStepXA;
//                
//            stepXA = initStepXA;
//            for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//            {
//                stepT[clusterIndex] = initStepT;
//            }
//        }
//
//        // double line search for updating {template, control points} from observations from the cluster and {momenta} from all observations - acceptation is based on the decrease of the sum of the likelihood over newly assigned cluster
////        bool foundMinInCluster = false;
//        for (unsigned int li = 0; li < Superclass::m_MaxLineSearchIterations; li++)
//        {
//            
//            std::vector< FunctionalValuesType*> valuesTest(m_MaximumNumberOfClasses);
//            VectorType<unsigned int> foundMinInCluster_T(m_MaximumNumberOfClasses);
//            for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//            {
//                std::cout << "cluster " << clusterIndex << " stepsizeX = " << stepX[clusterIndex] << "\tstepsizeT = " << stepT[clusterIndex] << std::endl;
//
//                this->GradientDescentStep(clusterIndex, Xnew, X, Anew, A, Tnew, T, stepX, stepA, stepT);
//                m_AtlasList[clusterIndex]->SetControlPoints(Xnew);
//                m_AtlasList[clusterIndex]->SetTemplateData(Tnew);
//                
//                valuesTest[clusterIndex] = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(Anew[clusterIndex], Superclass::m_TargetList); // parameters of the atlas should not be updated here, they will be with the new subjects' class assignment.
//                
////                // test whether the update of the control points and template decrease the sum of likelihood terms for subjects within the current cluster
////                TScalar values_T = 0.0;
////                for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
////                {
////                    if (SubjectClass[s] == clusterIndex)
////                    {
////                        values_T += valuesTest[clusterIndex]->GetLikelihoodWithoutPriorTermsForSubject(s);
////                    }
////                }
////                TScalar Q_T = lsqRef_T[clusterIndex] - values_T;
////                foundMinInCluster_T[clusterIndex] = (Q_T >=0);
//            }
//            std::cout << "stepsizeXA = " << stepsizeXA << std::endl;
//            
//            // test whether the update of the control points and template decrease the sum of likelihood terms for subjects within the current cluster
//            TScalar values_T = 0.0
//            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//            {
//                    values_T += valuesTest[SubjectClass[s]]->GetLikelihoodWithoutPriorTermsForSubject(s);
//            }
//            TScalar Q_T = lsqRef_T[clusterIndex] - values_T;
//            bool foundMinInCluster_T = (Q_T >=0);
//          
//           
//            // find the optimal class assignement for the updated momenta
//            VectorType newSubjectClass(Superclass::m_NumberOfSubjects);
//            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//            {
//                int newClass = 0;
//                
//                TScalar minllh = valuesTest[0]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(0));
//                for (int c = 1; c < m_MaximumNumberOfClasses; c++)
//                {
//                    TScalar value = valuesTest[c]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(c));
//                    if (value < minllh)
//                        newClass = c;
//                }
////                if (newClass != SubjectClass[s])
////                    std::cout << "Subject " << s << " class has changed from " << SubjectClass[s] << " to " << newClass << std::endl;
//                
//                newSubjectClass[s] = newClass;
//            }
//
//            // check whether the momenta update enables the decrease in the sum of the likelihood terms within the new cluster
//            TScalar values_XA = 0.0;
//            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//            {
//                values_XA += valuesTest[newSubjectClass[s]]->GetLikelihoodWithoutPriorTermsForSubject(s);
//            }
//            TScalar Q_XA = lsqRef_XA - values_XA;
//            foundMin_XA = (Q_XA >= 0.0);
//            
//            
//            // if foundmin: ne pas oublier de faire SubjectClass = newSubjectClass;
//            bool foundmin;
//            if (foundMinInCluster_T.sum())
//            {
//                foundmin = true;
//                for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//                {
//                    if (foundMinInCluster_T(clusterIndex))
//                    {
//                        T[clusterIndex] = Tnew[clusterIndex];
//                        stepsizeT[clusterIndex] *= Superclass::m_AdaptiveExpand;
//                    }
//                    else
//                    {
//                        stepsize[clusterIndex] *= Superclass::m_AdaptiveShrink;
//                    }
//                }
//            }
//            
//            
//            
//            
//
//            
//                // check whether the likelihood is decreased for the template update
//            
//                // check whether each subject's contribution to the criterion is decreased
// // Before LinAlg : vnl_vector<unsigned int> isCriterionDecreased(Superclass::m_NumberOfSubjects, 0);
//                LinearAlgebra<unsigned int>::Vector isCriterionDecreased(Superclass::m_NumberOfSubjects, 0);
//                for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//                {
//                    TScalar Q_As = lsqRef_A[clusterIndex][s] - valuesTest->GetLikelihoodWithoutPriorTermsForSubject(s);
//                    if (Q_As >= 0)
//                    {
//                        isCriterionDecreased(s) = 1;
//                    }
//                    else
//                    {
//                        // remove updates which did not decrease the criterion
//                        Anew[s] = A[clusterIndex][s];
//                    }
//                    
//                }
//                std::cout << isCriterionDecreased << std::endl;
//                bool foundMinInCluster_A = (isCriterionDecreased.sum() > 0);
//                
//                if (foundMinInCluster_A && foundMinInCluster_T)
//                {
//                    foundMinInCluster = true;
////                    LikelihoodWithoutPriors_Ref[clusterIndex] = valuesTest->Clone();
//                    break;
//                }
//                else
//                {
//                    if (!foundMinInCluster_T)
//                        stepT[clusterIndex] *= Superclass::m_AdaptiveShrink;
//                    
//                    if (!foundMinInCluster_A)
//                        stepXA[clusterIndex] *= Superclass::m_AdaptiveShrink;
//                }
//            } // end of line search
//            
//                
//            if (foundMinInCluster)
//            {
//                X[clusterIndex] = Xnew;
//                A[clusterIndex] = Anew;
//                T[clusterIndex] = Tnew;
//                
//                stepXA[clusterIndex] *= Superclass::m_AdaptiveExpand;
//                stepT[clusterIndex] *= Superclass::m_AdaptiveExpand;
//                
//                foundMin = true;
//            }
//            else
//            {
//                // Loop terminated without finding smaller functional in this cluster
//                std::cout << " number of loops exceeded in line search for cluster " << clusterIndex << std::endl;
//                m_AtlasList[clusterIndex]->SetControlPoints(X[clusterIndex]);
//                m_AtlasList[clusterIndex]->SetTemplateData(T[clusterIndex]);
////                LikelihoodWithoutPriorsHistory[iter+1][clusterIndex] = LikelihoodWithoutPriorsHistory[iter][clusterIndex]->Clone();
//            }
//            
//            // need to re-compute the matchings since Anew may have changed since last call of m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(Anew, Superclass::m_TargetList)
//            LikelihoodWithoutPriors_Ref[clusterIndex] = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(A[clusterIndex], Superclass::m_TargetList);
//            
//        } // end loop over clusters
//        
//        
//        // update the subjects' class by finding the class which minimizes the contribution of the given subject to the likelihood
//        for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//        {
//            int newClass = 0;
//            
//            TScalar minllh = LikelihoodWithoutPriors_Ref[0]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(0));
//            for (int c = 1; c < m_MaximumNumberOfClasses; c++)
//            {
//                TScalar value = LikelihoodWithoutPriors_Ref[c]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(c));
//                if (value < minllh)
//                    newClass = c;
//            }
//            if (newClass != SubjectClass[s])
//                std::cout << "Subject " << s << " class has changed from " << SubjectClass[s] << " to " << newClass << std::endl;
//            
//            SubjectClass[s] = newClass;
//        }
//        
//        // update class probabilities
//        ClassProba = this->ComputeClassProbabilities(SubjectClass);
//
//        // update parameters of atlases
//        totalLikelihood = 0.0;
//        for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
//        {
//            // isolate momenta and targets from the current cluster
//            std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
//            MatrixList MomentaInCluster, MomentaOther;
//            
//            bool isEmpty = this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
//            this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
//
//            if (isEmpty)
//                std::cout << "Cluster " << clusterIndex << " is empty" << std::endl;
//            
//            // compute the optimal covariance matrix
//            MatrixType OptimalCovMomInverse;
//            VectorType CovMomTermsInCluster;
//            TScalar CovMomPriorInCluster;
//            TScalar llh_alpha = m_AtlasList[clusterIndex]->ComputeOptimalCovMomInverse(MomentaInCluster, OptimalCovMomInverse, CovMomTermsInCluster, CovMomPriorInCluster);
//            
//            // isolate residual norms in the current cluster
//            std::vector< std::vector< TScalar> > RN = LikelihoodWithoutPriors_Ref[clusterIndex]->GetResiduals();
//            std::vector< std::vector< TScalar> > RNInCluster, RNOther;
//            this->Split(clusterIndex, SubjectClass, RN, RNInCluster, RNOther);
//            
//            // compute the optimal data sigma squared
//            VectorType  OptimalDataSigmaSquared;
//            std::vector< std::vector<TScalar > > DataTermsInCluster;
//            VectorType DataPriorsInCluster;
//            TScalar llh_dds = m_AtlasList[clusterIndex]->ComputeOptimalDataSigmaSquared(RNInCluster, OptimalDataSigmaSquared, DataTermsInCluster, DataPriorsInCluster);
//            
//            
//            // set the updated parameters (that are given to the atlas at the beginning of the loop)
//            CovMomInv[clusterIndex] = OptimalCovMomInverse;
//            DSS[clusterIndex] = OptimalDataSigmaSquared;
//            
//            int NumberOfSubjectInCluster = MomentaInCluster.size();
//            totalLikelihood += (llh_alpha + llh_dds + NumberOfSubjectInCluster * log(ClassProba[clusterIndex]));
//            
//            lsqRef_T[clusterIndex] = CovMomTermsInCluster.sum();
//            
//            m_AtlasList[clusterIndex]->SetControlPoints(X[clusterIndex]);
//            m_AtlasList[clusterIndex]->SetTemplateData(T[clusterIndex]);
//            m_AtlasList[clusterIndex]->SetDataSigmaSquared(DSS[clusterIndex]);
//            m_AtlasList[clusterIndex]->SetCovarianceMomentaInverse(CovMomInv[clusterIndex]);
//
//            // update the values of lsqRef_A
//            // compute CovMomTerms and DataTerms for matching the cluster's template to samples that are not in this cluster
//            VectorType CovMomTermsOther;
//            std::vector< TScalar > CovMomTerms;
//            std::vector< std::vector< TScalar > > DataTerms, DataTermsOther;
//            m_AtlasList[clusterIndex]->ComputeCovMomTerms(MomentaOther, CovMomTermsOther);
//            m_AtlasList[clusterIndex]->ComputeDataTerms(RNOther, DataTermsOther);
//            
//            this->Combine(clusterIndex, SubjectClass, CovMomTerms, this->VNLtoSTDVector(CovMomTermsInCluster), this->VNLtoSTDVector(CovMomTermsOther));
//            this->Combine(clusterIndex, SubjectClass, DataTerms, DataTermsInCluster, DataTermsOther);
//            
//            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//            {
//                TScalar aux = 0.0;
//                for (unsigned int i = 0; i < DataTerms[s].size(); i ++)
//                    aux += DataTerms[s][i];
//                
//                lsqRef_A[clusterIndex][s] = aux + CovMomTerms[s];
//            }
//            
//            // mimic the FunctionalValues->PrintIter()
//            VectorType ResidualsInClusterPerObjects(RNInCluster[0].size(), 0.0);
//            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
//            {
//                for (unsigned int i = 0; i < RNInCluster[0].size(); i++)
//                    ResidualsInClusterPerObjects[i] += RNInCluster[s][i];
//            }
//
//            std::cout << "Iter " << iter+1 << " Cluster = " << clusterIndex << "  >> - log_Likelihood = " << (llh_alpha + llh_dds) << std::endl;
//            std::cout << "\t DataSigmaSquared = " << OptimalDataSigmaSquared << std::endl;
//            std::cout << "\t Residual Norms Per Objects = " << ResidualsInClusterPerObjects << std::endl;
//            std::cout << "\t CovMomTerms = " << CovMomTermsInCluster.sum() << std::endl;
//        
//        }
//        m_ValuesHistory[iter+1] = totalLikelihood;
//        
//        if (!foundMin)
//        {
//            std::cout << "no cluster has been updated at this iteration" << std::endl;
//            break;
//        }
//        
//        std::cout << "Iter " << iter + 1 << " total likelihood = " << m_ValuesHistory[iter+1] <<std::endl << std::endl;
//
//        if ( (m_ValuesHistory[iter] - m_ValuesHistory[iter+1])
//            < Superclass::m_AdaptiveTolerance * (m_ValuesHistory[0] - m_ValuesHistory[iter+1]) )
//            break;
//        
//    } // for iter
//}



template <class TScalar, unsigned int Dimension>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::GradientDescent(MatrixList& X, std::vector<MatrixList>& A, std::vector<MatrixList>& T, VectorList& DSS, MatrixList& CovMomInv, VectorType& SubjectClass)
{
    
    MatrixType Xnew;
    MatrixList Anew( Superclass::m_NumberOfSubjects );
    MatrixList Tnew( Superclass::m_NumberOfObjects );
    
    std::vector<TScalar> stepXA( m_MaximumNumberOfClasses );
    std::vector<TScalar> stepT( m_MaximumNumberOfClasses );
    
    //    std::vector< std::vector< FunctionalValuesType* > > LikelihoodWithoutPriorsHistory(Superclass::m_MaxIterations +1 );
    
    for (unsigned int c = 0; c < m_MaximumNumberOfClasses; c++)
    {
        m_AtlasList[c]->SetControlPoints(X[c]);
        m_AtlasList[c]->SetTemplateData(T[c]);
        m_AtlasList[c]->SetDataSigmaSquared(DSS[c]);
        m_AtlasList[c]->SetCovarianceMomentaInverse(CovMomInv[c]);
    }
    
    // Update Class Proba according to the initial values in SubjectClass
    VectorType ClassProba = this->ComputeClassProbabilities(SubjectClass);
    
    
    // reference values
    TScalar totalLikelihood = 0.0;
    std::vector<TScalar> lsqRef_T(m_MaximumNumberOfClasses, 0.0);
    std::vector< std::vector<TScalar> > lsqRef_A(m_MaximumNumberOfClasses);
    for (int c = 0; c < m_MaximumNumberOfClasses; c++)
    {
        lsqRef_A[c].resize(Superclass::m_NumberOfSubjects, 0.0);
    }
    
    // Compute the optimal atlas parameters for each cluster (data sigma squared and momenta covariance matrix)
    for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
    {
        std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
        MatrixList MomentaInCluster, MomentaOther;
        
        bool isEmpty = this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
        this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
        
        if (isEmpty)
        {
            std::cout << "Cluster " << clusterIndex << " is empty" << std::endl;
        }
        
        // Compute optimal values of data sigma and Convariance matrices for each cluster (using observations in this cluster only)
        FunctionalValuesType* valuesInCluster = m_AtlasList[clusterIndex]->ComputeLikelihoodWithOptimalParameters(MomentaInCluster, targetsInCluster);
        
        // Save these optimal values
        DSS[clusterIndex] = valuesInCluster->GetDataSigmaSquared();
        CovMomInv[clusterIndex] = valuesInCluster->GetCovMomInverse();
        m_AtlasList[clusterIndex]->SetDataSigmaSquared(DSS[clusterIndex]);
        m_AtlasList[clusterIndex]->SetCovarianceMomentaInverse(CovMomInv[clusterIndex]);
        
        // Compute the reference values of the functional (lsqRef_T, lsqRef_A and totallikelihood)
        lsqRef_T[clusterIndex] = valuesInCluster->GetLikelihoodWithoutPriorTerms();
        
        int numSubjectsInCluster = MomentaInCluster.size();
        totalLikelihood += valuesInCluster->GetLikelihood() + numSubjectsInCluster * log(ClassProba[clusterIndex]);
        
        FunctionalValuesType* valuesOther = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(MomentaOther, targetsOther);
        
        std::vector<TScalar> LikelihoodInCluster = valuesInCluster->GetLikelihoodWithoutPriorTermsForAllSubjects();
        std::vector<TScalar> LikelihoodOther = valuesOther->GetLikelihoodWithoutPriorTermsForAllSubjects();
        std::vector<TScalar> Likelihood;
        this->Combine(clusterIndex, SubjectClass, lsqRef_A[clusterIndex], LikelihoodInCluster, LikelihoodOther);
        
        std::cout << "Cluster " << clusterIndex << " ";
        valuesInCluster->PrintIter(0);
        
        std::cout << "lsqRef_T = " << lsqRef_T[clusterIndex] << std::endl;
        std::cout << "lsqRef_A = " << std::endl;
        for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
            std::cout << lsqRef_A[clusterIndex][s] << " ";
        std::cout << std::endl;

    }
    m_ValuesHistory[0] = totalLikelihood;
    
    std::cout << "Iter 0 - total likelihood: " << totalLikelihood << std::endl << std::endl;
    

    
    
    // gradient descent begins
    for (unsigned int iter = 0; iter < Superclass::m_MaxIterations; iter++)
    {
        
        bool foundMin = false;
        
        if (!(iter % 1))
        {
            std::ostringstream oss;
            oss << "Iter" << std::setfill('0') << std::setw(4) << iter + 1 << "_" << Superclass::m_AtlasName;
            this->WriteOutput(oss.str());
        }
        
        std::vector< FunctionalValuesType* > LikelihoodWithoutPriors_Ref(m_MaximumNumberOfClasses);
        
        for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
        {
            
            // Compute likelihood gradient separating terms from observations assigned to the current cluster from others assigned to another cluster
            MatrixList MomentaInCluster, MomentaOther;
            MatrixList gradMomInCluster, gradMomOther, Aux2;
            MatrixType Aux;
            std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
            this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
            this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
            
            m_AtlasList[clusterIndex]->ComputeLikelihoodGradient(MomentaInCluster, targetsInCluster, gradMomInCluster, m_GradPos, m_GradTempL);
            m_AtlasList[clusterIndex]->ComputeLikelihoodGradient(MomentaOther, targetsOther, gradMomOther, Aux, Aux2);
            
            this->Combine(clusterIndex, SubjectClass, m_GradMom, gradMomInCluster, gradMomOther);
            
            // determine step size at first iteration
            if (iter == 0)
            {
                TScalar maxGrad = 1e-20;
                for (unsigned int i = 0; i < m_AtlasList[clusterIndex]->GetControlPoints().rows(); i++)
                {
                    TScalar g = 0;
                    g = m_GradPos.get_row(i).magnitude();
                    if (g > maxGrad)
                        maxGrad = g;
                    for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
                    {
                        g = m_GradMom[s].get_row(i).magnitude();
                        if (g > maxGrad)
                            maxGrad = g;
                    }
                }
                
                TScalar initStepXA = Superclass::m_InitialStepMultiplier * fabs(totalLikelihood) / maxGrad / maxGrad; // It is not clear whether totalLikelihood is here the good normalization for gradMom...
                TScalar initStepT = initStepXA;
                
                stepXA[clusterIndex] = initStepXA;
                stepT[clusterIndex] = initStepT;
            }
            
            // double line search for updating {template, control points} from observations from the cluster and {momenta} from all observations
            bool foundMinInCluster = false;
            for (unsigned int li = 0; li < Superclass::m_MaxLineSearchIterations; li++)
            {
                std::cout << "stepsizeXA = " << stepXA[clusterIndex] << "\tstepsizeT = " << stepT[clusterIndex] << std::endl;
                
                this->GradientDescentStep(Xnew, X[clusterIndex], Anew, A[clusterIndex], Tnew, T[clusterIndex], stepXA[clusterIndex], stepT[clusterIndex]);
                
                m_AtlasList[clusterIndex]->SetControlPoints(Xnew);
                m_AtlasList[clusterIndex]->SetTemplateData(Tnew);
                
                FunctionalValuesType* valuesTest = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(Anew, Superclass::m_TargetList); // parameters of the atlas should not be updated at this stage, they will be with the new subjects' class assignment.
                
                // check whether the likelihood is decreased for the template update
                TScalar values_T = 0.0;
                for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
                {
                    if (SubjectClass[s] == clusterIndex)
                    {
                        values_T += valuesTest->GetLikelihoodWithoutPriorTermsForSubject(s);
                    }
                }
                
                
                std::cout << "in line search: values_T = " << values_T << std::endl;
                
                TScalar Q_T = lsqRef_T[clusterIndex] - values_T;
                bool foundMinInCluster_T = (Q_T >=0);
                if (!foundMinInCluster_T)
                {
                    Tnew = T[clusterIndex];
                }
                
                // check whether each subject's contribution to the criterion is decreased
// Before LinAlg                vnl_vector<unsigned int> isCriterionDecreased(Superclass::m_NumberOfSubjects, 0);
                LinearAlgebra<unsigned int>::Vector isCriterionDecreased(Superclass::m_NumberOfSubjects, 0);
                for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
                {
                    TScalar Q_As = lsqRef_A[clusterIndex][s] - valuesTest->GetLikelihoodWithoutPriorTermsForSubject(s);
                    if (Q_As >= 0)
                    {
                        isCriterionDecreased(s) = 1;
                    }
                    else
                    {
                        // remove updates which did not decrease the criterion
                        Anew[s] = A[clusterIndex][s];
                    }
                    
                }
                std::cout << isCriterionDecreased << std::endl;
                bool foundMinInCluster_A = (isCriterionDecreased.sum() > 0);
                
                
                if (!foundMinInCluster_A)
                {
                    Xnew = X[clusterIndex];
                }

                stepT[clusterIndex] *= (foundMinInCluster_T)?(Superclass::m_AdaptiveExpand):(Superclass::m_AdaptiveShrink);
                stepXA[clusterIndex] *= (foundMinInCluster_A)?(Superclass::m_AdaptiveExpand):(Superclass::m_AdaptiveShrink);
                
                if (foundMinInCluster_A || foundMinInCluster_T)
                {
                    foundMinInCluster = true;
                    break;
                }
//                else
//                {
//                    if (!foundMinInCluster_T)
//                        stepT[clusterIndex] *= Superclass::m_AdaptiveShrink;
//                    
//                    if (!foundMinInCluster_A)
//                        stepXA[clusterIndex] *= Superclass::m_AdaptiveShrink;
//                }
            } // end of line search
            
            
            if (foundMinInCluster)
            {
                X[clusterIndex] = Xnew;
                A[clusterIndex] = Anew;
                T[clusterIndex] = Tnew;
                
                m_AtlasList[clusterIndex]->SetControlPoints(X[clusterIndex]);
                m_AtlasList[clusterIndex]->SetTemplateData(T[clusterIndex]);

//                stepXA[clusterIndex] *= Superclass::m_AdaptiveExpand;
//                stepT[clusterIndex] *= Superclass::m_AdaptiveExpand;
                
                foundMin = true;
            }
            else
            {
                // Loop terminated without finding smaller functional in this cluster
                std::cout << " number of loops exceeded in line search for cluster " << clusterIndex << std::endl;
                m_AtlasList[clusterIndex]->SetControlPoints(X[clusterIndex]);
                m_AtlasList[clusterIndex]->SetTemplateData(T[clusterIndex]);
            }
            
            // need to re-compute the matchings since Anew may have changed since last call of m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(Anew, Superclass::m_TargetList)
            LikelihoodWithoutPriors_Ref[clusterIndex] = m_AtlasList[clusterIndex]->ComputeLikelihoodWithoutPriors(A[clusterIndex], Superclass::m_TargetList);
            
            
            // SANITY CHECK: compute the total likelihood before class re-assignment
            TScalar values_T = 0.0;
            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
            {
                if (SubjectClass[s] == clusterIndex)
                {
                    values_T += LikelihoodWithoutPriors_Ref[clusterIndex]->GetLikelihoodWithoutPriorTermsForSubject(s);
                }
            }
            std::cout << "\t\t values_T = " << values_T << std::endl;
            
            
            
        } // end loop over clusters
        
        
        
        
        
        // update the subjects' class by finding the class which minimizes the contribution of the given subject to the likelihood
        for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
        {
            int newClass = 0;
            
            TScalar minllh = LikelihoodWithoutPriors_Ref[0]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(0));
            for (int c = 1; c < m_MaximumNumberOfClasses; c++)
            {
                TScalar value = LikelihoodWithoutPriors_Ref[c]->GetLikelihoodWithoutPriorTermsForSubject(s) + log(ClassProba(c));
                if (value < minllh)
                    newClass = c;
            }
            if (newClass != SubjectClass[s])
                std::cout << "Subject " << s << " class has changed from " << SubjectClass[s] << " to " << newClass << std::endl;
            
            SubjectClass[s] = newClass;
        }
        
        // update class probabilities
        ClassProba = this->ComputeClassProbabilities(SubjectClass);
        
        // update parameters of atlases
        totalLikelihood = 0.0;
        for (unsigned int clusterIndex = 0; clusterIndex < m_MaximumNumberOfClasses; clusterIndex++)
        {
            // isolate momenta and targets from the current cluster
            std::vector<DeformableMultiObjectType*> targetsInCluster, targetsOther;
            MatrixList MomentaInCluster, MomentaOther;
            
            bool isEmpty = this->Split(clusterIndex, SubjectClass, A[clusterIndex], MomentaInCluster, MomentaOther);
            this->Split(clusterIndex, SubjectClass, Superclass::m_TargetList, targetsInCluster, targetsOther);
            
            if (isEmpty)
                std::cout << "Cluster " << clusterIndex << " is empty" << std::endl;
            
            // compute the optimal covariance matrix
            MatrixType OptimalCovMomInverse;
            VectorType CovMomTermsInCluster;
            TScalar CovMomPriorInCluster;
            TScalar llh_alpha = m_AtlasList[clusterIndex]->ComputeOptimalCovMomInverse(MomentaInCluster, OptimalCovMomInverse, CovMomTermsInCluster, CovMomPriorInCluster);
            
            // isolate residual norms in the current cluster
            std::vector< std::vector< TScalar> > RN = LikelihoodWithoutPriors_Ref[clusterIndex]->GetResiduals();
            std::vector< std::vector< TScalar> > RNInCluster, RNOther;
            this->Split(clusterIndex, SubjectClass, RN, RNInCluster, RNOther);
            
            // compute the optimal data sigma squared
            VectorType  OptimalDataSigmaSquared;
            std::vector< std::vector<TScalar > > DataTermsInCluster;
            VectorType DataPriorsInCluster;
            TScalar llh_dds = m_AtlasList[clusterIndex]->ComputeOptimalDataSigmaSquared(RNInCluster, OptimalDataSigmaSquared, DataTermsInCluster, DataPriorsInCluster);
            
            
            // set the updated parameters (that are given to the atlas at the beginning of the loop)
            CovMomInv[clusterIndex] = OptimalCovMomInverse;
            DSS[clusterIndex] = OptimalDataSigmaSquared;
            
            int NumberOfSubjectInCluster = MomentaInCluster.size();
            totalLikelihood += (llh_alpha + llh_dds + NumberOfSubjectInCluster * log(ClassProba[clusterIndex]));
            
            TScalar lsqRef_T_c = CovMomTermsInCluster.sum();
            for (unsigned int s = 0; s < DataTermsInCluster.size(); s++)
            {
                for (unsigned int i = 0; i < DataTermsInCluster[0].size(); i++)
                {
                    lsqRef_T_c += DataTermsInCluster[s][i];
                }
            }
            lsqRef_T[clusterIndex] =  lsqRef_T_c;
            
//            m_AtlasList[clusterIndex]->SetControlPoints(X[clusterIndex]);
//            m_AtlasList[clusterIndex]->SetTemplateData(T[clusterIndex]);
            m_AtlasList[clusterIndex]->SetDataSigmaSquared(DSS[clusterIndex]);
            m_AtlasList[clusterIndex]->SetCovarianceMomentaInverse(CovMomInv[clusterIndex]);
            
            // update the values of lsqRef_A
            // compute CovMomTerms and DataTerms for matching the cluster's template to samples that are not in this cluster
            VectorType CovMomTermsOther;
            std::vector< TScalar > CovMomTerms;
            std::vector< std::vector< TScalar > > DataTerms, DataTermsOther;
            m_AtlasList[clusterIndex]->ComputeCovMomTerms(MomentaOther, CovMomTermsOther);
            m_AtlasList[clusterIndex]->ComputeDataTerms(RNOther, DataTermsOther);
            
            this->Combine(clusterIndex, SubjectClass, CovMomTerms, this->VNLtoSTDVector(CovMomTermsInCluster), this->VNLtoSTDVector(CovMomTermsOther));
            this->Combine(clusterIndex, SubjectClass, DataTerms, DataTermsInCluster, DataTermsOther);
            
            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
            {
                TScalar aux = 0.0;
                for (unsigned int i = 0; i < DataTerms[s].size(); i ++)
                    aux += DataTerms[s][i];
                
                lsqRef_A[clusterIndex][s] = aux + CovMomTerms[s];
            }
            
            
            
            std::cout << "lsqRef_T = " << lsqRef_T[clusterIndex] << std::endl;
            std::cout << "lsqRef_A = " << std::endl;
            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
                std::cout << lsqRef_A[clusterIndex][s] << " ";
            std::cout << std::endl;
            
            
            // mimic the FunctionalValues->PrintIter()
            VectorType ResidualsInClusterPerObjects(RNInCluster[0].size(), 0.0);
            for (unsigned int s = 0; s < Superclass::m_NumberOfSubjects; s++)
            {
                for (unsigned int i = 0; i < RNInCluster[0].size(); i++)
                    ResidualsInClusterPerObjects[i] += RNInCluster[s][i];
            }
            
            std::cout << "Iter " << iter+1 << " Cluster = " << clusterIndex << "  >> - log_Likelihood = " << (llh_alpha + llh_dds) << std::endl;
            std::cout << "\t DataSigmaSquared = " << OptimalDataSigmaSquared << std::endl;
            std::cout << "\t Residual Norms Per Objects = " << ResidualsInClusterPerObjects << std::endl;
            std::cout << "\t CovMomTerms = " << CovMomTermsInCluster.sum() << std::endl;
            
        }
        m_ValuesHistory[iter+1] = totalLikelihood;
        
        if (!foundMin)
        {
            std::cout << "no cluster has been updated at this iteration" << std::endl;
            break;
        }
        
        std::cout << "Iter " << iter + 1 << " total likelihood = " << m_ValuesHistory[iter+1] <<std::endl << std::endl;
        
        if ( (m_ValuesHistory[iter] - m_ValuesHistory[iter+1])
            < Superclass::m_AdaptiveTolerance * (m_ValuesHistory[0] - m_ValuesHistory[iter+1]) )
            break;
        
    } // for iter
}


template <class TScalar, unsigned int Dimension>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::GradientDescentStep(
                      MatrixType& Xtest, const MatrixType& X,
                      MatrixList& Atest, const MatrixList& A,
                      MatrixList& Ttest, const MatrixList& T,
                      TScalar stepXA, TScalar stepT)
{
    
    if (Superclass::m_UpdateCP)
        Xtest = X - m_GradPos*stepXA;
    else
        Xtest = X;
    
    for (int s = 0; s < Superclass::m_NumberOfSubjects; s++)
        Atest[s] = A[s] - m_GradMom[s] * stepXA;
    
    for (int i = 0; i < Superclass::m_NumberOfObjects; i++)
        Ttest[i] = T[i] - m_GradTempL[i] * stepT;
    
    return;
    
}



template <class TScalar, unsigned int Dimension> template<class ObjectType>
bool
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::Split(unsigned int c, const VectorType& SubjectClass, const std::vector<ObjectType>& vector, std::vector<ObjectType>& vectorInCluster, std::vector<ObjectType>& vectorOther)
{
    
    unsigned int NbInCluster = 0;
    unsigned int NbOther = 0;
    for (unsigned int s = 0; s < SubjectClass.size(); s++)
        (SubjectClass(s) == c)?(NbInCluster++):(NbOther++);

    
    vectorInCluster.resize(NbInCluster);
    vectorOther.resize(NbOther);
    unsigned int CInCluster = 0;
    unsigned int COther = 0;
    for (unsigned int s = 0; s < SubjectClass.size(); s++)
    {
        (SubjectClass(s) == c)?(vectorInCluster[CInCluster++] = vector[s]):(vectorOther[COther++] = vector[s]);
    }
    
    return (NbInCluster==0);
    
}

template <class TScalar, unsigned int Dimension> template<class ObjectType>
void
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::Combine(unsigned int c, const VectorType& SubjectClass, std::vector<ObjectType>& vector, const std::vector<ObjectType>& vectorInCluster, const std::vector<ObjectType>& vectorOther)
{
    
    vector.resize(SubjectClass.size());
    
    unsigned int CInCluster = 0;
    unsigned int COther = 0;
    for (unsigned int s = 0; s < SubjectClass.size(); s++)
    {
        (SubjectClass(s) == c)?(vector[s] = vectorInCluster[CInCluster++]):(vector[s] = vectorOther[COther++]);
    }
    
    return;
}

template <class TScalar, unsigned int Dimension>
std::vector< TScalar >
MultiClassBayesianAtlasEstimator<TScalar, Dimension>
::VNLtoSTDVector(const VectorType& in)
{
    std::vector< TScalar > out(in.size());

    for (unsigned int i = 0; i < in.size(); i++)
    {
        out[i] = in(i);
    }

    return out;
}





#endif /* _MultiClassBayesianEstimator_txx */

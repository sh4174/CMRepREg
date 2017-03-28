/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AbstractAtlasEstimator_h
#define _AbstractAtlasEstimator_h

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "LinearAlgebra.h"
#include <vector>

#include "readMatrixDLM.txx"

#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "Atlas.h"

/**
 *  \brief      Abstract Atlas estimation
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The AbstractAtlasEstimator class enables to compute the optimal parameters of a given atlas for a given set
 *              of shape data (i.e. the observations).
 */
template <class TScalar, unsigned int Dimension>
class AbstractAtlasEstimator
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Deformable multi-object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	/// List of deformable multi-objects type.
	typedef std::vector<DeformableMultiObjectType*> DeformableMultiObjectList;
	/// Atlas type type.
	typedef Atlas<TScalar, Dimension> AtlasType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	AbstractAtlasEstimator()
	{
		this->SetMaxIterations(10);
		m_MaxLineSearchIterations = 10;

		m_AdaptiveExpand = 1.2;
		m_AdaptiveShrink = 0.5;
		m_AdaptiveTolerance = 1e-4;

		m_InitialStepMultiplier = 1.0;

		m_UpdateCP = true;		
		m_UpdateTemplate = true;
	}


	virtual ~AbstractAtlasEstimator() {}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the list of targets to \e obj.
	inline void SetTargetList(DeformableMultiObjectList& obj)
	{
		m_TargetList = obj;
		m_NumberOfSubjects = obj.size();
		m_NumberOfObjects = obj[0]->GetNumberOfObjects();
	}

	/// Sets the atlas to \e obj.
	virtual void SetAtlas(AtlasType* obj) = 0;
	// {
	// 	m_Atlas = obj;
	// 	m_NumberOfCPs = m_Atlas->GetControlPoints().rows(); // this number of CP is not supposed to change during optimization
	// }
	/// Returns the output atlas.
	virtual AtlasType* GetAtlas() const = 0; // { return m_Atlas; }
	
	// // PIETRO -> Remark: Why not in AbstractAtlasEstimator?
	// /// Returns the list of the momentas.
	// inline MatrixList GetMomenta() const { return m_InitialMomentas; }
	// /// Initializes the momentas by reading the file \e fn.
	// inline void SetInitialMomentas(std::string fn)
	// {
	// 	if (strlen(fn.c_str())) // null file name means no file to be read
	// 		m_InitialMomentas = readMultipleMatrixDLM<TScalar>(fn.c_str());
	// }
	

	/// Disables optimization in control point positions.
	inline void SetFreezeCP() { m_UpdateCP = false; }
	/// Enables optimization in control point positions.
	inline void UnsetFreezeCP() { m_UpdateCP = true; }

	/// Disables optimization of the template shapes (for atlas-to-subject matching for instance)
	inline void SetFreezeTemplate() { m_UpdateTemplate = false; }
	/// Enables optimization of the template shape (for atlas estimation)
	inline void UnsetFreezeTemplate() { m_UpdateTemplate = true; }

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
	inline void SetAtlasName(std::string AtlasName) { m_AtlasName = AtlasName; }
	
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
	virtual void Update() = 0;
	
	/// Writes output (Atlas + subject-specific parameters).
	virtual void WriteOutput(std::string AtlasName)
	{		
		if ( ( AtlasName.size() == 0 ) || ( m_TemplateObjectsName.size() == 0 ) || ( m_TemplateObjectsNameExtension.size()) == 0 )
			throw std::runtime_error("No file names given");
		
		if ( ( m_TemplateObjectsName.size() != m_NumberOfObjects ) || ( m_TemplateObjectsNameExtension.size() != m_NumberOfObjects ) )
			throw std::runtime_error("Number of objects and objects'name mismatch in AbstractAtlasEstimator::WriteOutput"); 		
		
		if (m_UpdateTemplate)
			this->GetAtlas()->WriteTemplateData(AtlasName, m_TemplateObjectsName, m_TemplateObjectsNameExtension);
		
		this->GetAtlas()->WriteAtlasParameters(AtlasName);
	}
	
	/// Writes output (Atlas + subject-specific parameters).
	virtual void WriteOutput() { this->WriteOutput(m_AtlasName); }

	/// Writes template-to-subjects deformations
	virtual void WriteAtlasToSubjectDeformations() = 0;
	

protected:


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

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

	/// Determines if control points are optimized or not.
	bool m_UpdateCP;
	
	/// Determines if template shapes are optimized
	bool m_UpdateTemplate;

	/// \f$\nbsubj\f$ collections of deformable target objects.
	DeformableMultiObjectList m_TargetList;

	// /// Pointer on the atlas (template + control points + possibly covariance parameters for Bayesian Atlases)
	// AtlasType* m_Atlas;
	 
	/// Number of control points in the Atlas.
	int m_NumberOfCPs;

	/// Number of deformable objects of the target and the source.
	int m_NumberOfObjects;

	/// Number of subjects.
	int m_NumberOfSubjects;
	
	/// Atlas name.
	std::string m_AtlasName;
	
	/// Template Objects Name.
	std::vector< std::string > m_TemplateObjectsName;
	
	/// Extension of Objects Name.
	std::vector < std::string > m_TemplateObjectsNameExtension;
	

protected:


}; /* class AbstractAtlasEstimator */


#endif /* _AbstractAtlasEstimator_h */

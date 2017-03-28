/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Regression_h
#define _Regression_h

#include "DeformableObject.h"
#include "DeformableMultiObject.h"
#include "Diffeos.h"

#include "LinearAlgebra.h"
#include <vector>

#include "itkImage.h"

#include "itkSimpleFastMutexLock.h"

#include "readMatrixDLM.txx"

#include "DeformationFieldIO.h"

#include "RegressionFunctionalValues.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *	\brief      Regression object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A Regression is the combination of a template shape along with control points/initial momenta
 *		    that deform the template shape to match closely time-indexed shape observations. 
 */
template <class TScalar, unsigned int Dimension>
class Regression
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	/// List of matrices type.
	typedef std::vector<MatrixType> MatrixList;

	/// Multi-Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	typedef std::vector<DeformableObjectType*> DeformableObjectList;

	/// ITK image type.
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;
	
	/// Diffeos type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;

	/// Values of the different terms in the log-likelihood
	typedef RegressionFunctionalValues<TScalar> FunctionalValuesType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Default constructor
	Regression();
	/// Copy constructor.
	Regression(const Regression& other);

	/// Makes a copy of the object.
	Regression* Clone() { return new Regression(*this); }

	// Destructor
	~Regression();
	

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the sparsity prior.
	inline TScalar GetSparsityPrior() const { return m_SparsityPrior; }
	/// Sets the sparsity prior to \e d.
	inline void SetSparsityPrior(const TScalar d) { m_SparsityPrior = d; }

	/// Returns the deformable objects.
	inline DeformableMultiObjectType* GetTemplate() const { return m_Template; }
	/// Sets the deformable objetcs
	inline void SetTemplate(DeformableMultiObjectType* objects) { m_Template = objects; m_NumberOfObjects = m_Template->GetNumberOfObjects(); }
	
	/// Returns image intensity and landmark point coordinates of the template.
	inline MatrixList GetTemplateData() const { return m_Template->GetImageIntensityAndLandmarkPointCoordinates(); }
	/// Sets Image Intensity and landmark point coordinates of the template
	inline void SetTemplateData(const MatrixList& tempData) { m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(tempData); m_Template->Update(); }

	/// Returns data sigma squared of each object.
	inline VectorType GetDataSigmaSquared() const { return m_DataSigmaSquared; }
	/// Returns data sigma squared of object \e i.
	inline TScalar GetDataSigmaSquared(int i) const { return m_DataSigmaSquared[i]; }
	/// Sets data sigma squared of each object to \e d2.
	inline void SetDataSigmaSquared(VectorType d2) { m_DataSigmaSquared = d2; };
	/// Sets data sigma squared of object \e i to \e d2.
	inline void SetDataSigmaSquared(TScalar d2, int i) { m_DataSigmaSquared[i] = d2; };

	/// Sets control points stored in the file \e fn
	inline void SetControlPoints(std::string& fn);
	/// Updates positions of control points
	inline void SetControlPoints(MatrixType& CP){ m_ControlPoints = CP; }
	/// Returns control points positions.
	inline MatrixType GetControlPoints() const { return m_ControlPoints; }
		
	/// Sets control point spacing. Used in case no control points have been set to define a regular lattice of control points.
	inline void SetCPSpacing (TScalar s) { m_CPSpacing = s; }
	
	/// Sets the size of the smoothing kernel to \e d. See Regression::ConvolveGradTemplate().
	inline void SetSmoothingKernelWidth(TScalar d) { m_SmoothingKernelWidth = d; }
	
	/// Returns a tight bounding box including baseline shape objects and control points.
	inline MatrixType GetBoundingBox() { return m_BoundingBox; }
	
	/// Sets Diffeos to deform the template.
	inline void SetDiffeos(DiffeosType* def) { m_Def = def; }

				
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the bounding box and initializes the control points if needed.
	virtual void Update();

	/// Computes the value of the functional
	FunctionalValuesType* ComputeFunctional(const MatrixType& momenta, const std::vector< DeformableMultiObjectType* > target,
						std::vector<unsigned int> timeIndices);

	/// Computes the gradient of the functional
	void ComputeFunctionalGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
			MatrixType& gradMom, MatrixType& gradPos, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob);

	/// Computes the gradient of the functional
	void ComputeFunctionalGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
			MatrixType& gradMom, MatrixType& gradPos, MatrixList& gradTempL, bool do_SobolevGradient);

	/// Computes the gradient of regularity portion of the functional			
	void AddGradientRegularityTerm(MatrixType gradPos, MatrixType gradMom, const MatrixType& momenta);

	/// Saves the current template shape(s).
	void WriteTemplateFlow(const MatrixType& momenta, const std::string& regressionName,
		    	       const std::vector<std::string>& templateObjectsName, const std::vector<std::string>& templateObjectsNameExtension);
	
	/// Saves regression parameters.
	void WriteRegressionParameters(const MatrixType& momenta, const std::string& regressionName);
		

protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// If no control points have been set, this method generates a regular lattice of control points.
	void InitializeControlPoints();
	
	/**
	 *  \brief      Computes the residuals (i.e. the squared norm of the difference between each deformed baseline objects and target objects).
	 *
	 *  \details    For all \e timepoints, this method computes for each object k=1,..,NumberOfObjects the following residual :
	 *
	 *  \param[in]    Momentas      The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
	 *  \param[in]    target        The given target \f$S\f$.
	 *  \param[in]    timeIndices   List of time indices corresponding to targets
	 *  \param[out]   Residuals     List of the the residuals.
	 */
	bool ComputeResiduals(const MatrixType& momenta, const std::vector<DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices, 
			      std::vector< std::vector<TScalar> >& residuals);

	/**
	*  \brief      Computes the gradient ...
	*
	*  \details    For a series of \e targets over time, this method computes the gradient of ...
	*
	*   with respect to \f$\alpha_i\f$, \f$X_{0,k}\f$ and control points positions.
	*
	*  \param[in]   Momentas        The momentas parameterizing a series of deformations \f$\phi^{\alpha_i}\f$.
	*  \param[in]   target          The series of targets \f$S_i\f$.
	*  \param[in]   timeIndices     List of time indices corresponding to targets
	*  \param[out]  gradPos         The gradient of the residual w.r.t. the control points (summing contributions of each subject).
	*  \param[out]  gradMom         The gradient of the residual w.r.t. the series of momentas.
	*  \param[out]  gradTempL_L2    The gradient of the residual w.r.t. the template (summing contributions of each subject),
	*                               using the \f$L^2\f$ metric for template objects of landmark types (and children types).
	*  \param[out]  gradTempL_Sob   The gradient of the residual w.r.t. the template (summing contributions of each subject),
	*                               using a Sobolev metric for template objects of landmark types (and children types). It is
	*                               the convolution of \e gradTempL_L2 with a smoothing kernel, which guarantees that template meshes
	*                               do not self-intersect during optimization.
	*/
	void ComputeDataTermGradient(const MatrixType& momentas, const std::vector<DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
				     MatrixType& gradPos, MatrixType& gradMom, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob);
	/// Compute the data term gradient and return the Sobolev gradient of the template objects if \e do_Sobolev_Template_Gradient and the \f$L^2\f$ gradient otherwise
	void ComputeDataTermGradient(const MatrixType& momentas, const std::vector< DeformableMultiObjectType* > target, std::vector<unsigned int> timeIndices,
				     MatrixType& gradPos, MatrixType& gradMom, MatrixList& gradTempL, bool do_Sobolev_Template_Gradient);

	/// Computes the Sobolev gradient of template objects of landmark type (or children types) from the \f$L^2\f$ gradient
	MatrixList ConvolveGradTemplate(MatrixList& gradTemplate_L2);
	
	/// Updates the data domain of the deformation and the kernel.
	void UpdateDeformationAndKernelDataDomain(const std::vector< DeformableMultiObjectType* > target);
	
	/// Converts a VNLMatrix of size N x Dimension to a VNLVector \e V of length Dimension x N.
	inline VectorType Vectorize(const MatrixType& M) { return M.vectorise_row_wise(); }
	/// Converts a VNLVector \e V of length Dimension x N to a VNLMatrix of size N x Dimension. The inverse operation of Vectorize
	inline MatrixType VectorToMatrix(const VectorType& V)
	{
		//////// Before LinAlg ///////////////////////////////////////////////
		//		int nrow = V.size()/Dimension;
		//		return MatrixType(V.data_block(), nrow, Dimension);
		//////////////////////////////////////////////////////////////////////*
		int nrow = V.size()/Dimension;
		return V.convert_to_matrix_row_wise(nrow, Dimension);
	}
	
	

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Protected attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Template as a DeformableMultiObject.
	DeformableMultiObjectType* m_Template;
	
	/// Working deformation to compute regression deformation.
	DiffeosType* m_Def;

	/// Number of objects in the template.
	int m_NumberOfObjects;

	/// Control points.
	MatrixType m_ControlPoints;

	/// CP spacing used if no set of control points is given.
	TScalar m_CPSpacing;
		
	/// Template bounding box: the union of the template bounding box and the bounding box around control points.
	MatrixType m_BoundingBox;
	
	/// Data sigma squared for each object.
	VectorType m_DataSigmaSquared;
		
	/// Kernel width to compute the Sobolev gradient of the functional/likelihood w.r.t. template variable.
	TScalar m_SmoothingKernelWidth;
	
	/// Coefficient added to the \f$L^1\f$ penalty.
	TScalar m_SparsityPrior;


}; /* class Regression */


#ifndef MU_MANUAL_INSTANTIATION
#include "Regression.txx"
#endif


#endif /* _Regression_h */

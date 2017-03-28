/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Atlas_h
#define _Atlas_h

#include "DeformableMultiObject.h"
#include "Diffeos.h"

#include "LinearAlgebra.h"
#include <vector>

#include "itkImage.h"

#include "itkSimpleFastMutexLock.h"

#include "readMatrixDLM.txx"

// Added by ABGF@11aout14
#include "DeformationFieldIO.h"
// Added by ABGF@11aout14

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *	\brief      Atlas object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    An atlas is the combination of a template shape, control points, and parameters of variability
 *	            such as covariance matrices of the momentum vectors.
 */
template <class TScalar, unsigned int Dimension>
class Atlas
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Possible type of Atlas.
	typedef enum
	{
		null,				/*!< Null value. */
		DeterministicAtlas,	/*!< Simple deterministic atlas (See DeterministicAtlas). */
		BayesianAtlas,		/*!< Atlas with baysesian estimation of its parameter (See BayesianAtlas). */
	} AtlasType;
	

	/// Vector type
	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	/// Matrix type.
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	/// List of matrices type.
	typedef std::vector<MatrixType> MatrixList;

	/// Multi-Deformable object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

	/// ITK image type.
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;
	
	/// Diffeos type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Atlas();
	/// Copy constructor.
	Atlas(const Atlas& other);

	/// Makes a copy of the object.
	Atlas* Clone() { return new Atlas(*this); }

	virtual ~Atlas();
	


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns true if the atlas is of Deterministic Type.
	inline bool IsBayesian() { return (m_Type == BayesianAtlas); }
	/// Sets Bayesian Atlas Type.
	inline void SetBayesianAtlasType() { m_Type = BayesianAtlas; }
	/// Returns true if atlas is of Bayesian Type.
	inline bool IsDeterministic() { return (m_Type == DeterministicAtlas); }
	/// Sets simple Atlas type.
	inline void SetDeterministicAtlasType() { m_Type = DeterministicAtlas; }

	/// Returns the deformable objects.
	inline DeformableMultiObjectType* GetTemplate() const { return m_Template; }
	/// Sets the deformable objetcs to \e objects.
	inline void SetTemplate(DeformableMultiObjectType* objects) {
		m_Template = objects; m_NumberOfObjects = m_Template->GetNumberOfObjects(); }
	
	/// Returns image intensity and landmark point coordinates of the template.
	inline MatrixList GetTemplateData() const {	return m_Template->GetImageIntensityAndLandmarkPointCoordinates(); }
	/// Sets Image Intensity and landmark point coordinates of the template
	inline void SetTemplateData(const MatrixList& tempData) {
		m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(tempData); m_Template->Update(); }

	/// Returns data sigma squared of each object.
	inline VectorType GetDataSigmaSquared() const { return m_DataSigmaSquared; }
	/// Returns data sigma squared of object \e i.
	inline TScalar GetDataSigmaSquared(int i) const { return m_DataSigmaSquared[i]; }
	/// Sets data sigma squared of each object to \e d2.
	inline void SetDataSigmaSquared(VectorType d2) { m_DataSigmaSquared = d2; };
	// /// Sets data sigma squared of object \e i to \e d2.
	// inline void SetDataSigmaSquared(TScalar d2, int i) { m_DataSigmaSquared[i] = d2; };

	/// Sets control points stored in the file \e fn
	inline void SetControlPoints(std::string& fn);
	/// Updates positions of control points in the atlas.
	inline void SetControlPoints(MatrixType& CP){ m_ControlPoints = CP; }
	/// Returns control points positions.
	inline MatrixType GetControlPoints() const { return m_ControlPoints; }
	
	// PIETRO
	/// Sets inverse of covariance matrix of momentas from file fn..
	inline void SetCovarianceMomentaInverse(std::string& fn);
	/// Sets inverse of covariance matrix of momentas.
	inline void SetCovarianceMomentaInverse(MatrixType M) { m_CovMomInverse = M; }
	/// Returns the inverse of the covariance matrix of momentas.
	inline MatrixType GetCovarianceMomentaInverse() const { return m_CovMomInverse; }
	
	/// Sets control point spacing. Used in case no control points have been set to define a regular lattice of control points.
	inline void SetCPSpacing (TScalar s) { m_CPSpacing = s; }
	
	/// Sets the size of the smoothing kernel to \e d. See Atlas::ConvolveGradTemplate().
	inline void SetSmoothingKernelWidth(TScalar d) { m_SmoothingKernelWidth = d; }
	
	/// Returns a tight bounding box including template objects and control points.
	inline MatrixType GetBoundingBox() { return m_BoundingBox; }
	
	/// Returns the number of threads.
	inline unsigned int GetNumberOfThreads() const { return m_NumberOfThreads; }
	/// Sets the number of threads to \e n.
	void SetNumberOfThreads(unsigned int n) { m_NumberOfThreads = n; }
	
	/// Sets Diffeos to deform the template.
	inline void SetDiffeos(DiffeosType* def) { m_Def = def; }

		
		
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the bounding box and initializes the control points if needed.
	virtual void Update();

	/// Saves the current template.
	virtual void WriteTemplateData(const std::string& AtlasName,
			const std::vector<std::string>& TemplateObjectsNames, const std::vector<std::string>& TemplateObjectsNamesExtension);
	
	/// Saves the deformation of the template to every subject.
	virtual void WriteAtlasDeformation(const MatrixType& Momenta, const std::string& AtlasName,
			const std::vector<std::string>& TemplateObjectsNames, const std::vector<std::string>& TemplateObjectsNamesExtension);

	/// Saves atlas parameters (such as optimal control points, covariance matrices, etc..)
	virtual void WriteAtlasParameters(const std::string& AtlasName);
	
	// virtual TScalar ComputeLikelihood(const MatrixType& Momentas, const DeformableObjectType* target) = 0;
	// virtual TScalar ComputeLikelihood(const MatrixList& Momentas, const std::vector<DeformableObjectType*> target) = 0;
	

protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// If no control points have been set, this method generates a regular lattice of control points.
	void InitializeControlPoints();
	
	/**
	 *  \brief      Computes the residuals (i.e. the squared norm of the difference between each deformed template objects and target objects).
	 *
	 *  \details    For a given \e target, this method computes for each object k=1,..,NumberOfObjects the following residual :
	 *              \f[
	 *              \left\Vert\phi^{\alpha}(X_{0,k}) - S_k\right\Vert_{W}^2.
	 *              \f]
	 *
	 *  \param[in]    Momentas    The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
	 *  \param[in]    target      The given target \f$S\f$.
	 *  \param[out]   Residuals   The residuals.
	 */
	bool ComputeResidualsSubject(const MatrixType& Momentas, const DeformableMultiObjectType* target, std::vector< TScalar >& Residuals);
	/// Computes the residuals for a series of \e targets (representing a series of subjects).
	bool ComputeResiduals(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> target,
			std::vector< std::vector< TScalar > >& Residuals);
	
	/**
	 *  \brief      Computes the gradient of the data term with respect to momenta, control points and template objects. The data term is the sum of the residuals normalized by 2 times the data sigma squared.
	 *
	 *  \details    For a given \e target, this method computes the gradient
	 *              of :
	 *              \f[
	 *              \sum_{k=1}^{NumberOfObjects}\frac{1}{2\sigma_k^2}\left\Vert\phi^{\alpha}(X_{0,k}) - S_k\right\Vert_{W}^2.
	 *              \f]
     *   with respect to \f$\alpha\f$, \f$X_{0,k}\f$ and control points positions.
	 *
	 *  \param[in]    Momentas    The momentas parameterizing the deformation \f$\phi^{\alpha}\f$.
	 *  \param[in]    target      The given target \f$S\f$.
	 *  \param[out]   dPos        The gradient of the residual w.r.t. the control points.
	 *  \param[out]   dMom        The gradient of the residual w.r.t. the momentas.
	 *  \param[out]   dTempL      The gradient of the residual w.r.t. the template.
	 */
	void ComputeDataTermGradientSubject(const MatrixType& Momentas, const DeformableMultiObjectType* target,
			MatrixType& dPos, MatrixType& dMom, MatrixList& dTempL);
    /**
    *  \brief      Computes the gradient of the sum of subject's data terms with respect to momenta, control points and template objects. The data term is the sum of the residuals normalized by 2 times the data sigma squared.
    *
    *  \details    For a series of \e target, this method computes the gradient
    *              of :
    *              \f[
    * \sum_{i=1}^{NumberOfSubjects}\sum_{k=1}^{NumberOfObjects}\frac{1}{2\sigma_k^2}\left\Vert\phi^{\alpha_i}(X_{0,k}) - S_{i,k}\right\Vert_{W}^2.
    *              \f]
    *   with respect to \f$\alpha_i\f$, \f$X_{0,k}\f$ and control points positions.
    *
    *  \param[in]   Momentas        The momentas parameterizing a series of deformations \f$\phi^{\alpha_i}\f$.
    *  \param[in]   target          The series of targets \f$S_i\f$.
    *  \param[out]  gradPos         The gradient of the residual w.r.t. the control points (summing contributions of each subject).
    *  \param[out]  gradMom         The gradient of the residual w.r.t. the series of momentas.
    *  \param[out]  gradTempL_L2    The gradient of the residual w.r.t. the template (summing contributions of each subject),
    *                               using the \f$L^2\f$ metric for template objects of landmark types (and children types).
    *  \param[out]  gradTempL_Sob   The gradient of the residual w.r.t. the template (summing contributions of each subject),
    *                               using a Sobolev metric for template objects of landmark types (and children types). It is
    *                               the convolution of \e gradTempL_L2 with a smoothing kernel, which guarantees that template meshes
    *                               do not self-intersect during optimization.
    */
	void ComputeDataTermGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target,
			MatrixType& gradPos, MatrixList& gradMom, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob);
	/// Compute the data term gradient and return the Sobolev gradient of the template objects if \e do_Sobolev_Template_Gradient and the \f$L^2\f$ gradient otherwise
	void ComputeDataTermGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target,
			MatrixType& gradPos, MatrixList& gradMom, MatrixList& gradTempL, bool do_Sobolev_Template_Gradient);

	/// Computes the Sobolev gradient of template objects of landmark type (or children types) from the \f$L^2\f$ gradient
	MatrixList ConvolveGradTemplate(MatrixList& GradTemplate_L2);
	
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

	/// Type of atlas.
	AtlasType m_Type;

	/// Template as a DeformableMultiObject.
	DeformableMultiObjectType* m_Template;
	
	/// Working deformation to compute atlas deformation.
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
	
	/// Inverse of Covariance Matrix of the momentas.
	MatrixType m_CovMomInverse;
	
	/// Kernel width to compute the Sobolev gradient of the functional/likelihood w.r.t. template variable.
	TScalar m_SmoothingKernelWidth;
	
	/// Number of threads.
	unsigned int m_NumberOfThreads;
	


protected :

	/// \cond HIDE_FOR_DOXYGEN

	static ITK_THREAD_RETURN_TYPE _residualsThread(void* arg);

	static ITK_THREAD_RETURN_TYPE _gradientResidualsThread(void* arg);

	itk::SimpleFastMutexLock m_Mutex;

	unsigned int m_MT_SubjectCounter;
	unsigned int m_MT_NumberOfSubjects;

	bool m_MT_OutOfBox;

	MatrixList m_MT_Momentas;
	std::vector< DeformableMultiObjectType*> m_MT_Target;
	std::vector< std::vector< TScalar > > m_MT_Residuals;

	MatrixType m_MT_GradPos;
	MatrixList m_MT_GradMom;
	MatrixList m_MT_GradTempL_L2;

	/// \endcond


}; /* class Atlas */


#ifndef MU_MANUAL_INSTANTIATION
#include "Atlas.txx"
#endif


#endif /* _Atlas_h */

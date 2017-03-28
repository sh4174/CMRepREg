/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Diffeos_h
#define _Diffeos_h

#include "AbstractDeformations.h"

#include "itkImage.h"

/**
 *  \brief      Standard diffeomorphisms.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The Diffeos class inherited from AbstractDeformations represents the standard
 *              deformation which is usually employed in minimization problems such as registration
 *              or atlas construction. Deformations are encoded by control points and momentum vectors attached to them,
 *              which move in time according a to an Hamiltonian set of equations.\n \n
 */
template <class TScalar, unsigned int Dimension>
class Diffeos : public AbstractDeformations<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Abstract deformation type.
	typedef AbstractDeformations<TScalar, Dimension> Superclass;

	/// Vector type.
	typedef typename Superclass::VectorType VectorType;
	/// Matrix type.
	typedef typename Superclass::MatrixType MatrixType;
	/// List of matrices type.
	typedef typename Superclass::MatrixList MatrixList;

	/// Deformable multi-object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

	/// ITK image type.
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Diffeos();

	/// Copy constructor.
	Diffeos(const Diffeos& other);

	virtual Diffeos* Clone() { return new Diffeos(*this); };

	virtual ~Diffeos();
	
	// // Copy essential information, but not deformable multi-object, control points and momentas
	// virtual void CopyInformation(const Diffeos& other);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Return initial time.
	inline float GetT0() const { return m_T0; }
	/// Sets initial time to \e t0.
	inline void SetT0(float t0) { m_T0 = t0; this->SetModified(); }

	/// Return final time.
	inline float GetTN() const { return m_TN; }
	/// Sets final time to \e tn.
	inline void SetTN(float tn) { m_TN = tn; this->SetModified(); }

	/// Return the number of time points between \f$t_0\f$ and \f$t_n\f$.
	inline unsigned int GetNumberOfTimePoints() const { return m_NumberOfTimePoints; }
	/// Set the number of time points between \f$t_0\f$ and \f$t_n\f$.
	inline void SetNumberOfTimePoints(unsigned int n) { m_NumberOfTimePoints = n; this->SetModified();}

	///	Return the type of the kernel.
	inline KernelEnumType GetKernelType() const { return m_KernelType; }
	/// Set the type of the kernel to \e kernelType.
	inline void SetKernelType(KernelEnumType kernelType) { m_KernelType = kernelType; this->SetModified(); }

	///	Return the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Set the size of the kernel to \e kernelWidth.
	inline void SetKernelWidth(TScalar kernelWidth) { m_KernelWidth = kernelWidth; this->SetModified(); }

	/// Set the initial control points to \e X.
	inline void SetStartPositions(const MatrixType& X) { m_StartPositions = X; this->SetModified(); }

	/// Set the initial momenta to \e A.
	inline void SetStartMomentas(const MatrixType& A) { m_StartMomentas = A; this->SetModified(); }
	
	
	/// Get trajectories of control points
	inline MatrixList GetTrajectoryPositions() const { return m_PositionsT; }
	/// Get values of momentum vectors in time
	inline MatrixList GetTrajectoryMomentas() const { return m_MomentasT; }
	
	/// Returns the adjoint variable of CP positions (computed by solving adjoint equations).
	inline MatrixType GetAdjointPosAt0() const { return m_AdjointPosAt0; };
	/// Returns the adjoint variable of initial momenta (computed by solving adjoint equations).
	inline MatrixType GetAdjointMomAt0() const { return m_AdjointMomAt0; };
	/// Returns the adjoint variable of landmark points (computed by solving adjoint equations).
	inline MatrixType GetAdjointLandmarkPointsAt0() const { return m_AdjointLandmarkPointsAt0; };
	
	/// Return true if use improved Euler's method, false otherwise.
	inline bool ImprovedEuler() const { return m_UseImprovedEuler; }
	/// Set standard Euler's method.
	inline void UseStandardEuler() { m_UseImprovedEuler = false; this->SetModified(); }
	/// Set improved Euler's method.
	inline void UseImprovedEuler() { m_UseImprovedEuler = true; this->SetModified(); }

	/// Return the data domain.
	inline MatrixType GetDataDomain() const { return m_DataDomain; }
	/// Set the data domain to \e domain.
	inline void SetDataDomain(MatrixType& domain) { m_DataDomain = domain; }

	/// Return true if any point is out of the bounding box, false otherwise.
	inline bool OutOfBox() const { return m_OutOfBox; }

	/// Return the padding factor.
	inline TScalar GetPaddingFactor() const { return m_PaddingFactor; }
	/// Set the padding factor to \e paddingFactor.
	inline void SetPaddingFactor(TScalar paddingFactor) { m_PaddingFactor = paddingFactor; }

	/// Return true if the true inverse flow is used, false otherwise.
	inline bool ComputeTrueInverseFlow() const { return m_ComputeTrueInverseFlow; }
	/// Set true inverse flow (see Diffeos::m_ComputeTrueInverseFlow for details).
	inline void SetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = true; }
	/// Set direct flow to compute inverse deformation (see Diffeos::m_ComputeTrueInverseFlow for details).
	inline void UnsetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = false; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void Update();

	/// Computes the linear adjoint ODE of the flow equation (See AdjointEquationsIntegrator class for details).
	void IntegrateAdjointEquations(MatrixType& InitialConditionsLandmarkPoints, MatrixType& InitialConditionsImagePoints);
	void IntegrateAdjointEquations(MatrixList& InitialConditionsLandmarkPoints, MatrixList& InitialConditionsImagePoints, std::vector<unsigned int> jumpTimes);

	/**
	 *	\brief		Checks if a set of points is out of box or not.
	 *
	 *	\details	This function enables to see if the coordinates of X at time t is out of box or not.
	 *
	 *	\param[in]	X	List of matrices containing coordinates at different time steps.
	 *	\param[in]	t	Time index for X[t].
	 *	\return		True if X[t] is out of box, false otherwise.
	 */
	bool CheckBoundingBox(MatrixList& X, int t);

// /*
//  *	\brief		Implements the linear adjoint ODE of the flow equation
//  *
//  *	\details	The flow equations writes: \f$\dot X(t) = G(X(t),S(t))\f$. This function computes \f$\dot\theta(t) = \partial_1 G(X(t),S(t))^T\theta(t) \f$ with \f$ \theta(1) \f$ the gradient of the fidelity term.
//  *
//  *	\param[in]  InitialConditions \f$ \theta(1) \f$	
//  *	\param[out]	PointsT	Trajectories of data and image points \f$ X(t) \f$
//  *	\param[out]	VectorsT Trajectories of auxiliary variable \f$ \theta(t) \f$
//  *	\return		True if X[t] is out of box, false otherwise.
//  */
// 	// void TransportAlongGeodesic(MatrixList& InitialConditions, MatrixList& VectorsT, MatrixList& PointsT);

	virtual DeformableMultiObjectType* GetDeformedObject()
		{ return this->GetDeformedObjectAt( (this->GetNumberOfTimePoints()-1) ); }

	/// Returns the deformed objects at time \e t.
	DeformableMultiObjectType* GetDeformedObjectAt(unsigned int t);

	virtual void WriteFlow(const std::vector<std::string>& name, const std::vector<std::string>& extension);
	
	/// splat the residual (this deformed image - target image) defined on the final image map
	MatrixType SplatResidualImage(const DeformableMultiObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void FlowImagePointsTrajectory();
	virtual void FlowLandmarkPointsTrajectory();

	/**
	 *	\brief		Initializes the bounding box.
	 *
	 *	\details	This method initialize the Diffeos::m_BoundingBox attribute according to
	 *				the data domain, the padding factor and the size of the kernel.
	 */
	void InitBoundingBox();

	/// Solves the Hamiltonian system associated to the initial positions and momenta.
	void Shoot();



private :

	/// Compute voxels trajectories using the direct flow integrated backward with speed flipped (i.e. \f$\phi_t\circ\phi_1^{-1}\f$).
	void IntegrateImagePointsBackward();
	/// Compute voxels trajectories using the flow of inverse deformations \f$\phi_t^{-1}\f$.
	void IntegrateImagePointsWithTrueInverseFlow();


protected :

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Initial time \f$ t_0 \f$.
	TScalar m_T0;
	///	Final time \f$ t_n \f$.
	TScalar m_TN;
	///	Number of time points between \f$ t_0 \f$ and \f$ t_n \f$.
	unsigned int m_NumberOfTimePoints;

	///	Matrix containing the initial control points (Size : N x Dimension).
	MatrixType m_StartPositions;
	///	Matrix containing the initial momenta (Size : N x Dimension).
	MatrixType m_StartMomentas;

	///	List containing the position of the control points at different time points.
	MatrixList m_PositionsT;
	///	List containing the position of the momenta at different time points.
	MatrixList m_MomentasT;
	

	/// Type of the kernel.
	KernelEnumType m_KernelType;
	///	Size of the kernel associated to the deformation.
	TScalar m_KernelWidth;

	/// Boolean which indicates if we use improved Euler's method or not.
	bool m_UseImprovedEuler;

	/************************************ BEGIN CHANGE **********************************************/
	// MatrixType m_CPDomain;
	/************************************** END CHANGE ************************************************/

	///	Matrix containing the min (resp. the max) of the initial positions at first (resp. second) line.
	MatrixType m_DataDomain;
	///	Box where any trajectory must not exit (Size : 2 x Dimension with the min (resp. the max) at first (resp. second) line).
	MatrixType m_BoundingBox;
	///	Multiplier coefficient for the creation of the bounding box ("boundingBox = dataDomain +/- paddingFactor*kernelWidth/2").
	TScalar m_PaddingFactor;
	///	Boolean which prevents computations if any point (e.g. the coordinates of trajectory) is outside the bounding box.
	bool m_OutOfBox;

	/// This parameter is used if there is an image.
	/// If set, true inverse flow will be used (i.e. \f$\phi_t^{-1}\f$).
	/// If not, direct flow will be integrated backward (speed flipped) to compute inverse deformation (i.e. \f$\phi_t\circ\phi_1^{-1}\f$). 
	/// It moves particles backward from time t=1 back to time t=0 along the same trajectory as in the forward flow.
	/// In this case, results are different from inverse deformation if \f$t\f$ is different from \f$0\f$ and \f$1\f$.
	/// \warning For regression, always use true inverse flow.
	bool m_ComputeTrueInverseFlow;

	/// Flow of voxel positions (with backward integration i.e. \f$\phi_t^{-1}\f$)
	MatrixList m_MapsT;
	/// Flow of voxel positions (with true inverse flow i.e. \f$\phi_t\circ\phi_1^{-1}\f$).
	MatrixList m_InverseMapsT;

	/// Velocity of landmark points
	MatrixList m_LandmarkPointsVelocity;

	///	Trajectory of the whole vertices of Landmark type at different time steps.
	MatrixList m_LandmarkPointsT;

	/// Adjoint variable of control points at time 0 (computed by solving adjoint equations)
	MatrixType m_AdjointPosAt0;

	/// Adjoint variable of momentum vectors at time 0  (computed by solving adjoint equations)	
	MatrixType m_AdjointMomAt0;
	
	/// Adjoint variable of landmark points at time 0  (computed by solving adjoint equations)
	MatrixType m_AdjointLandmarkPointsAt0;

}; /* class Diffeos */


#ifndef MU_MANUAL_INSTANTIATION
#include "Diffeos.txx"
#endif


#endif /* _Diffeos_h */

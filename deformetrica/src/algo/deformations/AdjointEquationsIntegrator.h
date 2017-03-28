#ifndef _AdjointEquationsIntegrator_h
#define _AdjointEquationsIntegrator_h

#include "itkImage.h"
#include "itkVector.h"
#include "Diffeos.h"

#include "KernelFactory.h"

#include "LinearAlgebra.h"
#include <vector>

/**
 *  \brief      Computes adjoint equations.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The AdjointEquationsIntegrator class computes the adjoint equations of the geodesic shooting equations and
 *              flow equations, which moves the gradient of the data term back to time t = 0. Its value at time t=0 is used
 *              to update initial momenta, initial position of control points, and possibly vertices of template shapes.
 */
template <class TScalar, unsigned int Dimension>
class AdjointEquationsIntegrator
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

	/// Image type (itk).
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// Image pointer type (itk).
	typedef typename ImageType::Pointer ImageTypePointer;

	/// Vector type (itk).
	typedef itk::Vector<TScalar, Dimension> itkVectorType;
	/// Vector image type (itk).
	typedef itk::Image<itkVectorType, Dimension> VectorImageType;
	/// Vector image pointer type (itk).
	typedef typename VectorImageType::Pointer VectorImagePointer;

	/// Kernel factory type.
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	/// Exact kernel type.
	typedef typename KernelFactoryType::KernelBaseType KernelType;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	AdjointEquationsIntegrator();

	~AdjointEquationsIntegrator();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets initial time to \e t0
	inline void SetT0(float t0) { m_T0 = t0; }

	/// Sets final time to \e tn.
	inline void SetTN(float tn) { m_TN = tn; }

	inline void SetJumpTimes(std::vector<unsigned int> jumps) { m_JumpTimes = jumps; }

	/// Sets the number of time points between \f$t_0\f$ and \f$t_n\f$.
	inline void SetNumberOfTimePoints(unsigned int n) { m_NumberOfTimePoints = n; }

	/// Sets the type of the kernel to \e kernelType.
	inline void SetKernelType(KernelEnumType kernelType) { m_KernelType = kernelType; }

	/// Sets kernel width
	inline void SetKernelWidth (TScalar h) { m_KernelWidth = h; }
	
	/// Sets downsampled working image
	inline void SetFullResolutionImage(ImageTypePointer img) { m_FullResolutionImage = img; }

	/// Sets trajectories of control points
	inline void SetControlPointsTrajectory(MatrixList& PtsT) { m_PosT = PtsT; }
	
	/// Sets values of the momentum vectors in time
	inline void SetMomentaTrajectory(MatrixList& PtsT) { m_MomT = PtsT; }
	
	/// Sets trajectories of landmark points
	inline void SetLandmarkPointsTrajectory(MatrixList& PtsT) { m_LandmarkPointsT = PtsT; m_IsLandmarkPoints = true; }
	
	/// Sets trajectories of voxels in the image domain
	inline void SetImagePointsTrajectory(MatrixList& PtsT) { m_ImagePointsT = PtsT; m_IsImagePoints = true; }
	
	/// Sets initial conditions of adjoint equations located at final position of landmark points
	/// (initial conditions are actually final condition, i.e. at time t=1)
	inline void SetInitialConditionsLandmarkPoints(MatrixType& M)
	{ 
		m_ListInitialConditionsLandmarkPoints.resize(1); 
		m_ListInitialConditionsLandmarkPoints[0] =  M; 
		m_HasJumps = false;
	}
	inline void SetInitialConditionsLandmarkPoints(MatrixList& M) 
	{ 
		m_ListInitialConditionsLandmarkPoints = M; 
		if (m_ListInitialConditionsLandmarkPoints.size() > 1)
		{	
			m_HasJumps = true; 
		}
		else
		{
			m_HasJumps = false;
		}
	}
	
	/// Sets initial conditions of adjoint equations located at final position of voxels
	/// (initial conditions are actually final conditions, i.e. at time t=1, for TrueInverseFlow,
	/// otherwise they are initial condition at time t=0)
	void SetInitialConditionsImagePoints(MatrixType& M) 
	{ 
		m_ListInitialConditionsImagePoints.resize(1); 
		m_ListInitialConditionsImagePoints[0] =  M;
		m_HasJumps = false;
	}
	inline void SetInitialConditionsImagePoints(MatrixList& M) 
	{ 
		m_ListInitialConditionsImagePoints = M;
		if (m_ListInitialConditionsImagePoints.size() > 1)
		{	
			m_HasJumps = true; 
		}
		else
		{
			m_HasJumps = false;
		}
	}
	
	/// Sets standard Euler's method
	inline void UseStandardEuler() { m_UseImprovedEuler = false; }
	/// Sets improved Euler's method
	inline void UseImprovedEuler() { m_UseImprovedEuler = true; }
	
	/// Sets true inverse flow (see Diffeos::m_ComputeTrueInverseFlow for details).
	inline void SetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = true; }
	/// Sets direct flow to compute inverse deformation (see Diffeos::m_ComputeTrueInverseFlow for details).
	inline void UnsetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = false; }


	/// Returns the value of the adjoint variable of control points at time t = \e q.
	MatrixType GetAdjointPosAt(long q) const { return m_XiPosT[q]; }

	/// Returns the value of the adjoint variable of momentum vectors at time t = \e q.
	MatrixType GetAdjointMomAt(long q) const { return m_XiMomT[q]; }
	
	/// Returns the value of the adjoint variable of landmark points at time t = \e q.
	MatrixType GetAdjointLandmarkPointsAt(long q) const { return m_ThetaT[q]; }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Compute system of adjoint equations
	void Update();


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 *  \brief        Computes the time derivative of the ajoint equation.
	 *
	 *  \details      Computes the time derivative of the ajoint equation. It is used during the integration scheme.
	 *
	 *  \param[in]    s      Time index.
	 *  \param[out]   dPos   Output derivative of the auxiliary variable of size dimension times the number of control points.
	 *  \param[out]   dMom   Output derivative of the auxiliary variable of size dimension times the number of momentas.
	 */
	void ComputeUpdateAt(unsigned int s, MatrixType& dPos,  MatrixType& dMom);
	
	/// TODO .
	void IntegrateAdjointOfLandmarkPointsEquations();
	/// TODO .
	void IntegrateAdjointOfImagePointsBackward();
	/// TODO .
	void IntegrateAdjointOfImagePointsForward();
	/// TODO .
	void IntegrateAdjointOfDiffeoParametersEquations();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Initial time \f$ t_0 \f$.
	TScalar m_T0;
	///	Final time \f$ t_n \f$.
	TScalar m_TN;
	///	Number of time points between \f$ t_0 \f$ and \f$ t_n \f$.
	unsigned int m_NumberOfTimePoints;

	std::vector<unsigned int> m_JumpTimes;
	
	/// Type of the kernel.
	KernelEnumType m_KernelType;
	///	Size of the kernel associated to the deformation.
	TScalar m_KernelWidth;
	
	/// Working image at full resolution (at which the gradient of the data term is computed)
	ImageTypePointer m_FullResolutionImage;
		
	/// Boolean which indicates if image points trajectories were computed using a true inverse flow or not.
	bool m_ComputeTrueInverseFlow;
	
	/// Boolean which indicates if we use Heun's integration method or not.
	bool m_UseImprovedEuler;
	
	/// true if Landmark points are set
	bool m_IsLandmarkPoints;
	/// true if Image points are set
	bool m_IsImagePoints;
	
	bool m_HasJumps;

	
	/// Trajectories of control points.
	MatrixList m_PosT;
	
	/// Values of momentum vectors in time.
	MatrixList m_MomT;
	
	/// Trajectories of the landmark points.
	MatrixList m_LandmarkPointsT;
	
	/// Trajectories of voxels in the image domain.
	MatrixList m_ImagePointsT;
	
	/// Initial conditions of adjoint equations located at final position of landmark points
	/// (initial conditions are actually final condition, i.e. at time t=1)
	MatrixList m_ListInitialConditionsLandmarkPoints;
	
	/// Initial conditions of adjoint equations located at final position of voxels
	/// (initial conditions are actually final conditions, i.e. at time t=1,
	/// for TrueInverseFlow, otherwise they are initial condition at time t=0)
	MatrixList m_ListInitialConditionsImagePoints;

	// propagated from time Tn to time T0
	/// Adjoint variable of control points.
	MatrixList m_XiPosT;
	/// Adjoint variable of momentum vectors.
	MatrixList m_XiMomT;

	/// Adjoint variable of landmark points.
	MatrixList m_ThetaT;
	/// Adjoint variable of image points.
	MatrixList m_EtaT;

	/// \cond HIDE_FOR_DOXYGEN

	KernelType* m_KernelObj1;
	KernelType* m_KernelObj2;
	KernelType* m_KernelObj3;
	KernelType* m_KernelObj4;

	/// \endcond


}; /* class AdjointEquationsIntegrator */


#ifndef MU_MANUAL_INSTANTIATION
#include "AdjointEquationsIntegrator.txx"
#endif


#endif /* _AdjointEquationsIntegrator_h */

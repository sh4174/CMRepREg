/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _Atlas_txx
#define _Atlas_txx

#include "Atlas.h"
#include "DeformableMultiObject.h"

#include "readMatrixDLM.txx"
#include "writeMatrixDLM.txx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Atlas<TScalar, Dimension>
::Atlas()
 {
	m_Type = null;
	m_CPSpacing = 0.0;
	m_DataSigmaSquared.set_size(0);
	m_CovMomInverse.set_size(0,0);
	m_NumberOfThreads = 1;
	m_SmoothingKernelWidth = 0.0;

	m_NumberOfObjects = 0;
	m_Def = NULL;
	m_Template = NULL;
 }



template <class TScalar, unsigned int Dimension>
Atlas<TScalar, Dimension>
::~Atlas()
{}



template <class TScalar, unsigned int Dimension>
Atlas<TScalar, Dimension>
::Atlas(const Atlas& other)
{
	m_Type = other.m_Type;

	m_Template = other.m_Template->Clone();
	m_ControlPoints = other.m_ControlPoints;
	m_DataSigmaSquared.set_size(other.m_DataSigmaSquared.size());
	for (int i = 0; i < other.m_DataSigmaSquared.size(); i++)
		m_DataSigmaSquared(i) = other.m_DataSigmaSquared(i);

	m_CovMomInverse = other.m_CovMomInverse;
	m_Def = other.m_Def->Clone();

	m_BoundingBox = other.m_BoundingBox;
	m_CPSpacing = other.m_CPSpacing;
	m_SmoothingKernelWidth = other.m_SmoothingKernelWidth;

	m_NumberOfThreads = other.m_NumberOfThreads;
	m_NumberOfObjects = other.m_NumberOfObjects;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::SetControlPoints(std::string& fn)
{
	if (strlen(fn.c_str()))
	{
		m_ControlPoints = readMatrixDLM<TScalar>(fn.c_str());
		std::cout << "Using a set of " << m_ControlPoints.rows() << " control points in file " << fn << std::endl;
	}
}

template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::SetCovarianceMomentaInverse(std::string& fn)
{
	if (strlen(fn.c_str()))
	{
		m_CovMomInverse = readMatrixDLM<TScalar>(fn.c_str());
		std::cout << "Using an inverse of the momenta covariance matrix of size " << m_CovMomInverse.rows() << " x " << m_CovMomInverse.columns() << " from file: " << fn << std::endl;
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::Update()
 {

	m_Template->Update();
	m_NumberOfObjects = m_Template->GetNumberOfObjects();
	m_BoundingBox = m_Template->GetBoundingBox();

	if (m_Def == NULL)
		throw std::runtime_error("A deformation should be set to the atlas model");

	if (m_ControlPoints.rows() == 0)
	{
		if (m_CPSpacing == 0.0)
		{
			m_CPSpacing = m_Def->GetKernelWidth();
			std::cout << "No initial CP spacing given: using diffeo kernel width of " << m_CPSpacing << std::endl;
		}

		this->InitializeControlPoints();
	}

	for (unsigned int i = 0; i < m_ControlPoints.rows(); i++)
	{
		for (unsigned int dim = 0; dim < Dimension; dim++)
		{
			m_BoundingBox(dim,0) = (m_BoundingBox(dim,0)<m_ControlPoints(i,dim)?m_BoundingBox(dim,0):m_ControlPoints(i,dim));
			m_BoundingBox(dim,1) = (m_BoundingBox(dim,1)>m_ControlPoints(i,dim)?m_BoundingBox(dim,1):m_ControlPoints(i,dim));
		}
	}
		 
 }
	

template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::WriteTemplateData(const std::string& AtlasName,
		const std::vector<std::string>& TemplateObjectsName, const std::vector<std::string>& TemplateObjectsNameExtension)
{
	// write template objects
	std::vector<std::string> FullNames(TemplateObjectsName.size());
	for (int i = 0; i < TemplateObjectsName.size(); i++)
	{
		std::ostringstream oss;
		oss << AtlasName << "_template_" << TemplateObjectsName[i] << TemplateObjectsNameExtension[i] << std::ends;
		FullNames[i] = oss.str();
	}

	m_Template->WriteMultiObject(FullNames);
}


template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::WriteAtlasParameters(const std::string& AtlasName)
{
	// write control points
	std::ostringstream oss;
	oss << AtlasName << "_ControlPoints.txt" << std::ends;
	writeMatrixDLM<TScalar>(oss.str().c_str(), m_ControlPoints);
}


template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::WriteAtlasDeformation(const MatrixType& Momenta, const std::string& AtlasName,
		const std::vector<std::string>& TemplateObjectsName, const std::vector<std::string>& TemplateObjectsNameExtension)
 {
	if (Momenta.rows() != m_ControlPoints.rows())
		throw std::runtime_error("Number of Momentas and Control Points mismatch");


	m_Def->SetDeformableMultiObject(m_Template);
	m_Def->SetStartPositions(m_ControlPoints);
	m_Def->SetStartMomentas(Momenta);
	m_Def->Update();

	if (m_Def->OutOfBox())
	{
		throw std::runtime_error("Deformation out of box in Atlas::WriteAtlasToSubjectDeformation: this should not be");
	}

	std::vector< std::string > names(TemplateObjectsName.size());
	for (unsigned int i = 0; i < TemplateObjectsName.size(); i++)
	{
		std::ostringstream oss;
		oss << AtlasName << "_template_" << TemplateObjectsName[i];
		names[i] = oss.str();
	}

	m_Def->WriteFlow(names, TemplateObjectsNameExtension);

	// Added by ABGF@11aout14
	// output deformation field
	DeformationFieldIO<TScalar, Dimension> defFieldIO;
	defFieldIO.SetDiffeos(m_Def);
	defFieldIO.SetAnatomicalCoordinateSystemLabel("LPS");
	defFieldIO.WriteDeformationField(AtlasName,true);
	// Added by ABGF@11aout14

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::InitializeControlPoints()
 {
	VectorType Xmin = m_Template->GetBoundingBox().get_column(0);
	VectorType Xmax = m_Template->GetBoundingBox().get_column(1);

	//std::cout << "Xmin =\n " << Xmin << std::endl; // OK
	//std::cout << "Xmax =\n " << Xmax << std::endl; // OK

	std::vector<VectorType> pointList;
	VectorType v(Dimension);
	switch (Dimension)
	{
	case 2:
	{
		TScalar offsetX = 0.5*( Xmax[0] - Xmin[0] - m_CPSpacing * floor( ( Xmax[0] - Xmin[0] ) / m_CPSpacing ) );
		TScalar offsetY = 0.5*( Xmax[1] - Xmin[1] - m_CPSpacing * floor( ( Xmax[1] - Xmin[1] ) / m_CPSpacing ) );

		for (TScalar x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing)
			for (TScalar y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing)
			{
				v[0] = x;
				v[1] = y;
				pointList.push_back(v);
			}
		break;
	}
	case 3:
	{
		TScalar offsetX = 0.5*( Xmax[0] - Xmin[0] - m_CPSpacing * floor( ( Xmax[0] - Xmin[0] ) / m_CPSpacing ) );
		TScalar offsetY = 0.5*( Xmax[1] - Xmin[1] - m_CPSpacing * floor( ( Xmax[1] - Xmin[1] ) / m_CPSpacing ) );
		TScalar offsetZ = 0.5*( Xmax[2] - Xmin[2] - m_CPSpacing * floor( ( Xmax[2] - Xmin[2] ) / m_CPSpacing ) );

		for (TScalar x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing)
			for (TScalar y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing)
				for (TScalar z = Xmin[2] + offsetZ; z <= Xmax[2]; z += m_CPSpacing)
				{
					v[0] = x;
					v[1] = y;
					v[2] = z;
					pointList.push_back(v);
				}
		break;
	}
	default:
		throw std::runtime_error("GenerateInitialCP not implemented in Dimensions other than 2 and 3");
		break;
	}

	unsigned int NumCPs = pointList.size();

	m_ControlPoints.set_size(NumCPs, Dimension);
	for (unsigned int i = 0; i < NumCPs; i++)
		m_ControlPoints.set_row(i, pointList[i]);

	std::cout << "Set of " << NumCPs << " control points defined" << std::endl;

	// DEBUG
	// std::cout << "Control points =\n " << m_ControlPoints << std::endl; // OK

 }


template <class TScalar, unsigned int Dimension>
bool
Atlas<TScalar, Dimension>
::ComputeResidualsSubject(const MatrixType& Momentas, const DeformableMultiObjectType* target,
		std::vector< TScalar >& Residuals)
{

	if (Momentas.rows() != m_ControlPoints.rows())
		throw std::runtime_error("Number of Momentas and Control Points mismatch");

	// A copy of the deformation is needed since this method can be called by different threads at the same time
	DiffeosType* subjectDef = m_Def->Clone();

	subjectDef->SetDeformableMultiObject(m_Template);
	subjectDef->SetStartPositions(m_ControlPoints);
	subjectDef->SetStartMomentas(Momentas);
	subjectDef->Update();

	if (subjectDef->OutOfBox())
	{
		return true;
	}

	Residuals = subjectDef->GetDeformedObject()->ComputeMatch(target);

    // no normalization in the residuals
//	for (int i = 0; i < m_NumberOfObjects; i++)
//		Residuals[i] /= ( 2.0 * m_DataSigmaSquared(i) );
    
    delete subjectDef;

	return false;
}



template <class TScalar, unsigned int Dimension>
bool
Atlas<TScalar, Dimension>
::ComputeResiduals(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType*> target,
		std::vector< std::vector< TScalar > >& Residuals)
 {
	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of subjects in Momentas and target lists mismatch");

	int numberOfSubjects = Momentas.size();

	//
	// Version without multi-threading :
	//
	if (m_NumberOfThreads < 2)
	{
		Residuals.resize(numberOfSubjects);

		bool oob = false;
		for (int s = 0; s < numberOfSubjects; s++)
		{
			oob = this->ComputeResidualsSubject(Momentas[s], target[s], Residuals[s]);

			if (oob)
			{
				return true;
			}
		}

		return false;
	}

	//
	// Multi-threaded version :
	//
	m_MT_SubjectCounter = 0;
	m_MT_NumberOfSubjects = numberOfSubjects;

	m_MT_Momentas = Momentas;

	if ( m_MT_Target != target )
		m_MT_Target = target;

	m_MT_OutOfBox = false;

	m_MT_Residuals.resize(numberOfSubjects);
	for (int s = 0; s < numberOfSubjects; s++)
		m_MT_Residuals[s].resize(m_NumberOfObjects);

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_MT_NumberOfSubjects)
		numThreads = m_MT_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&Atlas::_residualsThread, (void*)this);
	threader->SingleMethodExecute();

//	Residuals = m_MT_Residuals;


	Residuals.resize(numberOfSubjects);
	for (int s = 0; s < numberOfSubjects; s++)
		Residuals[s] = m_MT_Residuals[s];

	return m_MT_OutOfBox;
 }



template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::ComputeDataTermGradientSubject(const MatrixType& Momentas, const DeformableMultiObjectType* target,
		MatrixType& dPos, MatrixType& dMom, MatrixList& dTempL)
 {

	if (Momentas.rows() != m_ControlPoints.rows())
		throw std::runtime_error("Number of Momentas and Control Points mismatch");

	// A copy of the deformation is needed since this method can be called by different threads at the same time
	DiffeosType* subjectDef = m_Def->Clone();
	subjectDef->SetDeformableMultiObject(m_Template);
	subjectDef->SetStartPositions(m_ControlPoints);
	subjectDef->SetStartMomentas(Momentas);
	subjectDef->Update();

	if (subjectDef->OutOfBox())
		throw std::runtime_error("Out of box in Atlas::ComputeResidualGradient (this should not be...)");


	/// Get the deformed template
	DeformableMultiObjectType* deformedTemplateObjects = subjectDef->GetDeformedObject();

	/// Get the gradient of the similarity metric between deformed template and target
	MatrixList GradientSimilarityMetric = deformedTemplateObjects->ComputeMatchGradient(target);

	/// Divide each gradient of the data term by 1/(2*DataSigmaSquared)
	for (int i = 0; i < m_NumberOfObjects; i++)
		GradientSimilarityMetric[i] /= ( 2.0 * m_DataSigmaSquared(i) );
	
	/// pool the gradients into big matrices for all objects of type landmarks (and children types) and image (and children types)
	MatrixType GradientDataTermOfLandmarkTypes;
	MatrixType GradientDataTermOfImageTypes;
	// subjectTemplate->GetTemplateObjects()->ListToMatrices(GradientSimilarityMetric, GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes);
	m_Template->ListToMatrices(GradientSimilarityMetric, GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes);

	//cout<<GradientDataTermOfLandmarkTypes<<"\n";

	/// Integrate the adjoint equations
	subjectDef->IntegrateAdjointEquations(GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes);

	/// Get the gradient w.r.t. deformation parameters and landmark points positions
	dPos = subjectDef->GetAdjointPosAt0();
	dMom = subjectDef->GetAdjointMomAt0();

	//cout<<"gradPos "<<dPos<<"\n";
	//cout<<"gradMom "<<dMom<<"\n";	

	MatrixType dTempLLandmarkTypes = subjectDef->GetAdjointLandmarkPointsAt0();

	/// Compute the gradient w.r.t. the intensities of the template image
	MatrixType dTempLImageTypes = subjectDef->SplatResidualImage(target);

	if ( m_Template->GetNumberOfImageKindObjects() )
		dTempLImageTypes /= m_DataSigmaSquared( m_Template->GetImageIndex() );

	/// turn matrices into list
	m_Template->MatricesToList(dTempL, dTempLLandmarkTypes, dTempLImageTypes);
     
     
     delete subjectDef;

 }



template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::ComputeDataTermGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target,
		MatrixType& gradPos, MatrixList& gradMom, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob)
 {

	if (Momentas.size() != target.size())
		throw std::runtime_error("Number of subjects in Momentas and target lists mismatch");

	int numberOfSubjects = Momentas.size();

	gradPos.set_size(m_ControlPoints.rows(), Dimension);
	gradPos.fill(0.0);
	gradMom.resize(numberOfSubjects);
	gradTempL_L2 = this->GetTemplateData(); // to get the correct matrices size
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradTempL_L2[i].fill(0.0);

	//cout<<"NUM SUBJECTS = "<<numberOfSubjects<<"\n";

	//
	// Version without multi-threading :
	//
	if (m_NumberOfThreads < 2)
	{
		for (int s = 0; s < numberOfSubjects; s++)
		{
			MatrixType dPos;
			MatrixType dMom;
			MatrixList dTempL(m_NumberOfObjects);
			this->ComputeDataTermGradientSubject(Momentas[s], target[s], dPos, dMom, dTempL);

			gradPos += dPos;
		// dMom /= m_NumberOfSubjects;
			gradMom[s] = dMom;

			for (int i = 0; i < m_NumberOfObjects; i++)
			{
			// gradTempL_L2[i] += (dTempL[i] / m_NumberOfSubjects);
				gradTempL_L2[i] += dTempL[i];
			}
		}

	// m_GradPos /= m_NumberOfSubjects;
		gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);

		return;
	}



	//
	// Multi-threading :
	//

	m_MT_SubjectCounter = 0;
	m_MT_NumberOfSubjects = numberOfSubjects;

	m_MT_Momentas = Momentas;
	if ( m_MT_Target != target )
		m_MT_Target = target;


	m_MT_GradPos = gradPos;
	m_MT_GradMom = gradMom;
	m_MT_GradTempL_L2 = gradTempL_L2;

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_MT_NumberOfSubjects)
		numThreads = m_MT_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&Atlas::_gradientResidualsThread, (void*)this);
	threader->SingleMethodExecute();

	// gradPos /= numberOfSubjects;


	gradPos = m_MT_GradPos;
	gradMom = m_MT_GradMom;
	gradTempL_L2 = m_MT_GradTempL_L2;
	
	
	

	gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);

 }



template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::ComputeDataTermGradient(const MatrixList& Momentas, const std::vector< DeformableMultiObjectType* > target,
		MatrixType& gradPos, MatrixList& gradMom, MatrixList& gradTempL, bool do_Sobolev_Template_Gradient)
 {
	 
	MatrixList Aux;
	if (do_Sobolev_Template_Gradient)
		ComputeDataTermGradient(Momentas, target, gradPos, gradMom, Aux, gradTempL);
	else
		ComputeDataTermGradient(Momentas, target, gradPos, gradMom, gradTempL, Aux);

	return;
 }



template <class TScalar, unsigned int Dimension>
typename Atlas<TScalar, Dimension>::MatrixList
Atlas<TScalar, Dimension>
::ConvolveGradTemplate(MatrixList& GradTemplate_L2)
 {

	if ( m_SmoothingKernelWidth < 1e-20 )
	{
		return GradTemplate_L2;
	}

	MatrixList GradTemplate_Sob(m_NumberOfObjects);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject(m_Def->GetKernelType());
	momKernelObj->SetKernelWidth(m_SmoothingKernelWidth);

	for (unsigned int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
		{
			MatrixType TempData = this->GetTemplateData()[i];

			momKernelObj->SetSources(TempData);
			momKernelObj->SetWeights(GradTemplate_L2[i]);

			GradTemplate_Sob[i] = momKernelObj->Convolve(TempData);
		}
		else
			GradTemplate_Sob[i] = GradTemplate_L2[i];

	}

	delete momKernelObj;

	return GradTemplate_Sob;

 }



template <class TScalar, unsigned int Dimension>
void
Atlas<TScalar, Dimension>
::UpdateDeformationAndKernelDataDomain(const std::vector< DeformableMultiObjectType* > target)
 {
	// There are 3 boxes:
	// m_BoundingBox: tighlty enclosing template data and control points
	// a padded version with padding equal to 0.5*PaddingFactor*DeformationKernelWidth: points should never move outside this box, it used in Diffeos->CheckBoundingBox
	// a padded version with padding equal to     PaddingFactor*DeformationKernelWidth: used to define p3MKernel grid to account for periodic boundary conditions

	MatrixType DataDomain = m_BoundingBox;
	for (unsigned int s = 0; s < target.size(); s++)
	{
		MatrixType BB = target[s]->GetBoundingBox();
		for (int d = 0; d < Dimension; d++)
		{
			DataDomain(d,0) = ( DataDomain(d,0)<BB(d,0)?DataDomain(d,0):BB(d,0) );
			DataDomain(d,1) = ( DataDomain(d,1)>BB(d,1)?DataDomain(d,1):BB(d,1) );
		}
	}

	m_Def->SetDataDomain( DataDomain ); // will be used to define bounding box with padding factor of 0.5*PaddingFactor*DeformationKernelWidth

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetDataDomain( DataDomain );  // will be used to define bounding box with padding factor of    PaddingFactor*DeformationKernelWidth

 }



//
// Multi-threaded components
//

template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
Atlas<TScalar, Dimension>
::_residualsThread(void* arg)
 {
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	Atlas* obj = static_cast<Atlas*>(infoStruct->UserData);

	unsigned int s = obj->m_MT_NumberOfSubjects;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_MT_NumberOfSubjects)
			break;

		std::vector<TScalar> values;

		bool oob = obj->ComputeResidualsSubject(obj->m_MT_Momentas[s], obj->m_MT_Target[s], values);

		obj->m_Mutex.Lock();
		obj->m_MT_OutOfBox = oob;
		obj->m_Mutex.Unlock();

		if (oob)
			break;

		for (int i = 0; i < obj->m_NumberOfObjects; i++)
		{
			obj->m_Mutex.Lock();
			obj->m_MT_Residuals[s][i] = values[i];
			obj->m_Mutex.Unlock();
		}

	}

	return ITK_THREAD_RETURN_VALUE;
 }



template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
Atlas<TScalar, Dimension>
::_gradientResidualsThread(void* arg)
 {
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	Atlas* obj = static_cast<Atlas*>(
			infoStruct->UserData);

	unsigned int s = obj->m_MT_NumberOfSubjects;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_MT_NumberOfSubjects)
			break;

		MatrixType dPos;
		MatrixType dMom;
		MatrixList dTempL(obj->m_NumberOfObjects);

		obj->ComputeDataTermGradientSubject(obj->m_MT_Momentas[s], obj->m_MT_Target[s], dPos, dMom, dTempL);

		obj->m_Mutex.Lock();
		obj->m_MT_GradPos += dPos;
		obj->m_Mutex.Unlock();

		obj->m_Mutex.Lock();
		// dMom /= m_NumberOfSubjects;
		obj->m_MT_GradMom[s] = dMom;
		obj->m_Mutex.Unlock();

		for (int i = 0; i < obj->m_NumberOfObjects; i++)
		{
			obj->m_Mutex.Lock();
			// gradTempL_L2[i] += (dTempL[i] / m_NumberOfSubjects);
			obj->m_MT_GradTempL_L2[i] += dTempL[i];
			obj->m_Mutex.Unlock();
		}
	}

	return ITK_THREAD_RETURN_VALUE;

 }



#endif /* _Atlas_txx */

/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _Regression_txx
#define _Regression_txx

#include "Regression.h"
#include "DeformableMultiObject.h"

#include "readMatrixDLM.txx"
#include "writeMatrixDLM.txx"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Regression<TScalar, Dimension>
::Regression()
{
	m_CPSpacing = 0.0;
	m_DataSigmaSquared.set_size(0);
	m_SmoothingKernelWidth = 0.0;

	m_NumberOfObjects = 0;
	m_Def = NULL;
	m_Template = NULL;
}


template <class TScalar, unsigned int Dimension>
Regression<TScalar, Dimension>
::~Regression()
{
}


template <class TScalar, unsigned int Dimension>
Regression<TScalar, Dimension>
::Regression(const Regression& other)
{
	m_Template = other.m_Template->Clone();
	m_ControlPoints = other.m_ControlPoints;
	m_DataSigmaSquared.set_size(other.m_DataSigmaSquared.size());
	for (int i = 0; i < other.m_DataSigmaSquared.size(); i++)
		m_DataSigmaSquared(i) = other.m_DataSigmaSquared(i);

	m_Def = other.m_Def->Clone();

	m_BoundingBox = other.m_BoundingBox;
	m_CPSpacing = other.m_CPSpacing;
	m_SmoothingKernelWidth = other.m_SmoothingKernelWidth;

	m_NumberOfObjects = other.m_NumberOfObjects;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::SetControlPoints(std::string& fn)
{
	if (strlen(fn.c_str()))
	{
		m_ControlPoints = readMatrixDLM<TScalar>(fn.c_str());
		std::cout << "Using a set of " << m_ControlPoints.rows() << " control points in file " << fn << std::endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::Update()
{
	m_Template->Update();
	m_NumberOfObjects = m_Template->GetNumberOfObjects();
	m_BoundingBox = m_Template->GetBoundingBox();

	if (m_Def == NULL)
		throw std::runtime_error("A deformation should be set to the regression model");

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
typename Regression<TScalar, Dimension>::FunctionalValuesType*
Regression<TScalar, Dimension>
::ComputeFunctional(const MatrixType& momenta, const std::vector< DeformableMultiObjectType* > target, std::vector<unsigned int> timeIndices)
{
	this->UpdateDeformationAndKernelDataDomain(target);

	int numObservations = target.size();
	FunctionalValuesType* out = new FunctionalValuesType();
	std::vector< std::vector< TScalar > > residuals(numObservations);

	bool oob = this->ComputeResiduals(momenta, target, timeIndices, residuals);
	if (oob)
	{
		out->SetOutOfBox();
		return out;
	}

	// dataTerm = sum of residuals
	TScalar dataTerm = 0.0;
	for (unsigned int s = 0; s < numObservations; s++)
	{
		for (unsigned int i = 0; i < this->m_NumberOfObjects; i++)
		{
			dataTerm += residuals[s][i] / (2.0 * this->m_DataSigmaSquared[i]);
		}
	}

	// Regularity term
	TScalar regTerm = 0.0;
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

    	KernelType* momKernelObj  = kfac->CreateKernelObject(this->m_Def->GetKernelType());
	momKernelObj->SetKernelWidth(this->m_Def->GetKernelWidth());
	momKernelObj->SetSources(this->m_ControlPoints);
	momKernelObj->SetWeights(momenta);
		
	MatrixType kMom = momKernelObj->Convolve(this->m_ControlPoints);
	for (unsigned int i = 0; i < this->m_ControlPoints.rows(); i++)
		regTerm += dot_product(kMom.get_row(i), momenta.get_row(i));
	
	delete momKernelObj;
	regTerm *= 0.5;

	// Sparsity Prior
	TScalar sparsityTerm = 0.0;
	for (unsigned int i = 0; i < this->m_ControlPoints.rows(); i++)
		sparsityTerm += momenta.get_row(i).magnitude();
	sparsityTerm *= m_SparsityPrior;

	// Set functional values
	out->SetDataTerm(dataTerm);
	out->SetRegularityTerm(regTerm);
	out->SetSparsityTerm(sparsityTerm);
	out->SetResiduals(residuals);

	return out;
}
	

template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::ComputeFunctionalGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
			MatrixType& gradMom, MatrixType& gradPos, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob)
{
	this->UpdateDeformationAndKernelDataDomain(target);

	this->ComputeDataTermGradient(momenta, target, timeIndices, gradPos, gradMom, gradTempL_L2, gradTempL_Sob);

	this->AddGradientRegularityTerm(gradPos, gradMom, momenta);
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::ComputeFunctionalGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
			    MatrixType& gradMom, MatrixType& gradPos, MatrixList& gradTempL, bool do_SobolevGradient)
{
	this->UpdateDeformationAndKernelDataDomain(target);

	this->ComputeDataTermGradient(momenta, target, timeIndices, gradPos, gradMom, gradTempL, do_SobolevGradient);

	this->AddGradientRegularityTerm(gradPos, gradMom, momenta);
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::AddGradientRegularityTerm(MatrixType gradPos, MatrixType gradMom, const MatrixType& momenta)
{
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject(this->m_Def->GetKernelType());
	momKernelObj->SetKernelWidth(this->m_Def->GetKernelWidth());
	momKernelObj->SetSources(this->m_ControlPoints);

	momKernelObj->SetWeights(momenta);
	gradMom += momKernelObj->Convolve(this->m_ControlPoints);
	gradPos += momKernelObj->ConvolveGradient(this->m_ControlPoints, momenta);

	delete momKernelObj;
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::WriteTemplateFlow(const MatrixType& momenta, const std::string& regressionName,
		    const std::vector<std::string>& templateObjectsName, const std::vector<std::string>& templateObjectsNameExtension)
{
	m_Def->SetDeformableMultiObject(m_Template);
	m_Def->SetStartPositions(m_ControlPoints);
	m_Def->SetStartMomentas(momenta);
	m_Def->Update();

	if (m_Def->OutOfBox())
	{
		throw std::runtime_error("Deformation out of box in Regression::WriteTemplateFlow: this should not be");
	}

	std::vector< std::string > names(templateObjectsName.size());
	for (unsigned int i = 0; i < templateObjectsName.size(); i++)
	{
		std::ostringstream oss;
		oss << regressionName << "_baseline_" << templateObjectsName[i];
		names[i] = oss.str();
	}

	m_Def->WriteFlow(names, templateObjectsNameExtension);

	// Output deformation field
	DeformationFieldIO<TScalar, Dimension> defFieldIO;
	defFieldIO.SetDiffeos(m_Def);
	defFieldIO.SetAnatomicalCoordinateSystemLabel("LPS");
	defFieldIO.WriteDeformationField(regressionName, true);
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::WriteRegressionParameters(const MatrixType& momenta, const std::string& regressionName)
{
	// Write control points
	std::ostringstream oss1;
	oss1 << regressionName << "_ControlPoints.txt" << std::ends;
	writeMatrixDLM<TScalar>(oss1.str().c_str(), m_ControlPoints);

	// Write momenta
	std::ostringstream oss2;
	oss2 << regressionName << "_Momenta.txt" << std::ends;
	writeMatrixDLM<TScalar>(oss2.str().c_str(), momenta);

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::InitializeControlPoints()
 {
	VectorType xMin = m_Template->GetBoundingBox().get_column(0);
	VectorType xMax = m_Template->GetBoundingBox().get_column(1);

	std::vector<VectorType> pointList;
	VectorType v(Dimension);
	switch (Dimension)
	{
		case 2:
		{
			TScalar offsetX = 0.5*( xMax[0] - xMin[0] - m_CPSpacing * floor( ( xMax[0] - xMin[0] ) / m_CPSpacing ) );
			TScalar offsetY = 0.5*( xMax[1] - xMin[1] - m_CPSpacing * floor( ( xMax[1] - xMin[1] ) / m_CPSpacing ) );

			for (TScalar x = xMin[0] + offsetX; x <= xMax[0]; x += m_CPSpacing)
				for (TScalar y = xMin[1] + offsetY; y <= xMax[1]; y += m_CPSpacing)
				{
					v[0] = x;
					v[1] = y;
					pointList.push_back(v);
				}
			break;
		}
		case 3:
		{
			TScalar offsetX = 0.5*( xMax[0] - xMin[0] - m_CPSpacing * floor( ( xMax[0] - xMin[0] ) / m_CPSpacing ) );
			TScalar offsetY = 0.5*( xMax[1] - xMin[1] - m_CPSpacing * floor( ( xMax[1] - xMin[1] ) / m_CPSpacing ) );
			TScalar offsetZ = 0.5*( xMax[2] - xMin[2] - m_CPSpacing * floor( ( xMax[2] - xMin[2] ) / m_CPSpacing ) );

			for (TScalar x = xMin[0] + offsetX; x <= xMax[0]; x += m_CPSpacing)
				for (TScalar y = xMin[1] + offsetY; y <= xMax[1]; y += m_CPSpacing)
					for (TScalar z = xMin[2] + offsetZ; z <= xMax[2]; z += m_CPSpacing)
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

	unsigned int numCPs = pointList.size();

	m_ControlPoints.set_size(numCPs, Dimension);
	for (unsigned int i = 0; i < numCPs; i++)
		m_ControlPoints.set_row(i, pointList[i]);

	std::cout << "Set of " << numCPs << " control points defined" << std::endl;
}


template <class TScalar, unsigned int Dimension>
bool
Regression<TScalar, Dimension>
::ComputeResiduals(const MatrixType& momenta, const std::vector< DeformableMultiObjectType*> target, std::vector<unsigned int> timeIndices,
		   std::vector< std::vector< TScalar > >& residuals)
{
	if (target.size() != timeIndices.size())
	{
		throw std::runtime_error("Mismatch between number of observations and number of time indicies.");
	}

	bool oob = false;

	// A copy of the deformation is needed since this method can be called by different threads (maybe?) at the same time
	DiffeosType* subjectDef = m_Def->Clone();

	subjectDef->SetDeformableMultiObject(m_Template);
	subjectDef->SetStartPositions(m_ControlPoints);
	subjectDef->SetStartMomentas(momenta);
	subjectDef->Update();

	if (subjectDef->OutOfBox())
	{
		return true;
	}

	unsigned int numberOfObservations = target.size();
	residuals.resize(numberOfObservations);


	// Loop over the number of observations to compute the value of the data matching term
	for (int s = 0; s < numberOfObservations; s++)
	{
		// Compute the residual
		residuals[s] = subjectDef->GetDeformedObjectAt(timeIndices[s])->ComputeMatch(target[s]);
	}

	delete subjectDef;
	return false;
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::ComputeDataTermGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType* > target, std::vector<unsigned int> timeIndices,
			  MatrixType& gradPos, MatrixType& gradMom, MatrixList& gradTempL_L2, MatrixList& gradTempL_Sob)
{
	if (target.size() != timeIndices.size())
	{
		throw std::runtime_error("Mismatch between number of observations and number of time indicies.");
	}

	gradTempL_L2 = this->GetTemplateData(); // to get the correct matrices size
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradTempL_L2[i].fill(0.0);

	// A copy of the deformation is needed since this method can be called by different threads (maybe?) at the same time
	DiffeosType* subjectDef = m_Def->Clone(); // shoot control points and momenta and flow baseline

	//std::cout<<"	Shooting control points and momenta..."<<std::endl;

	subjectDef->SetDeformableMultiObject(m_Template);
	subjectDef->SetStartPositions(m_ControlPoints);
	subjectDef->SetStartMomentas(momenta);
	subjectDef->Update();

	//std::cout<<"	DONE Shooting control points and momenta"<<std::endl;

	if (subjectDef->OutOfBox())
		throw std::runtime_error("Out of box in Regression::ComputeDataTermGradient (this should not be...)");

	// Set up the data structure for holding the data gradients
	unsigned int numberOfObservations = target.size();
	std::vector< MatrixList > allGradXTi;
	allGradXTi.resize(numberOfObservations);
	for (int i=0; i<numberOfObservations; i++)
	{
		MatrixList matrList;
		matrList.resize(m_NumberOfObjects);
		allGradXTi[i] = matrList;

		for (int s=0; s<m_NumberOfObjects; s++)
		{

			DeformableObjectType* obj = m_Template->GetObjectList()[s];
			unsigned int numPoints_i = obj->GetNumberOfPoints();

			allGradXTi[i][s] = MatrixType(numPoints_i, Dimension, 0.0);
		}
	}

	//std::cout<<"	Computing gradient of data term..."<<std::endl;

	// Compute gradient of data matching term
	for (int s = 0; s < numberOfObservations; s++)
	{
		DeformableMultiObjectType* deformedTemplateObjects = subjectDef->GetDeformedObjectAt(timeIndices[s]);

		/// Get the gradient of the similarity metric between deformed template and target
		MatrixList gradDataTi = deformedTemplateObjects->ComputeMatchGradient(target[s]);

		//std::cout<<gradDataTi[0]<<std::endl;

		// We repackage the gradient matrices in a vector of size m_NumObjects, as the data classes expect
		DeformableObjectList templateList = m_Template->GetObjectList();
		for (int i=0; i<m_NumberOfObjects; i++)
		{
			MatrixType curGradTi = gradDataTi[i];
			
			TScalar thisIsPointlessButNecessary = 1 / (2.0 * this->m_DataSigmaSquared[i]);
			allGradXTi[s][i] = curGradTi * thisIsPointlessButNecessary;
			//allGradXTi[s][i] = curGradTi / (2.0 * this->m_DataSigmaSquared[i]);
		}
	}

	MatrixList GradientDataTermOfLandmarkTypes;
	MatrixList GradientDataTermOfImageTypes;
	GradientDataTermOfLandmarkTypes.resize(numberOfObservations);
	GradientDataTermOfImageTypes.resize(numberOfObservations);
	// Convert each matrix list into a matrix, giving us a matrix list at the end!
	for (int s = 0; s < numberOfObservations; s++)
	{

		MatrixType CurGradLandmark;
		MatrixType CurGradImage;
		m_Template->ListToMatrices(allGradXTi[s], CurGradLandmark, CurGradImage);

		GradientDataTermOfLandmarkTypes[s] = CurGradLandmark;
		GradientDataTermOfImageTypes[s] = CurGradImage;
	}

	//std::cout<<"	DONE computing gradient of data term"<<std::endl;

	//std::cout<<"	Integrating adjoint equations..."<<std::endl;
	
	subjectDef->IntegrateAdjointEquations(GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes, timeIndices);

	//std::cout<<"	DONE integrating adjoint equations..."<<std::endl;

	/// Get the gradient w.r.t. deformation parameters and landmark points positions
	gradPos = subjectDef->GetAdjointPosAt0();
	gradMom = subjectDef->GetAdjointMomAt0();

	MatrixType dTempLLandmarkTypes = subjectDef->GetAdjointLandmarkPointsAt0();
	/// Compute the gradient w.r.t. the intensities of the template image
	MatrixType dTempLImageTypes = subjectDef->SplatResidualImage(target[0]);

	/// Turn matrices into list
	m_Template->MatricesToList(gradTempL_L2, dTempLLandmarkTypes, dTempLImageTypes);

	gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::ComputeDataTermGradient(const MatrixType& momenta, const std::vector< DeformableMultiObjectType* > target, std::vector<unsigned int> timeIndices,
			  MatrixType& gradPos, MatrixType& gradMom, MatrixList& gradTempL, bool do_Sobolev_Template_Gradient)
{
	MatrixList aux;
	if (do_Sobolev_Template_Gradient)
		ComputeDataTermGradient(momenta, target, timeIndices, gradPos, gradMom, aux, gradTempL);
	else
		ComputeDataTermGradient(momenta, target, timeIndices, gradPos, gradMom, gradTempL, aux);

	return;
}


template <class TScalar, unsigned int Dimension>
typename Regression<TScalar, Dimension>::MatrixList
Regression<TScalar, Dimension>
::ConvolveGradTemplate(MatrixList& gradTemplate_L2)
{
	if ( m_SmoothingKernelWidth < 1e-20 )
	{
		return gradTemplate_L2;
	}

	MatrixList gradTemplate_Sob(m_NumberOfObjects);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject(m_Def->GetKernelType());
	momKernelObj->SetKernelWidth(m_SmoothingKernelWidth);

	for (unsigned int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_Template->GetObjectList()[i]->IsOfLandmarkKind() )
		{
			MatrixType tempData = this->GetTemplateData()[i];

			momKernelObj->SetSources(tempData);
			momKernelObj->SetWeights(gradTemplate_L2[i]);

			gradTemplate_Sob[i] = momKernelObj->Convolve(tempData);
		}
		else
			gradTemplate_Sob[i] = gradTemplate_L2[i];

	}

	delete momKernelObj;

	return gradTemplate_Sob;
}


template <class TScalar, unsigned int Dimension>
void
Regression<TScalar, Dimension>
::UpdateDeformationAndKernelDataDomain(const std::vector< DeformableMultiObjectType* > target)
{
	// There are 3 boxes:
	// m_BoundingBox: tighlty enclosing template data and control points
	// a padded version with padding equal to 0.5*PaddingFactor*DeformationKernelWidth: points should never move outside this box, it used in Diffeos->CheckBoundingBox
	// a padded version with padding equal to     PaddingFactor*DeformationKernelWidth: used to define p3MKernel grid to account for periodic boundary conditions

	MatrixType dataDomain = m_BoundingBox;
	for (unsigned int s = 0; s < target.size(); s++)
	{
		MatrixType BB = target[s]->GetBoundingBox();
		for (int d = 0; d < Dimension; d++)
		{
			dataDomain(d,0) = ( dataDomain(d,0)<BB(d,0)?dataDomain(d,0):BB(d,0) );
			dataDomain(d,1) = ( dataDomain(d,1)>BB(d,1)?dataDomain(d,1):BB(d,1) );
		}
	}

	m_Def->SetDataDomain( dataDomain ); // will be used to define bounding box with padding factor of 0.5*PaddingFactor*DeformationKernelWidth

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetDataDomain( dataDomain );  // will be used to define bounding box with padding factor of    PaddingFactor*DeformationKernelWidth
}


#endif /* _Regression_txx */

#ifndef _deform_txx
#define _deform_txx

#include "KernelFactory.h"
#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "DeformableObject.h"

#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

#include "DeformableObjectReader.h"

#include "writeMatrixDLM.txx"
#include "readMatrixDLM.txx"

#if ITK_VERSION_MAJOR >= 4
#include "itkFFTWGlobalConfiguration.h"
#endif

#include "itksys/SystemTools.hxx"

#include <cstring>
#include <iostream>
#include <sstream>



template <unsigned int Dimension>
void deform(SparseDiffeoParameters* paramDiffeos,
bool useInverseFlow,
const char* CP_fn,
const char* Mom_fn,
int numObjects,
std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
const std::vector<char*> objectfnList,
int subjectId)
{
	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

#ifdef USE_DOUBLE_PRECISION
	typedef double TScalar;
	std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
	typedef float TScalar;
	std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

	std::cout << "Shoot and Flow\n===" << std::endl;
	std::cout << "Number of objects to flow: " << numObjects << std::endl;
	std::cout << "\n===\n" << std::endl;
	std::cout << "Deformations Parameters:" << std::endl;
	paramDiffeos->PrintSelf(std::cout);
	std::cout << "\n===\n" << std::endl;
	for (int i = 0; i < numObjects; i++)
	{
		std::cout << "Object " << i << ": " << std::endl;
		std::cout << "Source file: " << objectfnList[i] << std::endl;
		paramObjectsList[i]->PrintSelf(std::cout);
		std::cout << "\n===\n" << std::endl;
	}
	
#if ITK_VERSION_MAJOR >= 4
	itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

	typedef Diffeos<TScalar, Dimension> DiffeosType;
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	typedef std::vector<DeformableObjectType*> DeformableObjectList;
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	// typedef Template<TScalar, Dimension> TemplateType;

	typedef typename LinearAlgebra<TScalar>::Vector VectorType;
	typedef typename LinearAlgebra<TScalar>::Matrix MatrixType;
	typedef std::vector<MatrixType> MatrixList;


	// Read initial control point positions and momenta:
	MatrixType CP0 = readMatrixDLM<TScalar>(CP_fn);
	MatrixList MOM0 = readMultipleMatrixDLM<TScalar>(Mom_fn);
	MatrixType MOM0_i;
	
	if (subjectId < 1 || subjectId > MOM0.size())
		throw std::runtime_error("SubjectID is invalid!");

	//for (int i = 0; i < MOM0.size(); i++)
	//{
		MOM0_i = MOM0.at(subjectId - 1);
		std::cout << "number of CPs: " << CP0.rows() << std::endl;

		if (CP0.rows() != MOM0_i.rows())
			throw std::runtime_error("Number of CPs and Momentas mismatched for the given subject !");

	//}

	//std::cout << "number of CPs: " << CP0.rows() << std::endl;

	//if (CP0.rows() != MOM0.rows())
		//throw std::runtime_error("Number of CPs and Momentas mismatched");

	// Set up the kernel factory
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
	kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

	// create the deformation object
	DiffeosType* def = new DiffeosType();
	def->SetKernelWidth(paramDiffeos->GetKernelWidth());
	def->SetT0(paramDiffeos->GetT0());
	def->SetTN(paramDiffeos->GetTN());
	def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
	def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
	def->UseImprovedEuler();
	if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0)
	{
		def->SetKernelType(P3M);
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "fgt") == 0)
	{
		def->SetKernelType(FGT);
	}
#ifdef USE_CUDA
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0)
	{
		def->SetKernelType(CUDAExact);
	}
#endif
	else
	{
		if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
			std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
		def->SetKernelType(Exact);
	}
	

	// create source object
	DeformableObjectList objectList(numObjects);
	for (int i = 0; i < numObjects; i++)
	{
		typedef DeformableObjectReader<TScalar, Dimension> ObjectReaderType;
		ObjectReaderType* reader = new ObjectReaderType();
		reader->SetObjectParameters(paramObjectsList[i]);

		reader->SetFileName(objectfnList[i]);
		//reader->SetTemplateType();
		reader->Update();

		objectList[i] = reader->GetOutput();
	}

	DeformableMultiObjectType* object = new DeformableMultiObjectType();
	object->SetObjectList(objectList);
	object->Update();

	for (int i = 0; i < numObjects; i++)
		delete objectList[i];

	// Determine the bounding box
	MatrixType DataDomain = object->GetBoundingBox();
	for (int k = 0; k < CP0.rows(); k++)
	{
		for (int d = 0; d < Dimension; d++)
		{
			DataDomain(d,0) = ( CP0(k, d)<DataDomain(d,0)?CP0(k, d):DataDomain(d,0) );
			DataDomain(d,1) = ( CP0(k, d)>DataDomain(d,1)?CP0(k, d):DataDomain(d,1) );
		}
	}
	// Pass DataDomain to kernel factory (for p3m kernel in the default mode) and deformation (for the definition of the "out of box" exception)
	kfac->SetDataDomain(DataDomain);
	def->SetDataDomain(DataDomain);
	std::cout << "Working domain: origin =  [" << DataDomain.get_column(0) << "] length = [" << DataDomain.get_column(1) - DataDomain.get_column(0) << "]" << std::endl;

	// shoot and flow
	std::cout << "Shoot and Flow" << std::endl;
	def->SetStartPositions(CP0);
	def->SetStartMomentas(MOM0_i);
	def->SetDeformableMultiObject(object);
	def->Update();

	if (def->OutOfBox())
		throw std::runtime_error("Out of box");
			
	// Write trajectories of control points and momenta
	MatrixList TrajectoryPositions = def->GetTrajectoryPositions();
	MatrixList TrajectoryMomentas = def->GetTrajectoryMomentas();
	for (int t = 0; t < def->GetNumberOfTimePoints(); t++)
	{
		std::ostringstream oss;
		oss << "CP_t_" << t  << ".txt" << std::ends;
		writeMatrixDLM<TScalar>(oss.str().c_str(), TrajectoryPositions[t]);

		std::ostringstream oss2;
		oss2 << "MOM_t_" << t  << ".txt" << std::ends;
		writeMatrixDLM<TScalar>(oss2.str().c_str(), TrajectoryMomentas[t]);
	}
		
	// Write deformation of objects
	std::cout << "Write output files" << std::endl;

	std::vector<std::string> outname(numObjects);
	std::vector<std::string> outext(numObjects);
	std::vector<std::string> outfn(numObjects);
	for (int i = 0; i < numObjects; i++)
	{
		std::string sourcefn;
		sourcefn.assign(objectfnList[i]);
		int length = sourcefn.size();
		int index = length - 1;
		while ((sourcefn[index] != '.') && (index >= 0))
			index--;

		outname[i] = sourcefn.substr(0, index);
		outext[i] = sourcefn.substr(index, length - index);

		std::ostringstream oss;
		oss << outname[i] << "_flow";
		outname[i] = oss.str();
	}

	def->WriteFlow(outname, outext);

	delete object;

}

#endif

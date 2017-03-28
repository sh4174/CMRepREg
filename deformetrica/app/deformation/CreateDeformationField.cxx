#ifndef _deformField_cxx
#define _deformField_cxx

#include "DeformetricaConfig.h"

#include "LinearAlgebra.h"

#include "SparseDiffeoParametersXMLFile.h"
#include "KernelFactory.h"
// #include "SparseDiffeoMatcher.h"
#include "readMatrixDLM.txx"
#include "DeformationFieldIO.h"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
    std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

#ifdef USE_DOUBLE_PRECISION
	typedef double TScalar;
	std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
	typedef float TScalar;
	std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

	if (argc != 8)
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml Direction CP.txt MOM.txt AnatomicalOrientationCode save_option output_filename" << std::endl;
		std::cerr << "where Direction = -1 for using the inverse flow, +1 otherwise" << std::endl;
		std::cerr << "      Save_option = 0 (default) for saving the final deformation field only, 1 for saving all intermediate results as well" << std::endl;
		return -1;
	}

	// Read general parameters for diffeomorphic matching:
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	int dir = atoi(argv[2]);
	if (dir == 1)
		std::cout << "Use the forward flow of diffeomorphisms" << std::endl;
	else if (dir == -1)
		std::cout << "Use the backward flow of diffeomorphisms" << std::endl;
	else
	{
		std::cout << "Unknown value for direction: " << dir << " is neither +1 or -1" << std::endl;
		return -1;
	}

	//	bool useInverseFlow = (dir == -1);

	const unsigned int Dimension = 3;

	typedef LinearAlgebra<TScalar>::Matrix MatrixType;
	typedef std::vector<MatrixType> MatrixList;

	// Read initial CPs and Momentas:
	MatrixType CP0 = readMatrixDLM<TScalar>(argv[3]);

	MatrixList MOM0 = readMultipleMatrixDLM<TScalar>(argv[4]);
	MatrixType MOM0_i;

	for (int i = 0; i < MOM0.size(); i++)
	{
		MOM0_i = MOM0.at(i);
		std::cout << "number of CPs: " << CP0.rows() << std::endl;

		if (CP0.rows() != MOM0_i.rows())
			throw std::runtime_error("Number of CPs and Momentas mismatched for a given object !");


	}
	//std::cout << "number of CPs: " << CP0.rows() << std::endl;

	//if (CP0.rows() != MOM0.rows())
	//throw std::runtime_error("Number of CPs and Momentas mismatched");

	//string outFn(argv[5]);

	typedef Diffeos<TScalar, Dimension> DiffeosType;
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	typedef std::vector<DeformableObjectType*> DeformableObjectList;
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

	// Set up the kernel factory:
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
	kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

	// Create the deformation object:
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


	try
	{

		// Set bounding box:
		MatrixType DataDomain(Dimension, 2);

		// Hard coded DD for acoustic images @ Pitié
		DataDomain(0,0) = -1.0;
		DataDomain(0,1) = 40.875;
		DataDomain(1,0) = -50.2;
		DataDomain(1,1) = 27.925;
		DataDomain(2,0) = -51.567;
		DataDomain(2,1) = 7.7996;

		std::cout << "Working domain: origin =  [" << DataDomain.get_column(0) << "] length = [" <<
				DataDomain.get_column(1) - DataDomain.get_column(0) << "]" << std::endl;

		def->SetDataDomain(DataDomain);
		kfac->SetDataDomain(DataDomain);

		for (int i = 0; i < MOM0.size(); i++)
		{
			MOM0_i = MOM0.at(i);
			std::cout << "number of CPs: " << CP0.rows() << std::endl;


			// Shoot and flow
			std::cout << "Creating deformation field..." << std::endl;
			def->SetStartPositions(CP0);
			def->SetStartMomentas(MOM0_i);
			def->Update();

			if (def->OutOfBox())
				throw std::runtime_error("out of box");

			//std::cout << "... and Flow" << std::endl;
			DeformationFieldIO<TScalar, Dimension> defFieldIO;
			defFieldIO.SetDiffeos(def);
			defFieldIO.SetAnatomicalCoordinateSystemLabel(argv[5]);

			if (strcmp(argv[6], "1") == 0)
				defFieldIO.WriteDeformationField(argv[7], false); // outFn
			else
				defFieldIO.WriteDeformationField(argv[7], true); // outFn
		}
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
		return -1;
	}
	catch (std::exception& e)
	{
		std::cerr << "Exception: " << e.what() << std::endl;
		return -1;
	}
	catch (std::string& s)
	{
		std::cerr << "Exception: " << s << std::endl;
		return -1;
	}
	catch (...)
	{
		std::cerr << "Unknown exception" << std::endl;
		return -1;
	}


	return 0;
}

#endif // _deformField_cxx

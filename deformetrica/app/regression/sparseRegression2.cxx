#include "DeformetricaConfig.h"
#include "KernelFactory.h"

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "estimateGeodesicRegression.txx"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{

	if (argc < 7)
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml NumberOfObjects paramsObject1.xml InitialTemplate1 Observation1 Time1 Observation2 Time2 Observation3 Time3 ... paramsObject2.xml InitialTemplate2 Observation1 Time1 Observation2 Time2 Observation3 Time3 ... " << std::endl;
		return -1;
	}

	int numObjects = atoi(argv[2]);
	std::cout << "Number of Objects = " << numObjects << std::endl;

	if ((argc - 3) % numObjects != 0)
	{
		std::cerr << "Number of files mismatched with the number of objects" << std::endl;
		return -1;
	}

	// The first 3 arguments are required, so (argc-3) tells us how many arguments potentially contain observations.
	// The observation list repeats for each object, so we divide by the number of objects.
	// The first two parameters of list are a object parameter xml file and the initial template, so we subtract 2.
	// The observations come in tuples of (data, time) so we divide by 2.
	int numObservations = (int) ((argc - 3) / numObjects - 2) / 2;
	std::cout << "Number of observations = " << numObservations << std::endl;

	if (numObservations < 2)
	{
		std::cerr << "Regression requires at least 2 observations" << std::endl;
		return -1;
	}

	// Read general parameters for diffeomorphic matching
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	// Read the list of objects: params, source and target
	std::vector<char*> templatefn;
	std::vector<std::vector<char*> > observationfn;
	std::vector<std::vector<double> > observationTimes;
	std::vector<DeformableObjectParameters::Pointer> paramObjects;

	templatefn.resize(numObjects);
	paramObjects.resize(numObjects);
	observationfn.resize(numObjects);
	observationTimes.resize(numObjects);
	for (int i = 0; i < numObjects; i++)
	{
		observationfn[i].resize(numObservations);
		observationTimes[i].resize(numObservations);
	}

	unsigned int indx = 3;
	for (unsigned int i = 0; i < numObjects; i++)
	{
		paramObjects[i] = readDeformableObjectParametersXML(argv[indx++]);
		if (paramObjects[i].IsNull())
		{
			std::cerr << "Failed creating XML object, bad input file " << argv[indx] << "?" << std::endl;
			return -1;
		}
		templatefn[i] = argv[indx++];
		for (unsigned s = 0; s < numObservations; s++)
		{
			observationfn[i][s] = argv[indx++];
			observationTimes[i][s] = atof(argv[indx++]);
		}

	}

	try
	{
		estimateGeodesicRegression<2>(paramDiffeos, numObservations, numObjects, paramObjects, templatefn, observationfn, observationTimes);
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

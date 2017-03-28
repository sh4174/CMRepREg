#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkMassProperties.h"
#include "vtkDecimatePro.h"

#include "surfio.h"

int main(int argc, char** argv)
{
	if (argc != 4)
	{
    		std::cerr << "Usage: " << argv[0] << " <in-poly> <out-poly> <percent-reduction>" << std::endl;
    		return -1;
  	}
	
	float percentReduce = atof(argv[3]);

	if ((percentReduce <= 0) || (percentReduce >= 1.00))
	{
		std::cerr << "<percent-reduction> must be >0 and <1" << std::endl;
    		return -1;
	}


  	vtkSmartPointer<vtkPolyData> poly = readSurface(argv[1]);
	
	//std::cout << "Before decimation" << std::endl << "------------" << std::endl;
	//std::cout << "There are " << poly->GetNumberOfPoints() << " points." << std::endl;
	//std::cout << "There are " << poly->GetNumberOfPolys() << " polygons." << std::endl;	

  	vtkSmartPointer<vtkDecimatePro> decimate =  vtkSmartPointer<vtkDecimatePro>::New();
	decimate->SetInputData(poly);
	decimate->SetTargetReduction(percentReduce);
  	decimate->Update();

	vtkSmartPointer<vtkPolyData> decimated = vtkSmartPointer<vtkPolyData>::New();
  	decimated->ShallowCopy(decimate->GetOutput());
 
	//std::cout << "After decimation" << std::endl << "------------" << std::endl;
	//std::cout << "There are " << decimated->GetNumberOfPoints() << " points." << std::endl;
	//std::cout << "There are " << decimated->GetNumberOfPolys() << " polygons." << std::endl;

	// Write the output surface mesh
	writeSurface(argv[2], decimate->GetOutput());

  	return 0;
}

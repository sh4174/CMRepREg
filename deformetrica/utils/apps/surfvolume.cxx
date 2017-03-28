#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkMassProperties.h"

#include "surfio.h"

int main(int argc, char** argv)
{
	if (argc != 2)
	{
    		std::cerr << "Usage: " << argv[0] << " <in-poly>" << std::endl;
    		return -1;
  	}

  	vtkSmartPointer<vtkPolyData> poly = readSurface(argv[1]);
	vtkMassProperties *massProperty = vtkMassProperties::New(); 

        massProperty->SetInputData(poly);

	double vol = massProperty->GetVolume();
	std::cout << vol << std::endl;

  	return 0;
}

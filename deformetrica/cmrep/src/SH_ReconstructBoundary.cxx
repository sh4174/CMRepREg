#include "ScriptInterface.h"
//#include "BasisFunctions2D.h"
//#include "MedialAtom.h"
//#include "CartesianMedialModel.h"
//#include "OptimizationTerms.h"
//#include "CoefficientMapping.h"
//#include "MedialAtomGrid.h"
//#include "PrincipalComponents.h"
//#include "System.h"
//#include "TestSolver.h"
//#include "ITKImageWrapper.h"
//#include <itksys/SystemTools.hxx>
//#include "itk_to_nifti_xform.h"
//
//// ITK includes
//#include <itkOrientedRASImage.h>
//#include <itkVTKImageExport.h>
//
//// VTK includes
//#include <vtkImageData.h>
//#include <vtkImageImport.h>
//#include <vtkMarchingCubes.h>
//#include <vtkPolyDataWriter.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkTransform.h>
//#include <vtkMatrix4x4.h>
//
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataReader.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>

using namespace std;

int usage()
{
  cout << "boundary_recon: Reconstruct a shape boundary from output cmrep (.vtk)" << endl;
  cout << "usage: " << endl; 
  cout << "   boundary_recon -i [input cmrep vtk ] -o [outputpath]" << endl;
  cout << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Check the number of input parameters
  if(argc < 5)
    return usage();

  // vtkSmartPoint< vtkPolyData > cmrep_poly = vtkSmartPointer< vtkPolyData >::New();
  string cmrepPath = "";
  string outPath = "";

  size_t argoptmax = 6;
  for(size_t i = 1; i < argc; i++)
  {
    string arg = argv[i];
    if(arg == "-i" && i < argc-1)
    {
    	cmrepPath = argv[ i + 1 ];
		i++;
    }
    else if(arg == "-o" && i < argc-1)
    {
		outPath = argv[ i + 1 ];
    	i++;
    }
    else
    {
      cerr << "Unknown option " << arg << endl;
      return -1;
    }
  }

//  // Read the current model
//  SubdivisionMPDE smod(cmrepPath.c_str());
//  smod.SubdivideMeshes(0, 0 );
//
//  // Create a filename for the output model
//  cout << "Read SubDiv" << endl;
//  string subdivPath = cmrepPath.substr( 0, cmrepPath.find( '.' ) ) + "_subdiv.cmrep";
//  smod.SaveToParameterFile(subdivPath.c_str());
//
//  // Load the model
//  cout << "Read CMRep" << endl;

  MedialPDE mrep(cmrepPath.c_str());
  mrep.UpdateAtoms();

  // Save the model
  // mrep.SaveToParameterFile(outputPath);
  // fn_current = fn_result;

  // Write mesh file(s)

  // This is so that the float image is sampled
  mrep.SaveVTKMesh( "temp.vtk", outPath.c_str() );
}

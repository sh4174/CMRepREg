#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "CartesianMedialModel.h"
#include "OptimizationTerms.h"
#include "CoefficientMapping.h"
#include "MedialAtomGrid.h"
#include "PrincipalComponents.h"
#include "System.h"
#include "TestSolver.h"
#include "ITKImageWrapper.h"
#include <itksys/SystemTools.hxx>
#include "itk_to_nifti_xform.h"

// ITK includes 
#include <itkOrientedRASImage.h>
#include <itkVTKImageExport.h>

// VTK includes
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>

int usage()
{
  cout << "cmrep_fit: fit cm-rep to image" << endl;
  cout << "usage: " << endl; 
  cout << "  cmrep_fit [options] param.txt template.cmrep target.img output_dir" << endl;
  cout << "options: " << endl;
  cout << "  -s N   : Only run fitting stage N (first stage is stage 0)" << endl;
  cout << "  -g FN  : Supply 'gray' image filename for direct-to-image fitting" << endl;
  cout << "  -t     : Test gradient computation for each of the terms (debug)" << endl;
  cout << "  -d     : Dump out mesh with gradient vectors at each iteration (debug)" << endl;
  cout << "parameter file specification: " << endl;
  cout << "  http://alliance.seas.upenn.edu/~pauly2/wiki/index.php?n=Main.CM-RepFittingToolCmrFit" << endl;
  cout << endl;
  return -1;
}

// Generate contour and save to file
void GenerateContour(FloatImage *image, string file)
{
  // Get the contour
  vtkSmartPointer<vtkPolyData> mesh = GenerateContour(image);

  // Create a writer
  vtkSmartPointer<vtkPolyDataWriter> m_Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  m_Writer->SetFileName(file.c_str());
  m_Writer->SetInputData(mesh);
  m_Writer->Update();
}

// Data fitting modes
enum FitMode {ALIGN_MOMENTS, FIT_TO_BINARY, NONE};

struct SubdivisionRefinement
{
  size_t sub_atoms;
  size_t sub_control;
};

struct Refinement
{
  SubdivisionRefinement sub_ref;
};

struct StageInfo
{
  string name;
  FitMode mode;
  double blur;
  size_t max_iter;
  Registry param;
  Refinement ref;
};

int main(int argc, char *argv[])
{
  // Report errors
  SetupSignalHandlers();

  // Check the number of input parameters
  if(argc < 5)
    return usage();

  // Check optional parameters
  bool flag_one_stage = false;
  size_t selected_stage = 0;

  // Optimization flags
  OptimizationFlags flags_opt;

  size_t argoptmax = argc-4;
  string fn_gray = "";
  for(size_t i = 1; i < argoptmax; i++)
    {
    string arg = argv[i];
    if(arg == "-s" && i < argoptmax-1)
      {
      flag_one_stage = true;
      selected_stage = atoi(argv[i+1]);
      i++;
      }
    else if(arg == "-g" && i < argoptmax-1)
      {
      fn_gray = argv[i+1];
      i++;
      }
    else if(arg == "-t")
      {
      flags_opt.flagTestGradient = true;
      }
    else if(arg == "-d")
      {
      flags_opt.flagDumpGradientMesh = true;
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }

  // Read the registry
  Registry r;
  try
    {
    r.ReadFromFile(argv[argc-4]);
    }
  catch(...)
    {
    cerr << "Unable to read parameters from " << argv[argc-4] << endl;
    cerr << "Call without parameters for proper usage" << endl;
    return -1;
    }

  // Registry-enum thingy
  RegistryEnumMap<FitMode> xModeEnumMap;
  xModeEnumMap.AddPair(ALIGN_MOMENTS, "AlignMoments");
  xModeEnumMap.AddPair(FIT_TO_BINARY, "FitToBinary");
  xModeEnumMap.AddPair(NONE, "NONE");

  try 
    {
    // Read the number of stages
    size_t n_stages = r["Stage.ArraySize"][0];
    if(n_stages == 0)
      throw string("No optimization stages found in parameter file");

    // Create the vector of stages
    vector<StageInfo> stages;

    // Read the default parameters
    Registry pdef = r.Folder("DefaultParameters");

    // For each stage read the interesting information
    for(size_t i = 0; i < n_stages; i++)
      {
      stages.push_back(StageInfo());

      // Read the stage folder
      string key = r.Key("Stage.Element[%d]",i);
      if(!r.IsFolder(key))
        throw string("Missing folder ") + key;
      Registry rs = r.Folder(key);

      // Get the name
      ostringstream oss;
      oss << "stage " << i;
      stages[i].name = rs["Name"][oss.str()];

      // Get the mode
      stages[i].mode = rs["Mode"].GetEnum(xModeEnumMap, NONE);
      if(stages[i].mode == NONE)
        throw string("Missing 'Mode' for stage ") + stages[i].name;

      // Get the number of iterations
      stages[i].max_iter = rs["MaxIterations"][1000];
      stages[i].blur = rs["Blur"][1.2];

      // Get the refinement info (depends on model type)
      stages[i].ref.sub_ref.sub_atoms = rs["Refinement.Subdivision.Atoms"][0];
      stages[i].ref.sub_ref.sub_control = rs["Refinement.Subdivision.Controls"][0];

      // Read the overriding parameters
      stages[i].param = pdef;
      Registry rp = rs.Folder("Parameters");

      // Get the list of all overriding parameters
      Registry::StringListType klist;
      rp.CollectKeys(klist);
      for(Registry::StringListType::iterator it = klist.begin(); it != klist.end(); it++)
        {
        stages[i].param[*it] << rp.Entry(*it).GetInternalString();
        }

      }

    // Read the input template
    string fn_temp = argv[argc-3];
    try
      {
      cout << fn_temp.c_str() << endl;      
      MedialPDE cmr_template(fn_temp.c_str());
      }
    catch(...)
      {
      throw string("Unable to read template file ") + fn_temp;
      }

    // Read the input image
    string fn_image = argv[argc-2];
    BinaryImage img;
    try
      {
      img.LoadFromFile(fn_image.c_str());
      }
    catch(...)
      {
      throw string("Unable to read image ") + fn_image;
      }

    // Float image for future use
    FloatImage imgfloat;

    // Grayscale image (if supplied)
    FloatImage imggray;
    if(fn_gray.length())
      {
      imggray.LoadFromFile(fn_gray.c_str());
      }

    // Before starting the loop, we need a filename for the current model
    string fn_current = fn_temp;

    // Get the output directory
    string dir_out = argv[argc-1];

    // Create all the output directories
    itksys::SystemTools::MakeDirectory((dir_out + "/cmrep").c_str());
    itksys::SystemTools::MakeDirectory((dir_out + "/mesh").c_str());
    itksys::SystemTools::MakeDirectory((dir_out + "/image").c_str());

    // Process each of the stages in turn
    for(size_t i = 0; i < n_stages; i++)
      {
      // If we are processing only one stage, we skip the rest
      if(flag_one_stage)
        {
        if(selected_stage != i) continue;
        if(i > 0)
          {
          // Load the starting point for this stage
          fn_current = dir_out + "/cmrep/" + stages[i-1].name + ".cmrep";
          }
        }

      cout << "### STAGE " << i << " ###" << endl;

      // Blur the image at appropriate blur level
      if(i == 0 || flag_one_stage || stages[i].blur != stages[i-1].blur)
        {
        if(stages[i].blur == 0.0)
          {
          cout << "Treating input image as smooth speed image (bkg -1, fore: +1)" << endl;
          imgfloat.LoadFromFile(fn_image.c_str());
          imgfloat.SetOutsideValue(-1.0);
          }
        else
          {
          cout << "Treating input image as binary (background 0, foreground 1)" << endl;
          cout << "Gaussian smoothing with sigma = " << stages[i].blur << endl;
          try 
            {
            imgfloat.SetToBlurredBinary(&img, stages[i].blur);
            imgfloat.SetOutsideValue(-1.0);
            }
          catch(exception &exc)
            {
            cerr << exc.what() << endl;
            return -1;
            }
          cout << "Done blurring" << endl;
          }
        }

      // Subdivide the model if necessary
      if(stages[i].ref.sub_ref.sub_atoms > 0 || stages[i].ref.sub_ref.sub_control > 0)
        {
        // Read the current model
        SubdivisionMPDE smod(fn_current.c_str());
        smod.SubdivideMeshes(stages[i].ref.sub_ref.sub_control, stages[i].ref.sub_ref.sub_atoms);

        // Create a filename for the output model
        string fn_init = dir_out + "/cmrep/" + stages[i].name + "_init.cmrep";
        smod.SaveToParameterFile(fn_init.c_str());
        fn_current = fn_init;
        }

      // Create a filename for the output model
      cout << dir_out << endl;

      string fn_result = dir_out + "/cmrep/" + stages[i].name + ".cmrep";

      cout << fn_current.c_str() << endl;

      // Load the model
      MedialPDE mrep(fn_current.c_str());

      // Use the appropriate optimization approach depending on mode
      if(stages[i].mode == ALIGN_MOMENTS)
        {
        mrep.MatchImageByMoments(&imgfloat, 10);
        }
      else if(stages[i].mode == FIT_TO_BINARY)
        {
        mrep.RunOptimization(&imgfloat, stages[i].max_iter, stages[i].param, flags_opt, &imggray);
        }
      else throw string("Unknown optimization mode");

      // Save the model
      mrep.SaveToParameterFile(fn_result.c_str());
      fn_current = fn_result;

      // Write mesh file(s)
      string fn_mesh_med = dir_out + "/mesh/" + stages[i].name + ".med.vtk";
      string fn_mesh_bnd = dir_out + "/mesh/" + stages[i].name + ".bnd.vtk";
      
      // This is so that the float image is sampled
      mrep.SaveVTKMesh(fn_mesh_med.c_str(), fn_mesh_bnd.c_str());

      // Write the target image
      string fn_target_mesh = dir_out + "/mesh/" + stages[i].name + "_target.vtk";
      GenerateContour(&imgfloat, fn_target_mesh);

      }
    }
  catch(MedialModelException &exc)
    {
    cerr << "MedialModelException: " << exc.what() << endl;
    return -1;
    }
  catch(string &err)
    {
    cerr << "Exception: " << err << endl;
    return -1;
    }
  catch(...)
    {
    cerr << "Unknown exception" << endl;
    return -1;
    }
}

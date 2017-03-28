#include "itkImageFileReader.h"
#include "itkOrientedRASImage.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkMarchingCubes.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include <vtkMatrix4x4.h>
#include "vnl/vnl_matrix_fixed.h"
#include "ReadWriteVTK.h"
#include <iostream>
#include "itkLinearInterpolateImageFunction.h"
#include "vtkImplicitFunction.h"
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkCell.h"

#include "MeshTraversal.h"

#include "itk_to_nifti_xform.h"
#include <vnl/vnl_inverse.h>

using namespace std;

int usage()
{
  cout << "Usage: vtklevelset [options] input.img output.vtk threshold" << endl;
  cout << "Options: " << endl;
  cout << "  -c clipimage.img   : clip output (keep part where clipimage > 0" << endl;
  cout << "  -v                 : export mesh in voxel coordinates (not physical)" << endl;
  cout << "  -k                 : apply clean filter to the mesh" << endl;
  cout << "  -d                 : perform Delaunay edge flipping" << endl;
  return -1;
}

template<class TImage>
void ConnectITKToVTK(itk::VTKImageExport<TImage> *fltExport,vtkImageImport *fltImport)
{
  fltImport->SetUpdateInformationCallback( fltExport->GetUpdateInformationCallback());
  fltImport->SetPipelineModifiedCallback( fltExport->GetPipelineModifiedCallback());
  fltImport->SetWholeExtentCallback( fltExport->GetWholeExtentCallback());
  fltImport->SetSpacingCallback( fltExport->GetSpacingCallback());
  fltImport->SetOriginCallback( fltExport->GetOriginCallback());
  fltImport->SetScalarTypeCallback( fltExport->GetScalarTypeCallback());
  fltImport->SetNumberOfComponentsCallback( fltExport->GetNumberOfComponentsCallback());
  fltImport->SetPropagateUpdateExtentCallback( fltExport->GetPropagateUpdateExtentCallback());
  fltImport->SetUpdateDataCallback( fltExport->GetUpdateDataCallback());
  fltImport->SetDataExtentCallback( fltExport->GetDataExtentCallback());
  fltImport->SetBufferPointerCallback( fltExport->GetBufferPointerCallback());
  fltImport->SetCallbackUserData( fltExport->GetCallbackUserData());
}

template<class TImage>
class ClipFunction : public vtkImplicitFunction 
{
public:

  vtkTypeMacro(ClipFunction, vtkImplicitFunction);

  // This is used by the clip code
  double EvaluateFunction(double x[3]);

  // This is not really used in the clip code
  void EvaluateGradient(double x[3], double g[3]) {}

  // Set the image
  void SetImage(TImage *image);

private:

  typedef itk::LinearInterpolateImageFunction<TImage> func;



};

int main(int argc, char *argv[])
{
  if(argc < 4)
    return usage();

  // Clip image
  const char *imClip = NULL;
  bool voxelSpace = false;
  bool flag_clean = false, flag_delaunay = false;
  for(int i = 1; i < argc - 3; i++)
    {
    if(!strcmp(argv[i], "-c"))
      {
      imClip = argv[++i];
      }
    else if(!strcmp(argv[i], "-v"))
      {
      voxelSpace = true;
      }
    else if(!strcmp(argv[i], "-d"))
      {
      flag_delaunay = true;
      }
    else if(!strcmp(argv[i], "-k"))
      {
      flag_clean = true;
      }
    else 
      { 
      cerr << "Unknown option " << argv[i] << endl;
      return -1;
      }
    }

  // Read the input image
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[argc-3]);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  // Get the range of the input image
  float imax = imgInput->GetBufferPointer()[0];
  float imin = imax;
  for(size_t i = 0; i < imgInput->GetBufferedRegion().GetNumberOfPixels(); i++)
    {
    float x = imgInput->GetBufferPointer()[i];
    imax = std::max(imax, x);
    imin = std::min(imin, x);
    }

  float cut = atof(argv[argc-1]);
  cout << "Image Range: [" << imin << ", " << imax << "]" << endl;
  cout << "Taking level set at " << cut << endl;

  // Create an importer and an exporter in VTK
  typedef itk::VTKImageExport<ImageType> ExporterType;
  ExporterType::Pointer fltExport = ExporterType::New();
  fltExport->SetInput(imgInput);
  vtkImageImport *fltImport = vtkImageImport::New();
  ConnectITKToVTK(fltExport.GetPointer(), fltImport);

  // Run marching cubes on the input image
  vtkMarchingCubes *fltMarching = vtkMarchingCubes::New();
  fltMarching->SetInputConnection(fltImport->GetOutputPort());
  fltMarching->ComputeScalarsOff();
  fltMarching->ComputeGradientsOff();
  fltMarching->ComputeNormalsOn();
  fltMarching->SetNumberOfContours(1);
  fltMarching->SetValue(0,cut);
  fltMarching->Update();

  vtkPolyData *pipe_tail = fltMarching->GetOutput();

  // If the clean option is requested, use it
  if(flag_clean || flag_delaunay) 
    {
    vtkTriangleFilter *trifi = vtkTriangleFilter::New();
    trifi->SetInputData(pipe_tail);
    trifi->PassLinesOff();
    trifi->PassVertsOff();
    trifi->Update();

    vtkCleanPolyData *clean = vtkCleanPolyData::New();
    clean->SetInputConnection(trifi->GetOutputPort());
    clean->PointMergingOn();
    clean->SetTolerance(0.0);
    clean->Update();

    vtkTriangleFilter *trifi2 = vtkTriangleFilter::New();
    trifi2->SetInputConnection(clean->GetOutputPort());
    trifi2->PassLinesOff();
    trifi2->PassVertsOff();
    trifi2->Update();
    pipe_tail = trifi2->GetOutput();
    }

  if(flag_delaunay)
    {

    typedef vnl_vector_fixed<double, 3> Vec;
    Vec *X = new Vec[pipe_tail->GetNumberOfPoints()];
    for(size_t i = 0; i < pipe_tail->GetNumberOfPoints(); i++)
      for(size_t k = 0; k < 3; k++)
        X[i][k] = pipe_tail->GetPoint(i)[k];

    // Generate the mesh object
    TriangleMesh bmesh;
    TriangleMeshGenerator gen(&bmesh, pipe_tail->GetNumberOfPoints());
    for(int i = 0; i < pipe_tail->GetNumberOfCells(); i++)
      {
      vtkCell *c = pipe_tail->GetCell(i);
      gen.AddTriangle(c->GetPointId(0), c->GetPointId(1), c->GetPointId(2));
      }
    gen.GenerateMesh();

    // Perform Delaunay
    bmesh.MakeDelaunay(X);

    // Set new cell data on the mesh
    pipe_tail->DeleteCells();
    pipe_tail->Allocate(bmesh.triangles.size());
    for(int i = 0; i < bmesh.triangles.size(); i++)
      {
      Triangle &tri = bmesh.triangles[i];
      vtkIdType pts[3];
      pts[0] = tri.vertices[0];
      pts[1] = tri.vertices[1];
      pts[2] = tri.vertices[2];
      pipe_tail->InsertNextCell(VTK_TRIANGLE, 3, pts);
      }

    pipe_tail->BuildCells();
    pipe_tail->BuildLinks();
    }
  
  // Create the transform filter
  vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
  fltTransform->SetInputData(pipe_tail);
 
  // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
  typedef vnl_matrix_fixed<double, 4, 4> Mat44;
  Mat44 vtk2out;
  Mat44 vtk2nii = ConstructVTKtoNiftiTransform(
    imgInput->GetDirection().GetVnlMatrix(),
    imgInput->GetOrigin().GetVnlVector(),
    imgInput->GetSpacing().GetVnlVector());

  // If we actually asked for voxel coordinates, we need to fix that
  if(voxelSpace)
    {
    Mat44 vox2nii = ConstructNiftiSform(
      imgInput->GetDirection().GetVnlMatrix(),
      imgInput->GetOrigin().GetVnlVector(),
      imgInput->GetSpacing().GetVnlVector());
    
    Mat44 nii2vox = vnl_inverse(vox2nii);
    vtk2out = nii2vox * vtk2nii;
    }
  else
    {
    vtk2out = vtk2nii;
    }
  
  // Update the VTK transform to match
  vtkTransform *transform = vtkTransform::New();
  transform->SetMatrix(vtk2out.data_block());
  fltTransform->SetTransform(transform);
  fltTransform->Update();

  // Get final output
  vtkPolyData *mesh = fltTransform->GetOutput();

  // Flip normals if determinant of SFORM is negative
  if(transform->GetMatrix()->Determinant() < 0)
    {
    vtkPointData *pd = mesh->GetPointData();
    vtkDataArray *nrm = pd->GetNormals();
    for(size_t i = 0; i < (size_t)nrm->GetNumberOfTuples(); i++)
      for(size_t j = 0; j < (size_t)nrm->GetNumberOfComponents(); j++)
        nrm->SetComponent(i,j,-nrm->GetComponent(i,j));
    nrm->Modified();
    }

  // Write the output
  WriteVTKData(mesh, argv[argc-2]);
}

#include "PointSetHamiltonianSystem.h"
#include "lddmm_data.h"
#include <iostream>
#include "util/ReadWriteVTK.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "FastLinearInterpolator.h"

#include <cstdarg>

using namespace std;


int usage()
{
  cout << "lmtowarp: Generate warp field from landmark geodesic shooting" << endl;
  cout << "Usage:" << endl;
  cout << "  lmtowarp [options]" << endl;
  cout << "Required options:" << endl;
  cout << "  -r image.nii       : Reference image (defines warp space)" << endl;
  cout << "  -o warp.nii        : Output deformation field" << endl;
  cout << "  -m mesh.vtk        : Mesh with the InitialMomentum array defined" << endl;
  cout << "  -s sigma           : kernel standard deviation" << endl;
  cout << "Additional options:" << endl;
  cout << "  -d dim                     : problem dimension (3)" << endl;
  cout << "  -n N                       : number of time steps (100)" << endl;

  return -1;
}

void check(bool condition, char *format,...)
{
  if(!condition)
    {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsprintf (buffer,format, args);
    va_end (args);

    cerr << buffer << endl;
    exit(-1);
    }
}


struct WarpGenerationParameters
{
  string fnReference, fnMesh, fnOutWarp;
  double sigma;
  unsigned int dim;
  unsigned int N;

  WarpGenerationParameters():
    N(100), dim(3), sigma(0.0) {}
};

template <class TPixel, unsigned int VDim>
class PointSetGeodesicToWarp
{
public:
  typedef PointSetHamiltonianSystem<double, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;


  typedef LDDMMData<TPixel, VDim> LDDMM;
  typedef typename LDDMM::ImageType ImageType;
  typedef typename LDDMM::VectorImageType VectorImageType;
  typedef typename LDDMM::ImagePointer ImagePointer;
  typedef typename LDDMM::VectorImagePointer VectorImagePointer;

  typedef typename ImageType::PointType PointType;
  typedef FastLinearInterpolator<VectorImageType,TPixel,VDim> FastInterpolator;
  typedef itk::ContinuousIndex<TPixel,VDim> ContIndexType;

  static int run(const WarpGenerationParameters &param);
};

template <class TPixel, unsigned int VDim>
int
PointSetGeodesicToWarp<TPixel, VDim>
::run(const WarpGenerationParameters &param)
{
  // Read the VTK mesh containing the points
  vtkPolyData *mesh = ReadVTKData(param.fnMesh);
  if(!mesh)
    {
    cerr << "Failed to read mesh from " << param.fnMesh << endl;
    return -1;
    }

  // Read the momentum field
  vtkDataArray *arr_p0 = mesh->GetPointData()->GetArray("InitialMomentum");
  if(!arr_p0 || arr_p0->GetNumberOfComponents() != VDim)
    {
    cerr << "Failed to read initial momentum from " << param.fnMesh << endl;
    return -1;
    }

  // Count the number of non-null entries
  vector<unsigned int> index;
  for(unsigned int i = 0; i < arr_p0->GetNumberOfTuples(); i++)
    {
    bool has_value = 1;
    for(unsigned int a = 0; a < VDim; a++)
      {
      if(isnan(arr_p0->GetComponent(i,a)))
        {
        has_value = 0;
        break;
        }
      }
    if(has_value)
      index.push_back(i);
    }

  // Populate the q0 and p0 arrays
  unsigned int k = index.size();
  Matrix q0(k, VDim), p0(k,VDim), q1(k, VDim), p1(k,VDim);

  for(unsigned int i = 0; i < k; i++)
    {
    for(unsigned int a = 0; a < VDim; a++)
      {
      q0(i,a) = mesh->GetPoint(index[i])[a];
      p0(i,a) = arr_p0->GetComponent(i,a);
      }
    }

  // We have read the mesh successfully
  // Create the hamiltonian system
  HSystem hsys(q0, param.sigma, param.N);

  // Flow without gradients - we have streamlines
  hsys.FlowHamiltonian(p0, q1, p1);

  // Read the reference image
  ImagePointer imRef;
  LDDMM::img_read(param.fnReference.c_str(), imRef);

  // Image for splatting
  VectorImagePointer imSplat, imVelocity, imLagragean, imPhi;
  LDDMM::alloc_vimg(imSplat, imRef);
  LDDMM::alloc_vimg(imVelocity, imRef);
  LDDMM::alloc_vimg(imLagragean, imRef);
  LDDMM::alloc_vimg(imPhi, imRef);

  // Compute the scaling factor for the Gaussian smoothing
  double scaling = 1.0;
  for(unsigned int a = 0; a < VDim; a++)
    {
    scaling *= sqrt(2 * vnl_math::pi) * param.sigma / imRef->GetSpacing()[a];
    }

  // Direction of interpolation
  int tdir = 1;
  int tStart = (tdir > 0) ? 0 : param.N - 1;
  int tEnd = (tdir > 0) ? param.N : 0;

  // Iterate over time
  // for(unsigned int t = 0; t < param.N; t++)
  for(int t = tStart; t != tEnd; t+=tdir)
    {
    // Get all landmarks
    const Matrix &qt = hsys.GetQt(t);
    const Matrix &pt = hsys.GetPt(t);

    // Splatting approach
    // Splat each landmark onto the splat image
    imSplat->FillBuffer(0.0);

    // Create an interpolator for splatting
    FastInterpolator flint(imSplat);

    for(int i = 0; i < k; i++)
      {
      // Map landmark point to continuous index
      ContIndexType cix;
      PointType point;
      typename VectorImageType::InternalPixelType vec;
      for(unsigned int a = 0; a < VDim; a++)
        {
        // Here we have to correct for the fact that the landmark coordinates are 
        // assumed to be in RAS space, but point is in LPS space. We therefore have
        // to apply the LPS/RAS transform before splatting
        if(a < 2)
          {
          point[a] = -qt(i,a);
          vec[a] = -pt(i,a);
          }
        else
          {
          point[a] = qt(i,a);
          vec[a] = pt(i,a);
          }

        }
      imRef->TransformPhysicalPointToContinuousIndex(point, cix);

      // Splat landmark's momentum at the point location
      flint.Splat(cix.GetVnlVector().data_block(), &vec);
      }

    // Now all the momenta have been splatted. The next step is to smooth with a Gaussian
    typename LDDMM::Vec vec_sigma; vec_sigma.Fill(param.sigma);
    LDDMM::vimg_smooth(imSplat, imVelocity, vec_sigma);

    // Accumulate the velocity into phi - using imSplat as a temporary
    LDDMM::vimg_scale_in_place(imVelocity, scaling);


    /*
    // Direct approach - going to be slow!
    double f = -0.5 / (param.sigma * param.sigma);
    typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IterType;
    IterType iter(imVelocity,imVelocity->GetBufferedRegion());
    for(; !iter.IsAtEnd(); ++iter)
      {
      // Get coordinate of this pixel
      itk::Index<VDim> idx = iter.GetIndex();
      PointType X;
      imVelocity->TransformIndexToPhysicalPoint(idx, X);

      typename LDDMM::Vec vi(0.0);
      for(unsigned int i = 0; i < k; i++)
        {
        double dstsq = 0;
        for(unsigned int a = 0; a < VDim; a++)
          {
          double d = (X[a] - qt(i,a));
          dstsq += d * d;
          }
        double w = exp(dstsq * f);
        for(unsigned int a = 0; a < VDim; a++)
          {
          vi[a] += w * pt(i,a);
          }
        }

      iter.Set(vi);
      }
      */

    // Euler interpolation
    LDDMM::interp_vimg(imVelocity, imPhi, 1.0, imSplat, false, true);
    LDDMM::vimg_add_scaled_in_place(imPhi, imSplat, tdir * 1.0 / (param.N - 1));

    cout << "." << flush;
    }

  cout << endl;

  // Save the warp
  LDDMM::vimg_write(imPhi, param.fnOutWarp.c_str());

  return 0;
}



int main(int argc, char *argv[])
{
  WarpGenerationParameters param;

  if(argc < 2)
    return usage();

  // Process parameters
  for(int i = 1; i < argc; i++)
    {
    string arg = argv[i];
    if(arg == "-r" && i < argc-1)
      {
      param.fnReference = argv[++i];
      }
    else if(arg == "-m" && i < argc-1)
      {
      param.fnMesh = argv[++i];
      }
    else if(arg == "-o" && i < argc-1)
      {
      param.fnOutWarp = argv[++i];
      }
    else if(arg == "-s" && i < argc-1)
      {
      param.sigma = atof(argv[++i]);
      }
    else if(arg == "-d" && i < argc-1)
      {
      param.dim = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-n" && i < argc-1)
      {
      param.N = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-h")
      {
      return usage();
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }

  check(param.sigma > 0, "Missing or negative sigma parameter");
  check(param.N > 0 && param.N < 10000, "Incorrect N parameter");
  check(param.dim >= 2 && param.dim <= 3, "Incorrect N parameter");
  check(param.fnReference.length(), "Missing template filename");
  check(param.fnMesh.length(), "Missing target filename");
  check(param.fnOutWarp.length(), "Missing output filename");

  // Specialize by dimension
  if(param.dim == 2)
    return PointSetGeodesicToWarp<float,2>::run(param);
  else
    return PointSetGeodesicToWarp<float,3>::run(param);


  return 0;
}

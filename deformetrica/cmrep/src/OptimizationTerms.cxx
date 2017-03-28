#include "OptimizationTerms.h"
#include "MedialAtomGrid.h"
#include "MedialAtom.h"
#include "ITKImageWrapper.h"
#include "itkOrientedRASImage.h"
#include <iostream>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_random.h>
#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkQuadricClustering.h>
#include <vtkCell.h>

using namespace std;

namespace medialpde {

} // namespace
using namespace medialpde;

/*********************************************************************************
 * EUCLIDEAN FUNCTION JUNK
 * ******************************************************************************/
ImageSmoothSamplingEuclideanFunction
::ImageSmoothSamplingEuclideanFunction(
  FloatImage *image,
  double sigma, double cutoff)
{
  // Store the image
  this->image = image;

  // Get the underlying itk::Image
  typedef FloatImage::WrapperType::ImageType ImageType;
  ImageType::Pointer im = image->GetInternalImage()->GetInternalImage();

  // Get the origin and voxel size
  for(size_t d = 0; d < 3; d++)
    {
    xOrigin[d] = image->GetInternalImage()->GetImageOrigin(d);
    xInvVoxSize[d] = 1.0 / image->GetInternalImage()->GetImageSpacing(d);
    }

  // Compute the bounding box (in pixel coordinates)
  double bb_start[] = {0.0, 0.0, 0.0};
  double bb_end[3];
  for(size_t d = 0; d < 3; d++)
    bb_end[d] = im->GetBufferedRegion().GetSize()[d];

  // Create the smooth image sampler
  sis = new SmoothedImageSampler(
    sigma, cutoff, im->GetBufferPointer(), bb_start, bb_end);
}

/*********************************************************************************
 * SOLUTION DATA BASE CLASS
 ********************************************************************************/
SolutionDataBase::SolutionDataBase(
  MedialIterationContext *xAtomGrid)
{
  nAtoms = xAtomGrid->GetNumberOfAtoms();
  nBndPts = xAtomGrid->GetNumberOfBoundaryPoints();
  nMedTri = xAtomGrid->GetNumberOfTriangles();
  nBndTri = xAtomGrid->GetNumberOfBoundaryTriangles();
  xBoundaryWeights.resize(nBndPts, 0.0);
  xBoundaryTriangleUnitNormal.resize(nBndTri, SMLVec3d(0.0));
  xMedialWeights.resize(nAtoms, 0.0);
  xMedialTriangleUnitNormal.resize(nMedTri, SMLVec3d(0.0));
  xInteriorVolumeElement.resize(nBndPts, SMLVec3d(0.0));
  xBoundaryTriangleArea.resize(nBndTri, 0.0);
  xBoundaryAreaVector.resize(nBndPts, SMLVec3d(0.0));

  // Copy stuff
  this->xAtomGrid = xAtomGrid;
}

/*********************************************************************************
 * SOLUTION DATA                  
 ********************************************************************************/
SolutionData::SolutionData(
  MedialIterationContext *xGrid, MedialAtom *xAtoms) :
  SolutionDataBase(xGrid)
{
  this->xAtoms = xAtoms;
}

void SolutionData::ComputeIntegrationWeights()
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  const static double EIGHTEENTH = 1.0f / 18.0f;

  // Initialize the accumulators
  xMedialArea = 0.0;
  xBoundaryArea = 0.0;

  // Initialize the weight arrays to zero
  std::fill(xBoundaryWeights.begin(), xBoundaryWeights.end(), 0.0);
  std::fill(xMedialWeights.begin(), xMedialWeights.end(), 0.0);
  std::fill(xInteriorVolumeElement.begin(), xInteriorVolumeElement.end(), SMLVec3d(0.0));
  std::fill(xBoundaryAreaVector.begin(), xBoundaryAreaVector.end(), SMLVec3d(0.0));

  // Iterate over the medial triangles
  for(MedialTriangleIterator imt(xAtomGrid); !imt.IsAtEnd() ; ++imt)
    {
    // Access the four medial atoms
    size_t i0 = imt.GetAtomIndex(0);
    size_t i1 = imt.GetAtomIndex(1);
    size_t i2 = imt.GetAtomIndex(2);

    // Access the four medial points
    SMLVec3d X0 = xAtoms[i0].X; 
    SMLVec3d X1 = xAtoms[i1].X;
    SMLVec3d X2 = xAtoms[i2].X;

    // Compute the normal vector
    SMLVec3d N = THIRD * vnl_cross_3d(X1-X0,X2-X0);

    // Compute the area of triangle
    double Nmag = N.magnitude();
    double A = 0.5 * Nmag;

    // Store the quantity N / norm(N) for fast derivative computation
    assert(imt.GetIndex() < nMedTri);

    xMedialTriangleUnitNormal[imt.GetIndex()] = N / Nmag;
    
    // Add to the total area
    xMedialArea += A;

    // Assign a third of each weight to each corner
    assert(i0 < nAtoms && i1 < nAtoms && i2 < nAtoms);
    xMedialWeights[i0] += A; xMedialWeights[i1] += A; xMedialWeights[i2] += A;
    }

  // Scale the medial area by 3
  xMedialArea *= 3.0;

  // Iterate over the boundary triangles
  for(MedialBoundaryTriangleIterator ibt(xAtomGrid); !ibt.IsAtEnd() ; ++ibt)
    {
    // Access the four medial atoms
    size_t ia0 = ibt.GetAtomIndex(0);
    size_t ia1 = ibt.GetAtomIndex(1);
    size_t ia2 = ibt.GetAtomIndex(2);
    size_t ib0 = ibt.GetBoundaryIndex(0);
    size_t ib1 = ibt.GetBoundaryIndex(1);
    size_t ib2 = ibt.GetBoundaryIndex(2);

    // Access the four medial points
    SMLVec3d X0 = xAtoms[ia0].X;
    SMLVec3d X1 = xAtoms[ia1].X;
    SMLVec3d X2 = xAtoms[ia2].X;
    SMLVec3d Y0 = GetBoundaryPoint(ibt, xAtoms, 0).X;
    SMLVec3d Y1 = GetBoundaryPoint(ibt, xAtoms, 1).X;
    SMLVec3d Y2 = GetBoundaryPoint(ibt, xAtoms, 2).X;

    // Compute the area of the boundary triangle
    SMLVec3d NB = THIRD * vnl_cross_3d(Y1-Y0,Y2-Y0);
    double NBmag = NB.magnitude();
    double A = 0.5 * NBmag;

    // Add to the total area
    xBoundaryArea += A;

    // Store the quantity N / norm(N) for fast derivative computation
    xBoundaryTriangleArea[ibt.GetIndex()] = A;
    xBoundaryTriangleUnitNormal[ibt.GetIndex()] = NB / NBmag;
    assert(ibt.GetIndex() < nBndTri);
    
    // Compute the volume element coefficients
    SMLVec3d U0 = Y0 - X0, U1 = Y1 - X1, U2 = Y2 - X2;
    SMLVec3d W = EIGHTEENTH * (U0 + U1 + U2);
    SMLVec3d Za = vnl_cross_3d(X1-X0, X2-X0);
    SMLVec3d Zb = vnl_cross_3d(U1-U0, X2-X0) + vnl_cross_3d(X1-X0, U2-U0);
    SMLVec3d Zc = vnl_cross_3d(U1-U0, U2-U0);
    double va = dot_product(Za, W);
    double vb = dot_product(Zb, W);
    double vc = dot_product(Zc, W);
    SMLVec3d vvec(va,vb,vc);
    
    // Assign a third of each weight to each corner
    xBoundaryWeights[ib0] += A; xBoundaryWeights[ib1] += A; xBoundaryWeights[ib2] += A;

    xBoundaryAreaVector[ib0] += 0.5 * NB;
    xBoundaryAreaVector[ib1] += 0.5 * NB;
    xBoundaryAreaVector[ib2] += 0.5 * NB;

    // Update the volume element coeffs
    assert(ib0 < nBndPts && ib1 < nBndPts && ib2 < nBndPts);
    xInteriorVolumeElement[ib0] += vvec;
    xInteriorVolumeElement[ib1] += vvec;
    xInteriorVolumeElement[ib2] += vvec;
    }

  xBoundaryArea = 0.0;
  for(MedialBoundaryPointIterator ibp(xAtomGrid); !ibp.IsAtEnd() ; ++ibp)
    {
    xBoundaryWeights[ibp.GetIndex()] = dot_product(
      xAtoms[ibp.GetAtomIndex()].xBnd[ibp.GetBoundarySide()].N, 
      xBoundaryAreaVector[ibp.GetIndex()]);
    xBoundaryArea += xBoundaryWeights[ibp.GetIndex()];
    }
  

  // Scale the boundary area by 3
  // xBoundaryArea *= 3.0;
}

/*********************************************************************************
 * PARTIAL DERIVATIVE SOLUTION DATA
 ********************************************************************************/
PartialDerivativeSolutionData
::PartialDerivativeSolutionData(SolutionData *xReference, MedialAtom *dAtoms)
: SolutionDataBase(xReference->xAtomGrid)
{
  this->xReference = xReference;
  xAtoms = dAtoms;
}

void PartialDerivativeSolutionData
::ComputeIntegrationWeights()
{
  // A constant to hold 1/3
  const static double SIXTH = 1.0f / 6.0f;
  const static double EIGHTEENTH = 1.0f / 18.0f;

  // Initialize the accumulators
  xMedialArea = 0.0;
  xBoundaryArea = 0.0;

  // Initialize the weight arrays to zero
  std::fill(xBoundaryWeights.begin(), xBoundaryWeights.end(), 0.0);
  std::fill(xMedialWeights.begin(), xMedialWeights.end(), 0.0);
  std::fill(xInteriorVolumeElement.begin(), xInteriorVolumeElement.end(), SMLVec3d(0.0));
  std::fill(xBoundaryAreaVector.begin(), xBoundaryAreaVector.end(), SMLVec3d(0.0));

  // Iterate over the medial triangles
  for(MedialTriangleIterator imt(xAtomGrid); !imt.IsAtEnd() ; ++imt)
    {
    // Access the four medial atoms
    size_t i0 = imt.GetAtomIndex(0);
    size_t i1 = imt.GetAtomIndex(1);
    size_t i2 = imt.GetAtomIndex(2);

    // Check dependency
    if(xAtoms[i0].order == 0 || xAtoms[i1].order == 0 || xAtoms[i2].order == 0)
      {
      // Access the four medial points
      SMLVec3d DX0 = xAtoms[i0].X; 
      SMLVec3d DX1 = xAtoms[i1].X;
      SMLVec3d DX2 = xAtoms[i2].X;
      SMLVec3d X0 = xReference->xAtoms[i0].X; 
      SMLVec3d X1 = xReference->xAtoms[i1].X;
      SMLVec3d X2 = xReference->xAtoms[i2].X;

      // Compute the normal vector
      SMLVec3d DN_over_2 = SIXTH * (
        vnl_cross_3d(DX1-DX0,X2-X0) + vnl_cross_3d(X1-X0,DX2-DX0));

      // Compute the area of triangle
      double dA = dot_product(DN_over_2, xReference->xMedialTriangleUnitNormal[imt.GetIndex()]);

      // Add to the total area
      xMedialArea += dA;

      // Assign a third of each weight to each corner
      xMedialWeights[i0] += dA; xMedialWeights[i1] += dA; xMedialWeights[i2] += dA;
      }
    }

  // Scale the medial area by 3
  xMedialArea *= 3.0;

  // Iterate over the boundary triangles
  for(MedialBoundaryTriangleIterator ibt(xAtomGrid); !ibt.IsAtEnd() ; ++ibt)
    {
    // Access the four medial atoms
    size_t ia0 = ibt.GetAtomIndex(0);
    size_t ia1 = ibt.GetAtomIndex(1);
    size_t ia2 = ibt.GetAtomIndex(2);

    // Check dependency of the triangle
    if(xAtoms[ia0].order <= 1 || xAtoms[ia1].order <= 1 || xAtoms[ia2].order <= 1)
      {
      size_t ib0 = ibt.GetBoundaryIndex(0);
      size_t ib1 = ibt.GetBoundaryIndex(1);
      size_t ib2 = ibt.GetBoundaryIndex(2);

      // Access the four medial points
      SMLVec3d DX0 = xAtoms[ia0].X;
      SMLVec3d DX1 = xAtoms[ia1].X;
      SMLVec3d DX2 = xAtoms[ia2].X;
      SMLVec3d DY0 = GetBoundaryPoint(ibt, xAtoms, 0).X;
      SMLVec3d DY1 = GetBoundaryPoint(ibt, xAtoms, 1).X;
      SMLVec3d DY2 = GetBoundaryPoint(ibt, xAtoms, 2).X;
      SMLVec3d X0 = xReference->xAtoms[ia0].X;
      SMLVec3d X1 = xReference->xAtoms[ia1].X;
      SMLVec3d X2 = xReference->xAtoms[ia2].X;
      SMLVec3d Y0 = GetBoundaryPoint(ibt, xReference->xAtoms, 0).X;
      SMLVec3d Y1 = GetBoundaryPoint(ibt, xReference->xAtoms, 1).X;
      SMLVec3d Y2 = GetBoundaryPoint(ibt, xReference->xAtoms, 2).X;

      // Compute the area of the boundary triangle
      SMLVec3d dNB_over_two = SIXTH * (
        vnl_cross_3d(DY1-DY0,Y2-Y0) + vnl_cross_3d(Y1-Y0,DY2-DY0));
      double dA = dot_product(
        dNB_over_two, xReference->xBoundaryTriangleUnitNormal[ibt.GetIndex()]);

      xBoundaryTriangleArea[ibt.GetIndex()] = dA;


      // Add to the total area
      xBoundaryArea += dA;

      // Assign a third of each weight to each corner
      xBoundaryWeights[ib0] += dA; 
      xBoundaryWeights[ib1] += dA; 
      xBoundaryWeights[ib2] += dA;

      xBoundaryAreaVector[ib0] += dNB_over_two;
      xBoundaryAreaVector[ib1] += dNB_over_two;
      xBoundaryAreaVector[ib2] += dNB_over_two;

      // Compute the volume element coefficients
      SMLVec3d U0 = Y0 - X0, U1 = Y1 - X1, U2 = Y2 - X2;
      SMLVec3d DU0 = DY0 - DX0, DU1 = DY1 - DX1, DU2 = DY2 - DX2;
      SMLVec3d W = EIGHTEENTH * (U0 + U1 + U2);
      SMLVec3d DW = EIGHTEENTH * (DU0 + DU1 + DU2);
      SMLVec3d Za = vnl_cross_3d(X1-X0, X2-X0);
      SMLVec3d DZa = vnl_cross_3d(DX1-DX0, X2-X0) + vnl_cross_3d(X1-X0, DX2-DX0);
      SMLVec3d Zb = vnl_cross_3d(U1-U0, X2-X0) + vnl_cross_3d(X1-X0, U2-U0);
      SMLVec3d DZb =
        vnl_cross_3d(DU1-DU0, X2-X0) + vnl_cross_3d(DX1-DX0, U2-U0) +
        vnl_cross_3d(U1-U0, DX2-DX0) + vnl_cross_3d(X1-X0, DU2-DU0);
      SMLVec3d Zc = vnl_cross_3d(U1-U0, U2-U0);
      SMLVec3d DZc = vnl_cross_3d(DU1-DU0, U2-U0) + vnl_cross_3d(U1-U0, DU2-DU0);
      double dva = dot_product(DZa, W) + dot_product(Za, DW);
      double dvb = dot_product(DZb, W) + dot_product(Zb, DW);
      double dvc = dot_product(DZc, W) + dot_product(Zc, DW);
      SMLVec3d dvvec(dva,dvb,dvc);

      // Update the volume element coeffs
      xInteriorVolumeElement[ib0] += dvvec;
      xInteriorVolumeElement[ib1] += dvvec;
      xInteriorVolumeElement[ib2] += dvvec;
      }
    }

  xBoundaryArea = 0.0;
  for(MedialBoundaryPointIterator ibp(xAtomGrid); !ibp.IsAtEnd() ; ++ibp)
    {
    xBoundaryWeights[ibp.GetIndex()] = 
      dot_product(
        xAtoms[ibp.GetAtomIndex()].xBnd[ibp.GetBoundarySide()].N, 
        xReference->xBoundaryAreaVector[ibp.GetIndex()]) +
      dot_product(
        xReference->xAtoms[ibp.GetAtomIndex()].xBnd[ibp.GetBoundarySide()].N, 
        xBoundaryAreaVector[ibp.GetIndex()]);
    xBoundaryArea += xBoundaryWeights[ibp.GetIndex()];
    }  

  // Scale the medial area by 3
  // xBoundaryArea *= 3.0;
}

/*********************************************************************************
 * ENERGY TERM
 ********************************************************************************/

/*********************************************************************************
 * MEDIAL INTEGRATION ENERGY TERM
 ********************************************************************************/
MedialIntegrationEnergyTerm
::MedialIntegrationEnergyTerm(GenericMedialModel *model)
{
  // Generate an array of weights for irregular grids
  xDomainWeights.set_size(model->GetNumberOfAtoms());
  xDomainArea = ComputeMedialDomainAreaWeights(
    model->GetIterationContext(), model->GetAtomArray(), 
    xDomainWeights.data_block());
}

/*********************************************************************************
 * BOUNDARY IMAGE MATCH TERM
 ********************************************************************************/

// Print a verbose report
void 
BoundaryImageMatchTerm
::PrintReport(ostream &sout)
{
  sout << "  Boundary Image Match Term " << endl;
  sout << "    image match   : " << xImageMatch << endl;
  sout << "    boundary area : " << xBoundaryArea << endl;
  sout << "    ratio         : " << xFinalMatch << endl;
}

BoundaryImageMatchTerm
::BoundaryImageMatchTerm(
  GenericMedialModel *model, FloatImage *image)
{
  this->xImage = image; 
  xGradI = new SMLVec3d[model->GetNumberOfBoundaryPoints()];
  xImageVal = new double[model->GetNumberOfBoundaryPoints()];
}

double
BoundaryImageMatchTerm::UnifiedComputeEnergy(SolutionData *S, bool gradient_mode)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageSquareValueFunctionAdapter fImage(xImage);
  // ImageSmoothSamplingEuclideanFunction fImage(xImage, 2.0, 5.5);

  // Compute the image and image gradient at each point in the image
  xImageMatch = 0.0;

  // Loop over all boundary points
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Compute the image gradient
    SMLVec3d &X = GetBoundaryPoint(it, S->xAtoms).X;
    if(gradient_mode)
      xImageVal[it.GetIndex()] = fImage.ComputeFunctionAndGradient(X, xGradI[it.GetIndex()]);
    else
      xImageVal[it.GetIndex()] = fImage.Evaluate(X);

    // Accumulate to get weighted match
    xImageMatch += xImageVal[it.GetIndex()] * S->xBoundaryWeights[it.GetIndex()];
    }
  
  // We will need the area in many calculations
  xBoundaryArea = S->xBoundaryArea;

  // Compute the final match
  xFinalMatch = xImageMatch / xBoundaryArea;

  // Return the solution
  return xFinalMatch;
}

double
BoundaryImageMatchTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *DS)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageSquareValueFunctionAdapter fImage(xImage);
  // ImageSmoothSamplingEuclideanFunction fImage(xImage, 2.0, 5.5);

  // Accumulator for the partial derivative of the weighted match function
  double dMatchdC = 0.0;

  // Compute the partial derivative for this coefficient
  
  for(MedialBoundaryPointIterator it(S->xAtomGrid) ; !it.IsAtEnd(); ++it)
    {
    // Get the index of this boundary point
    size_t iPoint = it.GetIndex();

    // Get the change in the boundary point
    SMLVec3d dX = GetBoundaryPoint(it, DS->xAtoms).X;

    // Get the area weights for this point
    double w = S->xBoundaryWeights[ iPoint ];
    double dw = DS->xBoundaryWeights[ iPoint ];

    // Compute the change in intensity per change in coefficient
    double I = xImageVal[iPoint];
    double dIdC = dot_product( xGradI[iPoint] , dX);
    
    // Increment the partial derivative of the weighted match
    dMatchdC += dw * I + w * dIdC;
    }
  
  // Compute the derivative
  double dFinaldC = 
    ( dMatchdC * S->xBoundaryArea - xImageMatch * DS->xBoundaryArea ) / 
    ( S->xBoundaryArea * S->xBoundaryArea );

  // Return the final match derivative
  return dFinaldC;
}

BoundaryImageMatchTerm::~BoundaryImageMatchTerm()
{
  // Clean up
  delete xGradI;
  delete xImageVal;
}

/*********************************************************************************
 * Iterative Closest Point style term
 ********************************************************************************/

// TODO: put this into some header file
extern vtkSmartPointer<vtkPolyData> GenerateContour(FloatImage *image);
#include <vtkPolyDataWriter.h>

SymmetricClosestPointMatchTerm
::SymmetricClosestPointMatchTerm(
  GenericMedialModel *model, FloatImage *image, int nClusterDivisions)
{
  xIterCount = 0;

  nClusterDivisions = 96;
  // Make a copy of the image
  this->xImage = image;
  this->xModel = model;

  // Use the VTK contour code to extract a surface mesh from the image
  xMesh = GenerateContour(image);

  // Create a locator
  xTargetLocator = vtkSmartPointer<vtkCellLocator>::New();
  xTargetLocator->SetDataSet(xMesh.GetPointer());
  xTargetLocator->CacheCellBoundsOn();
  xTargetLocator->BuildLocator();

  // Need the bounds on the target
  xMesh->ComputeBounds();

  // Set up the clustering
  vtkSmartPointer<vtkQuadricClustering> clu =
      vtkSmartPointer<vtkQuadricClustering>::New();
  clu->SetInputData(xMesh);
  clu->SetDivisionOrigin(xMesh->GetCenter());
  double spacing = xMesh->GetLength() / nClusterDivisions;
  clu->SetDivisionSpacing(spacing, spacing, spacing);
  clu->Update();

  // Save the samples
  vtkSmartPointer<vtkPolyDataWriter> wr =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  wr->SetInputConnection(clu->GetOutputPort());
  wr->SetFileName("clusty.vtk");
  wr->Update();

  // Get the reduced target
  xMeshReduced = clu->GetOutput()->GetPoints();
  // xMeshReduced = xMesh->GetPoints();

  // Perform the closest point computation
  this->FindClosestPoints();
}

void 
SymmetricClosestPointMatchTerm
::FindClosestPoints()
{
  MedialIterationContext *context = xModel->GetIterationContext();

  // Perform the closest to model computation

  // Output vector
  xClosestToModel.resize(xModel->GetNumberOfBoundaryPoints());

  // Create a VTK points object used for target-model location
  vtkSmartPointer<vtkPoints> out_pts = vtkSmartPointer<vtkPoints>::New();
  out_pts->Allocate(xModel->GetNumberOfBoundaryPoints());

  for(MedialBoundaryPointIterator bip(context); !bip.IsAtEnd(); ++bip)
    {
    // Get the boundary point
    SMLVec3d xBnd = GetBoundaryPoint(bip, xModel->GetAtomArray()).X;

    // Find the closest point on the other mesh
    double xs[3], d2, d;
    int subid;
    vtkIdType cellid;

    xTargetLocator->FindClosestPoint(xBnd.data_block(), xs, cellid, subid, d2);
    xClosestToModel[bip.GetIndex()] = SMLVec3d(xs);

    // Store the point 
    out_pts->InsertNextPoint(xBnd[0], xBnd[1], xBnd[2]);
    }

  // Perform the closest to target computation

  // Output vector
  xClosestToTarget.reserve(xMeshReduced->GetNumberOfPoints());

  // Create the polydata rep of the boundary
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  poly->SetPoints(out_pts);
  poly->Allocate(xModel->GetNumberOfBoundaryTriangles());

  for(MedialBoundaryTriangleIterator bt(context); !bt.IsAtEnd(); ++bt)
    {
    vtkIdType v[] = {bt.GetBoundaryIndex(0), bt.GetBoundaryIndex(1), bt.GetBoundaryIndex(2) } ;
    poly->InsertNextCell(VTK_TRIANGLE, 3, v);
    }

  // Build everything
  poly->BuildCells();

  // Create locator for finding closest points
  vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(poly);
  locator->BuildLocator();

  // Sample points from the target mesh
  for(int i = 0; i < xMeshReduced->GetNumberOfPoints(); i++)
    {
    MatchLocation loc;
    loc.xTarget.set(xMeshReduced->GetPoint(i));

    double xs[3], d2, d;
    int subid;
    vtkIdType cellid;

    // Find the closest point
    locator->FindClosestPoint(loc.xTarget.data_block(), xs, cellid, subid, d2);

    // Solve a system for the barycentric coordinates
    vtkCell *c = poly->GetCell(cellid);
    if(c->GetNumberOfPoints() != 3)
      throw MedialModelException("Bad cell in input");

    SMLVec3d A = SMLVec3d(poly->GetPoint(c->GetPointId(0)));
    SMLVec3d B = SMLVec3d(poly->GetPoint(c->GetPointId(1)));
    SMLVec3d C = SMLVec3d(poly->GetPoint(c->GetPointId(2)));

    vnl_matrix<double> W(3, 2);
    W.set_column(0, B-A);
    W.set_column(1, C-A);
    SMLVec3d v = SMLVec3d(xs) - A;

    vnl_matrix<double> WtW = W.transpose() * W;
    vnl_vector<double> Wtv = W.transpose() * v;
    vnl_vector<double> q = vnl_inverse(WtW) * Wtv;

    // The barycentric coordinates are (1-q1-q2, q1, q2)
    loc.iAtom[0] = context->GetBoundaryPointAtomIndex(c->GetPointId(0));
    loc.iAtom[1] = context->GetBoundaryPointAtomIndex(c->GetPointId(1));
    loc.iAtom[2] = context->GetBoundaryPointAtomIndex(c->GetPointId(2));
    loc.iSide[0] = context->GetBoundaryPointSide(c->GetPointId(0));
    loc.iSide[1] = context->GetBoundaryPointSide(c->GetPointId(1));
    loc.iSide[2] = context->GetBoundaryPointSide(c->GetPointId(2));
    loc.xBary[0] = 1 - (q[0] + q[1]);
    loc.xBary[1] = q[0];
    loc.xBary[2] = q[1];

    xClosestToTarget.push_back(loc);
    }
}


double
SymmetricClosestPointMatchTerm
::UnifiedComputeEnergy(SolutionData *S, bool gradient)
{
  xIterCount++;
  if(xIterCount % 10 == 0)
    FindClosestPoints();

  saDistSqToTarget.Reset();
  saDistSqToModel.Reset();

  double a = 0, b = 0;

  // Compute the distance from the model to the target
  for(MedialBoundaryPointIterator bip(S->xAtomGrid); !bip.IsAtEnd(); ++bip)
    {
    // Compute the distance to closest point
    // Get the medial atom in question
    SMLVec3d &X = GetBoundaryPoint(bip, S->xAtoms).X;

    // Compute the distance squared
    const SMLVec3d &Y = xClosestToModel[bip.GetIndex()];

    // Compute the distance squared
    double d2 = (X-Y).squared_magnitude();
    a += d2;

    // Update the distance squared
    saDistSqToTarget.Update(d2);

    // TODO: scale by BoundaryWeight
    }

  // Compute the distance from the target to the model
  for(int j = 0; j < xClosestToTarget.size(); j++)
    {
    // Compute the weighted point on the triangle
    const MatchLocation &match = xClosestToTarget[j];

    // Atom and side for each locatoin
    const SMLVec3d &X0 = S->xAtoms[match.iAtom[0]].xBnd[match.iSide[0]].X;
    const SMLVec3d &X1 = S->xAtoms[match.iAtom[1]].xBnd[match.iSide[1]].X;
    const SMLVec3d &X2 = S->xAtoms[match.iAtom[2]].xBnd[match.iSide[2]].X;
    SMLVec3d Xw = X0 * match.xBary[0] + X1 * match.xBary[1] + X2 * match.xBary[2];

    // Compute the distance
    double d2 = (Xw - match.xTarget).squared_magnitude();
    b += d2;

    // Update the distance squared
    saDistSqToModel.Update(d2);
    }

  // TODO: modulate properly
  // return saDistSqToTarget.GetMean() + saDistSqToTarget.GetMean();
  //
  return a / xModel->GetNumberOfBoundaryPoints() + 0.1 * b / xClosestToTarget.size();
}


double 
SymmetricClosestPointMatchTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  StatisticsAccumulator saDistSqToTargetDeriv, saDistSqToModelDeriv;
  double da, db;

  // Compute the distance from the model to the target
  for(MedialBoundaryPointIterator bip(xModel->GetIterationContext()); !bip.IsAtEnd(); ++bip)
    {
    // Compute the distance to closest point
    // Get the medial atom in question
    SMLVec3d &X = GetBoundaryPoint(bip, S->xAtoms).X;
    SMLVec3d &dX = GetBoundaryPoint(bip, dS->xAtoms).X;

    // Compute the distance squared
    const SMLVec3d &Y = xClosestToModel[bip.GetIndex()];

    // Compute the distance squared
    // diff: double d2 = (X-Y).squared_magnitude();
    double d_d2 = 2 * dot_product(X-Y, dX);

    // Update the distance squared
    saDistSqToTargetDeriv.Update(d_d2);

    // TODO: scale by BoundaryWeight
    da += d_d2;
    }

  // Compute the distance from the target to the model
  for(int j = 0; j < xClosestToTarget.size(); j++)
    {
    // Compute the weighted point on the triangle
    const MatchLocation &match = xClosestToTarget[j];

    // Atom and side for each locatoin
    const SMLVec3d &X0 = S->xAtoms[match.iAtom[0]].xBnd[match.iSide[0]].X;
    const SMLVec3d &X1 = S->xAtoms[match.iAtom[1]].xBnd[match.iSide[1]].X;
    const SMLVec3d &X2 = S->xAtoms[match.iAtom[2]].xBnd[match.iSide[2]].X;
    const SMLVec3d &dX0 = dS->xAtoms[match.iAtom[0]].xBnd[match.iSide[0]].X;
    const SMLVec3d &dX1 = dS->xAtoms[match.iAtom[1]].xBnd[match.iSide[1]].X;
    const SMLVec3d &dX2 = dS->xAtoms[match.iAtom[2]].xBnd[match.iSide[2]].X;
    SMLVec3d Xw = X0 * match.xBary[0] + X1 * match.xBary[1] + X2 * match.xBary[2];
    SMLVec3d dXw = dX0 * match.xBary[0] + dX1 * match.xBary[1] + dX2 * match.xBary[2];

    // Compute the distance
    // diff: double d2 = (Xw - match.xTarget).squared_magnitude();
    double d_d2 = 2 * dot_product(Xw - match.xTarget, dXw);

    // Update the distance squared
    saDistSqToModelDeriv.Update(d_d2);
    db+=d_d2;
    }

  return da / xModel->GetNumberOfBoundaryPoints() + 0.1 * db / xClosestToTarget.size();
  // return saDistSqToTargetDeriv.GetMean() + saDistSqToModelDeriv.GetMean();
}

void
SymmetricClosestPointMatchTerm
::PrintReport(ostream &sout)
{
  sout << " Symmetric Closest Point Match term " << endl;
  sout << "    total penalty  : " << saDistSqToModel.GetMean() + saDistSqToTarget.GetMean() << endl;
  sout << "    RMS dist model to target: " << sqrt(saDistSqToTarget.GetMean()) << endl;  
  sout << "    RMS dist target to model: " << sqrt(saDistSqToModel.GetMean()) << endl;  
}

SymmetricClosestPointMatchTerm
::~SymmetricClosestPointMatchTerm()
{
}

/*********************************************************************************
 * Distance preservation term
 ********************************************************************************/
LocalDistanceDifferenceEnergyTerm
::LocalDistanceDifferenceEnergyTerm(
  GenericMedialModel *model)
{
  // Store the model
  this->xModel = model;

  // Set the size of arrays
  xRefDistData.resize(model->GetNumberOfTriangles());
  xDistData.resize(model->GetNumberOfTriangles());
}

void
LocalDistanceDifferenceEnergyTerm
::SetParameters(Registry &r)
{
  // Get the array of reference models
  Registry &mdlarray = r.Folder("ReferenceModel");
  if(!mdlarray.GetArraySize())
    throw MedialModelException("Missing reference model array in LocalDistancePenaltyTerm term");

  // Load the reference models
  vector<string> fnRef = mdlarray.GetArray(string(""));
  vector<MedialPDE *> xReference;
  vector<GenericMedialModel *> xReferenceModel;
  for(size_t q = 0; q < fnRef.size(); q++)
    {
    MedialPDE *refModel = new MedialPDE(fnRef[q].c_str());
    
    // Check the model
    if(refModel->GetMedialModel()->GetNumberOfAtoms() != xModel->GetNumberOfAtoms())
      throw MedialModelException("Reference Model Incompatibility in LocalDistanceDET");
    
    xReference.push_back(refModel);
    xReferenceModel.push_back(refModel->GetMedialModel());    
    }

  // Initialize what we need from the reference data
  InitializeReferenceData(xReferenceModel);

  // Delete the reference models
  for(size_t q = 0; q < fnRef.size(); q++)
    delete xReference[q];
}

void
LocalDistanceDifferenceEnergyTerm
::InitializeReferenceData(vector<GenericMedialModel *> &mdlReference)
{
  // Iterate over all boundary triangles
  size_t n = mdlReference.size();
  for(MedialTriangleIterator mit(mdlReference[0]->GetIterationContext()); 
      !mit.IsAtEnd(); ++mit)
    {
    size_t k = mit.GetIndex();
    for(size_t i = 0; i < 3; i++)
      {
      size_t j1 = (i+1) % 3;
      size_t j2 = (i+2) % 3;

      double dm = 0.0;
      double db0 = 0.0;
      double db1 = 0.0;
      double r = 0.0;

      for(size_t z = 0; z < mdlReference.size(); z++)
        {
        MedialAtom *refatoms = mdlReference[z]->GetAtomArray();
        MedialAtom &A0 = refatoms[mit.GetAtomIndex(i)];
        MedialAtom &A1 = refatoms[mit.GetAtomIndex(j1)];
        MedialAtom &A2 = refatoms[mit.GetAtomIndex(j2)];

        dm += (A1.X-A2.X).magnitude();
        db0 += (A1.xBnd[0].X-A2.xBnd[0].X).magnitude();
        db1 += (A1.xBnd[1].X-A2.xBnd[1].X).magnitude();
        r += A0.R;
        }

      xRefDistData[k].dm[i] = dm / n;
      xRefDistData[k].db0[i] = db0 / n;
      xRefDistData[k].db1[i] = db1 / n;
      xRefDistData[k].r[i] = r / n;
      }
    }
}

double
LocalDistanceDifferenceEnergyTerm
::UnifiedComputeEnergy(SolutionData *S, bool gradient_mode)
{
  // Reset the penalty counter
  saDist.Reset();

  // Iterate over all boundary triangles
  for(MedialTriangleIterator mit(S->xAtomGrid); !mit.IsAtEnd(); ++mit)
    {
    size_t k = mit.GetIndex();
    for(size_t i = 0; i < 3; i++)
      {
      size_t j1 = (i+1) % 3;
      size_t j2 = (i+2) % 3;
      MedialAtom &A0 = S->xAtoms[mit.GetAtomIndex(i)];
      MedialAtom &A1 = S->xAtoms[mit.GetAtomIndex(j1)];
      MedialAtom &A2 = S->xAtoms[mit.GetAtomIndex(j2)];

      xDistData[k].dm[i] = (A1.X-A2.X).magnitude();
      xDistData[k].db0[i] = (A1.xBnd[0].X-A2.xBnd[0].X).magnitude();
      xDistData[k].db1[i] = (A1.xBnd[1].X-A2.xBnd[1].X).magnitude();
      xDistData[k].r[i] = A0.R;

      double delm =  xDistData[k].dm[i]  - xRefDistData[k].dm[i];
      double delb0 = xDistData[k].db0[i] - xRefDistData[k].db0[i];
      double delb1 = xDistData[k].db1[i] - xRefDistData[k].db1[i];
      double delr =  xDistData[k].r[i]   - xRefDistData[k].r[i];

      saDist.Update(delm * delm);
      saDist.Update(delb0 * delb0);
      saDist.Update(delb1 * delb1);
      saDist.Update(delr * delr);
      }
    }

  // Report mean square difference
  xTotalPenalty = sqrt(saDist.GetMean());
  return xTotalPenalty;
}

double
LocalDistanceDifferenceEnergyTerm::
ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Initialize the derivative accumulator
  saDistDeriv.Reset();

  // Iterate over all boundary triangles
  for(MedialTriangleIterator mit(S->xAtomGrid); !mit.IsAtEnd(); ++mit)
    {
    size_t k = mit.GetIndex();

    // If all three triangles are non-affected, we can safely set the
    // derivative to zero and contunue
    MedialAtom &D0 = dS->xAtoms[mit.GetAtomIndex(0)];
    MedialAtom &D1 = dS->xAtoms[mit.GetAtomIndex(1)];
    MedialAtom &D2 = dS->xAtoms[mit.GetAtomIndex(2)];
    if(D0.order <= 1 || D1.order <= 1 || D2.order <= 1)
    {
      for(size_t i = 0; i < 3; i++)
        {
        size_t j1 = (i+1) % 3;
        size_t j2 = (i+2) % 3;
        MedialAtom &A1 = S->xAtoms[mit.GetAtomIndex(j1)];
        MedialAtom &A2 = S->xAtoms[mit.GetAtomIndex(j2)];

        MedialAtom &dA0 = dS->xAtoms[mit.GetAtomIndex(i)];
        MedialAtom &dA1 = dS->xAtoms[mit.GetAtomIndex(j1)];
        MedialAtom &dA2 = dS->xAtoms[mit.GetAtomIndex(j2)];

        double d_dm = dot_product(A1.X - A2.X, dA1.X - dA2.X) / xDistData[k].dm[i];
        double d_db0 = dot_product(A1.xBnd[0].X - A2.xBnd[0].X, dA1.xBnd[0].X - dA2.xBnd[0].X) 
          / xDistData[k].db0[i];
        double d_db1 = dot_product(A1.xBnd[1].X - A2.xBnd[1].X, dA1.xBnd[1].X - dA2.xBnd[1].X) 
          / xDistData[k].db1[i];
        double d_r = dA0.R;

        double delm =  xDistData[k].dm[i]  - xRefDistData[k].dm[i];
        double delb0 = xDistData[k].db0[i] - xRefDistData[k].db0[i];
        double delb1 = xDistData[k].db1[i] - xRefDistData[k].db1[i];
        double delr =  xDistData[k].r[i]   - xRefDistData[k].r[i];

        saDistDeriv.Update(delm * d_dm);
        saDistDeriv.Update(delb0 * d_db0);
        saDistDeriv.Update(delb1 * d_db1);
        saDistDeriv.Update(delr * d_r);
        }
      }
    }

  return saDistDeriv.GetSum() / (saDist.GetCount() * xTotalPenalty);
}


void
LocalDistanceDifferenceEnergyTerm
::PrintReport(ostream &sout)
{
  sout << " Local Distance Distortion Image Match Term " << endl;
  sout << "    total penalty  : " << xTotalPenalty << endl;
  sout << "    min sq dist diff : " << saDist.GetMin() << endl;  
  sout << "    max sq dist diff : " << saDist.GetMax() << endl;  
  sout << "    avg sq dist diff : " << saDist.GetMean() << endl;  
}


/*********************************************************************************
 * CROSS-CORRELATION ENERGY TERM
 ********************************************************************************/
CrossCorrelationImageMatchTerm
::CrossCorrelationImageMatchTerm(GenericMedialModel *model, FloatImage *xGrayImage)
{
  this->fTarget = new FloatImageEuclideanFunctionAdapter(xGrayImage);
  this->fReference = NULL;
  this->xModel = model;
}

CrossCorrelationImageMatchTerm
::~CrossCorrelationImageMatchTerm()
{
  delete fTarget;
  if(fReference)
    delete fReference;
}

void 
CrossCorrelationImageMatchTerm::
SetParameters(Registry &r)
{
  // Get the reference image filename
  std::string fnRefImage = r["ReferenceImage"][""];
  std::string fnRefModel = r["ReferenceModel"][""];

  // If the references don't exist, throw an exception
  if(!fnRefModel.length())
    throw MedialModelException("Missing reference model in CrossCorrelation term");

  if(!fnRefImage.length())
    throw MedialModelException("Missing reference image in CrossCorrelation term");

  // Load the reference image
  FloatImage refImage;
  refImage.LoadFromFile(fnRefImage.c_str());

  // Load the reference model
  MedialPDE ref(fnRefModel.c_str());
  GenericMedialModel *xRefModel = ref.GetMedialModel();

  // Make sure the number of points matches
  if(xRefModel->GetNumberOfAtoms() != xModel->GetNumberOfAtoms())
    throw MedialModelException("Model mismatch in CrossCorrelationPenaltyTerm");

  // Get the number of cuts and the extent of the cross-correlation model
  size_t nCuts = r["NumberOfCuts"][16];
  double xiMax = r["MaximumXi"][2.0];

  // Create the sampling functions
  fReference = new FloatImageEuclideanFunctionAdapter(&refImage);

  // Initialize the profile array
  nSamplesPerAtom = nCuts + 2;
  xProfile.resize(
    xModel->GetNumberOfBoundaryPoints(), 
    ProfileData(nSamplesPerAtom));

  // Compute the sample point values along xi
  xSamples.resize(nSamplesPerAtom, 0.0);
  for(size_t i = 0; i < nSamplesPerAtom; i++)
    {
    xSamples[i] = i * xiMax / (nCuts + 1.0);
    }

  // Init the reference data
  InitializeReferenceData(fReference, xRefModel);
}

void
CrossCorrelationImageMatchTerm::
InitializeReferenceData(
  EuclideanFunction *fReference, GenericMedialModel *mdlReference)
{
  // Sample the reference image along the reference model. It is assumed that
  // the reference and target models have exactly the same mesh characteristics.
  for(MedialBoundaryPointIterator bip(mdlReference->GetIterationContext());
    !bip.IsAtEnd(); ++bip)
    {
    size_t ibnd = bip.GetIndex(), iatom = bip.GetAtomIndex();
    ProfileData &p = xProfile[ibnd];

    // Get the medial atom in question
    MedialAtom &a = mdlReference->GetAtomArray()[iatom];

    // Get the vector from medial to the boundary in the reference model
    SMLVec3d U = a.xBnd[bip.GetBoundarySide()].X - a.X;

    // Accumulators for mean and standard deviation
    double xSumSq = 0.0, xSum = 0.0;

    // Compute the volume element and image value for the intermediate points
    for(size_t j = 0; j < nSamplesPerAtom; j++)
      {
      // Compute the sample points
      SMLVec3d Xj = a.X + xSamples[j] * U;

      // Sample the image intentisy
      double v = fReference->Evaluate(Xj);
      p.xRefImgVal[j] = v;

      // Accumulate
      xSumSq += v * v; xSum += v;
      }

    // Compute the mean 
    p.xRefMean = xSum / nSamplesPerAtom;

    // Compute the standard deviation. We don't scale by (n-1) for simplicity
    // since correlation is just a ratio of quantities scaled by (n-1)
    p.xRefSD = sqrt((xSumSq - xSum * p.xRefMean));
    }
}

double
CrossCorrelationImageMatchTerm::
UnifiedComputeEnergy(SolutionData *S, bool gradient_mode)
{
  // Reset the statistical accumulators
  saPenalty.Reset();

  // Sample the target image in the space of the target model.  
  for(MedialBoundaryPointIterator bip(S->xAtomGrid); !bip.IsAtEnd(); ++bip)
    {
    size_t ibnd = bip.GetIndex(), iatom = bip.GetAtomIndex();
    ProfileData &p = xProfile[ibnd];

    // Get the medial atom in question
    MedialAtom &a = S->xAtoms[iatom];

    // Get the vector from medial to the boundary in the reference model
    SMLVec3d U = a.xBnd[bip.GetBoundarySide()].X - a.X;

    // Accumulators for mean and standard deviation
    double xCovAcc = 0.0, xSumSq = 0.0, xSum = 0.0;

    // Compute the volume element and image value for the intermediate points
    for(size_t j = 0; j < nSamplesPerAtom; j++)
      {
      // Compute the sample points
      SMLVec3d Xj = a.X + xSamples[j] * U;

      // Sample the image intentisy at the sample point
      double v;

      // For gradient computations, store the image value and the gradient
      if(gradient_mode)
        {
        v = fTarget->ComputeFunctionAndGradient(Xj, p.xImageGrad[j]);
        p.xImageVal[j] = v;
        }
      else
        {
        v = fTarget->Evaluate(Xj);
        }

      // Compute the contribution to the covariance
      xCovAcc += v * (p.xRefImgVal[j] - p.xRefMean);

      // Compute the contribution to the sum and sum of squares
      xSumSq += v * v;
      xSum += v;

      }

    // Compute the covariance of X and Y
    p.xCov = xCovAcc;
    p.xMean = xSum / nSamplesPerAtom;
    p.xSD = sqrt((xSumSq - xSum * p.xMean));

    // Compute the correlation of X and Y
    p.xCorr = p.xCov / (p.xSD * p.xRefSD);
    saCorr.Update(p.xCorr);

    // Scale the correlation by the area element
    double xContrib = (1.0 - p.xCorr) * S->xBoundaryWeights[ibnd];
    saPenalty.Update(xContrib);
    }

  // Return the mean penalty scaled by the total area
  xTotalPenalty = saPenalty.GetSum() / S->xBoundaryArea;

  if(vnl_math_isnan(xTotalPenalty))
    throw MedialModelException("NAN");

  return xTotalPenalty;
}


double
CrossCorrelationImageMatchTerm::
ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Create an accumulator for total contibution
  double dAccumPenalty = 0.0;

  // Sample the target image in the space of the target model.  
  for(MedialBoundaryPointIterator bip(S->xAtomGrid); !bip.IsAtEnd(); ++bip)
    {
    size_t ibnd = bip.GetIndex(), iatom = bip.GetAtomIndex();
    ProfileData &p = xProfile[ibnd];

    // Get the medial atom in question
    MedialAtom &a = S->xAtoms[iatom];
    MedialAtom &da = dS->xAtoms[iatom];

    // Get the vector from medial to the boundary in the reference model
    SMLVec3d U = a.xBnd[bip.GetBoundarySide()].X - a.X;
    SMLVec3d dU = da.xBnd[bip.GetBoundarySide()].X - da.X;

    // Accumulators for mean and standard deviation
    double dCovAcc = 0.0, dSumSqHalf = 0.0, dSum = 0.0;

    // Compute the volume element and image value for the intermediate points
    for(size_t j = 0; j < nSamplesPerAtom; j++)
      {
      // Compute the sample points
      SMLVec3d Xj = a.X + xSamples[j] * U;
      SMLVec3d Dj = da.X + xSamples[j] * dU;

      // Get the image value and the gradient of the image
      double v = p.xImageVal[j];
      double dv = dot_product(p.xImageGrad[j], Dj);

      // Compute the contribution to the covariance
      dCovAcc += dv * (p.xRefImgVal[j] - p.xRefMean);

      // Compute the contribution to the sum and sum of squares
      dSumSqHalf += v * dv;
      dSum += dv;
      }

    // Compute the covariance of X and Y
    double dCov = dCovAcc;
    double dSD = (dSumSqHalf - dSum * p.xMean) / p.xSD;

    // Compute the correlation of X and Y
    double dCorr = (dCov * p.xSD - p.xCov * dSD) / (p.xSD * p.xSD * p.xRefSD);

    // Scale the correlation by the area element
    double dContrib = 
      (1.0 - p.xCorr) * dS->xBoundaryWeights[ibnd] -
      dCorr * S->xBoundaryWeights[ibnd];

    dAccumPenalty += dContrib;
    }

  // Return the mean penalty scaled by the total area
  double dTotalPenalty = 
    (dAccumPenalty - xTotalPenalty * dS->xBoundaryArea) / S->xBoundaryArea;

  return dTotalPenalty;
}


void
CrossCorrelationImageMatchTerm
::PrintReport(ostream &sout)
{
  sout << " Cross Correlation Image Match Term " << endl;
  sout << "    total penalty  : " << xTotalPenalty << endl;
  sout << "    min corr coeff : " << saCorr.GetMin() << endl;  
  sout << "    max corr coeff : " << saCorr.GetMax() << endl;  
  sout << "    avg corr coeff : " << saCorr.GetMean() << endl;  
}

/*********************************************************************************
 * BOUNDARY JACOBIAN TERM
 ********************************************************************************/

inline double ComputeJacobian(const SMLVec3d &Xu, const SMLVec3d &Xv, 
  const SMLVec3d &Yu, const SMLVec3d &Yv)
{
  return 
    ( dot_product(Yu,Xu) * dot_product(Yv,Xv) - 
      dot_product(Yu,Xv) * dot_product(Yv,Xu) ) / 
    ( dot_product(Xu,Xu) * dot_product(Xv,Xv) - 
      dot_product(Xu,Xv) * dot_product(Xu,Xv) );
}

const double BoundaryJacobianEnergyTerm::xDefaultPenaltyA = 10;
const double BoundaryJacobianEnergyTerm::xDefaultPenaltyB = 20;

BoundaryJacobianEnergyTerm::BoundaryJacobianEnergyTerm()
{
  xPenaltyA = xDefaultPenaltyA;
  xPenaltyB = xDefaultPenaltyB;
}

/**
 * This is a penalty function of the form exp(alpha * (x-x0)) that 
 * has a 'cap' C, such that if exp(alpha * (x-x0)) > C, the function
 * becomes linear rather than exponential.
 */
class ExponentialBarrierFunction
{
public:
  ExponentialBarrierFunction(double alpha, double x0, double C)
    {
    this->a = alpha;
    this->x0 = x0;
    this->C = C;
    this->aC = alpha * C;
    this->b = C * (1 - log(C));
    this->t = x0 + log(C) / alpha;
    }

  double f(double x)
    {
    return x < t ? exp(a * (x - x0)) : aC * (x - x0) + b;
    }

  double df(double x, double f)
    {
    return x < t ? a * f : aC;
    }
private:
  double a, x0, C, aC, b, t;
};

double BoundaryJacobianEnergyTerm::ComputeEnergy(SolutionData *S)
{
  // Place to store the Jacobian
  saJacobian.Reset();

  // Reset the penalty accumulator
  saLower.Reset();
  saUpper.Reset();
  saPenalty.Reset();

  // Create a barrier function
  ExponentialBarrierFunction ebfA(xPenaltyA, 0.0, 100);
  ExponentialBarrierFunction ebfB(1.0, xPenaltyB, 100);

  // Reset the TriangleEntry vector
  if(xTriangleEntries.size() != S->xAtomGrid->GetNumberOfTriangles())
    xTriangleEntries.resize(S->xAtomGrid->GetNumberOfTriangles());

  // Keep track of the entry index
  TriangleVector::iterator eit = xTriangleEntries.begin();
  
  // Create a Triangle-iterator through the atoms
  for(MedialTriangleIterator itt(S->xAtomGrid); !itt.IsAtEnd(); ++itt, ++eit)
    {
    // Get the four atoms in this quad
    MedialAtom &A0 = S->xAtoms[itt.GetAtomIndex(0)];
    MedialAtom &A1 = S->xAtoms[itt.GetAtomIndex(1)];
    MedialAtom &A2 = S->xAtoms[itt.GetAtomIndex(2)];

    // Compute the Xu and Xv vectors and the normal
    eit->XU = A1.X - A0.X; eit->XV = A2.X - A0.X;
    eit->NX = vnl_cross_3d(eit->XU, eit->XV);
    eit->gX2 = dot_product(eit->NX, eit->NX);

    // Compute side-wise entries
    for(size_t z = 0; z < 2; z++)
      {
      // Compute the same for the upper and lower boundaries
      eit->YU[z] = A1.xBnd[z].X - A0.xBnd[z].X;
      eit->YV[z] = A2.xBnd[z].X - A0.xBnd[z].X;
      eit->NY[z] = vnl_cross_3d(eit->YU[z], eit->YV[z]);
      
      // Compute the Jacobian
      eit->J[z] = dot_product(eit->NY[z], eit->NX) / eit->gX2;

      // Add to the average Jacobian
      saJacobian.Update(eit->J[z]);
      
      // Compute the penalty terms
      // return exp(-a * x) + exp(x - b); 
      eit->PenA[z] = ebfA.f(-eit->J[z]);
      eit->PenB[z] = ebfB.f( eit->J[z]);
      
      saLower.Update(eit->PenA[z]);
      saUpper.Update(eit->PenB[z]);
      saPenalty.Update(eit->PenA[z] + eit->PenB[z]);
      }
    }

  // Return the total value
  return saPenalty.GetSum();
}

double BoundaryJacobianEnergyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Place to store the Jacobian
  double dTotalPenalty = 0.0;

  // Make sure that the quad entry array has been initialized
  assert(xTriangleEntries.size() == S->xAtomGrid->GetNumberOfTriangles());
  
  // Create a barrier function
  ExponentialBarrierFunction ebfA(xPenaltyA, 0.0, 100);
  ExponentialBarrierFunction ebfB(1.0, xPenaltyB, 100);

  // Keep track of the entry index
  TriangleVector::iterator eit = xTriangleEntries.begin();
  
  // Create a Triangle-iterator through the atoms
  for(MedialTriangleIterator itt(S->xAtomGrid); !itt.IsAtEnd(); ++itt, ++eit)
    {
    // Get the derivative atoms too
    MedialAtom &dA0 = dS->xAtoms[itt.GetAtomIndex(0)];
    MedialAtom &dA1 = dS->xAtoms[itt.GetAtomIndex(1)];
    MedialAtom &dA2 = dS->xAtoms[itt.GetAtomIndex(2)];

    // If all three triangles are non-affected, we can safely set the
    // derivative to zero and contunue
    if(dA0.order <= 1 || dA1.order <= 1 || dA2.order <= 1)
      {
      // Compute the average Xu and Xv vectors and derivatives
      SMLVec3d dXU = dA1.X - dA0.X; SMLVec3d dXV = dA2.X - dA0.X;
      SMLVec3d dNX = vnl_cross_3d(dXU,  eit->XV) + vnl_cross_3d(eit->XU,  dXV);

      // Compute G and its derivative
      double dgX2 = 2.0 * dot_product(dNX, eit->NX);

      // Compute side-wise derivatives
      for(size_t z = 0; z < 2; z++)
        {
        // Compute boundary vector derivatives
        SMLVec3d dYU = dA1.xBnd[z].X - dA0.xBnd[z].X;
        SMLVec3d dYV = dA2.xBnd[z].X - dA0.xBnd[z].X;
        SMLVec3d dNY = vnl_cross_3d(dYU, eit->YV[z]) + vnl_cross_3d(eit->YU[z], dYV);

        // Compute the Jacobian derivative
        double dJ = (
          dot_product(dNY, eit->NX) + 
          dot_product(eit->NY[z], dNX) - eit->J[z] * dgX2) / eit->gX2;

        // Compute the penalty terms
        dTotalPenalty += dJ * (
          - ebfA.df(-eit->J[z], eit->PenA[z]) 
          + ebfB.df( eit->J[z], eit->PenB[z]));
        }
      }
    }

  // Return the total value
  return dTotalPenalty;
}
    
// Print a verbose report
void 
BoundaryJacobianEnergyTerm
::PrintReport(ostream &sout)
{
  sout << "  Boundary Jacobian Term " << endl;
  sout << "    total penalty      : " << saPenalty.GetSum() << endl;
  sout << "    min jacobian       : " << saJacobian.GetMin() << endl;  
  sout << "    max jacobian       : " << saJacobian.GetMax() << endl;  
  sout << "    avg jacobian       : " << saJacobian.GetMean() << endl; 
  sout << "    min lower penalty  : " << saLower.GetMin() << endl;
  sout << "    max lower penalty  : " << saLower.GetMax() << endl;
  sout << "    avg lower penalty  : " << saLower.GetMean() << endl;
  sout << "    min upper penalty  : " << saUpper.GetMin() << endl;
  sout << "    max upper penalty  : " << saUpper.GetMax() << endl;
  sout << "    avg upper penalty  : " << saUpper.GetMean() << endl;
}

VolumeIntegralEnergyTerm
::VolumeIntegralEnergyTerm(
  GenericMedialModel *model, EuclideanFunction *function, size_t nCuts)
{
  // Store the inputs
  this->nCuts = nCuts;
  this->nAtoms = model->GetNumberOfAtoms();
  this->function = function;

  // Allocate an array to store sampled image values
  nSamplesPerAtom = nCuts + 2;
  xProfile.resize(
    model->GetNumberOfBoundaryPoints(), 
    ProfileData(nSamplesPerAtom));

  double delta = 1.0 / (1.0 + nCuts);
  double hd = 0.5 * delta;

  // Compute the xi values associated with the samples
  xSamples.resize(nSamplesPerAtom);
  for(size_t i = 0; i < nSamplesPerAtom; i++)
    xSamples[i] = i / (nCuts + 1.0);

  // Compute the coefficients for volume element computation at each
  // depth level (xi)
  xSampleCoeff.resize(nSamplesPerAtom);
  xSampleCoeff[0] = SMLVec3d(hd, hd*hd, hd*hd*hd);
  xSampleCoeff[nCuts+1] = SMLVec3d(hd, hd * (1-hd), hd*(1-hd)*(1-hd));
  for(size_t i = 1; i <= nCuts; i++)
    xSampleCoeff[i] = SMLVec3d(delta, delta * xSamples[i], delta * 
      (xSamples[i] * xSamples[i] + hd * hd));
}

double VolumeIntegralEnergyTerm
::UnifiedComputeEnergy(SolutionData *S, bool gradient_mode)
{
  // Clear the intersection volume accumulator
  xObjectIntegral = 0;
  xVolumeIntegral = 0;

  saImage.Reset();

  // Iterate over the medial atoms in the model
  for(MedialBoundaryPointIterator bip(S->xAtomGrid); !bip.IsAtEnd(); ++bip)
    {
    size_t ibnd = bip.GetIndex(), iatom = bip.GetAtomIndex();
    SMLVec3d &vvec = S->xInteriorVolumeElement[ibnd];
    MedialAtom &a = S->xAtoms[iatom];
    ProfileData &p = xProfile[ibnd];

    // Get the vector from medial to the boundary
    SMLVec3d U = a.xBnd[bip.GetBoundarySide()].X - a.X;

    // Compute the volume element and image value for the intermediate points
    for(size_t j = 0; j < nSamplesPerAtom; j++)
      {
      // Compute the volume element
      p.xVolumeElt[j] = dot_product(vvec, xSampleCoeff[j]);

      // Compute the sample points
      SMLVec3d Xj = a.X + xSamples[j] * U;

      // Sample the image intentisy
      p.xImageVal[j] = function->Evaluate(Xj);
      saImage.Update(p.xImageVal[j]);

      // Compute the image gradient for these sites
      if(gradient_mode)
        function->ComputeGradient(Xj, p.xImageGrad[j]);


      // printf("Spoke %04d:%04d, Sample %02d: VE=%8.4f\t IM=%8.4f\t X=%f, %f, %f\n",
      //        ibnd, iatom, j, p.xVolumeElt[j], p.xImageVal[j], Xj[0], Xj[1], Xj[2]);

      // Compute the contribution to the total
      xVolumeIntegral += p.xVolumeElt[j];
      xObjectIntegral += p.xVolumeElt[j] * p.xImageVal[j];
      }
    }

  // Compute an estimate of volume overlap
  return xObjectIntegral;
}

double VolumeIntegralEnergyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Clear the intersection volume accumulator
  dObjectIntegral = 0;
  dVolumeIntegral = 0;

  // Iterate over the medial atoms in the model
  for(MedialBoundaryPointIterator bip(S->xAtomGrid); !bip.IsAtEnd(); ++bip)
    {
    size_t iatom = bip.GetAtomIndex();
    MedialAtom &da = dS->xAtoms[iatom];

    // Check dependency
    if(da.order <= 2)
      {
      size_t ibnd = bip.GetIndex();

      SMLVec3d &dvvec = dS->xInteriorVolumeElement[ibnd];
      ProfileData &p = xProfile[ibnd];

      // Get the vector from medial to the boundary
      SMLVec3d dU = da.xBnd[bip.GetBoundarySide()].X - da.X;

      // Compute the volume element and image value for the intermediate points
      for(size_t j = 0; j < nSamplesPerAtom; j++)
        {
        // Compute the volume element
        double dVolumeElt = dot_product(dvvec, xSampleCoeff[j]);

        // Compute the sample points
        SMLVec3d DXj = da.X + xSamples[j] * dU;

        // Sample the image intentisty
        dVolumeIntegral += dVolumeElt;
        dObjectIntegral +=
          dVolumeElt * p.xImageVal[j] + 
          p.xVolumeElt[j] * dot_product(DXj, p.xImageGrad[j]);
        }
      }
    }

  // Compute an estimate of volume overlap
  return dObjectIntegral;
}

void VolumeIntegralEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Volume Integral Energy Term: " << endl;
  sout << "    object integral: " << xObjectIntegral << endl;
  sout << "    object volume  : " << xVolumeIntegral << endl;
  sout << "    integral/vol   : " << xObjectIntegral/xVolumeIntegral << endl;
  sout << "    sample pointsx : " << saImage.GetCount() << endl;
  sout << "    sample min/max : " << saImage.GetMin() << "; " << saImage.GetMax() << endl;
  sout << "    sample mean/sd : " << saImage.GetMean() << " +- " << saImage.GetStdDev() << endl;
}

/*********************************************************************************
 * VOLUME OVERLAP IMAGE MATCH TERM
 ********************************************************************************/
VolumeOverlapEnergyTerm
::VolumeOverlapEnergyTerm(
  GenericMedialModel *model, FloatImage *xImage, size_t nCuts)
{
  // Initialize the worker
  function = new FloatImageEuclideanFunctionAdapter(xImage);
  worker = new VolumeIntegralEnergyTerm(model, function, nCuts);

  // Compute the image volume (it remains a constant throughout)
  xImageIntegral = xImage->IntegratePositiveVoxels();
}

double VolumeOverlapEnergyTerm
::ComputeEnergy(SolutionData *S)
{
  double xObjectIntegral = worker->ComputeEnergy(S);
  return xRatio = 1.0 - xObjectIntegral / xImageIntegral;
}

double VolumeOverlapEnergyTerm
::BeginGradientComputation(SolutionData *S)
{
  double xObjectIntegral = worker->BeginGradientComputation(S);
  return xRatio = 1.0 - xObjectIntegral / xImageIntegral;
}

double VolumeOverlapEnergyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dObjectIntegral = worker->ComputePartialDerivative(S, dS);
  return - dObjectIntegral / xImageIntegral;
}

void VolumeOverlapEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Volume Overlap Energy Term: " << endl;
  sout << "    object integral: " << worker->xObjectIntegral << endl;
  sout << "    object volume  : " << worker->xVolumeIntegral << endl;
  sout << "    image integral : " << xImageIntegral << endl;
  sout << "    ratio          : " << worker->xObjectIntegral / xImageIntegral << endl;
  sout << "    final value    : " << xRatio << endl;
  worker->PrintReport(sout);
}

/*********************************************************************************
 * ProbabilityIntegralEnergyTerm
 ********************************************************************************/
ProbabilityIntegralEnergyTerm
::ProbabilityIntegralEnergyTerm(
  GenericMedialModel *model, FloatImage *xImage, size_t nCuts)
{
  // Initialize the worker
  function = new FloatImageEuclideanFunctionAdapter(xImage);
  worker = new VolumeIntegralEnergyTerm(model, function, nCuts);

  // Compute the image volume (it remains a constant throughout)
  xIntPOverImage = xImage->IntegratePositiveVoxels();
  xInt1OverImage = xImage->ComputeImageVolume();
}

double ProbabilityIntegralEnergyTerm
::ComputeEnergy(SolutionData *S)
{
  worker->ComputeEnergy(S);
  double xIntPOverModel = worker->xObjectIntegral;
  double xInt1OverModel = worker->xVolumeIntegral;
  xIntP = 2*xIntPOverModel + xInt1OverImage - (xIntPOverImage + xInt1OverModel);
  xResult = 1.0 - xIntP / xInt1OverImage;
  return xResult;
}

double ProbabilityIntegralEnergyTerm
::BeginGradientComputation(SolutionData *S)
{
  worker->BeginGradientComputation(S);
  double xIntPOverModel = worker->xObjectIntegral;
  double xInt1OverModel = worker->xVolumeIntegral;
  xIntP = 2*xIntPOverModel + xInt1OverImage - (xIntPOverImage + xInt1OverModel);
  xResult = 1.0 - xIntP / xInt1OverImage;
  return xResult;
}

double ProbabilityIntegralEnergyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  worker->ComputePartialDerivative(S, dS);
  double dIntPOverModel = worker->dObjectIntegral;
  double dInt1OverModel = worker->dVolumeIntegral;
  double dIntP = 2*dIntPOverModel  - dInt1OverModel;
  double dResult = -dIntP / xInt1OverImage;
  return dResult;
}

void ProbabilityIntegralEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Probability Integral Energy Term: " << endl;
  sout << "    integral of p(x) over model: " << worker->xObjectIntegral << endl;
  sout << "    integral of 1    over model: " << worker->xVolumeIntegral << endl;
  sout << "    integral of p(x) over image: " << xIntPOverImage << endl;
  sout << "    integral of 1    over image: " << xInt1OverImage << endl;
  sout << "    function to be maximized   : " << xIntP << endl;
  sout << "    final (normalized) value   : " << xResult << endl;
}

/*********************************************************************************
 * BoundaryGradRPenaltyTerm
 ********************************************************************************/
const double BoundaryGradRPenaltyTerm::xScale = 10.0;

double 
BoundaryGradRPenaltyTerm
::ComputeEnergy(SolutionData *S)
{
  // Reset the stats arrays
  saGradR.Reset();
  saGradRInt.Reset();
  saPenalty.Reset();
  
  // Iterate over all crest atoms
  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    MedialAtom &a = S->xAtoms[i];

    // First the penalty that penalizes atoms at the edge for having a gradR that
    // is not equal to 1. This is only in effect for BruteForce models. For the 
    // PDE model this should always give a zero penalty
    if(a.flagCrest)
      {
      // Get the badness of the atom. At boundary atoms, the badness
      // is a function of how far |gradR| is from zero. We can set the
      // penalty to have the form alpha * (1 - |gradR|^2)^2
      // double devn = (1.0 - S->xAtoms[i].xGradRMagSqr);
      // double penalty = (xScale * devn * devn);
      double penalty = 0.0001 * (1.0 / a.xGradRMagSqrOrig);

      // Register the gradR
      saGradR.Update(a.xGradRMagSqrOrig);
      saPenalty.Update(penalty); 
      }

    // Next the penalty term that penalizes gradR for being close to one at internal
    // atoms, and that penalizes dRdS for being close to one at edge atoms. The same
    // penalty structure applies in both cases.
    double x = (a.flagCrest) ? a.Rs2 : a.xGradRMagSqrOrig;

    // For atoms that aren't at the medial edge, we want to keep gradR
    // from being anywhere near 1, since that will make the derivatives
    // of the boundary nodes go to infinity, killing the optimizer. So 
    // we let the penalty be of the form alpha / (1-gradR^2)^2
    //
    // Why this particular shape of the penalty? Because the penalty for
    // |gradR|^2 = 0.95 is 60 times greater than for |gradR|^2 = 0.5 and
    // |gradR|^2 = 0.99 is 1400 times greater than for |gradR|^2 = 0.5
    // which I feel is an acceptably steep polynomial penalty
    saGradRInt.Update(x);

    // p_ref is the reference penalty for |gradR|^2 = 0.5
    const double p_ref = 0.5625;
    double t = (1.0 - x);
    double t2 = t * t;
    double penalty = 0.0001 * (p_ref / t2);
    saPenalty.Update(penalty);
    }

  // Return the mean penalty
  return saPenalty.GetMean();
}  

double 
BoundaryGradRPenaltyTerm::
ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Initialize the accumulators
  double dTotalPenalty = 0.0;
  
  // Iterate over all crest atoms
  size_t nAtoms = S->xAtomGrid->GetNumberOfAtoms();
  for(size_t i = 0; i < nAtoms; i++)
    {
    MedialAtom &a = S->xAtoms[i], &da = dS->xAtoms[i];

    if(da.order <= 1)
      {
      // The part for not being 1 at edge
      if(a.flagCrest)
        {
        // double devn = (1.0 - S->xAtoms[i].xGradRMagSqr);
        // double ddevn = - dS->xAtoms[i].xGradRMagSqr;
        // double d_penalty = 2 * xScale * devn * ddevn; 
        // dTotalPenalty += d_penalty;
        double d_penalty = 0.0001 * 
          (-da.xGradRMagSqrOrig / (a.xGradRMagSqrOrig * a.xGradRMagSqrOrig));
        dTotalPenalty += d_penalty;
        }

      // The part for being close to 1 on the interior or for the along-edge
      // component being close to 1 at edge
      double x = (a.flagCrest) ? a.Rs2 : a.xGradRMagSqrOrig;
      double dx = (a.flagCrest) ? da.Rs2 : da.xGradRMagSqrOrig;

      // p_ref is the reference penalty for |gradR|^2 = 0.5
      const double p_ref = 0.5625;
      double t = (1.0 - x);
      double dt = - dx;
      double t2 = t * t;
      double dt2 = 2.0 * t * dt;
      double penalty = 0.0001 * (p_ref / t2);
      double d_penalty = penalty * (- dt2 / t2); 
      dTotalPenalty += d_penalty;
      }
    }

  // dTotalPenalty /= nAtoms;
  return dTotalPenalty / saPenalty.GetCount();
}

void BoundaryGradRPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Boundary |gradR|^2 Penalty: " << endl;
  sout << "    number of boundary atoms          : " << saGradR.GetCount() << endl; 
  sout << "    range of boundary |gradR|^2       : " <<
    saGradR.GetMin() << " to " << saGradR.GetMax() << endl;
  sout << "    mean (std) of boundary |gradR|^2  : " <<
    saGradR.GetMean() << " (" << saGradR.GetStdDev() << ")" << endl;
  sout << "    number of internal atoms          : " << saGradRInt.GetCount() << endl; 
  sout << "    range of internal |gradR|^2       : " <<
    saGradRInt.GetMin() << " to " << saGradRInt.GetMax() << endl;
  sout << "    mean (std) of internal |gradR|^2  : " <<
    saGradRInt.GetMean() << " (" << saGradRInt.GetStdDev() << ")" << endl;

  sout << "    total penalty            : " << saPenalty.GetMean() << endl;
}



/*********************************************************************************
 * Atom Badness Penalty term
 ********************************************************************************/
double AtomBadnessTerm::ComputeEnergy(SolutionData *S)
{
  // This term computes the irregularity of medial atoms. The following are examples
  // of irregular atom conditions that must be penalized. These penalties are quite
  // relative and should be set on case-by-case basis
  // 
  //   -  Internal atom has |gradR| > 1 - eps1
  //   -  Boundary atom has |gradR| < 1 - eps2

  // Initialize the penalty array
  nAtoms = S->xAtomGrid->GetNumberOfAtoms();
  xPenalty.set_size(nAtoms);
  
  // Initialize the accumulators
  xMinBadness = 1.0;
  xTotalPenalty = xAvgBadness = 0.0;
  nBadAtoms = 0;
  
  // Iterate over all crest atoms
  for(size_t i = 0; i < nAtoms; i++)
    {
    // Penalize internal atoms where gradR is excessive
    if(!S->xAtoms[i].flagCrest)
      {
      double badness = 1.0 - S->xAtoms[i].xGradRMagSqr;
      xPenalty[i] = exp(300 * - (badness - 0.05));
      xTotalPenalty += xPenalty[i];
      xMinBadness = std::min(badness, xMinBadness);
      }

    // Penalize boundary atoms where gradR is far from 1
    else
      {
      // Get the badmess of the atom. At boundary atoms, the badness
      // is a function of how far |gradR| is from zero. We can set the
      // penalty to have the form alpha * (1 - |gradR|^2)^2
      double devn = (1.0 - S->xAtoms[i].xGradRMagSqr);
      xPenalty[i] = (10 * devn * devn);
      xTotalPenalty += xPenalty[i];
      }
    }

  // Finish up
  // xTotalPenalty /= nAtoms;

  // Return the total penalty
  return xTotalPenalty;
}  

double AtomBadnessTerm::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Initialize the accumulators
  double dTotalPenalty = 0.0;
  
  // Iterate over all crest atoms
  nAtoms = S->xAtomGrid->GetNumberOfAtoms();
  for(size_t i = 0; i < nAtoms; i++)
    {
    if(!S->xAtoms[i].flagCrest)
      {
      double d_badness = -dS->xAtoms[i].xGradRMagSqr;
      double d_penalty = xPenalty[i] * (-100 * d_badness);
      dTotalPenalty += d_penalty;
      }
    else
      {
      double devn = (1.0 - S->xAtoms[i].xGradRMagSqr);
      double ddevn = - dS->xAtoms[i].xGradRMagSqr;
      double d_penalty = 2 * 10 * devn * ddevn; 
      dTotalPenalty += d_penalty;
      }
    }

  // dTotalPenalty /= nAtoms;
  return dTotalPenalty;
}

void AtomBadnessTerm::PrintReport(ostream &sout)
{
  sout << "  Atom Penalty Term:           " << endl;
  sout << "    number of atoms          : " << nAtoms << endl; 
  sout << "    number of bad atoms      : " << nBadAtoms << endl; 
  sout << "    largest badness value    : " << xMinBadness << endl; 
  sout << "    average badness value    : " << xAvgBadness << endl; 
  sout << "    total penalty            : " << xTotalPenalty << endl;
}

/********************************************************************************
 * Penalty based on the minimum angle of the triangles
 *******************************************************************************/
//
// Compute squared cosines of four angles in a quad
double CosineSquareTuple(MedialAtom *A, MedialTriangleIterator &it)
{
  // Get the four vectors
  const SMLVec3d &X0 = A[it.GetAtomIndex(0)].X;
  const SMLVec3d &X1 = A[it.GetAtomIndex(1)].X;
  const SMLVec3d &X2 = A[it.GetAtomIndex(2)].X;

  // Compute the differences
  SMLVec3d D10 = X1 - X0;
  SMLVec3d D21 = X2 - X1;
  SMLVec3d D02 = X0 - X2;

  // Compute the lengths squared
  double L10 = dot_product(D10, D10);
  double L21 = dot_product(D21, D21);
  double L02 = dot_product(D02, D02);

  // Compute the cosines squared
  double A0 = dot_product(D10, D02);
  double A1 = dot_product(D21, D10);
  double A2 = dot_product(D02, D21);

  // Compute the actual values
  double C0 = (A0 * A0) / (L10 * L02);
  double C1 = (A1 * A1) / (L21 * L10);
  double C2 = (A2 * A2) / (L02 * L21);

  // Compute the weighted sum
  return C0 + C1 + C2;
}

// Compute squared cosines of four angles in a quad
double CosineSquareTupleDerivative(MedialAtom *A, MedialAtom *dA, MedialTriangleIterator &it)
{
  size_t i0 = it.GetAtomIndex(0);
  size_t i1 = it.GetAtomIndex(1);
  size_t i2 = it.GetAtomIndex(2);

  // Compute the differences and their derivatives
  SMLVec3d D10 = A[i1].X - A[i0].X;
  SMLVec3d D21 = A[i2].X - A[i1].X;
  SMLVec3d D02 = A[i0].X - A[i2].X;

  SMLVec3d dD10 = dA[i1].X - dA[i0].X;
  SMLVec3d dD21 = dA[i2].X - dA[i1].X;
  SMLVec3d dD02 = dA[i0].X - dA[i2].X;

  // Compute the lengths squared
  double L10 = dot_product(D10, D10);
  double L21 = dot_product(D21, D21);
  double L02 = dot_product(D02, D02);

  double dL10 = 2.0 * dot_product(D10, dD10);
  double dL21 = 2.0 * dot_product(D21, dD21);
  double dL02 = 2.0 * dot_product(D02, dD02);

  // Compute the cosines squared
  double A0 = dot_product(D10, D02);
  double A1 = dot_product(D21, D10);
  double A2 = dot_product(D02, D21);

  double dA0 = dot_product(D10, dD02) + dot_product(dD10, D02);
  double dA1 = dot_product(D21, dD10) + dot_product(dD21, D10);
  double dA2 = dot_product(D02, dD21) + dot_product(dD02, D21);

  // Compute the derivatives of the actual values
  double C0 = (A0 * A0) / (L10 * L02);
  double C1 = (A1 * A1) / (L21 * L10);
  double C2 = (A2 * A2) / (L02 * L21);

  // (a/b)' = (a'b - ab')/b^2 = (a' - (a/b)b')/b
  double dC0 = (2.0 * A0 * dA0 - C0 * (L10 * dL02 + dL10 * L02)) / (L10 * L02);
  double dC1 = (2.0 * A1 * dA1 - C1 * (L21 * dL10 + dL21 * L10)) / (L21 * L10);
  double dC2 = (2.0 * A2 * dA2 - C2 * (L02 * dL21 + dL02 * L21)) / (L02 * L21);

  // Compute the weighted sum
  return dC0 + dC1 + dC2;
}

MedialTriangleAnglePenaltyTerm
::MedialTriangleAnglePenaltyTerm(GenericMedialModel *model)
{

}

double 
MedialTriangleAnglePenaltyTerm
::UnifiedComputeEnergy(SolutionData *S, bool flagGradient)
{
  // Reset the penalty terms
  sCosSquare.Reset(); sPenalty.Reset();

  // Loop over all medial triangles
  for(MedialTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    for(size_t v = 0; v < 3; v++)
      {
      // Get the vertices around angle
      const SMLVec3d &X0 = S->xAtoms[it.GetAtomIndex(v)].X;
      const SMLVec3d &X1 = S->xAtoms[it.GetAtomIndex((v+1) % 3)].X;
      const SMLVec3d &X2 = S->xAtoms[it.GetAtomIndex((v+2) % 3)].X;

      // Get the edges
      SMLVec3d D1 = X1 - X0, D2 = X2 - X0;

      // Get the squared lengths
      double L1 = dot_product(D1, D1);
      double L2 = dot_product(D2, D2);

      // Get the dot product
      double Z = dot_product(D1, D2);

      // Get the squared cosine
      double csq = (Z / L1) * (Z / L2);

      // Generate penalty
      double penalty = 0.001 / (1 - csq);

      // Update penalties
      sCosSquare.Update(csq);
      sPenalty.Update(penalty);
      }
    }

  return sPenalty.GetMean();
}

double
MedialTriangleAnglePenaltyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Reset derivative accumulator
  sDPenalty.Reset();

  // Loop over all medial triangles
  for(MedialTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    for(size_t v = 0; v < 3; v++)
      {
      // Get the vertices around angle
      const SMLVec3d &X0 = S->xAtoms[it.GetAtomIndex(v)].X;
      const SMLVec3d &X1 = S->xAtoms[it.GetAtomIndex((v+1) % 3)].X;
      const SMLVec3d &X2 = S->xAtoms[it.GetAtomIndex((v+2) % 3)].X;
      const SMLVec3d &dX0 = dS->xAtoms[it.GetAtomIndex(v)].X;
      const SMLVec3d &dX1 = dS->xAtoms[it.GetAtomIndex((v+1) % 3)].X;
      const SMLVec3d &dX2 = dS->xAtoms[it.GetAtomIndex((v+2) % 3)].X;

      // Get the edges
      SMLVec3d D1 = X1 - X0, D2 = X2 - X0;
      SMLVec3d dD1 = dX1 - dX0, dD2 = dX2 - dX0;

      // Get the squared lengths
      double L1 = dot_product(D1, D1);
      double L2 = dot_product(D2, D2);
      double dL1 = 2.0 * dot_product(D1, dD1);
      double dL2 = 2.0 * dot_product(D2, dD2);

      // Get the dot product
      double Z = dot_product(D1, D2);
      double dZ = dot_product(dD1, D2) + dot_product(D1, dD2);

      // Get the squared cosine
      double csq = (Z / L1) * (Z / L2);
      double dcsq = 
        (Z / L1) * ((dZ * L2 - Z * dL2) / (L2 * L2)) +
        ((dZ * L1 - Z * dL1) / (L1 * L1)) * (Z / L2);

      // Generate penalty
      double penalty = 0.001 / (1 - csq);
      double dpenalty = penalty * (dcsq / (1. - csq));

      // Update penalties
      sDPenalty.Update(dpenalty);
      }
    }

  return sDPenalty.GetMean();
}

void 
MedialTriangleAnglePenaltyTerm
::PrintReport(ostream &sout)
{
  sout << "  Medial Triangle Angle Penalty Term : " << endl;
  sout << "    min angle                : " << acos(sqrt(sCosSquare.GetMin())) << endl;
  sout << "    max angle                : " << acos(sqrt(sCosSquare.GetMax())) << endl;
  sout << "    avg cos sq.              : " << sCosSquare.GetMean() << endl;
  sout << "    max penalty              : " << sPenalty.GetMax() << endl;
  sout << "    mean penalty (r.v.)      : " << sPenalty.GetMean() << endl;
}

/*********************************************************************************
 * Boundary Triangle Angle penalty term
*********************************************************************************/

BoundaryTriangleAnglePenaltyTerm
::BoundaryTriangleAnglePenaltyTerm(GenericMedialModel *model)
{

}

double 
BoundaryTriangleAnglePenaltyTerm
::UnifiedComputeEnergy(SolutionData *S, bool flagGradient)
{
  // Reset the penalty terms
  sCosSquare.Reset(); sPenalty.Reset();

  // Loop over all medial triangles
  for(MedialBoundaryTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    for(size_t v = 0; v < 3; v++)
      {
      // Get the vertices around angle
      const SMLVec3d &X0 = GetBoundaryPoint(it, S->xAtoms, (v+0) % 3).X; 
      const SMLVec3d &X1 = GetBoundaryPoint(it, S->xAtoms, (v+1) % 3).X; 
      const SMLVec3d &X2 = GetBoundaryPoint(it, S->xAtoms, (v+2) % 3).X; 

      // Get the edges
      SMLVec3d D1 = X1 - X0, D2 = X2 - X0;

      // Get the squared lengths
      double L1 = dot_product(D1, D1);
      double L2 = dot_product(D2, D2);

      // Get the dot product
      double Z = dot_product(D1, D2);

      // Get the squared cosine
      double csq = (Z / L1) * (Z / L2);

      // Generate penalty
      double penalty = 0.001 / (1 - csq);

      // Update penalties
      sCosSquare.Update(csq);
      sPenalty.Update(penalty);
      }
    }

  return sPenalty.GetMean();
}

double
BoundaryTriangleAnglePenaltyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Reset derivative accumulator
  sDPenalty.Reset();

  // Loop over all medial triangles
  for(MedialBoundaryTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    for(size_t v = 0; v < 3; v++)
      {
      // Get the vertices around angle
      const SMLVec3d &X0 = GetBoundaryPoint(it, S->xAtoms, (v+0) % 3).X; 
      const SMLVec3d &X1 = GetBoundaryPoint(it, S->xAtoms, (v+1) % 3).X; 
      const SMLVec3d &X2 = GetBoundaryPoint(it, S->xAtoms, (v+2) % 3).X; 

      const SMLVec3d &dX0 = GetBoundaryPoint(it, dS->xAtoms, (v+0) % 3).X; 
      const SMLVec3d &dX1 = GetBoundaryPoint(it, dS->xAtoms, (v+1) % 3).X; 
      const SMLVec3d &dX2 = GetBoundaryPoint(it, dS->xAtoms, (v+2) % 3).X; 

      // Get the edges
      SMLVec3d D1 = X1 - X0, D2 = X2 - X0;
      SMLVec3d dD1 = dX1 - dX0, dD2 = dX2 - dX0;

      // Get the squared lengths
      double L1 = dot_product(D1, D1);
      double L2 = dot_product(D2, D2);
      double dL1 = 2.0 * dot_product(D1, dD1);
      double dL2 = 2.0 * dot_product(D2, dD2);

      // Get the dot product
      double Z = dot_product(D1, D2);
      double dZ = dot_product(dD1, D2) + dot_product(D1, dD2);

      // Get the squared cosine
      double csq = (Z / L1) * (Z / L2);
      double dcsq = 
        (Z / L1) * ((dZ * L2 - Z * dL2) / (L2 * L2)) +
        ((dZ * L1 - Z * dL1) / (L1 * L1)) * (Z / L2);

      // Generate penalty
      double penalty = 0.001 / (1 - csq);
      double dpenalty = penalty * (dcsq / (1. - csq));

      // Update penalties
      sDPenalty.Update(dpenalty);
      }
    }

  return sDPenalty.GetMean();
}

void 
BoundaryTriangleAnglePenaltyTerm
::PrintReport(ostream &sout)
{
  sout << "  Boundary Triangle Angle Penalty Term : " << endl;
  sout << "    min angle                : " << acos(sqrt(sCosSquare.GetMin())) << endl;
  sout << "    max angle                : " << acos(sqrt(sCosSquare.GetMax())) << endl;
  sout << "    avg cos sq.              : " << sCosSquare.GetMean() << endl;
  sout << "    max penalty              : " << sPenalty.GetMax() << endl;
  sout << "    mean penalty (r.v.)      : " << sPenalty.GetMean() << endl;
}


/*********************************************************************************
 * LoopTangentSchemeValidityPenaltyTerm
 ********************************************************************************/

LoopTangentSchemeValidityPenaltyTerm
::LoopTangentSchemeValidityPenaltyTerm(GenericMedialModel *model)
: MedialIntegrationEnergyTerm(model)
{
}

/** 
 * The penalty is of the form c / sin(alpha)^2
 * where c = sin(20)^2, and alpha is the angle between the Xu and Xv vectors
 * This way, the penalty for angle of 20 degrees is 1, for 1 degree is 384,
 * for 0.1 degree is 34800. The idea is that everything in the ballpark of
 * 20 degrees is acceptable, but we don't want smaller angles.
 *
 * sin(alpha)^2 = 1 - G12^2 / (G11 * G22)
 *
 * We need this penalty because the LoopTangentScheme, which is based on the
 * 1994 Hoppe paper "Piecewise Smooth Surface Reconstruction" can create some
 * degenerate (parallel) tangent vectors for boundary vertices with valence > 4
 */

double LoopTangentSchemeValidityPenaltyTerm::ComputeEnergy(SolutionData *S)
{
  // Reset the statistics accumulator
  saPenalty.Reset();

  // We use 20 degrees as the baseline (at which the per-vertex penalty is 2.0)
  static const double c = pow(sin(vnl_math::pi / 9.0), 2);

  // We will simply look at the ratio |Xu x Xv| to |Xu| |Xv| and
  // take its inverse as the penalty, with a very small weight
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    double g11 = a.G.xCovariantTensor[0][0];
    double g22 = a.G.xCovariantTensor[1][1];
    double g12 = a.G.xCovariantTensor[0][1];

    double g11_times_g12 = g11 * g22;
    double sin_alpha_sqr = (g11_times_g12 - g12 * g12) / g11_times_g12;
    double xpen = c / sin_alpha_sqr;

    saPenalty.Update(xpen);
    }

  // return 0.25 * xTotalPenalty / S->xAtomGrid->GetNumberOfAtoms();
  return saPenalty.GetMean();
}

double LoopTangentSchemeValidityPenaltyTerm::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  // We use 20 degrees as the baseline (at which the per-vertex penalty is 2.0)
  static const double c = pow(sin(vnl_math::pi / 9.0), 2);

  // Compute the square of the cosine of the angle between Xu and Xv
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &da = dS->xAtoms[it.GetIndex()];
    if(da.order <= 2)
      {
      MedialAtom &a = S->xAtoms[it.GetIndex()];
      double g11 = a.G.xCovariantTensor[0][0];
      double g22 = a.G.xCovariantTensor[1][1];
      double g12 = a.G.xCovariantTensor[0][1];
      double dg11 = da.G.xCovariantTensor[0][0];
      double dg22 = da.G.xCovariantTensor[1][1];
      double dg12 = da.G.xCovariantTensor[0][1];
    
      double g11_times_g12 = g11 * g22;
      double d_g11_times_g12 = g11 * dg22 + dg11 * g22;

      double sin_alpha_sqr = (g11_times_g12 - g12 * g12) / g11_times_g12;
      double d_sin_alpha_sqr = 
        ((d_g11_times_g12 - 2 * g12 * dg12) - sin_alpha_sqr * d_g11_times_g12) / g11_times_g12;

      double dpen = - c * d_sin_alpha_sqr / (sin_alpha_sqr * sin_alpha_sqr);

      dTotalPenalty += dpen;
      }
    }

  return dTotalPenalty / saPenalty.GetCount();
}

void LoopTangentSchemeValidityPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  LoopTangentSchemeValidity Penalty Term : " << endl;
  sout << "    max penalty              : " << saPenalty.GetMax() << endl;
  sout << "    mean penalty             : " << saPenalty.GetMean() << endl;
}



/*********************************************************************************
 * Angles penalty term
 ********************************************************************************/

MedialAnglesPenaltyTerm
::MedialAnglesPenaltyTerm(GenericMedialModel *model)
: MedialIntegrationEnergyTerm(model)
{
}

double MedialAnglesPenaltyTerm::ComputeEnergy(SolutionData *S)
{
  xTotalPenalty = xMaxPenalty = 0.0;

  // We will simply look at the ratio |Xu x Xv| to |Xu| |Xv| and
  // take its inverse as the penalty, with a very small weight
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    double g11 = a.G.xCovariantTensor[0][0];
    double g22 = a.G.xCovariantTensor[1][1];

    double xpen = 0.01 * g11 * g22 * a.G.gInv;

    xTotalPenalty += xpen;
    xMaxPenalty = std::max(xpen, xMaxPenalty);
    }

  // return 0.25 * xTotalPenalty / S->xAtomGrid->GetNumberOfAtoms();
  xTotalPenalty /= S->xAtomGrid->GetNumberOfAtoms();
  return xTotalPenalty;
}

double MedialAnglesPenaltyTerm::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  // Compute the square of the cosine of the angle between Xu and Xv
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &da = dS->xAtoms[it.GetIndex()];
    //if(da.order <= 2)
    //  {
      MedialAtom &a = S->xAtoms[it.GetIndex()];
      double g11 = a.G.xCovariantTensor[0][0];
      double g22 = a.G.xCovariantTensor[1][1];
      double dg11 = da.G.xCovariantTensor[0][0];
      double dg22 = da.G.xCovariantTensor[1][1];
    
      double dxpen = 0.01 * (
        dg11 *  g22 *  a.G.gInv + 
         g11 * dg22 *  a.G.gInv +
         g11 *  g22 * da.G.gInv);

      dTotalPenalty += dxpen;
    //  }
    }

  return dTotalPenalty / S->xAtomGrid->GetNumberOfAtoms();
}

void MedialAnglesPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Square Cosine Penalty Term : " << endl;
  sout << "    max penalty              : " << xMaxPenalty << endl;
  sout << "    total penalty            : " << xTotalPenalty << endl;
}

/*********************************************************************************
 * Boundary Elasticity Prior
 ********************************************************************************/

/*********************************************************************************
 * Boundary Curvature Penalty Term
 ********************************************************************************/
// const double BoundaryCurvaturePenalty::xPower = 1;
// const double BoundaryCurvaturePenalty::xScale = 1;

BoundaryCurvaturePenalty
::BoundaryCurvaturePenalty(GenericMedialModel *model)
{
  // Define the regularization term as a spring energy term
}

double
BoundaryCurvaturePenalty
::ComputeEnergy(SolutionData *S)
{  
  // Reset the penalty
  saSqrEdgeLen.Reset();

  // Compute the squared distances between all neighboring nodes
  TriangleMesh *mesh = S->xAtomGrid->GetBoundaryMesh();

  // Loop once over every edge
  for(size_t t = 0; t < mesh->triangles.size(); t++)
    {
    for(size_t q = 0; q < 3; q++)
      {
      if(t < mesh->triangles[t].neighbors[q])
        {
        size_t ib1 = mesh->triangles[t].vertices[ror(q)];
        size_t ib2 = mesh->triangles[t].vertices[rol(q)];
        size_t ia1 = S->xAtomGrid->GetBoundaryPointAtomIndex(ib1);
        size_t ia2 = S->xAtomGrid->GetBoundaryPointAtomIndex(ib2);
        size_t is1 = S->xAtomGrid->GetBoundaryPointSide(ib1);
        size_t is2 = S->xAtomGrid->GetBoundaryPointSide(ib2);
        SMLVec3d X1 = S->xAtoms[ia1].xBnd[is1].X;
        SMLVec3d X2 = S->xAtoms[ia2].xBnd[is2].X;

        // Compute the energy
        SMLVec3d edge = X2 - X1;
        saSqrEdgeLen.Update(dot_product(edge, edge));
        }
      }
    }

  return saSqrEdgeLen.GetSum();
}

double 
BoundaryCurvaturePenalty
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Reset the penalty
  double dPenalty = 0.0;

  // Compute the squared distances between all neighboring nodes
  TriangleMesh *mesh = S->xAtomGrid->GetBoundaryMesh();

  // Loop once over every edge
  for(size_t t = 0; t < mesh->triangles.size(); t++)
    {
    for(size_t q = 0; q < 3; q++)
      {
      if(t < mesh->triangles[t].neighbors[q])
        {
        size_t ib1 = mesh->triangles[t].vertices[ror(q)];
        size_t ib2 = mesh->triangles[t].vertices[rol(q)];
        size_t ia1 = S->xAtomGrid->GetBoundaryPointAtomIndex(ib1);
        size_t ia2 = S->xAtomGrid->GetBoundaryPointAtomIndex(ib2);
        size_t is1 = S->xAtomGrid->GetBoundaryPointSide(ib1);
        size_t is2 = S->xAtomGrid->GetBoundaryPointSide(ib2);
        
        SMLVec3d X1 = S->xAtoms[ia1].xBnd[is1].X;
        SMLVec3d X2 = S->xAtoms[ia2].xBnd[is2].X;
        SMLVec3d dX1 = dS->xAtoms[ia1].xBnd[is1].X;
        SMLVec3d dX2 = dS->xAtoms[ia2].xBnd[is2].X;

        // Compute the energy
        SMLVec3d edge = X2 - X1;
        SMLVec3d d_edge = dX2 - dX1;

        dPenalty += 2 * dot_product(d_edge, edge);
        }
      }
    }

  return dPenalty;
}

void BoundaryCurvaturePenalty
::PrintReport(ostream &sout)
{
  sout << "  Boundary Spring Penalty: " << endl;
  sout << "    rms edge length   : " << sqrt(saSqrEdgeLen.GetMean()) << endl;
  sout << "    total penalty     : " << saSqrEdgeLen.GetSum() << endl;
}

/**
BoundaryCurvaturePenalty
::BoundaryCurvaturePenalty(GenericMedialModel *model)
{
  // Initialize the curvature vector array
  nBnd = model->GetNumberOfBoundaryPoints();
  xMeanCurvVec.resize(nBnd, SMLVec3d(0.0));
  dMeanCurvVec.resize(nBnd, SMLVec3d(0.0));
}

double
BoundaryCurvaturePenalty
::ComputeEnergy(SolutionData *S)
{  
  double xTest = 0.0;

  saDenom.Reset();
  saCurv.Reset();

  // Set all the curvature vectors to zero
  std::fill(xMeanCurvVec.begin(), xMeanCurvVec.end(), SMLVec3d(0.0));

  // Iterate over the triangles on the boundary
  for(MedialBoundaryTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Compute the cotangent of each angle in the triangle
    for(size_t j = 0; j < 3; j++)
      {

      size_t j1 = (j+1) % 3, j2 = (j+2) % 3;
      SMLVec3d &Xj = GetBoundaryPoint(it, S->xAtoms, j).X;
      SMLVec3d &Xj1 = GetBoundaryPoint(it, S->xAtoms, j1).X;
      SMLVec3d &Xj2 = GetBoundaryPoint(it, S->xAtoms, j2).X;
      SMLVec3d A = Xj2 - Xj, B = Xj1 - Xj;

      // Get the normal vector at this node (is this fair?)
      SMLVec3d &Nj = GetBoundaryPoint(it, S->xAtoms, j).N;

      // Compute the cotangent of the angle at j
      double AB = dot_product(A, B);
      SMLVec3d AxB = vnl_cross_3d(A, B);
      
      // double magAxB = dot_product(AxB.magnitude();
      double magAxB = - dot_product(Nj, AxB);
      double cotAlphaJ = AB / magAxB;
      
      xTest += cotAlphaJ;

      saDenom.Update(magAxB);      

      // Compute the contribution to the normal vector
      SMLVec3d xcontr = (Xj2 - Xj1) * cotAlphaJ;
      xMeanCurvVec[it.GetBoundaryIndex(j2)] += xcontr;
      xMeanCurvVec[it.GetBoundaryIndex(j1)] -= xcontr;
      }
    }

  // Integrate the sqaured mean curvature
  xIntegralSqrMeanCrv = 0.0;
  xCrestCurvatureTerm = 0.0;
  for(MedialBoundaryPointIterator ip(S->xAtomGrid); !ip.IsAtEnd(); ++ip)
    {
    size_t ib = ip.GetIndex();
    double w = S->xBoundaryWeights[ib];
    double mc2 = xMeanCurvVec[ib].squared_magnitude();
    xIntegralSqrMeanCrv += mc2 / (3.0 * w);
    // cout << mc2 << " " << w << " " << mc2 / (w * w) << endl;
    saCurv.Update(mc2 / (w * w));

    // We must also subtract one of the principal curvatures along the boundary
    // which is just 1/r
    if(ip.IsEdgeAtom())
      {
      double F = S->xAtoms[ip.GetAtomIndex()].F;
      xCrestCurvatureTerm += w / F;
      }
    }

  // TURN OFF CREST RADIUS TERM
  xCrestCurvatureTerm = 0;

  xPenalty = (4 * xIntegralSqrMeanCrv - xCrestCurvatureTerm) 
    / S->xBoundaryArea;
  
  return xTest;
}

double 
BoundaryCurvaturePenalty
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTest = 0.0;

  // Set all the curvature vectors to zero
  std::fill(dMeanCurvVec.begin(), dMeanCurvVec.end(), SMLVec3d(0.0));

  // Iterate over the triangles on the boundary
  for(MedialBoundaryTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Compute the cotangent of each angle in the triangle
    for(size_t j = 0; j < 3; j++)
      {
      // Only consider the point if it's order 2 or less
      if(dS->xAtoms[it.GetAtomIndex(j)].order <= 2)
        {
        size_t j1 = (j+1) % 3, j2 = (j+2) % 3;
        SMLVec3d &Xj = GetBoundaryPoint(it, S->xAtoms, j).X;
        SMLVec3d &Xj1 = GetBoundaryPoint(it, S->xAtoms, j1).X;
        SMLVec3d &Xj2 = GetBoundaryPoint(it, S->xAtoms, j2).X;
        SMLVec3d &DXj = GetBoundaryPoint(it, dS->xAtoms, j).X;
        SMLVec3d &DXj1 = GetBoundaryPoint(it, dS->xAtoms, j1).X;
        SMLVec3d &DXj2 = GetBoundaryPoint(it, dS->xAtoms, j2).X;
        SMLVec3d A = Xj2 - Xj, B = Xj1 - Xj;
        SMLVec3d dA = DXj2 - DXj, dB = DXj1 - DXj;

        // Get the normal vector at this node (is this fair?)
        SMLVec3d &Nj = GetBoundaryPoint(it, S->xAtoms, j).N;
        SMLVec3d &dNj = GetBoundaryPoint(it, dS->xAtoms, j).N;

        // Compute the cotangent of the angle at j
        double AB = dot_product(A, B);
        double dAB = dot_product(dA, B) + dot_product(A, dB);

        SMLVec3d AxB = vnl_cross_3d(A, B);
        SMLVec3d d_AxB = vnl_cross_3d(dA, B) + vnl_cross_3d(A, dB);

        // double magAxB = AxB.magnitude();
        // double d_magAxB = dot_product(AxB, d_AxB) / magAxB;
        
        double magAxB = - dot_product(Nj, AxB);
        double d_magAxB = - dot_product(Nj, d_AxB) - dot_product(dNj, AxB);

        

        double cotAlphaJ = AB / magAxB;
        double d_cotAlphaJ = (dAB * magAxB - AB * d_magAxB) / (magAxB * magAxB);
        cout << d_cotAlphaJ << endl;

        dTest += d_cotAlphaJ;

        // Compute the contribution to the normal vector
        SMLVec3d xcontr = (Xj2 - Xj1) * cotAlphaJ;
        SMLVec3d d_xcontr = 
          (Xj2 - Xj1) * d_cotAlphaJ + (DXj2 - DXj1) * cotAlphaJ;

        dMeanCurvVec[it.GetBoundaryIndex(j2)] += d_xcontr;
        dMeanCurvVec[it.GetBoundaryIndex(j1)] -= d_xcontr;
        }
      }
    }

  // Integrate the sqaured mean curvature
  double dIntegralSqrMeanCrv = 0.0;
  double dCrestCurvatureTerm = 0.0;
  for(MedialBoundaryPointIterator ip(S->xAtomGrid); !ip.IsAtEnd(); ++ip)
    {
    if(dS->xAtoms[ip.GetAtomIndex()].order <= 2)
      {
      size_t ib = ip.GetIndex();

      double w = S->xBoundaryWeights[ib];
      double dw = dS->xBoundaryWeights[ib];
      dIntegralSqrMeanCrv +=
        (2 * dot_product(xMeanCurvVec[ib], dMeanCurvVec[ib]) * w -
         xMeanCurvVec[ib].squared_magnitude() * dw) / (3.0 * w * w);

      // We must also subtract one of the principal curvatures along the boundary
      // which is just 1/r
      if(ip.IsEdgeAtom())
        {
        double F = S->xAtoms[ip.GetAtomIndex()].F;
        double dF = dS->xAtoms[ip.GetAtomIndex()].F;
        dCrestCurvatureTerm += (dw * F - w * dF) / (F * F);
        }
      }
    }

  // TURN OFF CREST RADIUS TERM
  dCrestCurvatureTerm = 0;

  
  double d_penalty = 
    ((4 * dIntegralSqrMeanCrv - dCrestCurvatureTerm) - xPenalty * dS->xBoundaryArea)
    / S->xBoundaryArea; 

  return dTest;
}

void BoundaryCurvaturePenalty
::PrintReport(ostream &sout)
{
  sout << "  Boundary Curvature Penalty: " << endl;
  sout << "    int sqr mean curv           : " << xIntegralSqrMeanCrv << endl;
  sout << "    int crest kappa1            : " << xCrestCurvatureTerm << endl;
  sout << "    mean curv sq (min/mean/max) : " << saCurv.GetMin() 
    << " / " << saCurv.GetMean() << " / " << saCurv.GetMax() << endl;
  sout << "    magAxB (min/mean/max) : " << saDenom.GetMin() 
    << " / " << saDenom.GetMean() << " / " << saDenom.GetMax() << endl;
  sout << "    total penalty     : " << xPenalty << endl;
}
*/

/*********************************************************************************
 * Medial Curvature Penalty Term
 ********************************************************************************/
const double MedialCurvaturePenalty::xPower = 5;
const double MedialCurvaturePenalty::xScale = 2;

double
MedialCurvaturePenalty
::ComputeEnergy(SolutionData *S)
{
  // Reset the accumulators
  saMeanCurv.Reset();
  saGaussCurv.Reset();
  saPenalty.Reset();
  saSumSqKappa.Reset();
  saFeature.Reset();
  saRad.Reset();

  // Iterate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    saMeanCurv.Update(a.xMeanCurv);
    saGaussCurv.Update(a.xGaussCurv);
    saRad.Update(a.R);

    // Compute sum of the squares of principal curvatures
    double k2 = 4 * a.xMeanCurv * a.xMeanCurv - 2 * a.xGaussCurv;
    saSumSqKappa.Update(k2);

    // Compute the feature that gets penalized
    double f = k2;
    saFeature.Update(f);

    // Compute the penalty. The penalty should be such that the 
    double p = pow(f / xScale, xPower);
    saPenalty.Update(p);
    }

  // Return the sum of the penalty terms
  return saPenalty.GetSum();
}

double 
MedialCurvaturePenalty
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  // Iterate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    MedialAtom &da = dS->xAtoms[it.GetIndex()];

    // Compute sum of the squares of principal curvatures
    double k2 = 4 * a.xMeanCurv * a.xMeanCurv - 2 * a.xGaussCurv;
    double dk2 = 8 * a.xMeanCurv * da.xMeanCurv - 
      4 * a.xGaussCurv * da.xGaussCurv;

    // Compute the feature that gets penalized
    double f = k2;
    double df = dk2;

    // Compute the penalty. The penalty should be such that the 
    double dp = xPower * pow(f / xScale, xPower-1) * df / xScale;
    dTotalPenalty += dp;
    }

  return dTotalPenalty;
}

void MedialCurvaturePenalty
::PrintReport(ostream &sout)
{
  sout << "  Medial Curvature Penalty: " << endl;
  sout << "    mean curvature stats: " << saMeanCurv << endl;
  sout << "    gaussuan curv. stats: " << saGaussCurv << endl;
  sout << "    sum sqr. kappa stats: " << saSumSqKappa << endl;
  sout << "    radius stats: " << saRad << endl;
  sout << "    feautre stats: " << saFeature << endl;
  sout << "    total penalty: " << saPenalty.GetSum() << endl;
}

/*********************************************************************************
 * Bending Energy Term
 ********************************************************************************/
MedialBendingEnergyTerm::MedialBendingEnergyTerm(GenericMedialModel *model)
  : MedialIntegrationEnergyTerm(model)
{
}

double MedialBendingEnergyTerm::ComputeEnergy(SolutionData *S)
{
  // Bending energy defined as Xuu * Xuu + Xvv * Xvv + 2 Xuv * Xuv
  xMaxBending = 0;
  xTotalBending = 0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Get the medial atom
    MedialAtom &a = S->xAtoms[it.GetIndex()];

    // Compute the regularity term here
    double be = 
      dot_product(a.Xuu, a.Xuu) + 
      dot_product(a.Xvv, a.Xvv) + 
      2.0 * dot_product(a.Xuv, a.Xuv);

    // Update the maximum
    xMaxBending = std::max(xMaxBending, be);

    // Update the integral
    xTotalBending += be * xDomainWeights[it.GetIndex()];
    }

  xMeanBending = xTotalBending / xDomainArea;
  return xMeanBending;
}

double MedialBendingEnergyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalBending = 0.0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    MedialAtom &da = dS->xAtoms[it.GetIndex()];

    // Compute the regularity term here
    double d_be = 
      2.0 * dot_product(a.Xuu, da.Xuu) + 
      2.0 * dot_product(a.Xvv, da.Xvv) + 
      4.0 * dot_product(a.Xuv, da.Xuv);

    // Update the integral
    dTotalBending += xDomainWeights[it.GetIndex()] * d_be;
    }

  // Return the error
  return dTotalBending / xDomainArea; 
}

void MedialBendingEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Medial Bending Energy Term: " << endl;
  sout << "    maximum bending energy:     " << xMaxBending << endl;
  sout << "    total bending energy: " << xTotalBending << endl;
  sout << "    average bending energy: " << xMeanBending << endl;
}




/*********************************************************************************
 * Regularity term (used to maintain correspondence)
 ********************************************************************************/
MedialRegularityTerm::MedialRegularityTerm(GenericMedialModel *model)
  : MedialIntegrationEnergyTerm(model)
{
}

double MedialRegularityTerm::ComputeEnergy(SolutionData *S)
{
  // Integral of (G2(1,11) + G(2,12))^2 + (G(1,21) + G(2,22)) du dv 
  // (no area elt. scaling)
  xGradMagIntegral = 0;
  xGradMagMaximum = 0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Get the medial atom
    MedialAtom &a = S->xAtoms[it.GetIndex()];

    // Only compute the penalty if the atom is internal
    // if(a.flagCrest || !a.flagValid) 
    //  continue;

    // Compute the regularity term here
    double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
    double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
    double reg = reg1 * reg1 + reg2 * reg2;

    // Update the sum and the maximum
    if(reg > xGradMagMaximum)
      xGradMagMaximum = reg;

    // Update the integral
    xGradMagIntegral += xDomainWeights[it.GetIndex()] * reg;
    }

  return xGradMagIntegral / xDomainArea;
}

double MedialRegularityTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dIntegral = 0.0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &da = dS->xAtoms[it.GetIndex()];
    if(da.order <= 2)
      {
      MedialAtom &a = S->xAtoms[it.GetIndex()];

      // Only compute the penalty if the atom is internal
      // if(a.flagCrest || !a.flagValid) 
      //  continue;

      // Compute the regularity term here
      double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
      double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
      double dreg1 = da.G.xChristoffelSecond[0][0][0] + da.G.xChristoffelSecond[1][0][1];
      double dreg2 = da.G.xChristoffelSecond[0][1][0] + da.G.xChristoffelSecond[1][1][1];
      double dreg = 2.0 * (reg1 * dreg1 + reg2 * dreg2);

      // Update the integral
      dIntegral += xDomainWeights[it.GetIndex()] * dreg;
      }
    }

  // Return the error
  return dIntegral / xDomainArea; 
}

void MedialRegularityTerm::PrintReport(ostream &sout)
{
  sout << "  Medial Regularization Term: " << endl;
  sout << "    maximum grad mag sqr:     " << xGradMagMaximum << endl;
  sout << "    grad mag sqr integral: " << xGradMagIntegral << endl;
  sout << "    domain area: " << xDomainArea << endl;
}

/*
double RadiusPenaltyTerm::ComputeEnergy(SolutionData *S)
{ 
  xTotalPenalty = 0.0;
  xMinR2 = 1e100;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    // Get the square of R
    double phi = S->xAtoms[i].F;
    if(phi < xMinR2)
      xMinR2 = phi;

    // Apply the penalty function
    double x = phi / xScale;
    double x2 = x * x;
    double x4 = x2 * x2;
    
    xTotalPenalty += 1.0 / x4;
    }

  xTotalPenalty /= S->xAtomGrid->GetNumberOfAtoms();

  return xTotalPenalty;
}

double RadiusPenaltyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    if(dS->xAtoms[i].order == 0)
      {
      // Get the square of R
      double phi = S->xAtoms[i].F;
      double dphi = dS->xAtoms[i].F;

      // Apply the penalty function
      double x = phi / xScale, dx = dphi / xScale;
      double x2 = x * x;
      double x4 = x2 * x2;

      dTotalPenalty += -4.0 * dx / (x * x4);
      }
    }

  dTotalPenalty /= S->xAtomGrid->GetNumberOfAtoms();

  return dTotalPenalty;
}  

void RadiusPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Radius Penalty Term : " << endl;
  sout << "    total penalty            : " << xTotalPenalty << endl;
  sout << "    smallest R^2             : " << xMinR2 << endl;
}
*/

/*********************************************************************************
 * MEDIAL RADIUS PENALTY TERM
 ********************************************************************************/
double RadiusPenaltyTerm::ComputeEnergy(SolutionData *S)
{ 
  saRadius.Reset();
  saPenHigh.Reset();
  saPenLow.Reset();

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    // Get the radius
    double r = S->xAtoms[i].R;
    saRadius.Update(r);

    // Set the inverse radius penalty
    saPenLow.Update(1.0 / r);

    // Compare to min and max
    if(r > rMax)
      {
      double dr = r - rMax;
      saPenHigh.Update(dr * dr);
      }
    }

  return saPenLow.GetSum() * sLower + saPenHigh.GetSum() * sUpper;
}

double RadiusPenaltyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dUpper = 0.0, dLower = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    if(dS->xAtoms[i].order == 0)
      {
      // Get the radius
      double r = S->xAtoms[i].R;
      double d_r = dS->xAtoms[i].R;

      // Compare to min and max
      dLower += - d_r / (r * r);
      
      if(r > rMax)
        {
        dUpper += (r - rMax) * d_r;
        }
      }
    }

  return 2 * (sUpper * dUpper) + sLower * dLower;
}  

void RadiusPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Radius Penalty Term : " << endl;
  sout << "    min R                    : " << saRadius.GetMin() << endl;
  sout << "    max R                    : " << saRadius.GetMax() << endl;
  sout << "    mean R                   : " << saRadius.GetMean() << endl;
  sout << "    lower bound              : " << rMin << endl;
  sout << "    upper bound              : " << rMax << endl;
  sout << "    lower bound violations   : " << saPenLow.GetCount() << endl;
  sout << "    lower penalty            : " << saPenLow.GetSum() << endl;
  sout << "    upper bound violations   : " << saPenHigh.GetCount()<< endl;
  sout << "    upper penalty            : " << saPenHigh.GetSum() << endl;
  sout << "    total penalty            : " 
    << saPenLow.GetSum() * sLower + saPenHigh.GetSum() * sUpper << endl;
}

void RadiusPenaltyTerm::SetParameters(Registry &R)
{
  rMin = R["LowerBound"][0.0];
  rMax = R["UpperBound"][1e100];
  sLower = R["LowerScale"][1.0];
  sUpper = R["UpperScale"][1.0];
}

/*********************************************************************************
 * FIT TO POINT SET ENERGY TERM
 ********************************************************************************/
DistanceToPointSetEnergyTerm
::DistanceToPointSetEnergyTerm(GenericMedialModel *model,
  double *x, double *y, double *z)
: MedialIntegrationEnergyTerm(model)
{
  // Store the XYZ values
  for(size_t i = 0; i < model->GetNumberOfAtoms(); i++)
    target.push_back(SMLVec3d(x[i], y[i], z[i]));
}

double
DistanceToPointSetEnergyTerm
::ComputeEnergy(SolutionData *S)
{
  xTotalDist = 0.0;
  xTotalArea = 0.0;
  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); ++i)
    {
    SMLVec3d delta = S->xAtoms[i].X - target[i];
    double weight = S->xAtoms[i].aelt * xDomainWeights[i];
    xTotalDist += delta.squared_magnitude() * weight;
    xTotalArea += weight;
    }
  xTotalMatch = xTotalDist / xTotalArea;
  return xTotalMatch;
}

double
DistanceToPointSetEnergyTerm
::ComputePartialDerivative(
  SolutionData *S,
  PartialDerivativeSolutionData *dS)
{
  double dTotalDist = 0.0;
  double dTotalArea = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); ++i)
    {
    SMLVec3d delta = S->xAtoms[i].X - target[i];
    SMLVec3d d_delta = dS->xAtoms[i].X;
    double weight = S->xAtoms[i].aelt * xDomainWeights[i];
    double d_weight = dS->xAtoms[i].aelt * xDomainWeights[i];
    dTotalDist += 
      (2.0 * dot_product(delta, d_delta) * weight + 
       delta.squared_magnitude() * d_weight);
    dTotalArea += d_weight;
    }

  double dTotalMatch = (dTotalDist - dTotalArea * xTotalMatch) / xTotalArea;
  return dTotalMatch;
}

void
DistanceToRadiusFieldEnergyTerm
::PrintReport(ostream &sout)
{
  sout << "  DistanceToRadiusFieldEnergyTerm:" << endl;
  sout << "    Mean squared phi-distance     : " << xTotalMatch << endl;
}



/*********************************************************************************
 * FIT TO RADIUS FIELD ENERGY TERM
 ********************************************************************************/
DistanceToRadiusFieldEnergyTerm
::DistanceToRadiusFieldEnergyTerm(GenericMedialModel *model, double *rad)
: MedialIntegrationEnergyTerm(model),
  xTargetPhi(rad, model->GetNumberOfAtoms()), 
  xLastDelta(model->GetNumberOfAtoms(), 0.0)
{
  // Put r^2 in the target field
  for(size_t i = 0; i < model->GetNumberOfAtoms(); i++)
    if(xTargetPhi[i] < 0)
      xTargetPhi[i] = 0.0;
    else
      xTargetPhi[i] = xTargetPhi[i] * xTargetPhi[i];
}

double
DistanceToRadiusFieldEnergyTerm
::ComputeEnergy(SolutionData *S)
{
  xTotalDiff = 0.0;
  xTotalArea = 0.0;
  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); ++i)
    {
    double delta = S->xAtoms[i].F - xTargetPhi[i];
    double weight = S->xAtoms[i].aelt * xDomainWeights[i];

    xLastDelta[i] = delta;
    xTotalDiff += delta * delta * weight;
    xTotalArea += weight;
    }
  xTotalMatch = xTotalDiff / xTotalArea;
  return xTotalMatch;
}

double
DistanceToRadiusFieldEnergyTerm
::ComputePartialDerivative(
  SolutionData *S,
  PartialDerivativeSolutionData *dS)
{
  double dTotalDiff = 0.0;
  double dTotalArea = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); ++i)
    {
    double delta = xLastDelta[i];
    double d_delta = dS->xAtoms[i].F;
    double weight = S->xAtoms[i].aelt * xDomainWeights[i];
    double d_weight = dS->xAtoms[i].aelt * xDomainWeights[i];
    dTotalDiff += delta * (2.0 * d_delta * weight + delta * d_weight);
    dTotalArea += d_weight;
    }

  double dTotalMatch =( dTotalDiff - xTotalMatch * dTotalArea) / xTotalArea;
  return dTotalMatch;
}

void
DistanceToPointSetEnergyTerm
::PrintReport(ostream &sout)
{
  sout << "  DistanceToRadiusFieldEnergyTerm:" << endl;
  sout << "    Mean squared distance     : " << xTotalMatch << endl;
}

/*********************************************************************************
 * MESH REGULARIZATION PENALTY TERM
 ********************************************************************************/
MeshRegularizationPenaltyTerm
::MeshRegularizationPenaltyTerm(
  GenericMedialModel *model, size_t nu, size_t nv)
{
  // We first compute the weight matrix W. For the given basis with coeffs c
  // the matrix gives x = Wc. Then the regularization prior for x with respect to
  // the basis is given by |x|^2 - x'W.inv(W'W).W'x, or |x|^2 - |Qx|^2, where
  // (U, S, V) = svd(W) and Q = U[:,1:m]', i.e, the first m rows of the U matrix
  
  // Set m, the number of basis functions, and n, the number of atoms
  size_t m = nu * nv;
  size_t n = model->GetNumberOfAtoms();
  size_t i, j;

  // Compute the range of u and v in the atoms for the basis functions
  StatisticsAccumulator sau, sav;
  for(size_t i = 0; i < n; i++)
    { 
    sau.Update(model->GetAtomArray()[i].u);
    sav.Update(model->GetAtomArray()[i].v);
    }

  // Fill out the weight matrix
  vnl_matrix<double> W(n, m);
  for(i = 0; i < n; i++)
    {
    // Get the medial atom
    const MedialAtom &a = model->GetAtomArray()[i];

    // Scale u and v to the unit square (isn't this a problem for funky shapes?)
    double u = (a.u - sau.GetMin()) / (sau.GetMax() - sau.GetMin());
    double v = (a.v - sav.GetMin()) / (sav.GetMax() - sav.GetMin());

    // Shift u so no point is on the boundary
    u = 0.9 * u + 0.05; v = 0.9 * v + 0.05;

    // Compute the weight matrix
    j = 0;
    for(size_t iu = 0; iu < nu; iu++) for(size_t iv = 0; iv < nv; iv++)
      {
      // Compute the basis functions for u and v
      W(i, j++) = cos(vnl_math::pi * iu * u) * cos(vnl_math::pi * iv * v);
      }
    }

  // Compute the SVD of the matrix W and the matrix Q
  vnl_svd<double> svd(W);
  Qt = svd.U().get_n_columns(0,m);
  Q = Qt.transpose();

  // TEST code
  Z = W * vnl_svd<double>(W.transpose() * W).inverse() * W.transpose();
}

double
MeshRegularizationPenaltyTerm
::ComputeEnergy(SolutionData *S)
{
  // The penalty term
  saPenalty.Reset();

  // Compute the distortion in X, Y and Z (I guess we can let R alone)
  vnl_vector<double> x(S->xAtomGrid->GetNumberOfAtoms());

  // Repeat for each component
  for(size_t c = 0; c < 3; c++)
    {
    size_t i;

    // Populate the x-vector
    for(i = 0; i < x.size(); i++)
      {
      x[i] = S->xAtoms[i].X[c];
      }

    vnl_vector<double> delta = x - Z * x;
    saDelta[c].Reset();
    for(i = 0; i < x.size(); i++)
      saDelta[c].Update(delta[i] * delta[i]);

    // Compute the penalty
    vnl_vector<double> Qx = Q * x;
    saPenalty.Update(dot_product(x,x) - dot_product(Qx,Qx));

    // Update the vector b
    b[c] = 2.0 * (x - Qt * Qx);
    }



  // Return the total penalty
  return saPenalty.GetSum() / S->xAtomGrid->GetNumberOfAtoms();
}

// Compute the partial derivative
double 
MeshRegularizationPenaltyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Initialize the derivative
  double dPenalty = 0.0;

  // Generate the derivative vectors
  vnl_vector<double> dx(S->xAtomGrid->GetNumberOfAtoms());

  // Repeat for each component
  for(size_t c = 0; c < 3; c++)
    {
    // Populate the x-vector
    for(size_t i = 0; i < dx.size(); i++)
      dx[i] = dS->xAtoms[i].X[c];

    // Compute the gradient
    dPenalty += dot_product(b[c], dx);
    }

  return dPenalty / S->xAtomGrid->GetNumberOfAtoms();
}

// Describe the terms of the penalty
void 
MeshRegularizationPenaltyTerm
::PrintReport(ostream &sout)
{
  sout << "  MeshRegularizationPenaltyTerm:" << endl;
  sout << "    Total penalty     : " << saPenalty.GetSum() << endl;
  sout << "    X distance mean   : " << saDelta[0].GetMean() << endl;
  sout << "    Y distance mean   : " << saDelta[1].GetMean() << endl;
  sout << "    Z distance mean   : " << saDelta[2].GetMean() << endl;
}

/*********************************************************************************
 * MEDIAL OPTIMIZATION PROBLEM
 ********************************************************************************/
const double MedialOptimizationProblem::xPrecision = 1.0e-10;
const double MedialOptimizationProblem::xEpsilon = 1.0e-6;

MedialOptimizationProblem
::MedialOptimizationProblem(
  GenericMedialModel *xMedialModel, CoefficientMapping *xCoeff)
{
  nGradCalls = nEvalCalls = 0;

  this->xCoeff = xCoeff;
  this->xMedialModel = xMedialModel;
  this->nCoeff = xCoeff->GetNumberOfParameters();

  flagLastEvalAvailable = false;
  flagPhiGuessAvailable = false;
  flagGradientComputed = false;
  flagDumpGradientMesh = false;

  // Save the initial coefficients currently in the model
  xInitialCoefficients = xMedialModel->GetCoefficientArray();
  
  // Allocate the array of medial atoms
  dAtoms = new MedialAtom[xMedialModel->GetNumberOfAtoms()];

  // Create the solution data and its 

  // Initialize the basis array
  xBasis.set_size(nCoeff, xCoeff->GetNumberOfCoefficients());

  // Prepare medial model for gradient computation by specifying the set of 
  // variations for which variational derivatives are to be computed. This
  // may only be done for coefficient mappings that are linear. For non-linear
  // coefficient mappings, the basis has to be specified at each iteration
  // because the variations change depending on where in parameter space we are
  if(xCoeff->IsLinear()) 
    {
    // Compute the all the variations in the space of model's coefficients
    for(size_t i = 0; i < nCoeff; i++)
      xBasis.set_row(i, xCoeff->GetVariationForParameter(
        xInitialCoefficients, Vec(xCoeff->GetNumberOfParameters(), 0.0), i));

    // Pass these variations to the model
    xMedialModel->SetVariationalBasis(xBasis);
    }

  // Initialize the statistical arrays
  xGradSum.set_size(nCoeff);
  xGradSumSqr.set_size(nCoeff);

  // Initialize the solution and derivative solution
  S = new SolutionData(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());
  dS = new PartialDerivativeSolutionData(S, dAtoms);

  flagQuiet = false;
}

MedialOptimizationProblem
::~MedialOptimizationProblem()
{
  delete S;
  delete dS;
  delete dAtoms;
}

// Add the energy term
void MedialOptimizationProblem::AddEnergyTerm(EnergyTerm *term, double xWeight)
{
  xWeights.push_back(xWeight);
  xTerms.push_back(term);
  xTimers.push_back(CodeTimer());
  xGradTimers.push_back(CodeTimer());
  xLastGradientPerTerm.push_back(vnl_vector<double>(nCoeff, 0.0));
}


bool MedialOptimizationProblem::SolvePDE(double *xEvalPoint)
{
  // Make a vector from the eval point
  vnl_vector<double> X(xEvalPoint, nCoeff);

  // Check if we are re-evaluating at a previously tried location
  if(flagLastEvalAvailable && X == xLastEvalPoint)
    return false;

  // Turn of the last eval flag in case that exception is thrown below
  flagLastEvalAvailable = false;

  // Solve the equation
  xSolveTimer.Start();

  // Update the medial model with the new coefficients
  xMedialModel->SetCoefficientArray(xCoeff->Apply(xInitialCoefficients, X));

  // Solve the PDE using the last phi field if we have it
  if(flagGradientComputed)
    xMedialModel->ComputeAtoms(false, xLastGradHint.data_block());
  else
    xMedialModel->ComputeAtoms(false);

  // Stop timing
  xSolveTimer.Stop();

  // Set the last evaluation point (this is to ensure we never double-evaluate)
  xLastEvalPoint = X;
  flagLastEvalAvailable = true;

  // Return true : solution updated
  return true;
}

double MedialOptimizationProblem::Evaluate(double *xEvalPoint)
{
  // Solve the PDE - if there is no update, return the last solution value
  double t0 = clock();

  bool flagUpdate = SolvePDE(xEvalPoint);

  // cout << "[SLV: " << (clock() - t0)/CLOCKS_PER_SEC << " s] " << flush;

  // If the solution changed compute the terms
  if(flagUpdate)
    {
    t0 = clock();

    // Compute the solution for each term
    xLastSolutionValue = 0.0;
    xLastTermValues.set_size(xTerms.size());

    // Compute the medial integration terms
    S->ComputeIntegrationWeights();

    // Evaluate each of the terms
    ++nEvalCalls;
    // printf("%4d   %4d   ",nGradCalls,++nEvalCalls);
    for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
      { 
      xTimers[iTerm].Start();
      xLastTermValues[iTerm] = xTerms[iTerm]->ComputeEnergy(S);
      xLastSolutionValue += xLastTermValues[iTerm] * xWeights[iTerm]; 
      xTimers[iTerm].Stop();
      // printf("%7.3le  ",xLastTermValues[iTerm]);
      }
    // printf("  |  %7.3le\n", xLastSolutionValue);

    // cout << "[MAP: " << (clock() - t0)/CLOCKS_PER_SEC << " s] " << endl;
    }

  // Return the result
  return xLastSolutionValue;
}

double
MedialOptimizationProblem
::TestGradientComputation(double *x, double eps)
{
  // Allocate the analytic and finite difference gradients
  vnl_vector<double> gcd(nCoeff), ga(nCoeff), acc(nCoeff);

  // Compute the gradient at the current point
  ComputeGradient(x, ga.data_block());

  // Compute the central difference approximation
  for(size_t i = 0; i < nCoeff; i++)
    { 
    // Compute central difference
    vnl_vector<double> x1(x, nCoeff), x2(x, nCoeff);
    x1[i] -= eps;
    x2[i] += eps;
    gcd[i] = (Evaluate(x2.data_block()) - Evaluate(x1.data_block())) / (2 * eps);

    // Compute the accuracy of the approximation
    acc[i] = fabs(gcd[i] - ga[i]) / (eps + 0.5 * (fabs(gcd[i]) + fabs(ga[i])));
    }

  // Construct a test value
  return acc.inf_norm();
}

#include "PDESubdivisionMedialModel.h"

double 
MedialOptimizationProblem
::ComputeGradient(double *xEvalPoint, double *xGradient)
{
  size_t iTerm;

  // TODO: REMOVE THIS!!!
  // 
  /* 
  PDESubdivisionMedialModel *smm = dynamic_cast<PDESubdivisionMedialModel *>(xMedialModel);
  if(smm)
    {
    SubdivisionSurface::MeshLevel mlCoeffOld = *smm->GetCoefficientMesh(); 

    size_t nc = smm->GetNumberOfComponents();

    // We need to get a list of coordinates for remeshing
    typedef vnl_vector_fixed<double, 3> Vec;
    Vec *X = new Vec[mlCoeffOld.nVertices];
    for(size_t i = 0; i < mlCoeffOld.nVertices; i++)
      for(size_t k = 0; k < 3; k++)
        X[i][k] = smm->GetCoefficient(i * nc + k);

    mlCoeffOld.MakeDelaunay(X);

    smm->SetMesh(mlCoeffOld, 
      smm->GetCoefficientArray(), 
      smm->GetCoefficientU(),
      smm->GetCoefficientV(),
      smm->GetSubdivisionLevel(), 0);

    flagLastEvalAvailable = false;
    }
  */

  // Solve the PDE
  SolvePDE(xEvalPoint);

  // Begin the gradient computation timer
  xSolveGradTimer.Start();
  
  // If the coefficient mapping is non-linear, we have to specify the basis
  // for gradient computation at each iteration, because the set of variations
  // changes depending on where we are in search space. 
  if(!xCoeff->IsLinear()) 
    {
    // Compute the all the variations in the space of model's coefficients
    for(size_t i = 0; i < nCoeff; i++)
      xBasis.set_row(i, xCoeff->GetVariationForParameter(
        xInitialCoefficients, Vec(xEvalPoint, nCoeff), i));

    // Pass these variations to the model
    xMedialModel->SetVariationalBasis(xBasis);
    }

  // Begin the gradient computation for the model
  xMedialModel->BeginGradientComputation();

  // Pause the solver gradient timer
  xSolveGradTimer.Stop();

  // Compute the integration weights
  S->ComputeIntegrationWeights();


  // Verbose reporting
  if(!flagQuiet)
    {
    // Print header line
    if((nGradCalls % 10) == 0)
      {
      printf("\nGRAD   EVAL  ");
      for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
        printf("   %7s ",xTerms[iTerm]->GetShortName().c_str());
      printf("   |  TOTAL\n");
      }

    printf("%4lu   %4lu   ", (long unsigned) ++nGradCalls, (long unsigned) nEvalCalls);
    }
    
  // Begin the gradient computation for each of the energy terms
  xLastSolutionValue = 0.0;
  if(xLastGradEvalTermValues.size() == 0)
    xLastGradEvalTermValues.set_size(xTerms.size());
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    xTimers[iTerm].Start();
    double lval = xLastGradEvalTermValues[iTerm];
    xLastGradEvalTermValues[iTerm] = xTerms[iTerm]->BeginGradientComputation(S);
    xLastSolutionValue += xWeights[iTerm] * xLastGradEvalTermValues[iTerm];
    xTimers[iTerm].Stop();
    
    if(!flagQuiet)
      printf("%7.3le%s ",
        xWeights[iTerm] * xLastGradEvalTermValues[iTerm],
        (xLastGradEvalTermValues[iTerm] < lval) ? "*" : " ");
    }

  if(!flagQuiet)
    printf("  |  %7.3le\n", xLastSolutionValue);

  // Iterate variation by variation to compute the gradient
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Compute the variational derivative
    xSolveGradTimer.Start();
    xMedialModel->ComputeAtomVariationalDerivative(iCoeff, dAtoms);
    xSolveGradTimer.Stop();

    // Compute integration weights
    dS->ComputeIntegrationWeights();
    
    // Compute the partial derivatives for each term
    xGradient[iCoeff] = 0.0;
    for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
      {
      xGradTimers[iTerm].Start();
      xLastGradientPerTerm[iTerm][iCoeff] = xTerms[iTerm]->ComputePartialDerivative(S, dS);
      xGradient[iCoeff] += xWeights[iTerm] * xLastGradientPerTerm[iTerm][iCoeff];
      xGradTimers[iTerm].Stop();
      }
    }

  // Clear up gradient computation
  xMedialModel->EndGradientComputation();
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xTerms[iTerm]->EndGradientComputation();

  

  // Store the information about the gradient
  xLastGradPoint = vnl_vector<double>(xEvalPoint, nCoeff);
  xLastGradient = vnl_vector<double>(xGradient, nCoeff);
  xLastGradHint = xMedialModel->GetHintArray();
  flagGradientComputed = true;

  // cout << "Last gradient : " << xLastGradient << endl;

  // Random quality control check
  /*
  vnl_random randy;
  if(randy.lrand32(20) == 0)
    {
    Vec xVar(nCoeff,0.0);
    double eps = 0.0001;
    for(size_t i=0; i < nCoeff;i++)
      xVar[i] = randy.drand32(-1.0, 1.0);
    Vec x1 = Vec(xEvalPoint, nCoeff) + eps * xVar;
    Vec x2 = Vec(xEvalPoint, nCoeff) - eps * xVar;
    double dfn = 0.5 * 
      (Evaluate(x1.data_block()) - Evaluate(x2.data_block())) / eps;
    double dfa = dot_product(xVar, Vec(xGradient, nCoeff));
    printf(
      "QA GRAD CHECK: ANDRV = %4.2le  "
      "FDDRV = %4.2le  ABSER = %4.2le  RELER = %4.2le\n",
      dfa, dfn, fabs(dfa-dfn), fabs(dfa-dfn) / fabs(eps+dfa+dfn));
    Evaluate(xEvalPoint);
    }
    */

  // If gradient mesh dumping is on, do it
  if(flagDumpGradientMesh)
    DumpGradientMesh();

  // Return the solution value
  return xLastSolutionValue;
}

void MedialOptimizationProblem::PrintReport(ostream &sout)
{
  sout << "Optimization Summary: " << endl;
  sout << "  # variables           : " << xLastEvalPoint.size() << endl; 
  sout << "  energy value          : " << xLastSolutionValue << endl; 
  sout << "  solver time           : " << xSolveTimer.Read() << endl;
  sout << "  solver gradient time  : " << xSolveGradTimer.Read() << endl;
  sout << "  weights time          : " << xWeightsTimer.Read() << endl;
  sout << "  weights gradient time : " << xWeightsGradTimer.Read() << endl;
  sout << "Per-Term Report:" << endl;

  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    sout << "Term " << iTerm << endl;
    xTerms[iTerm]->PrintReport(sout);
    sout << "  Contribution: " << 
      xWeights[iTerm] << " * " << xLastTermValues[iTerm] <<
      " = " << xWeights[iTerm] * xLastTermValues[iTerm] << endl;
    sout << "  Elapsed time: " << xTimers[iTerm].Read() << endl;
    sout << "  Gradient time: " << xGradTimers[iTerm].Read() << endl;
    }
}

#include "SubdivisionMedialModel.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "MedialModelIO.h"

void MedialOptimizationProblem::DumpGradientMesh()
{
  // Compute the partial derivative of each term with respect to each
  // coefficient. We are assuming that the coefficients are mapping one
  // to one to the coefficient mesh of a medial model.
  SubdivisionMedialModel *model 
    = dynamic_cast<SubdivisionMedialModel *>(xMedialModel);
  if(!model)
    throw string("DumpGradientMesh operation applied only to subdivision models");

  // Get the control mesh
  const SubdivisionSurface::MeshLevel *mctl = model->GetCoefficientMesh();

  // Check that the number of parameters matches
  size_t nc = model->GetNumberOfComponents();
  if(xCoeff->GetNumberOfCoefficients() != mctl->nVertices * nc)
    throw string("Coefficient mapping does not support DumpGradientMesh");

  // Save the cm-rep file for this model
  char *fnout = new char[100];
  sprintf(fnout, "dumpgrad%04lu.cmrep", (unsigned long) this->nGradCalls);
  SubdivisionMedialModelIO::WriteModel(model, fnout);
  
  // Now load the VTK file and append with some gradient fields
  sprintf(fnout, "dumpgrad%04lu.vtk", (unsigned long) this->nGradCalls);
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fnout);
  reader->Update();
  vtkPolyData *poly = reader->GetOutput();

  /*
  // Generate a VTK mesh
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(mctl.nVertices);
  for(size_t i = 0; i < mctl.nVertices; i++)
    lPoints->InsertNextPoint(model->GetCoefficientArray().data_block() + i * 4);

  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(lPoints);
  SubdivisionSurface::ExportLevelToVTK(mctl, poly);
  */

  // Generate the net force arrays
  vtkFloatArray *gradXnet = vtkFloatArray::New();
  string nameX = string("X_NET");
  gradXnet->SetName(nameX.c_str());
  gradXnet->SetNumberOfComponents(3);
  gradXnet->SetNumberOfTuples(nCoeff / nc);
  poly->GetPointData()->AddArray(gradXnet);
  
  vtkFloatArray *gradRnet = vtkFloatArray::New();
  string nameR = string("R_NET");
  gradRnet->SetName(nameR.c_str());
  gradRnet->SetNumberOfComponents(1);
  gradRnet->SetNumberOfTuples(nCoeff / nc);
  poly->GetPointData()->AddArray(gradRnet);

  // For each term, generate a new data array
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    vtkFloatArray *gradX = vtkFloatArray::New();
    string nameX = string("X_") + xTerms[iTerm]->GetShortName();
    gradX->SetName(nameX.c_str());
    gradX->SetNumberOfComponents(3);
    gradX->SetNumberOfTuples(nCoeff / nc);
    poly->GetPointData()->AddArray(gradX);
    
    vtkFloatArray *gradR = vtkFloatArray::New();
    string nameR = string("R_") + xTerms[iTerm]->GetShortName();
    gradR->SetName(nameR.c_str());
    gradR->SetNumberOfComponents(1);
    gradR->SetNumberOfTuples(nCoeff / nc);
    poly->GetPointData()->AddArray(gradR);

    for(size_t i = 0; i < nCoeff / nc; i++)
      {
      SMLVec3d dX(
        xLastGradientPerTerm[iTerm][i * nc + 0],
        xLastGradientPerTerm[iTerm][i * nc + 1],
        xLastGradientPerTerm[iTerm][i * nc + 2]);
      SMLVec3d wdX = xWeights[iTerm] * dX;
      double dR = xLastGradientPerTerm[iTerm][i * nc + 3];
      double wdR = xWeights[iTerm] * dR;

      gradX->SetTuple3(i, wdX[0], wdX[1], wdX[2]);
      gradR->SetTuple1(i, wdR);

      if(iTerm == 0)
        {
        gradXnet->SetTuple3(i, wdX[0], wdX[1], wdX[2]);
        gradRnet->SetTuple1(i, wdR);
        }
      else
        {
        SMLVec3d XX = SMLVec3d(gradXnet->GetTuple3(i)) + wdX;
        gradXnet->SetTuple3(i, XX[0], XX[1], XX[2]);
        gradRnet->SetTuple1(i, gradRnet->GetTuple1(i) + wdR);
        }      }
    }

  vtkPolyDataWriter *write = vtkPolyDataWriter::New();
  write->SetInputData(poly);
  write->SetFileName(fnout);
  write->Update();

  delete fnout;  
}


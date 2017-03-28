#include "CoefficientMapping.h"
#include "PrincipalComponents.h"
#include "SubdivisionMedialModel.h"

PCACoefficientMapping
::PCACoefficientMapping(PrincipalComponents *pca, size_t nModes)
{ 
  this->pca = pca;
  this->n = nModes;
  this->m = pca->GetMean().size(); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::Apply(const Vec &C, const Vec &P)
{
  return C + pca->MapToFeatureSpaceZeroMean(P);
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
{ 
  return pca->MapToFeatureSpaceZeroMean(varP); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
{ 
  return varC; 
}

PCAPlusAffineCoefficientMapping
::PCAPlusAffineCoefficientMapping(
  GenericMedialModel *model, PrincipalComponents *pca, size_t nModes) :
  CompositionCoefficientMapping(
    new AffineTransformCoefficientMapping(model),
    new PCACoefficientMapping(pca, nModes))
{
}

PCAPlusAffineCoefficientMapping
::~PCAPlusAffineCoefficientMapping()
{
  delete this->f;
  delete this->g;
}
  

ReflectionCoefficientMapping
::ReflectionCoefficientMapping(
  SubdivisionMedialModel *model,
  SMLVec3d &p, double b)
{
  this->p = p;
  this->b = b;
  this->n = model->GetNumberOfCoefficients();
  this->m = model->GetNumberOfCoefficients();
  this->nc = model->GetNumberOfComponents();

  // Assuming 4 coefficients per node (standard for this code)
  this->k = m * nc;

  // Resize opp array
  opp.resize(k);

  // Compute the opposites map. This means for each point (u,v) in the model
  // to find the point closest to (u,-v). We do this by brute force here
  for(size_t i = 0; i < k; i++)
    {
    double ui = model->GetCoefficientU()[i];
    double vi = model->GetCoefficientV()[i];
    double dist = 0;
    size_t iopp = NOID;
    for(size_t j = 0; j < k; j++)
      {
      double uj = model->GetCoefficientU()[j];
      double vj = model->GetCoefficientV()[j];
      double d = (ui-uj)*(ui-uj) + (vi+vj)*(vi+vj);
      if(j == 0 || dist > d)
        { iopp = j; dist = d; }
      }
    opp[i] = iopp;
    // printf("Map (%f,%f) to (%f,%f)\n",ui,vi,
    //   model->GetCoefficientU()[opp[i]], model->GetCoefficientV()[opp[i]]);
    }
}

ReflectionCoefficientMapping::Vec
ReflectionCoefficientMapping
::Apply(const Vec &C, const Vec &P)
{
  Vec X = C+P, Y(n,0.0);
  for(size_t i = 0; i < k; i++)
    {
    size_t ci = i * nc, co = opp[i] * nc;

    // Get the vector component
    SMLVec3d Xi = X.extract(3, ci);
    SMLVec3d Xj = X.extract(3, co);
    SMLVec3d Yi = 0.5 * (Xi + Xj) - (dot_product(p,Xj) - b) * p;
    double Ri = X[ci + 3], Rj = X[co + 3];
    double Qi = 0.5 * (Ri + Rj);
    
    // Set the output
    Y.update(Yi, ci);
    Y[ci + 3] = Qi;
    }
  return Y;
}

/** Compute the variation J_T(P) * v_P given C, P and v_P */
ReflectionCoefficientMapping::Vec 
ReflectionCoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
{
  Vec dY(n, 0.0);
  for(size_t i = 0; i < k; i++)
    {
    size_t ci = i * nc, co = opp[i] * nc;

    SMLVec3d Vi = varC.extract(3,ci);
    SMLVec3d Vj = varC.extract(3,co);
    SMLVec3d dYi = 0.5 * (Vi + Vj) - dot_product(Vj,p) * p;
    double dQi = 0.5 * (varC[ci+3] + varC[co+3]);
    dY.update(dYi,ci);
    dY[ci+3] = dQi;
    }
  return dY;
}

ReflectionCoefficientMapping::Vec 
ReflectionCoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP)
{
  Vec dY(n, 0.0);
  for(size_t i = 0; i < k; i++)
    {
    size_t ci = i * nc, co = opp[i] * nc;

    SMLVec3d Vi = varP.extract(3,ci);
    SMLVec3d Vj = varP.extract(3,co);
    SMLVec3d dYi = 0.5 * (Vi + Vj) - dot_product(Vj,p) * p;
    double dQi = 0.5 * (varP[ci+3] + varP[co+3]);
    dY.update(dYi,ci);
    dY[ci+3] = dQi;
    }
  return dY;
}

/* ============================================================ */

#include <vnl/algo/vnl_generalized_eigensystem.h>

MeshBasisCoefficientMapping
::MeshBasisCoefficientMapping(const TriangleMesh *mesh, size_t basisSize, size_t nComponents)
  {
  // Get the characeteristics of the mesh
  this->nb = basisSize;
  this->nv = mesh->nVertices;
  this->nc = nComponents;

  // Number of parameters
  this->n = nb * nComponents;

  // Number of coefficients
  this->m = nv * nComponents;

  // Solve the problem A x = l B x
  // A is the laplace matrix (for now trivially implemented)
  // B is a diagonal weight matrix (hmm)

  // Create the adjacency matrix
  // TODO: use sparse code
  typedef vnl_matrix<double> Mat;
  Mat A(nv, nv, 0.0), B(nv, nv, 0.0);

  // Initialize the matrices
  for(size_t i = 0; i < nv; i++)
    {
    for(EdgeWalkAroundVertex it(mesh, i); !it.IsAtEnd(); ++it)
      {
      size_t j = it.MovingVertexId();
      A[i][i] += 1.0;
      A[i][j] -= 1.0;
      B[i][i] += 1.0;
      }
    }
  
  // Solve the eigenvalue problem
  vnl_generalized_eigensystem gev(A, B);


  // Get the first nb rows (this is the worst possible way to compute eigenvalues,
  // but it does not require linking to other libraries)
  this->V = gev.V.get_n_columns(0, nb).transpose();
  Vec test = gev.V.get_row(1);
  }

/** Apply the coefficient mapping C' = T(C, P) */
MeshBasisCoefficientMapping::Vec 
MeshBasisCoefficientMapping
::Apply(const Vec &C, const Vec &p)
  {
  Vec X = C;
  
  for(size_t iv = 0, j = 0; iv < nv; iv++)
    for(size_t ic = 0; ic < nc; ic++, j++)
      for(size_t ib = 0, ip = ic; ib < nb; ib++, ip+=nc)
        X[j] += p[ip] * V[ib][iv];
  
  return X;
  }

MeshBasisCoefficientMapping::Vec 
MeshBasisCoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP)
  {
  Vec dX(C.size(), 0.);

  for(size_t iv = 0, j = 0; iv < nv; iv++)
    for(size_t ic = 0; ic < nc; ic++, j++)
      for(size_t ib = 0, ip = ic; ib < nb; ib++, ip+=nc)
        dX[j] += varP[ip] * V[ib][iv];

  return dX;
  }

MeshBasisCoefficientMapping::Vec 
MeshBasisCoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
  {
  return varC;
  }


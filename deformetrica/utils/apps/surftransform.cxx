#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkPerspectiveTransform.h"

#include "itkTransformFileReader.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkVersorRigid3DTransform.h"

#include "surfio.h"

int main(int argc, char** argv)
{
	if (argc != 4)
	{
    		std::cerr << "Usage: " << argv[0] << " <in-poly> <out-poly> <itk-transformation-file>" << std::endl;
    		return -1;
  	}

	// Read the input surface mesh
  	vtkSmartPointer<vtkPolyData> poly = readSurface(argv[1]);
	
	//-------------------------------------------------------------------------------------
	// Read the versor transform and store matrix and offset
	//-------------------------------------------------------------------------------------

	typedef itk::TransformFileReader TransformReaderType;
	TransformReaderType::Pointer reader = TransformReaderType::New();
	reader->SetFileName(argv[3]);
	reader->Update();

	typedef itk::MatrixOffsetTransformBase< double, 3, 3 > TransformType;
	typedef itk::VersorRigid3DTransform<double> VersorTransformType;
		
	typedef TransformReaderType::TransformListType * TransformListType;
	TransformListType transforms = reader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator tit = transforms->begin();
    
	VersorTransformType::Pointer versor_read = static_cast<VersorTransformType*>((*tit).GetPointer());
    VersorTransformType::Pointer versor_transform =  dynamic_cast< VersorTransformType * >( versor_read.GetPointer() );

	VersorTransformType::MatrixType matrix = versor_transform->GetMatrix();
	VersorTransformType::OffsetType offset = versor_transform->GetOffset();

	//-------------------------------------------------------------------------------------
	// 4x4 matrices needed to compute the transform in RAS space
	//-------------------------------------------------------------------------------------

	vtkSmartPointer<vtkMatrix4x4> vtkmat = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkmat->Identity();

	vtkSmartPointer<vtkMatrix4x4> lps2ras = vtkSmartPointer<vtkMatrix4x4>::New();
	lps2ras->Identity();
	(*lps2ras)[0][0] = (*lps2ras)[1][1] = -1.0;

	vtkSmartPointer<vtkMatrix4x4> ras2lps = vtkSmartPointer<vtkMatrix4x4>::New();
	ras2lps->Identity();
	(*ras2lps)[0][0] = (*ras2lps)[1][1] = -1.0;

	// Load the versor transform into the 4x4 matrix
	for (int i=0; i < 3; i++)
	{
		for (int j=0; j < 3; j++)
		{
			(*vtkmat)[i][j] = matrix[i][j];
		}
		(*vtkmat)[i][3] = offset[i];
	}

	vtkMatrix4x4::Multiply4x4(lps2ras, vtkmat, vtkmat);
    vtkMatrix4x4::Multiply4x4(vtkmat, ras2lps, vtkmat);
	vtkmat->Invert();

	//-------------------------------------------------------------------------------------
	// Apply the 4x4 matrix transformation using a perspective transform
	//-------------------------------------------------------------------------------------

	vtkSmartPointer<vtkPerspectiveTransform> perspectiveTransform = vtkSmartPointer<vtkPerspectiveTransform>::New();
  	perspectiveTransform->SetMatrix(vtkmat);

	vtkSmartPointer<vtkPolyData> outPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Loop over the points and transform them
	for(vtkIdType i = 0; i < poly->GetNumberOfPoints(); i++)
	{
		double p[3], newp[3];
		poly->GetPoint(i,p);
		
		perspectiveTransform->TransformPoint(p, newp);
		points->InsertPoint(i, newp);
	}
	outPoly->SetPoints(points);
	outPoly->SetPolys(poly->GetPolys());

	// Write the output surface mesh
	writeSurface(argv[2], outPoly);
}

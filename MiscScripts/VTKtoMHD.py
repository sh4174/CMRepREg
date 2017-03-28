# VTK to MHD

import vtk

def vtk2mhd( vtkPath, outImgPath, refImgPath ):
		# Read Poly Data
	polyReader = vtk.vtkPolyDataReader()
	polyReader.SetFileName( vtkPath )
	polyReader.Update()
	
	# Poly Data
	polyData = polyReader.GetOutput()

	# Flip Poly Data 
	nPt = polyData.GetNumberOfPoints()
	flippedPts = vtk.vtkPoints()

	for i in range( nPt ):
		pt = polyData.GetPoint( i )
		flippedPts.InsertNextPoint( pt[0], pt[1], pt[2] )

	polyData.SetPoints( flippedPts )

	# Output Label Image (MHD)
	mhdImg = vtk.vtkImageData()

	if refImgPath == 0:
		mhdImg.SetOrigin( [ -100, -100, -100 ] )
		mhdImg.SetSpacing( [ 1, 1, 1 ] )
		mhdImg.SetDimensions( [ 200, 200, 200 ] )
		mhdImg.SetExtent( 0, mhdImg.GetDimensions()[0] - 1, 0, mhdImg.GetDimensions()[1] - 1, 0, mhdImg.GetDimensions()[2] - 1 )
	else:		
		# Reference Image Data with Spacing/Origin/Dimension
		refReader = vtk.vtkMetaImageReader()
		refReader.SetFileName( refImgPath )
		refReader.Update()

		# Reference Image 
		refImg = refReader.GetOutput()

		# Copy Meta Information from Reference Image
		mhdImg.SetOrigin( refImg.GetOrigin() )
		mhdImg.SetSpacing( refImg.GetSpacing() )
		mhdImg.SetDimensions( refImg.GetDimensions() )
		mhdImg.SetExtent( 0, refImg.GetDimensions()[0] - 1, 0, refImg.GetDimensions()[1] - 1, 0, refImg.GetDimensions()[2] - 1 )

	# Allocate Out Image
	mhdImg.SetScalarTypeToUnsignedChar()
	mhdImg.AllocateScalars()

	count = mhdImg.GetNumberOfPoints()
	for i in range( count ):
		mhdImg.GetPointData().GetScalars().SetTuple1( i, 1 )
	mhdImg.Update()


	# Poly Data to Image Stencil
	stenciler = vtk.vtkPolyDataToImageStencil()
	stenciler.SetInput( polyData )
	stenciler.SetOutputOrigin( mhdImg.GetOrigin() )
	stenciler.SetOutputSpacing( mhdImg.GetSpacing() )
	stenciler.SetOutputWholeExtent( mhdImg.GetExtent() )
	stenciler.Update()

	# Cut the corresponding white image and set the background
	stencil = vtk.vtkImageStencil()
	stencil.SetInput( mhdImg )
	stencil.SetStencil( stenciler.GetOutput() )
	stencil.ReverseStencilOff()
	stencil.SetBackgroundValue( 0 )
	stencil.Update()

	# Write Label Image 
	writer = vtk.vtkMetaImageWriter()
	writer.SetInput( stencil.GetOutput() )
	writer.SetFileName( outImgPath )
	writer.Write()


def main():
	folderPath = "/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepTest/hd_example/GroupRegression/Cont/"
	filePrefix = "Regression_baseline_Atlas_template_init_basline_shape_left_caudate_ctrl_to_subject_0__t_9_trajectory___t_"
	for i in range( 66 ):
		vtkPath = folderPath + filePrefix + str( i ) +".vtk"
		outImgPath = folderPath + filePrefix + str( i ) +".mha"
				
		vtk2mhd( vtkPath, outImgPath, 0 )

	# for i in range( 1, 4 ):
	# 	vtkPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/CMRepTest/hd_example/Synthetic2/Data/extrinsic_ellipse" + str( i ) +".vtk"
	# 	outImgPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/CMRepTest/hd_example/Synthetic2/Data/extrinsic_ellipse" + str( i ) +".mha"
				
	# 	vtk2mhd( vtkPath, outImgPath, 0 )




if __name__ == '__main__':
	main()
 

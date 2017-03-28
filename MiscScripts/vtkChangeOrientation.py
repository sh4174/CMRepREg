import vtk


iPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/left_putamen_00.vtk"
oPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/left_putamen_00.vtk"

reader = vtk.vtkPolyDataReader()
reader.SetFileName( iPath )
reader.Update()

data = reader.GetOutput()

nPt = data.GetNumberOfPoints()


oData = vtk.vtkPolyData()
oPts = vtk.vtkPoints()

for i in range( nPt ):
	pt = data.GetPoint( i )
	ptR = [ -pt[0], -pt[1], pt[2] ]

	oPts.InsertNextPoint( ptR ) 

oData.SetPoints( oPts )
oData.SetPolys( data.GetPolys() )

writer = vtk.vtkPolyDataWriter()
writer.SetFileName( oPath )
writer.SetInput( oData )
writer.Write()

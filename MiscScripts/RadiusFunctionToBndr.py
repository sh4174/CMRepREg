# Change Radius, Rho Function Name
import vtk

vtkBndPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/left_putamen_00.vtk"
vtkPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/lp_med_00.vtk"

reader = vtk.vtkPolyDataReader()
reader.SetFileName( vtkPath )
reader.Update()
polyData = reader.GetOutput()

radiusFuncArr = polyData.GetPointData().GetArray( "Radius" )
rhoFuncArr = polyData.GetPointData().GetArray( "Rho" )

radiusArr = vtk.vtkFloatArray()
radiusArr.DeepCopy( radiusFuncArr )
radiusArr.SetName( "Radius" )

rhoArr = vtk.vtkFloatArray()
rhoArr.DeepCopy( rhoFuncArr )
rhoArr.SetName( "Rho" )

bndReader = vtk.vtkPolyDataReader()
bndReader.SetFileName( vtkBndPath )
bndReader.Update()

bndData = bndReader.GetOutput()

bndData.GetFieldData().AddArray( radiusArr )
bndData.GetFieldData().AddArray( rhoArr )
bndData.Update()

writer = vtk.vtkPolyDataWriter()
writer.SetFileName( vtkBndPath )
writer.SetInput( bndData )

writer.Update()
writer.Write()

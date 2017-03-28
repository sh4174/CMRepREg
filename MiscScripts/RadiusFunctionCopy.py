# Change Radius, Rho Function Name
import vtk

vtkPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/lp_med_02.vtk"

reader = vtk.vtkPolyDataReader()
reader.SetFileName( vtkPath )
reader.Update()

polyData = reader.GetOutput()

radiusFuncArr = polyData.GetPointData().GetArray( "Radius Function" )
rhoFuncArr = polyData.GetPointData().GetArray( "Rho Function" )
textCoordArr = polyData.GetPointData().GetArray( "Texture Coordinates" )


radiusArr = vtk.vtkFloatArray()
radiusArr.DeepCopy( radiusFuncArr )
radiusArr.SetName( "Radius" )

rhoArr = vtk.vtkFloatArray()
rhoArr.DeepCopy( rhoFuncArr )
rhoArr.SetName( "Rho" )

tCoordArr = vtk.vtkFloatArray()
tCoordArr.DeepCopy( textCoordArr )
tCoordArr.SetName( "tcoords" )

polyData.GetPointData().RemoveArray( "Area Element" )
polyData.GetPointData().RemoveArray( "Atom Normal" )
polyData.GetPointData().RemoveArray( "Bending Energy" )
polyData.GetPointData().RemoveArray( "Covariant Tensor Determinant" )
polyData.GetPointData().RemoveArray( "Curvature Penalty Feature" )
polyData.GetPointData().RemoveArray( "Dummy1" )
polyData.GetPointData().RemoveArray( "Gauss Curvature" )
polyData.GetPointData().RemoveArray( "Grad R Magnitude (original)" )
polyData.GetPointData().RemoveArray( "GradR" )
polyData.GetPointData().RemoveArray( "Kappa1" )
polyData.GetPointData().RemoveArray( "Kappa2" )
polyData.GetPointData().RemoveArray( "LaplaceBasis" )
polyData.GetPointData().RemoveArray( "Mean Curvature" )
polyData.GetPointData().RemoveArray( "Metric Angle" )
polyData.GetPointData().RemoveArray( "normals" )
polyData.GetPointData().RemoveArray( "Off Diagonal Term of Contravariant MT" )
polyData.GetPointData().RemoveArray( "Phi" )
polyData.GetPointData().RemoveArray( "Radius Function" )
polyData.GetPointData().RemoveArray( "Regularity Penalty" )
polyData.GetPointData().RemoveArray( "Rho Function" )
polyData.GetPointData().RemoveArray( "Rs2" )
polyData.GetPointData().RemoveArray( "Spoke1" )
polyData.GetPointData().RemoveArray( "Spoke2" )
polyData.GetPointData().RemoveArray( "Stretch" )
polyData.GetPointData().RemoveArray( "Texture Coordinates" )
polyData.GetPointData().RemoveArray( "U Coordinate" )
polyData.GetPointData().RemoveArray( "V Coordinate" )
polyData.GetPointData().RemoveArray( "Xu" )
polyData.GetPointData().RemoveArray( "Xv" )

polyData.GetPointData().RemoveArray( "Radius" )
polyData.GetPointData().RemoveArray( "Rho" )
polyData.GetPointData().RemoveArray( "tcoords" )


polyData.GetPointData().AddArray( radiusArr )
polyData.GetPointData().AddArray( rhoArr )
polyData.GetPointData().AddArray( tCoordArr )

polyData.Update()

writer = vtk.vtkPolyDataWriter()
writer.SetFileName( vtkPath )
writer.SetInput( polyData )

writer.Update()
writer.Write()







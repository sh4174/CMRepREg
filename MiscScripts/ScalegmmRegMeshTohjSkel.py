import os
import vtk
import math
import numpy as np
import matplotlib.pyplot as plt

def main():
	# HamiltonJacobi Dimensions
	gmmPolyDataPath = "/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepTest/InitializationMatching/def1_org.vtk"
	gmmReader = vtk.vtkPolyDataReader()
	gmmReader.SetFileName( gmmPolyDataPath )
	gmmReader.Update()
	gmmEllipse = gmmReader.GetOutput()

	# Scale-Origin Matched HJ Skeleton
	regHJSkelPath = "/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepTest/InitializationMatching/reg_HJSkel_03.vtk"
	hjReader = vtk.vtkPolyDataReader()
	hjReader.SetFileName( regHJSkelPath )
	hjReader.Update()
	hjSkel = hjReader.GetOutput()

	# Output
	gmmPolyDataPathS = "/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepTest/InitializationMatching/def1_org_a.vtk"
	gmmPoints = vtk.vtkPoints()

	# Compute HJ poly surf bound
	hjSkel.ComputeBounds()
	hjSkel.Update()
	hjBounds = hjSkel.GetBounds()

	print hjBounds
	hjOrg = [ hjBounds[ 0 ], hjBounds[ 2 ], hjBounds[ 4 ] ]
	hjDim = [ hjBounds[ 1 ] - hjBounds[ 0 ], hjBounds[ 3 ] - hjBounds[ 2 ], hjBounds[ 5 ] - hjBounds[ 4 ] ]

	hjCM = [ 0, 0, 0 ]

	for i in range( hjSkel.GetNumberOfPoints() ):
		hjPt = hjSkel.GetPoint( i )
		hjCM[ 0 ] = hjCM[ 0 ] + hjPt[ 0 ]
		hjCM[ 1 ] = hjCM[ 1 ] + hjPt[ 1 ]
		hjCM[ 2 ] = hjCM[ 2 ] + hjPt[ 2 ]

	hjCM[ 0 ] = hjCM[ 0 ] / hjSkel.GetNumberOfPoints()
	hjCM[ 1 ] = hjCM[ 1 ] / hjSkel.GetNumberOfPoints()
	hjCM[ 2 ] = hjCM[ 2 ] / hjSkel.GetNumberOfPoints()

	# Compute Ellipse Bound
	gmmEllipse.ComputeBounds()
	gmmEllipse.Update()
	gmmBounds = gmmEllipse.GetBounds()

	print gmmBounds
	gmmOrg = [ gmmBounds[ 0 ], gmmBounds[ 2 ], gmmBounds[ 4 ] ]
	gmmDim = [ gmmBounds[ 1 ] - gmmBounds[ 0 ], gmmBounds[ 3 ] - gmmBounds[ 2 ], gmmBounds[ 5 ] - gmmBounds[ 4 ] ]

	for i in range( gmmEllipse.GetNumberOfPoints() ):
		pt = gmmEllipse.GetPoint( i )
		ptReg = [ 0, 0, 0 ]
		ptReg[ 0 ] = ( pt[ 1 ] - gmmOrg[ 1 ] ) / gmmDim[ 1 ] * hjDim[ 0 ] + hjOrg[ 0 ]
		ptReg[ 1 ] = ( pt[ 0 ] - gmmOrg[ 0 ] ) / gmmDim[ 0 ] * hjDim[ 1 ] + hjOrg[ 1 ]
		ptReg[ 2 ] = ( pt[ 2 ] - gmmOrg[ 2 ] ) / gmmDim[ 0 ] * hjDim[ 2 ] + hjOrg[ 2 ]
		gmmPoints.InsertNextPoint( ptReg[ 0 ], ptReg[ 1 ], ptReg[ 2 ] )

	gmmEllipse.SetPoints( gmmPoints )
	geCM = [ 0, 0, 0 ]

	for i in range( gmmEllipse.GetNumberOfPoints() ):
		pt = gmmEllipse.GetPoint( i )
		geCM[ 0 ] = geCM[ 0 ] + pt[ 0 ]		
		geCM[ 1 ] = geCM[ 1 ] + pt[ 1 ]		
		geCM[ 2 ] = geCM[ 2 ] + pt[ 2 ]		

	geCM[ 0 ] = geCM[ 0 ] / gmmEllipse.GetNumberOfPoints()
	geCM[ 1 ] = geCM[ 1 ] / gmmEllipse.GetNumberOfPoints()
	geCM[ 2 ] = geCM[ 2 ] / gmmEllipse.GetNumberOfPoints()


	# Match Rotation
	gmmPointsMC = vtk.vtkPoints()

	for i in range( gmmEllipse.GetNumberOfPoints() ):
		pt = gmmEllipse.GetPoint( i )
		pt2 = [ 0, 0, 0 ]
		pt2[ 0 ] = ( pt[ 0 ] - geCM[ 0 ] )
		pt2[ 1 ] = -( pt[ 1 ] - geCM[ 1 ] )
		pt2[ 2 ] = ( pt[ 2 ] - geCM[ 2 ] )

		gmmPointsMC.InsertNextPoint( pt2 )
	gmmEllipse.SetPoints( gmmPointsMC )

	transform = vtk.vtkTransform()
	transform.RotateWXYZ( 45, 1, 0, 0 )

	trFilter = vtk.vtkTransformPolyDataFilter()
	trFilter.SetTransform( transform )
	trFilter.SetInput( gmmEllipse )
	trFilter.Update()

	# Match Center
	gmmPointsMC2 = vtk.vtkPoints()
	for i in range( gmmEllipse.GetNumberOfPoints() ):
		pt = trFilter.GetOutput().GetPoint( i )
		pt2 = [ 0, 0, 0 ]
		pt2[ 0 ] = pt[ 0 ] + hjCM[ 0 ]
		pt2[ 1 ] = pt[ 1 ] + hjCM[ 1 ]
		pt2[ 2 ] = pt[ 2 ] + hjCM[ 2 ]

		gmmPointsMC2.InsertNextPoint( pt2 )
	gmmEllipse.SetPoints( gmmPointsMC2 )

	writer = vtk.vtkPolyDataWriter()
	writer.SetInput( gmmEllipse )
	writer.SetFileName( gmmPolyDataPathS )
	writer.Update()
	writer.Write()

if __name__ == '__main__':
	main()
 
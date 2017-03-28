# Change Radius, Rho Function Name
import vtk
import numpy as np
from scipy import interpolate

# #Cont
# folderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/50368_Cont/"

# Synth
extfolderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/"
nofolderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/" 

folderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/"

nObs = 3

# High
timeList = [ 42.23682409, 43.00342231, 46.01505818 ]

# # Cont
# timeList = [ 43.61670089, 45.76317591, 46.68309377 ]

# # Synth
# timeList = [ 1.0, 2.0 ]

obsFileNameList = [ 'left_putamen_00.vtk', 'left_putamen_01.vtk', 'left_putamen_02.vtk' ]

# Calculte Time Step
nTimePt = 30
timeRange = timeList[ len( timeList ) - 1 ] - timeList[ 0 ] 
timeStep = timeRange / ( nTimePt - 1 )

# Create Index List
idxList = []
for i in range( len( timeList ) ):
	idxList.append( round( ( timeList[ i ] - timeList[ 0 ] ) / timeStep ) )

timePtList = np.arange( 0, nTimePt )
print timePtList

obsList = []

# Read Radius Function Data
for fileName in obsFileNameList:
	filePath = extfolderPath + fileName

	print filePath

	reader = vtk.vtkPolyDataReader()
	reader.SetFileName( filePath )
	reader.Update()
	polyData = reader.GetOutput()
	obsList.append( polyData )

radArrList = []
rhoArrList = []

# Read Radius Array Function
for i in range( nTimePt ):
	radArr = vtk.vtkFloatArray()
	radArrList.append( radArr )

	rhoArr = vtk.vtkFloatArray()
	rhoArrList.append( rhoArr )

# Get Number of Points
nPtData = obsList[ 0 ].GetFieldData().GetArray( "Radius" ).GetNumberOfTuples()
print nPtData 

radArrList_reg = []
rhoArrList_reg = []

# Create Radius/Rho Arrays for Regression Results
for i in range( nTimePt ):
	radArr_i = vtk.vtkFloatArray()
	radArr_i.SetName( "Radius" )
	radArr_i.SetNumberOfComponents( 1 )
	radArrList_reg.append( radArr_i )

	rhoArr_i = vtk.vtkFloatArray()
	rhoArr_i.SetName( "Rho" )
	rhoArr_i.SetNumberOfComponents( 1 )
	rhoArrList_reg.append( rhoArr_i )

# Interpolate Rad/Rho 
for i in range( nPtData ):
	radList_i = []
	rhoList_i = []

	for j in range( nObs ):
		radList_i.append( obsList[ j ].GetFieldData().GetArray("Radius").GetValue( i ) )
		rhoList_i.append( obsList[ j ].GetFieldData().GetArray("Rho").GetValue( i ) )

	f_rad = interpolate.interp1d( idxList, radList_i, 'linear' )
	f_rho = interpolate.interp1d( idxList, rhoList_i, 'linear' )

	radList_i_interp = f_rad( timePtList )
	rhoList_i_interp = f_rho( timePtList )

	for j in range( nTimePt ):
		radArrList_reg[ j ].InsertNextValue( radList_i_interp[ j ] )	
		rhoArrList_reg[ j ].InsertNextValue( rhoList_i_interp[ j ] )	

# Assign Rad/Rho Array to Regression Results
regPrefix = "Regression_baseline_lc_med_00_trajectory___t_"
for i in range( nTimePt ):
	filePath_i = nofolderPath + regPrefix + str( i ) + ".vtk"
	outfilePath_i = folderPath + regPrefix + str( i ) + ".vtk"

	reader_i = vtk.vtkPolyDataReader()
	reader_i.SetFileName( filePath_i )
	reader_i.Update()

	regData = reader_i.GetOutput()

	regData.GetPointData().AddArray( radArrList_reg[ i ] )
	regData.GetPointData().AddArray( rhoArrList_reg[ i ] )
		

	writer_i = vtk.vtkPolyDataWriter()
	writer_i.SetFileName( outfilePath_i )
	writer_i.SetInput( regData )
	writer_i.Write()



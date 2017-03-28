# Reconstruct Boundary 
import subprocess
import os

# # Cont
# folderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/50368_Cont/"

# Synth
folderPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/DeformetricaTest/CMRep_Regression_Test/TestData/52351_High_Multi/"

regPrefix = "Regression_baseline_lp_med_00_trajectory___t_"
reconBinPath = "/media/shong/IntHard1/4DAnalysis/Code/SPTSkeleton/CMRepTest/cmrep_recon_bin/ReconstructBoundary"

nTimePt = 30;

for i in range( nTimePt ):
	cmrep_filename = "recon" + str( i ) + ".cmrep"
	cmrep_filePath = folderPath + cmrep_filename

	cmrep_file_i = open( cmrep_filePath, "w" )
	filePath_i = folderPath + regPrefix + str( i ) + ".vtk"
	fileReconPath_i = folderPath + regPrefix +"_recon_" + str( i ) + ".vtk"

	cmrep_file_i.write( "Grid.Type = LoopSubdivision\n" )
	cmrep_file_i.write( "Grid.Model.SolverType = PDE\n" )	
	cmrep_file_i.write( "Grid.Model.Atom.SubdivisionLevel = 0\n" )	
	cmrep_file_i.write( "Grid.Model.Coefficient.FileName = " + filePath_i + "\n" )	
	cmrep_file_i.write( "Grid.Model.Coefficient.FileType = VTK\n" )	
	cmrep_file_i.close()

	reconCmd = [ reconBinPath, "-i", cmrep_filePath, "-o", fileReconPath_i ]
	p = subprocess.Popen( reconCmd ).wait()

	os.remove( cmrep_filePath )


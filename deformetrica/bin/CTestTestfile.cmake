# CMake generated Testfile for 
# Source directory: /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/app
# Build directory: /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(surface_matching "sparseMatching3" "paramDiffeos.xml" "paramSurface.xml" "sourceSurfaceMesh.vtk" "targetSurfaceMesh.vtk")
set_tests_properties(surface_matching PROPERTIES  WORKING_DIRECTORY "../../examples/surface_matching/")
subdirs(/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611)

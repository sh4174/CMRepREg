#!/bin/bash

CMREPREG_PATH="/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin"

$CMREPREG_PATH/sparseGeodesicRegression3CMRep paramDiffeos.xml 1 paramSurface.xml cmrep_med_00.vtk cmrep_pair_00.dat 43.61670089 cmrep_pair_01.dat 45.76317591 cmrep_pair_02.dat 46.68309377

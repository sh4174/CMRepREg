#!/bin/bash


/media/shong/IntHard1/Projects/4DShapeAnalysis/Code/deformetrica-2.1/deformetrica/bin/sparseMatching3 paramDiffeos.xml paramSurface.xml def1_org_a.vtk reg_HJSkel_03.vtk

mkdir Results

cp def1_org_a__t_19.vtk Results/def1_org.vtk


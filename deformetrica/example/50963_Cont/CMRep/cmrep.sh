#!/bin/bash

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_caudate_00.nii.gz lc_00

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_caudate_01.nii.gz lc_01

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_caudate_02.nii.gz lc_02

mkdir TimeSeries

cp lc_00/mesh/def2.med.vtk TimeSeries/cmrep_med_00.vtk

cp lc_01/mesh/def2.med.vtk TimeSeries/cmrep_med_01.vtk

cp lc_02/mesh/def2.med.vtk TimeSeries/cmrep_med_02.vtk

cp lc_00/mesh/def2.bnd.vtk TimeSeries/cmrep_bnd_00.vtk

cp lc_01/mesh/def2.bnd.vtk TimeSeries/cmrep_bnd_01.vtk

cp lc_02/mesh/def2.bnd.vtk TimeSeries/cmrep_bnd_02.vtk




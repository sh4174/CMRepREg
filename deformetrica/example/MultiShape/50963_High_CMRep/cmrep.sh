#!/bin/bash

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_putamen_00.nii.gz left_putamen_00

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_putamen_01.nii.gz left_putamen_01

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep left_putamen_02.nii.gz left_putamen_02

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_putamen_00.nii.gz right_putamen_00

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_putamen_01.nii.gz right_putamen_01

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_putamen_02.nii.gz right_putamen_02

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_caudate_00.nii.gz right_caudate_00

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_caudate_01.nii.gz right_caudate_01

../../cmrep_eclipse_bin/cmrep_fit cmrparam_init.txt def1.cmrep right_caudate_02.nii.gz right_caudate_02

mkdir multishapes

cp left_putamen_00/mesh/def2.med.vtk multishapes/lp_med_00.vtk

cp left_putamen_01/mesh/def2.med.vtk multishapes/lp_med_01.vtk

cp left_putamen_02/mesh/def2.med.vtk multishapes/lp_med_02.vtk

cp left_putamen_00/mesh/def2.bnd.vtk multishapes/lp_bnd_00.vtk

cp left_putamen_01/mesh/def2.bnd.vtk multishapes/lp_bnd_01.vtk

cp left_putamen_02/mesh/def2.bnd.vtk multishapes/lp_bnd_02.vtk


cp right_putamen_00/mesh/def2.med.vtk multishapes/rp_med_00.vtk

cp right_putamen_01/mesh/def2.med.vtk multishapes/rp_med_01.vtk

cp right_putamen_02/mesh/def2.med.vtk multishapes/rp_med_02.vtk

cp right_putamen_00/mesh/def2.bnd.vtk multishapes/rp_bnd_00.vtk

cp right_putamen_01/mesh/def2.bnd.vtk multishapes/rp_bnd_01.vtk

cp right_putamen_02/mesh/def2.bnd.vtk multishapes/rp_bnd_02.vtk


cp right_caudate_00/mesh/def2.med.vtk multishapes/rc_med_00.vtk

cp right_caudate_01/mesh/def2.med.vtk multishapes/rc_med_01.vtk

cp right_caudate_02/mesh/def2.med.vtk multishapes/rc_med_02.vtk

cp right_caudate_00/mesh/def2.bnd.vtk multishapes/rc_bnd_00.vtk

cp right_caudate_01/mesh/def2.bnd.vtk multishapes/rc_bnd_01.vtk

cp right_caudate_02/mesh/def2.bnd.vtk multishapes/rc_bnd_02.vtk




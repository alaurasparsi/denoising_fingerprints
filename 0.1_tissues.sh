#!/usr/bin/env bash
# Author: Laura Belli
# segment tissues
# make sure you previously registered the segmented T1w to the functional sbref
ddir=/mnt/data
pdir=/mnt/derivatives
polort=4
for sub in $(seq -f %03g 1 10); do
    anat=${ddir}/sub-${sub}/ses-01/anat/sub-${sub}_ses-01_acq-uni_T1w
    func=${ddir}/sub-${sub}/ses-01/func/sub-${sub}_ses-01_task-rest_run-01_optcom_bold
    ddir=/mnt/StefHDD/EuskalIBUR/derivatives/physiodenoise/data_correct
    3dcalc -a ${anat}_seg2sbref.nii.gz -expr 'equals(a,1)' -prefix ${anat}_CSF.nii.gz -overwrite
    3dcalc -a ${anat}_seg2sbref.nii.gz -expr 'equals(a,3)' -prefix ${anat}_WM.nii.gz -overwrite
    3dcalc -a ${anat}_seg2sbref.nii.gz -expr 'equals(a,2)' -prefix ${anat}_GM.nii.gz -overwrite
    
    # use CSF and WM reference masks to get average values
    for ses in $(seq -f %02g 1 10); do
        func_in=${ddir}/sub-${sub}/ses-${ses}/func/00.sub-${sub}_ses-${ses}_task-rest_run-01_optcom_bold_native_preprocessed
        func=${ddir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-rest_run-01_optcom_bold
        fslmeants -i ${func_in}.nii.gz -m ${anat}_CSF.nii.gz -o ${func}_avg_CSF.1D
        fslmeants -i ${func_in}.nii.gz -m ${anat}_WM.nii.gz -o ${func}_avg_WM.1D
    done
done
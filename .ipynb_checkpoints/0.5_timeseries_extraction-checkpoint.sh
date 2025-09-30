#!/usr/bin/env bash
#Author: Laura Belli
# timeseries extraction
models=("md" "md_bp" "pd" "pd_bp" "Sara_denoised")
ddir=/mnt/StefHDD/data
for sub in $(seq -f %03g 1 10); do

    echo "processing sub-${sub}"

    for ses in $(seq -f %02g 1 10); do
        fdir=/mnt/data/sub-${sub}/ses-${ses}/func

         for model in "${models[@]}"; do
            3dROIstats -quiet -mask ${ddir}/sub-${sub}_Shendil2sbref_masked.nii.gz ${fdir}/sub-${sub}_ses-${ses}_task-rest_run-01_optcom_bold_${model}.nii.gz > ${fdir}/sub-${sub}_ses-${ses}_task-rest_run-01_Shen_dil_timeseries_${model}.1D
        done
    done
done
# !/usr/bin/env bash
# to register the atlas to the functional single band reference
wdr=/mnt/StefHDD/EuskalIBUR/derivatives/physiodenoise/data_correct
cwd=/mnt/StefHDD/EuskalIBUR/derivatives/physiodenoise
# Schaefer 400 atlas
for sub in $(seq -f "%03g" 1 10 ); do

    echo "Registering Schaefer 400 atlas in sub-${sub} functional space"

    antsApplyTransforms -d 3 \
        -i "${cwd}/ref_atlas/schaefer400MNI.nii" \
        -r "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2sbref.nii.gz" \
        -o "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_Schaefer400MNI2sbref.nii.gz" \
        -n MultiLabel \
        -t "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_acq-uni_T1w2std1InverseWarp.nii.gz" \
        -t "[${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_acq-uni_T1w2std0GenericAffine.mat,1]" \
        -t "[${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2acq-uni_T1w0GenericAffine.mat,1]" \
        -t "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2sbref0GenericAffine.mat" \
        -v
done
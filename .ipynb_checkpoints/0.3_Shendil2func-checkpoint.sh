# !/usr/bin/env bash
# Author: Laura Belli, Stefano Moia
# register the atlas to the functional single band reference
wdr=/mnt/data
cwd=/mnt/derivatives
# Shen dilated atlas
for sub in $(seq -f %03g 1 10); do

    echo "Registering Shen dilated atlas in sub-${sub} functional space"

    antsApplyTransforms -d 3 -i "${cwd}/ref_atlas/shen_dil.nii" \
        -r "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_sbref.nii.gz" \
        -o "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_Shendil2sbref.nii.gz" \
        -n MultiLabel \
        -t "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2sbref0GenericAffine.mat" \
        -t "[${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2acq-uni_T1w0GenericAffine.mat,1]" \
        -t "[${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_acq-uni_T1w2std0GenericAffine.mat,1]" \
        -t "${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_ses-01_acq-uni_T1w2std1InverseWarp.nii.gz" \
        -v
        
# mask atlas
    fslmaths ${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_Shendil2sbref.nii.gz -mas ${wdr}/sub-${sub}/ses-01/anat/sub-${sub}_ses-01_acq-uni_T1w_GM.nii.gz ${wdr}/sub-${sub}/ses-01/reg/sub-${sub}_Shendil2sbref_masked.nii.gz
    
done

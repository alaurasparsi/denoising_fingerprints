# !/usr/bin/env bash
# Author: Laura Belli
# to register the segmented T1w in the functional single band reference
asegfx=acq-uni_T1w

for sub in $(seq -f %03g 1 10); do
    aseg=sub-${sub}/ses-01/anat/sub-${sub}_ses-01_acq-uni_T1w
	anat=sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w
	mref=sub-${sub}/ses-01/reg/sub-${sub}_sbref
	anat2mref=sub-${sub}/ses-01/reg/sub-${sub}_ses-01_T2w2sbref0GenericAffine

	antsApplyTransforms -d 3 -i ${aseg}_seg.nii.gz \
					-r ${mref}.nii.gz -o ${aseg}_seg2sbref.nii.gz \
					-n Multilabel -v \
					-t ${anat2mref}.mat \
					-t [${anat}2${asegfx}0GenericAffine.mat,1]
done

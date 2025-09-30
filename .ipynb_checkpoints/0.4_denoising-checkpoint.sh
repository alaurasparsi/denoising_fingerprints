#!/usr/bin/env bash
#Author: Laura Belli, Stefano Moia
# physiological denoising: 5 models
ddir=/mnt/data
pdir=/mnt/derivatives
polort=4

for sub in $(seq -f %03g 1 10); do
    for ses in $(seq -f %02g 1 10); do
        mref=${ddir}/sub-${sub}/ses-${ses}/reg/sub-${sub}_sbref
        RVT=${pdir}/RVT/correct/sub-${sub}_ses-${ses}_task-rest_run-01_RVT_resampled_lag-0
        HRV=${pdir}/HRV/correct/sub-${sub}_ses-${ses}_task-rest_run-01_HRV_resampled_convolved
        func=${ddir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-rest_run-01_optcom_bold
        func_in=${ddir}/sub-${sub}/ses-${ses}/func/00.sub-${sub}_ses-${ses}_task-rest_run-01_optcom_bold_native_preprocessed
        RETROICOR=${pdir}/RETROICOR/sub-${sub}_ses-${ses}_task-rest_run-01_RETROICOR_regressors
        # prepare "physiologically denoised" WM and CSF regressors
        1dtranspose ${func}_avg_CSF.1D > ${tissue}_avg_CSF_f.1D
        3dTproject -polort ${polort} -input ${tissue}_avg_CSF_f.1D \
        -ort ${RVT}.1D \
        -ort ${HRV}.1D \
        -prefix ${tissue}_CSF_reg.1D

        1dtranspose ${tissue}_avg_WM.1D > ${tissue}_avg_WM_f.1D
        3dTproject -polort ${polort} -input ${tissue}_avg_WM_f.1D \
        -ort ${RVT}.1D \
        -ort ${HRV}.1D \
        -prefix ${tissue}_WM_reg.1D
        # 1.
        # minimally denoised: motor regressors (6), Legendre polynomials, WM and CSF
        3dTproject -polort ${polort} -input ${func_in}.nii.gz  -mask ${mref}_brain_mask.nii.gz \
        -ort ${func}_nuisreg_uncensored_mat.1D \
        -ort ${func}_CSF_reg_f.1D \
        -ort ${func}_WM_reg_f.1D \
        -prefix ${func}_md.nii.gz -overwrite
        # 2.
        # minimally denoised with bandpass
        3dTproject -polort ${polort} -input ${func_in}.nii.gz  -mask ${mref}_brain_mask.nii.gz \
        -ort ${func}_nuisreg_uncensored_mat.1D \
        -ort ${func}_CSF_reg_f.1D \
        -ort ${func}_WM_reg_f.1D \
        -bandpass 0.01 0.15 \
        -prefix ${func}_md_bp.nii.gz -overwrite
        # 3.
        # physiologically denoised: motor regressors (6), Legendre polynomials, WM and CSF, HRV, RVT, RETROICOR
        3dTproject -polort ${polort} -input ${func_in}.nii.gz  -mask ${mref}_brain_mask.nii.gz \
        -ort ${func}_nuisreg_uncensored_mat.1D \
        -ort ${func}_CSF_reg_f.1D \
        -ort ${func}_WM_reg_f.1D \
        -ort ${HRV}.1D \
        -ort ${RVT}.1D \
        -ort ${RETROICOR}.1D \
        -prefix ${func}_pd.nii.gz
        # 4.
        # physiologically denoised with bandpass
        3dTproject -polort ${polort} -input ${func_in}.nii.gz  -mask ${mref}_brain_mask.nii.gz \
        -ort ${func}_nuisreg_uncensored_mat.1D \
        -ort ${func}_CSF_reg_f.1D \
        -ort ${func}_WM_reg_f.1D \
        -ort ${HRV}.1D \
        -ort ${RVT}.1D \
        -ort ${RETROICOR}.1D \
        -bandpass 0.01 0.15 \
        -prefix ${func}_pd_bp.nii.gz
        # 5.
        # Stampacchia denoised
        3dTproject -polort ${polort} -input ${func_in}.nii.gz  -mask ${mref}_brain_mask.nii.gz \
        -ort ${func}_nuisreg_uncensored_mat.1D \
        -ort ${func}_avg_CSF.1D \
        -ort ${func}_avg_WM.1D \
        -bandpass 0.01 0.15 \
        -prefix ${func}_Sara_denoised.nii.gz
    done
done
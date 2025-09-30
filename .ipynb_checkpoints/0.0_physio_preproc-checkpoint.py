##!/usr/bin/env python
# Author: Laura Belli
#physiological preprocessing
# a complete script to run breathing and cardiac automatic peaks detection, 
# and to compute associated metrics (HRV, RVT and RV) as well as RETROICOR regressors
# Based on peakdet and phys2denoise packages developed by the Physiopy community https://physiopy.github.io/ 
# most of the code structure comes from a tutorial written by Stefano Moia https://github.com/smoia


#0. import libraries and functions
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import peakdet as pk
from phys2denoise.metrics import (
    heart_rate_variability,
    respiratory_variance,
    respiratory_variance_time,
)
from phys2denoise.metrics.utils import export_metric
import pathlib


#1. breathing_edit
def breathing_edit(i_sub, i_session, task_name):

    print(f"Subject {i_sub}, session {i_session} for task {task_name}")
    #select channels
    belt = 1 if i_session <= 6 else 2
    col_co2 = 4 if i_session <= 6 else 5
    #import data
    data = np.genfromtxt(f"./{i_sub:02d}/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_physio.tsv.gz", usecols=[0, belt, col_co2])
    idx_0 = np.argmax(data[:,0] >= 0)
    #create object physio
    phys = pk.Physio(data[idx_0:, 1], fs=10000, suppdata=data[idx_0:, 2])
    #filtering
    ph = pk.operations.filter_physio(phys, 2, method="lowpass", order=3)
    #downsample the signal
    ph = pk.operations.interpolate_physio(ph, 40.0, kind="linear")
    thr = 0.8
    dist = 100  # 1.5 s * 40 Hz 
    ph = pk.operations.peakfind_physio(ph, thresh=thr, dist=dist)
    #manually edit
    ph = pk.edit_physio(ph)
    phys.history
    #save breathing edited
    pk.save_physio(f"./checked/{i_sub:02d}/resp/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_resp_peaks.phys", phys)

    return (phys, i_sub, i_session, task_name)


#2. cardiac_edit
def cardio_edit(i_sub, i_session, task_name):
    print(f"Subject {i_sub:02d}, session {i_session:02d} for task {task_name}")
    #import data
    data = np.genfromtxt(f"./{i_sub:02d}/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_physio.tsv.gz")
    idx_0 = np.argmax(data[:, 0] >= 0)
    #select channels
    pulsometer = 3 if i_session <= 6 else 4
    phys = pk.Physio(data[idx_0:, pulsometer], fs=10000)
    #filtering
    ph = pk.operations.filter_physio(phys, 2, method="lowpass", order=7)
    #downsample the signal
    ph = pk.operations.interpolate_physio(phys, 40.0, kind="linear")
    thr = 0.3
    dist = 30 # 1s * 40 Hz
    ph = pk.operations.peakfind_physio(ph, thresh=thr, dist=dist)
    #manually edit
    ph = pk.edit_physio(ph)
    ph.history
    #save breathing edited
    pk.save_physio(f"./checked/{i_sub:02d}/cardio/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_card_peaks.phys", ph)

    return (ph, i_sub, i_session, task_name)


#3. physio metrics
def physio_model(i_sub, i_session, task_name):
    card = pk.load_physio(f"./checked/{i_sub:02d}/cardio/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_card_peaks.phys", allow_pickle=True)
    card.history
    #compute HRV 
    HRV = heart_rate_variability(card.data, card.peaks, card.fs)
    # export HRV
    export_metric(
    HRV, card.fs, tr=1.5, fileprefix=f"./checked/{i_sub:02d}/metrics/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_HRV", ntp=400)
    # Load the saved respiratory physiological file → allow pickle!!!
    resp = pk.load_physio(f"./checked/{i_sub:02d}/resp/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_resp_peaks.phys", allow_pickle=True)
    #check file
    resp.history
    # Compute RVT
    RVT = respiratory_variance_time(resp.data, resp.peaks, resp.troughs, resp.fs, lags=(0, 4, 8, 12))
    # export RVT
    export_metric(RVT, resp.fs, tr=1.5, fileprefix=f"./checked/{i_sub:02d}/metrics/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_RVT",
    ntp=400, is_convolved=False, has_lags=True)
    # Compute RV
    RV = respiratory_variance(resp.data, resp.fs)
    # export RV
    export_metric(
    RV, resp.fs, tr=1.5, fileprefix=f"./checked/{i_sub:02d}/metrics/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_RV", ntp=400)


#loop with breathing, cardiac cycles and physiological modelling
for i_sub in range(1, 11):
    for i_session in range(1, 11):
        for task_name in ["motor", "pinel", "rest_run-01", "rest_run-02", "rest_run-03", "rest_run-04", "simon"]:
            if not pathlib.Path(f"./checked/{i_sub:02d}/resp/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_resp_peaks.phys").exists():
                breathing_edit(i_session=i_session, task_name=task_name)
            
            if not pathlib.Path(f"./checked/{i_sub:02d}/cardio/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_card_peaks.phys").exists():
                cardio_edit(i_session=i_session, task_name=task_name)
            
            if not pathlib.Path(f"./checked/{i_sub:02d}/metrics/sub-{i_sub:02d}_ses-{i_session:02d}_task-{task_name}_HRV.phys").exists():
                physio_model(i_session=i_session, task_name=task_name)

# RETROICOR regressors for denoising
import peakdet as pk
for sub in range(1, 11):
    for ses in range(1, 11):
        data_resp = pk.load_physio(f"./breathing/sub-{sub:03d}/sub-{sub:03d}_ses-{ses:02d}_task-rest_run-01_resp_peaks.phys", allow_pickle=True)
        data_card = pk.load_physio(f"./cardiac/sub-{sub:03d}/sub-{sub:03d}_ses-{ses:02d}_task-rest_run-01_physio_peaks_ch-ppg.phys", allow_pickle=True)
        retroicor_regressors = retroicor(slice_timings= [0], 
                  resp_data=data_resp, 
                  card_data=data_card, 
                  resp_fs=40.0, 
                  card_fs=40.0,
                  n_harmonics=2,
                  base_offset=0,
                  nbins="p2d",
                  compute_interaction=True)
        results = retroicor_regressors.items()
        data = list(results)
        as_array = np.array(data)
        print(type(as_array), as_array.shape)
        print(as_array)
        np.save(f"./sub-{sub:03d}_ses-{ses:02d}_retroicor_reg", as_array)
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:59 2023

@author: Nate Whitsett
"""

import aflare as af
import Flare
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os

def amp_log_normal():
    value = 1000
    while value > 1 or value < 0.01:
        value = np.random.lognormal(0, sigma=2)
        
    return value

def FWHM_uniform():
    return np.random.uniform(0.001388888, 0.041)



# flare_baseline = dict()
# for T in range(2500, 6000):
#     flare_baseline[T] = pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000)

TESS_Folder_ID = [x[1] for x in os.walk('/data/whitsett.n/TESS_Light_Curves/All_TOI')]
TOI_Catalog = pd.read_csv('/data/whitsett.n/Reference_Files/All_TOI_12_17_23.csv')


TOI_ID_list = []
test_flare_list_T1 = []
test_flare_list_T2 = []
lc_num = 0
item_count = 0
random_sample = []
check = 0
    

for stars in range(200):
    sample = np.random.randint(0, len(TESS_Folder_ID[0]))
    while TESS_Folder_ID[0][sample] in random_sample:
        sample = np.random.randint(0, len(TESS_Folder_ID[0]))
    random_sample.append(TESS_Folder_ID[0][sample])
    
for M_dwarves in random_sample:
    ##Iteration Scheme
    a = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves)
    print(item_count, M_dwarves)
    
    ##Relevant parameters from TOI catalog
    flare_count = 1
    file = 0
    ##Trackable values per star
    for folders in a:
        if file == 0:
            inject_location_index = []

            flare_inject_dict_T1 = dict()
            flare_inject_dict_T2 = dict()

            b, pdcsap_flux, pdcsap_error = Flare.TESS_data_extract('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves + '/' + folders, PDCSAP_ERR=True)
            time, flux, pdcsap_error = Flare.delete_nans(b, pdcsap_flux, pdcsap_error)
            for flares in range(15):
                location = np.random.randint(50, len(flux)-50)
                inject_location_index.append(location)
                sample_baseline = flux[location-50:location+50]
                normalized_sample = sample_baseline/np.median(sample_baseline)
                FWHM = FWHM_uniform()
                amp = amp_log_normal()
                flare_inject = af.aflare1(time[location-50:location+50], time[location], FWHM, amp)
                normalized_sample_inject = normalized_sample + flare_inject
                flux[location-50:location+50] = np.median(sample_baseline)*normalized_sample_inject
                baseline_error = np.std(normalized_sample)
                flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, False]
                flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, False]

            detrend_flux = Flare.SMA_detrend(time, flux, 80, LS_Iterations=5)
            flares, lengths = Flare.flare_ID(detrend_flux, 3)
            index = 0
            for flare_events in flares:
                if flare_events >= 100 and len(flux) - flare_events > 100:
                    new_time = time[flare_events-100:flare_events+100]
                    new_data = flux[flare_events-100:flare_events+100]
                elif flare_events < 100:
                    new_time = time[0+flare_events:flare_events+100]
                    new_data = flux[0+flare_events:flare_events+100]
                elif len(flux) - flare_events < 100:
                    new_time = time[flare_events:]
                    new_data = flux[flare_events:]
                new_data = new_data/np.median(new_data)
                new_error = pdcsap_error[flare_events-100:flare_events+100]
                recenter = np.max(new_data[int(len(new_data)/2-15):int(len(new_data)/2+15)])
                norm_time = time[flare_events]
                events = np.where(new_data == recenter)[0][0]
                peak_centered_event = flare_events + (events-100)
                if peak_centered_event in inject_location_index:
                    flare_inject_dict_T1[peak_centered_event][3] = True
                criteria1 = False
                if recenter > np.mean(new_data)+2*(np.std(new_data)):
                    criteria1 = True
                if criteria1 == True and new_data[events+1] > np.mean(new_data)+(np.std(new_data)):
                    new_time = (new_time - new_time[events])*24*60
                    if lengths[index] >= 25:
                        # new_time = np.array(new_time[events-10:events+50])*24*60
                        # new_data = np.array(new_data[events-10:events+50])
                        alles_data = new_data/np.median(new_data)
                        BJD_time = time[events-10:events+50]
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
                        squares = (alles_data[events:events+30] - Flare.exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
                        chi2_cutoff = 20.514
                    elif lengths[index] >= 15 and lengths[index] < 25:
                        # new_time = np.array(new_time[events-10:events+30])*24*60
                        # new_data = np.array(new_data[events-10:events+30])
                        alles_data = new_data/np.median(new_data)
                        BJD_time = time[events-10:events+30]
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
                        squares = (alles_data[events:events+20] - Flare.exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
                        chi2_cutoff = 12.554
                    elif lengths[index] > 5 and lengths[index] < 15:
                        # new_time = np.array(new_time[events-10:events+20])*24*60
                        # new_data = np.array(new_data[events-10:events+20])
                        alles_data = new_data/np.median(new_data)
                        BJD_time = time[events-10:events+20]
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
                        squares = (alles_data[events:events+10] - Flare.exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
                        chi2_cutoff = 3.55
                    elif lengths[index] <= 5:
                        # new_time = np.array(new_time[events-10:events+14])*24*60
                        # new_data = np.array(new_data[(events-10):(events+14)])
                        alles_data = new_data/np.median(new_data)
                        BJD_time = time[events-10:events+14]
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
                        squares = (alles_data[events:events+7] - Flare.exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
                        chi2_cutoff = 1.455
                    chi_squared = np.sum(squares)
                    if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
                        if peak_centered_event in inject_location_index:
                            flare_inject_dict_T2[peak_centered_event][3] = True
                        index += 1

            test_flare_list_T1.append(flare_inject_dict_T1)
            test_flare_list_T2.append(flare_inject_dict_T2)
            lc_num += 1
            print('LC #: ' + str(lc_num))
            file += 1
        else:
            continue
        
amp_list_T1=[]
FWHM_list_T1 = []
result_T1 = []
error_T1 = []

amp_list_T2=[]
FWHM_list_T2 = []
result_T2 = []
error_T2 = []
for d in test_flare_list_T1:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T1.append(value[0])
        FWHM_list_T1.append(value[1])
        error_T1.append(value[2])
        result_T1.append(value[3])

for d in test_flare_list_T2:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T2.append(value[0])
        FWHM_list_T2.append(value[1])
        error_T2.append(value[2])
        result_T2.append(value[3])
        
ZZ = np.stack((amp_list_T1, FWHM_list_T1, error_T1, result_T1))
ZZ = np.transpose(ZZ)
np.savetxt('/data/whitsett.n/Pipeline_Validation/Injection_Recovery_All_TOI_T1.csv', ZZ, delimiter=',')

Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2, result_T2))
Z = np.transpose(Z)
np.savetxt('/data/whitsett.n/Pipeline_Validation/Injection_Recovery_All_TOI_T2.csv', Z, delimiter=',')
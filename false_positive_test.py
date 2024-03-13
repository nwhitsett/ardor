# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 21:47:16 2023

@author: Nate Whitsett
"""

import aflare as af
import Flare
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os
import time as t

def amp_log_normal():
    value = 1000
    while value > 0.8 or value < 0.2:
        value = np.random.lognormal(0, sigma=2)
        
    return value

def FWHM_uniform():
    return np.random.uniform(0.002388888, 0.041)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# flare_baseline = dict()
# for T in range(2500, 6000):
#     flare_baseline[T] = pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000)

TESS_Folder_ID = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data')

TOI_ID_list = []
test_flare_list = []
lc_num = 0
item_count = 0
random_sample = []
check = 0
test_time = 0
files = 0

total_flares_detected_T1 = 0
total_flares_detected_T2 = 0
inject_flare_count = 0

inject_flare_recovered_T1 = 0
inject_flare_recovered_T2 = 0

false_positive_T1 = 0
false_positive_T2 = 0

manual_inspect_flare = 0

tier_1_time = 0
tier_2_time = 0
for G_Stars in TESS_Folder_ID:
    ##Iteration Scheme
    a = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data/' + G_Stars)
    
    ##Trackable values per star
    location_list = []
    for folders in a:
        inject_location_index = []
        flare_inject_dict = dict()
        b, pdcsap_flux, pdcsap_error = Flare.TESS_data_extract('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data/' + G_Stars + '/' + folders, PDCSAP_ERR=True)
        time, flux = Flare.delete_nans(b, pdcsap_flux)
        if folders.endswith('a_fast-lc.fits') == True:
            test_time += len(time)*(0.33333333)
            print(0)
        elif folders.endswith('a_fast-lc.fits') == False:  
            test_time += len(time)*2
            print(1)
        
        for flares in range(10):
            location = np.random.randint(50, len(flux)-50)
            for locations in location_list:
                if locations + 300 > location or locations - 300 < location:
                    location = np.random.randint(50, len(flux)-50)
            inject_location_index.append(location)
            sample_baseline = flux[location-300:location+300]
            normalized_sample = sample_baseline/np.median(sample_baseline)
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            flux[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            baseline_error = np.std(normalized_sample)
            flare_inject_dict[location] = [amp, FWHM, baseline_error, False]
            inject_flare_count += 1
            location_list.append(location)
        t0 = t.time()
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
            recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
            c, d = Flare.flare_ID(np.array(new_data), 3)
            norm_time = time[flare_events]
            events = np.where(new_data == recenter)[0][0]
            peak_centered_event = flare_events + (events-100)
            criteria1 = False
            
            
            flare_location = find_nearest(inject_location_index, flare_events)
            if np.abs(flare_location - flare_events) < 50:
                inject_flare_recovered_T1 += 1
            elif np.abs(flare_location - flare_events) > 50:
                false_positive_T1 += 1
            total_flares_detected_T1 += 1
                    
                    
            if recenter > np.mean(new_data)+3*(np.std(new_data)):
                criteria1 = True
            t1 = t.time()
            tier_1_time += t1-t0
            if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)) and len(c) > 0:
                t2 = t.time()
                new_time = (new_time - new_time[events])*24*60
                if lengths[index] >= 25:
                    # new_time = np.array(new_time[events-10:events+50])*24*60
                    # new_data = np.array(new_data[events-10:events+50])
                    alles_data = new_data/np.median(new_data)
                    BJD_time = time[events-10:events+50]
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
                    squares = (alles_data[events:events+30] - Flare.exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
                    chi2_cutoff = 18
                elif lengths[index] >= 15 and lengths[index] < 25:
                    # new_time = np.array(new_time[events-10:events+30])*24*60
                    # new_data = np.array(new_data[events-10:events+30])
                    alles_data = new_data/np.median(new_data)
                    BJD_time = time[events-10:events+30]
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
                    squares = (alles_data[events:events+20] - Flare.exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
                    chi2_cutoff = 9.5
                elif lengths[index] > 5 and lengths[index] < 15:
                    # new_time = np.array(new_time[events-10:events+20])*24*60
                    # new_data = np.array(new_data[events-10:events+20])
                    alles_data = new_data/np.median(new_data)
                    BJD_time = time[events-10:events+20]
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
                    squares = (alles_data[events:events+10] - Flare.exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
                    chi2_cutoff = 3
                elif lengths[index] <= 5:
                    # new_time = np.array(new_time[events-10:events+14])*24*60
                    # new_data = np.array(new_data[(events-10):(events+14)])
                    alles_data = new_data/np.median(new_data)
                    BJD_time = time[events-10:events+14]
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
                    squares = (alles_data[events:events+7] - Flare.exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
                    chi2_cutoff = 2
                chi_squared = np.sum(squares)
                t3 = t.time()
                tier_2_time += t3 - t2
                if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
                    flare_location = find_nearest(inject_location_index, flare_events)
                    if np.abs(flare_location - flare_events) < 50:
                        inject_flare_recovered_T2 += 1
                    elif np.abs(flare_location - flare_events) > 50:
                        plt.scatter(new_time, alles_data)
                        plt.show()
                        a = input('Continue?')
                        plt.clf()
                        if a == 'y':
                            manual_inspect_flare += 1
                        else:
                            false_positive_T2 += 1
                    total_flares_detected_T2 += 1
                    index += 1
        
        test_flare_list.append(flare_inject_dict)
        lc_num += 1
        print('LC #: ' + str(lc_num))
        files += 1

print('Result of test:')
print('Files Tested: ' + str(files))
print('Total Time Analyzed: ' + str(test_time) + 'minutes')

print('Tier 1 Detections: ' + str(total_flares_detected_T1))
print('Tier 2 Detections: ' + str(total_flares_detected_T2))

print('Tier 1 False Detections: '+ str(false_positive_T1))
print('Tier 2 False Detections: '+ str(false_positive_T2))

print('Tier 1 True Flares Recovered: ' + str(inject_flare_recovered_T1))
print('Tier 2 True Flares Recovered: ' + str(inject_flare_recovered_T2))

print('Tier 1 True Flares Missed: ' + str(inject_flare_count - inject_flare_recovered_T1))
print('Tier 2 True Flares Missed: ' + str(inject_flare_count - inject_flare_recovered_T2))

print('Tier 1 %-True Flares Recovered: ' + str(inject_flare_recovered_T1/inject_flare_count))
print('Tier 2 %-True Flares Recovered: ' + str(inject_flare_recovered_T2/inject_flare_count))

print('Tier 1 False Positive Rate: ' + str(false_positive_T1/test_time))
print('Tier 2 False Positive Rate: ' + str(false_positive_T2/test_time))

print('Tier 2 Manually Inspected, Non Injected Flares: ' + str(manual_inspect_flare))
amp_list=[]
FWHM_list = []
result = []
error = []
for d in test_flare_list:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list.append(value[0])
        FWHM_list.append(value[1])
        error.append(value[2])
        result.append(value[3])
        
Z = np.stack((amp_list, FWHM_list, error, result))
Z = np.transpose(Z)
np.savetxt('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/False_Positive Test2.csv', Z, delimiter=',')
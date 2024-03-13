# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:20:10 2023

@author: Nate Whitsett
"""
from astropy.io import fits
from astropy.timeseries import LombScargle as LS
# from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
import pandas as pd
import statistics as st
import os
import warnings
import planck_law as pl
import planck_law as pl
import aflare
warnings.filterwarnings("ignore")
def exp_decay(x, a, b, c):
    '''
    

    Parameters
    ----------
    x : numpy array
        Time axis data
    a : float
        multiplicative parameter of the exponential decay function
    b : float
        time constant parameter of the exponential decay function
    c : float
        Time-offset for the exponential decay function

    Returns
    -------
    numpy array
        Gives exponential decay model output with the given parameters

    '''
    return a * np.exp(-b * x) + c


def TESS_FITS_csv(input_file, csv_directory, csv_name=None):
    '''

    Parameters
    ----------
    input_file : string 
        Directory to the TESS light curve fits (...lc.fits) file
    csv_directory : string
        Directory for output file.
    csv_name : string, optional
        Name for output csv file

    Returns
    -------
    .cvs
        Outputs the SAP_FLUX, PDCSAP_FLUX, and time parameters from the .fits file 
        to an easily readable .csv file

    '''
    rev_file = input_file[::-1]
    index = rev_file.index('/')
    file = ((rev_file[:index])[::-1])[:-5]
    if csv_name == None:
        directory = csv_directory + '/' + file + '.csv'
    elif csv_name != None:
        directory = csv_directory + '/' + csv_name + '.csv'
    hdul = fits.open(input_file)
    time = hdul[1].data['TIME']
    sap_flux = hdul[1].data['SAP_FLUX']
    pdcsap_flux = hdul[1].data['PDCSAP_FLUX']


    
    grand_list = pd.DataFrame({'time': time, 'sap_flux': sap_flux, 'pdcsap_flux': pdcsap_flux})
    grand_list = pd.DataFrame({'time': time, 'sap_flux': sap_flux, 'pdcsap_flux': pdcsap_flux})
    return grand_list.to_csv(directory)

def TESS_data_extract(fits_lc_file, SAP_ERR=False, PDCSAP_ERR=False):
    '''

    Parameters
    ----------
    fits_lc_file : string
        Directory of the TESS light curve fits file
    SAP_ERR : bool, optional
        True will return SAP_FLUX error. The default is False.
    PDCSAP_ERR : bool, optional
        True will return PDCSAP_FLUX error. The default is False.

    Returns
    -------
    ND np.array
        Returns an Nd array of the time, PDCSAP_FLUX, SAP_FLUX, and/or the
        SAP_FLUX and PDCSAP_FLUX error. Min 3D array, max 5D array.

    '''
    hdul = fits.open(fits_lc_file)
    time = hdul[1].data['TIME']
    sap_flux = hdul[1].data['SAP_FLUX']
    sap_flux_error = hdul[1].data['SAP_FLUX_ERR']
    pdcsap_flux = hdul[1].data['PDCSAP_FLUX']
    pdcsap_flux_error = hdul[1].data['PDCSAP_FLUX_ERR']
    #pdcsap_flux_method = hdul[1].header['PDCMETHD']
    if SAP_ERR == False and PDCSAP_ERR == False:
        return time, pdcsap_flux
    if SAP_ERR == False and PDCSAP_ERR == True:
        return time, pdcsap_flux, pdcsap_flux_error
    if SAP_ERR == True and PDCSAP_ERR == False:
        return time, sap_flux, sap_flux_error, pdcsap_flux
    if SAP_ERR == True and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux, sap_flux_error, pdcsap_flux


def phase_folder(time, period, epoch):
    '''

    Parameters
    ----------
    time : numpy array
        The time array to be phase folded
    period : float
        The period of the planet, in the same units as the time array
    epoch : float
        The time of the first transit (i.e. transit epoch), same unit as time

    Returns
    -------
    phase : numpy array
        Numpy array of the same dimension as the input time array, but phase
        folded to the period of the planet. Transit centered at 0.


    '''
    phase = (time - (epoch+period/2)) % period
    phase = phase - period/2
    return phase

def flare_ID(data, sigma):
    '''
    

    Parameters
    ----------
    data : numpy array
        Flux data for potential flares to be identified
    sigma : float
        The detection sensitivity for flares, in standard deviations. 
        For example, a sigma value of 1.0 will count any data one standard
        deviation away from the mean as a potential flare

    Returns
    -------
    flare_indices : numpy array 
        Outputs a list of potential flares as lists of indices in the provided
        data. The index begins at the triggering data point.

    '''
    count = 0
    flare_indices = []
    flare_length = 0
    flare_length_list = []
    flare = False
    begin = 0
    end = 250
    shift = 0
    peak_index = 0
    for flux in data:
        if end >= len(data):
            sig = sigma*np.std(data[begin:len(data)-1])
            mu = np.mean(data[begin:len(data)-1])
        else:
            sig = sigma*np.std(data[begin:end])
            mu = np.mean(data[begin:end])
        if flux > (mu + sig) and flare == False:
            flare = True
            flare_length += 1
            peak_index = count
        try:
            if flare == True and data[count+1] > (mu + sig/3):
                flare_length += 1
            if flare == True and data[count+1] < (mu + sig/3) and flare_length >= 3:
                flare = False
                flare_indices.append(peak_index)
                flare_length_list.append(flare_length)
                flare_length = 0
                peak_index = 0
            count += 1
            shift += 1
            if shift == 10:
                begin += 10
                end += 10
                shift = 0
        except:
            continue
    return flare_indices, flare_length_list

def delete_nans(time, data):
    '''
    
    Parameters
    ----------
    time : numpy array
        The time array to be cleared of NANs. Must be the same dimensionality
        and 1-to-1 with the data array.
    data : numpy array
        The data array to be cleared of NANs. Must be the same dimensionality
        and 1-to-1 with the time array

    Returns
    -------
    time1 : numpy array
        Returns the original time array, but with any NANs in the data or
        time array removed for both arrays, at the same indices, such that
        both arrays are still 1-to-1.
    data1 : numpy array
        Returns the original data array, but with any NANs in the data or
        time array removed for both arrays, at the same indices, such that
        both arrays are still 1-to-1

    '''
    time = time.tolist()
    data = data.tolist()
    nan_set = set()
    count_data = 0
    count_time = 0
    for indices in data:
        if np.isnan(indices) == True:
            nan_set.add(count_data)
        count_data += 1
    for indices in time:
        if np.isnan(indices) == True:
            nan_set.add(count_time)
        count_time += 1
    time1 = [i for j, i in enumerate(time) if j not in nan_set]
    data1 = [i for j, i in enumerate(data) if j not in nan_set]
    time1 = np.array(time1)
    data1 = np.array(data1)
    return time1, data1

def SMA_detrend(time, data, time_scale, LS_Iterations=3):
    '''
    
    This applies a windowed, Single Moving Average to detrend potentially
    periodic data.
    Parameters
    ----------
    data : numpy array
        The flux data to be detrended
    time_scale : int
        The time-scale to apply the detrending.

    Returns
    -------
    numpy array
        Returns the detrended data array.

    '''
    mov_average = []
    j = 0
    i = 0
    for a in range(time_scale - 1):
        window = data[a : 1 + a + j]
        mov_average.append(sum(window)/(j+1))
        j += 1
    while i < len(data) - time_scale + 1:
        window = data[i : i + time_scale]
        window_average = round(sum(window) / time_scale, 2)
        mov_average.append(window_average)
        i += 1
    SMA = data - np.array(mov_average)
    count = 0
    ls = LS(time, SMA)
    freq, power = ls.autopower(minimum_frequency=1, maximum_frequency=1000, method='fast')
    cutoff = ls.false_alarm_probability(power.max())
    while cutoff < 0.99 and count <= LS_Iterations:
        best_frequency = freq[np.argmax(power)]
        ls = LS(time, SMA)
        theta = ls.model_parameters(best_frequency)
        offset = ls.offset()
        design_matrix = ls.design_matrix(best_frequency, time)
        SMA = SMA -( offset +design_matrix.dot(theta))
        cutoff = ls.false_alarm_probability(power.max())
        ls = LS(time, SMA)
        freq, power = ls.autopower(minimum_frequency=1, maximum_frequency=1000, method='fast')
        count += 1
    return SMA
        
def flare_phase_folded_ID(phase, flare_array, period, epoch):
    new_ID_list = []
    for indices in flare_array:
        new_ID_list.append(((phase[indices] - (epoch+period/2)) % period)-period/2)
    return np.array(new_ID_list)

def bolo_flare_energy(parameters, Teff, R_stellar, planck_ratio, t_unit='days', function=exp_decay):
    a, b, c = parameters
    if t_unit == 'days':
        multiplier = 86400
        length_cap = 0.08333
    if t_unit == 'minutes':
        multiplier = 60
        length_cap = 120
    if function == exp_decay:
        integral, err = quad(function, 0, length_cap, args=(a, b, (c-1)))
    elif function == aflare.aflare1:
        integral, err =quad(function, -length_cap, length_cap, args=(a, b, c))
    energy = (5.67e-8)*(9000**4)*(integral)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*planck_ratio*(1e7)*multiplier
    return energy

def analyze_flares(data_directory, flare_csv_directory, period, epoch, host_T, host_radius, host_ID):
    """
    

    Parameters
    ----------
    data_directory : string
        The directory of the .fits file you wish to analyze flares for
    flare_csv_directory : string
        The directory of the folder you wish to save the flare output .csv's
    period : float
        The period of the planet, in days.
    epoch : float
        The epoch (time of center transit) of the planet, in BJD
    host_T : float
        The host's stellar effective temperature, in K
    host_radius : float
        The hosts's radius, in solar radii
    host_ID : string
        The name of the host star

    Returns
    -------
    host_ID_list : list
        A list of the host ID associated with the flares
    accepted_flare_number : list
        The flare number in the processing
    peak_time : list
        The peak time of the flare, in BJD
    flare_amplitude : list
        The relative amplitude of the flare, in relative flux
    flare_time_scale : list
        The time constant of the flare associated with the exponential decay function fit
    flare_energy : list
        The approximate bolometric flare energy, in ergs
    Teff : list
        The stellar effective temperature of the host
    radius : list
        The radius of the host
    flare_phase : list
        The phase of the flare relative to the planet's transit
    flare_count : int
        The number of flares found in the file (Useful if iterating over many files)

    """
    flare_baseline = dict()
    for T in range(2500, 6000):
        flare_baseline[T] = pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000)
    observation_time = 0
    possible_flares = 0
    flare_count = 0
    total_flares = 0
    
    flare_time_scale = []
    flare_amplitude = []
    peak_time = []
    flare_number = []
    host_ID_list = []
    flare_energy = []
    Teff = []
    radius = []
    flare_phase = []
    accepted_flare_index = []
    accepted_flare_number = []
    
    
    
    file = data_directory
    b, pdcsap_flux, pdcsap_error = TESS_data_extract(file, PDCSAP_ERR=True)
    time, flux = delete_nans(b, pdcsap_flux)
    detrend_flux = SMA_detrend(time, flux, 80, LS_Iterations=5)
    flares, lengths = flare_ID(detrend_flux, 3)
    if file.endswith('a_fast-lc.fits') == True:
        observation_time += len(time)*(0.33333333)
    elif file.endswith('a_fast-lc.fits') == False:  
        observation_time += len(time)*2
    phase = phase_folder(time, period, epoch)
    index = 0
    possible_flares += len(flares)
    for flare_events in flares:
        print(1)
        if flare_events >= 100 and len(flux) - flare_events > 100:
            new_time = time[flare_events-100:flare_events+100]
            new_data = flux[flare_events-100:flare_events+100]
        elif flare_events < 100:
            new_time = time[0+flare_events:flare_events+100]
            new_data = flux[0+flare_events:flare_events+100]
        elif len(flux) - flare_events < 100:
            new_time = time[flare_events:]
            new_data = flux[flare_events:]
        if period != '':
            phase = phase_folder(new_time, period, epoch)
        elif period == '':
            phase = ''
        new_error = pdcsap_error[flare_events-100:flare_events+100]
        recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
        c, d =flare_ID(np.array(new_data), 3)
        norm_time = time[flare_events]
        events = np.where(new_data == recenter)[0][0]
        flare_phase_value = phase[events]
        criteria1 = False
        if recenter > np.mean(new_data)+3*(np.std(new_data)):
            criteria1 = True
        if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)) and len(c) > 0:
            new_time = (new_time - new_time[events])*24*60
            if lengths[index] >= 25:
                alles_data = new_data/np.median(new_data)
                error = new_error/np.median(new_data)
                popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
                squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
                chi2_cutoff = 18
            elif lengths[index] >= 15 and lengths[index] < 25:
                alles_data = new_data/np.median(new_data)
                error = new_error/np.median(new_data)
                popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
                squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
                chi2_cutoff = 9.5
            elif lengths[index] > 5 and lengths[index] < 15:
                alles_data = new_data/np.median(new_data)
                error = new_error/np.median(new_data)
                popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
                squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
                chi2_cutoff = 2.167
            elif lengths[index] <= 5:
                alles_data = new_data/np.median(new_data)
                error = new_error/np.median(new_data)
                popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
                squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
                chi2_cutoff = 1.2
            chi_squared = np.sum(squares)
            if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
                flare_time_scale.append(popt[1])
                flare_amplitude.append(popt[0])
                peak_time.append(norm_time)               
                host_ID_list.append(host_ID)
                flare_number.append(flare_count)
                try:
                    energy = bolo_flare_energy(popt, host_T, host_radius, flare_baseline[host_T], t_unit='days')
                except:
                    energy = np.NaN
                flare_energy.append(energy)
                Teff.append(host_T)
                radius.append(host_radius)
                X = np.column_stack((new_time[events-30:events+40], alles_data[events-30:events+40], error[events-30:events+40]))
                flare_phase.append(flare_phase_value)
                accepted_flare_index.append(flares[index])
                accepted_flare_number.append(flare_count)
                np.savetxt(flare_csv_directory + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
                flare_count += 1
                total_flares += 1
        index += 1
    return host_ID_list, accepted_flare_number, peak_time, flare_amplitude, flare_time_scale, flare_energy, Teff, radius, flare_phase, flare_count

# TESS_Folder_ID = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/')
# TOI_Catalog = pd.read_csv('/data/whitsett.n/Reference_Files/All_TOI_12_17_23.csv')
# total_flares = 0
# total_possible_flares = 0
# total_observation_time = 0
# TOI_ID_list = []
# flare_number = []
# peak_time = []
# amplitude = []
# time_scale = []
# Teff = []
# radius = []
# flare_phase = []
# TOI_period = []
# total_flare_energies = []
# total_flare_phases = []
# item_count = 0
# total_periastron_list = []
# total_periastron_epoch_list = []
# total_epoch_list = []
# list_index = 0
# for M_dwarves in TESS_Folder_ID[1270:]:
#     ##Iteration Scheme
#     TOI_ID = M_dwarves
#     print(TOI_ID)

#     a = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves)
#     os.mkdir('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/')
#     print(item_count, M_dwarves)
#     ##Relevant parameters from TOI catalog
#     period = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_orbper'])[0]
#     epoch = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_tranmid'])[0]
#     stellar_radius = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_rad'])[0]
#     T = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_teff'])[0]
#     # periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbper'])[0]
#     # epoch_periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbtper'])[0]
#     if T == '':
#         T = np.NAN
#     if stellar_radius == '':
#         stellar_radius = np.NAN
#     # if epoch_periastron == '':
#     #     epoch_periastron = np.NAN
#     # if periastron == '':
#     #     periastron = np.NAN
#     if epoch == '':
#         epoch = np.NAN
#     flare_count = 1
    
#     ##Trackable values per star
#     phase_folded_flare_list = []
#     flare_amplitude = []
#     flare_time_scale = []
#     flare_energy = []
#     accepted_flare_index = []
#     flare_phase = []
#     accepted_flare_number = []
#     observation_time = 0
#     possible_flares = 0
#     # periastron_list = []
#     # periastron_epoch_list = []
#     epoch_list = []
    
#     ##Total Trackable Lists for all Data
#     TOI_ID_list = []
#     flare_number = []
#     peak_time = []
#     amplitude = []
#     time_scale = []
#     Teff = []
#     radius = []
#     flare_phase = []
#     TOI_period = []
#     total_flare_energies = []
#     total_flare_phases = []
#     item_count = 0
#     # total_periastron_list = []
#     # total_periastron_epoch_list = []
#     total_epoch_list = []
#     for folders in a:
#         try:
#             b, pdcsap_flux, pdcsap_error = TESS_data_extract('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves + '/' + folders, PDCSAP_ERR=True)
#         except:
#             continue
#         time, flux = delete_nans(b, pdcsap_flux)
#         detrend_flux = SMA_detrend(time, flux, 80, LS_Iterations=5)
#         flares, lengths = flare_ID(detrend_flux, 3)
#         if folders.endswith('a_fast-lc.fits') == True:
#             observation_time += len(time)*(0.33333333)
#             total_observation_time += len(time)*(0.33333333)
#         elif folders.endswith('a_fast-lc.fits') == False:  
#             observation_time += len(time)*2
#             total_observation_time += len(time)*2
#         index = 0
#         possible_flares += len(flares)
#         total_possible_flares += len(flares)

#         for flare_events in flares:
#             if flare_events >= 100 and len(flux) - flare_events > 100:
#                 new_time = time[flare_events-100:flare_events+100]
#                 new_data = flux[flare_events-100:flare_events+100]
#             elif flare_events < 100:
#                 new_time = time[0+flare_events:flare_events+100]
#                 new_data = flux[0+flare_events:flare_events+100]
#             elif len(flux) - flare_events < 100:
#                 new_time = time[flare_events:]
#                 new_data = flux[flare_events:]
#             new_error = pdcsap_error[flare_events-100:flare_events+100]
#             recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
#             c, d =flare_ID(np.array(new_data), 3)
#             norm_time = time[flare_events]
#             events = np.where(new_data == recenter)[0][0]
#             if period != '':
#                 phase = phase_folder(new_time, period, epoch)
#                 flare_phase_value = phase[events]
#             elif period == '':
#                 phase = []
#                 flare_phase_value = np.NAN
#                 period = np.NAN
                
#             criteria1 = False
#             if recenter > np.mean(new_data)+3*(np.std(new_data)):
#                 criteria1 = True
#             if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)) and len(c) > 0:
#                 new_time = (new_time - new_time[events])*24*60
#                 if lengths[index] >= 25:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
#                     squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
#                     chi2_cutoff = 20.843
#                 elif lengths[index] >= 15 and lengths[index] < 25:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
#                     squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
#                     chi2_cutoff = 11.912
#                 elif lengths[index] > 5 and lengths[index] < 15:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
#                     squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
#                     chi2_cutoff = 3.455
#                 elif lengths[index] <= 5:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
#                     squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
#                     chi2_cutoff = 1.3
#                 chi_squared = np.sum(squares)
#                 if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
#                     half_max = (alles_data[8:15].max()-np.median(alles_data[0:8]))/2
#                     time_scale.append(popt[1])
#                     flare_time_scale.append(popt[1])
#                     amplitude.append(popt[0])
#                     flare_amplitude.append(popt[0])
#                     peak_time.append(norm_time)
                    
#                     TOI_ID_list.append(TOI_ID)
#                     flare_number.append(flare_count)
#                     try:
#                         energy = bolo_flare_energy(popt, T, stellar_radius, pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000), t_unit='days')
#                     except:
#                         energy = np.NaN
#                     flare_energy.append(energy)
#                     total_flare_energies.append(energy)
#                     Teff.append(T)
#                     radius.append(stellar_radius)
#                     try:
#                         X = np.column_stack((new_time[events-30:events+40], alles_data[events-30:events+40], error[events-30:events+40]))
#                     except:
#                         X = np.column_stack(([0],[0],[0]))
#                     baseline = st.median(new_data)*(lengths[index])*2
#                     median = st.median(new_data)
#                     flare_phase.append(flare_phase_value)
#                     accepted_flare_index.append(flares[index])
#                     accepted_flare_number.append(flare_count)
#                     total_flare_phases.append(flare_phase_value)
#                     TOI_period.append(period)
#                     # total_periastron_list.append(periastron)
#                     # total_periastron_epoch_list.append(epoch_periastron)
#                     # periastron_list.append(periastron)
#                     # periastron_epoch_list.append(epoch_periastron)
#                     total_epoch_list.append(epoch)
#                     epoch_list.append(epoch)
#                     print(flare_count)
#                     if lengths[index] > 5:
#                         print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
#                     np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
#                     flare_count += 1
#                     total_flares += 1
                    
                        
#             index += 1
#     if observation_time == 0:
#         observation_time = -1
#     if possible_flares == 0 :
#         possible_flares = -1
#     Y = np.column_stack((flare_amplitude, flare_time_scale, flare_energy))
#     np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/All_Flares.csv', Y, delimiter=',')
#     f = open('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Host_Statistics.txt', 'w')
#     f.write('Flares/Day: ' + str(flare_count*60*24/observation_time) + '\n' + 'Possible Flares: ' + str(possible_flares) + '\n' +
#             'Accepted Flares: ' + str(flare_count) + '\n' + 'Accepted/Possible Flares: ' + str(flare_count/possible_flares) + '\n' +
#             'Observation Time (min): ' + str(observation_time))
#     f.close()
#     item_count += 1
#     ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, total_flare_energies, Teff, radius, TOI_period))
#     with open("/data/whitsett.n/Tier_3/All_TOI/All_TOI_Flares2.csv", "a") as f:
#         np.savetxt(f, ZZ, delimiter=",", fmt='%s')
#         f.close()
#     print(len(TOI_ID_list), len(flare_number), len(peak_time), len(amplitude), len(time_scale), len(total_flare_energies), len(Teff), len(radius), len(TOI_period))
#     print(list_index)
# g = open('/data/whitsett.n/Tier_3/All_TOI/All_statistic_TOI.txt', 'w')
# g.write('Total Flares: ' + str(total_flares) + '\n' + 'Net Observation Time (Days): ' + str(total_observation_time/(60*24)))
# g.close()
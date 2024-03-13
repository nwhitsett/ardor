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
import aflare
import time as timer
import copy
import allesfitter_priors
import shutil
import allesfitter
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

def flare_ID(data, sigma, fast = False):
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
    if fast == False:
        for flux in data:
            try:
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
                if flare == True and data[count+1] > (mu + sig/3):
                    flare_length += 1
                if flare == True and data[count+1] < (mu + sig/3) and flare_length >= 3:
                    flare = False
                    flare_indices.append(peak_index)
                    flare_length_list.append(flare_length)
                    flare_length = 0
                    peak_index = 0
                elif flare == True and data[count+1] < (mu + sig/3) and flare_length < 3:
                    flare = False
                    flare_length = 0
                count += 1
                shift += 1
                if shift == 250:
                    begin += 250
                    end += 250
                    shift = 0
            except:
                continue
    elif fast == True:
        for flux in data:
            try:
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
                if flare == True and data[count+1] > (mu + sig/3):
                    flare_length += 1
                if flare == True and data[count+1] < (mu + sig/3) and flare_length >= 3:
                    flare = False
                    flare_indices.append(peak_index)
                    flare_length_list.append(flare_length)
                    flare_length = 0
                    peak_index = 0
                elif flare == True and data[count+1] < (mu + sig/3) and flare_length < 3:
                    flare = False
                    flare_length = 0
                count += 1
                shift += 1
                if shift == 1000:
                    begin += 1000
                    end += 1000
                    shift = 0
            except:
                continue
    return flare_indices, flare_length_list

def delete_nans(time, data, error):
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
    error = error.tolist()
    nan_set = set()
    count_data = 0
    count_time = 0
    count_error = 0
    for indices in data:
        if np.isnan(indices) == True:
            nan_set.add(count_data)
        count_data += 1
    for indices in time:
        if np.isnan(indices) == True:
            nan_set.add(count_time)
        count_time += 1
    for indices in error:
        if np.isnan(indices) == True:
            nan_set.add(count_error)
        count_time += 1
    time1 = [i for j, i in enumerate(time) if j not in nan_set]
    data1 = [i for j, i in enumerate(data) if j not in nan_set]
    error1 = [i for j, i in enumerate(error) if j not in nan_set]
    time1 = np.array(time1)
    data1 = np.array(data1)
    error1 = np.array(error1)
    return time1, data1, error1

def SMA_detrend(time, data, time_scale, LS_Iterations=25):
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
    ls = LS(time,data)
    frequency,power = ls.autopower(minimum_frequency=0.1, maximum_frequency=100)
    cutoff = ls.false_alarm_probability(power.max()) 
    if cutoff < 0.10:
        ls = LS(time,data, nterms=3)
        frequency, power = ls.autopower(minimum_frequency=0.1, maximum_frequency=100)
        best_frequency = frequency[np.argmax(power)]
        theta = ls.model_parameters(best_frequency)
        offset = ls.offset()
        design_matrix = ls.design_matrix(best_frequency, time)
        data = data - (offset + design_matrix.dot(theta))
        LS_model = (offset + design_matrix.dot(theta))
    tmp_data = copy.deepcopy(data)
    sigma = np.std(tmp_data)
    mean = np.median(tmp_data)
    for index, num in enumerate(tmp_data):
        if num >= 3*sigma + mean:
            tmp_data[index] = mean
    mov_average = []
    j = 0
    i = 0
    for a in range(time_scale - 1):
        window = tmp_data[a : 1 + a + j]
        mov_average.append(sum(window)/(j+1))
        j += 1
    while i < len(data) - time_scale + 1:
        window = tmp_data[i : i + time_scale]
        window_average = round(sum(window) / time_scale, 2)
        mov_average.append(window_average)
        i += 1
    SMA = data - np.array(mov_average)
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
    time, flux, pdcsap_error = delete_nans(b, pdcsap_flux, pdcsap_error)
    detrend_flux = SMA_detrend(time, flux, 80, LS_Iterations=5)
    if file.endswith('a_fast-lc.fits') == True:
        observation_time += len(time)*(0.33333333)
        cadence = 'sec'
    elif file.endswith('a_fast-lc.fits') == False:  
        observation_time += len(time)*2
        cadence = 'min'
    flares, lengths = flare_ID(detrend_flux, 3, cadence)
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

def tier0(TESS_fits_file):
    '''
    

    Parameters
    ----------
    TESS_fits_file : string
        The TESS light curve you wish to detrend and clean up.
    Returns
    -------
    time : numpy array
        The time axis, given in BJD - 2457000 (Days). Used next in tier 2
    flux : numpy array
        The raw, cleaned up pdcsap flux from the TESS file. To be used in tier 2.
    detrend_flux : numpy array
        Median centered pdcsap flux, given in electrons/second. Used next in tier 1
    pdcsap_error : numpy array
        Error in the detrended flux. Used next in tier 2

    '''
    b, pdcsap_flux, pdcsap_error = TESS_data_extract(TESS_fits_file, PDCSAP_ERR=True)
    time, flux, pdcsap_error = delete_nans(b, pdcsap_flux, pdcsap_error)
    detrend_flux = SMA_detrend(time, flux, 80, LS_Iterations=5)
    if TESS_fits_file.endswith('a_fast-lc.fits') == True:
        fast = True
    elif TESS_fits_file.endswith('a_fast-lc.fits') == False:  
        fast = False
    return time, flux, detrend_flux, pdcsap_error, fast

def tier1(detrend_flux, sigma, fast=False):
    '''
    

    Parameters
    ----------
    detrend_flux : numpy array
        The median centered, detrended pdcsap flux data of the TESS light curve file, in units of electrons/second. 
        Intended to be from the second output of the tier 0 function
    cadence : string
        Either 'min' or 'sec'. Depends on what TESS fits file is used: a_fast-lc = 'sec', else 'min'
    sigma : float
        The sensitivity cutoff in which you wish to search for flares. Typically, this is 3 sigma, though some pipelines
        go as low as 2.5. Will NOT work well below 1 sigma.
    Returns
    -------
    flares : numpy array
        Returns an array of the **indices** of the flares in the time array passed. To get the BJD time of a given flare, 
        pass 'time[flares[0]]', etc. Used in tier 2.
    lengths : numpy array
        A 1:1 numpy array with flares, giving the approximate duration of the flare in units of **indices** of the time 
        axis. This will either be two minutes/index or 20 seconds an index, depending on what type of light curve is used.
        Used in tier 2.

    '''
    flares, lengths = flare_ID(detrend_flux, sigma, fast=fast)
    return flares, lengths

def tier2(time, flux, pdcsap_error, flares, lengths, output_dir, host_name = 'My_Host', T = 4000, host_radius = 1):
    '''
    
    Parameters
    ----------
    time : numpy array
        The time data of the TESS light curve file. Intended to be from the first output of the tier 0 function.
        Units of BJD - 255700 (Days)
    detrend_flux : numpy array
        The median centered, detrended pdcsap flux data of the TESS light curve file, in units of electrons/second. 
        Intended to be from the second output of the tier 0 function
    pdcsap_error : numpy array
        The uncertainty in the detrended flux. Intended to be from the third output of the tier 0 function.
    flares : numpy array
        Returns an array of the **indices** of the flares in the time array passed. To get the BJD time of a given flare, 
        pass 'time[flares[0]]', etc. Intended as the first output from tier 2.
    lengths : numpy array
        A 1:1 numpy array with flares, giving the approximate duration of the flare in units of **indices** of the time 
        axis. This will either be two minutes/index or 20 seconds an index, depending on what type of light curve is used.
        Intended as the first output from tier 2.
    output_dir : string
        The directory in which you wish to save the flare csv files to.
    Teff : float
        The effective stellar temperature of the host star. If not passed, 4000 K will be assumed. Used to estimate 
        flare energies.
    host_radiu : float
        The radius of the host star. If not passed, a value of 1 solar radii is assumed. Used to estimate the flare energies.

    Returns
    -------
    Snippets of the flares found from 

    '''
    TOI_ID_list = []
    flare_number = []
    peak_time = []
    amplitude = []
    time_scale = []
    Teff = []
    radius = []
    total_flare_energies = []
    flare_amplitude = []
    flare_time_scale = []
    flare_energy = []
    accepted_flare_index = []
    accepted_flare_number = []
    index = 0
    flare_count = 1
    total_flares = 0
    os.makedirs(output_dir + '/' + str(host_name), exist_ok=True)
    for flare_events in flares:
        if flare_events >= 50 and len(flux) - flare_events > 50:
            new_time = time[flare_events-50:flare_events+50]
            new_data = flux[flare_events-50:flare_events+50]
        elif flare_events < 50:
            new_time = time[0+flare_events:flare_events+50]
            new_data = flux[0+flare_events:flare_events+50]
        elif len(flux) - flare_events < 50:
            new_time = time[flare_events:]
            new_data = flux[flare_events:]
        new_error = pdcsap_error[flare_events-50:flare_events+50]
        recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
        norm_time = time[flare_events]
        events = np.where(new_data == recenter)[0][0]
        new_time = (new_time - new_time[events])*24*60
        if lengths[index] >= 25:
            alles_data = new_data/np.median(new_data)
            error = new_error/np.median(new_data)
            popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
            squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
            chi2_cutoff = 20.843
        elif lengths[index] >= 15 and lengths[index] < 25:
            alles_data = new_data/np.median(new_data)
            error = new_error/np.median(new_data)
            popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
            squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
            chi2_cutoff = 11.912
        elif lengths[index] > 5 and lengths[index] < 15:
            alles_data = new_data/np.median(new_data)
            error = new_error/np.median(new_data)
            popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
            squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
            chi2_cutoff = 3.455
        elif lengths[index] <= 5:
            alles_data = new_data/np.median(new_data)
            error = new_error/np.median(new_data)
            popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
            squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
            chi2_cutoff = 1.3
        chi_squared = np.sum(squares)
        if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
            time_scale.append(popt[1])
            flare_time_scale.append(popt[1])
            amplitude.append(popt[0])
            flare_amplitude.append(popt[0])
            peak_time.append(norm_time)
            TOI_ID_list.append(host_name)
            flare_number.append(flare_count)
            try:
                energy = bolo_flare_energy(popt, T, host_radius, pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000), t_unit='days')
            except:
                energy = np.NaN
            flare_energy.append(energy)
            total_flare_energies.append(energy)
            Teff.append(T)
            radius.append(host_radius)
            try:
                X = np.column_stack((new_time, alles_data, error))
            except:
                X = np.column_stack(([0],[0],[0]))
            accepted_flare_index.append(flares[index])
            accepted_flare_number.append(flare_count)
            print(flare_count)
            if lengths[index] > 5:
                print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
            np.savetxt(output_dir + '/' + str(host_name) + '/Flare_' + str(index) + '.csv', X, delimiter=',')
            flare_count += 1
            total_flares += 1
            ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, Teff, radius, flare_energy))
            with open(output_dir + '/' + str(host_name) + '/All_Flare_Parameters.csv', "a") as f:
                np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                f.close()
            index += 1
                    
def tier3(tier_2_output_dir, tier_3_working_dir, tier_3_output_dir, settings_template_dir, params_template_dir, host_name = 'My_Host', T = 4000, host_radius = 1, MCMC_CPUS = 1):
    #Check each folder
    flare_files = os.listdir(tier_2_output_dir)
    for csvs in flare_files:
        parameter_list_list = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
        # try:
        #Make output directory to store plots and flare data for each target
        #Omit the statistics files
        if csvs == 'All_Flare_Parameters.csv' or csvs == 'Host_Statistics.txt' or csvs == 'Flare_Phase.csv':
            continue
        # Copy relevant file to allesfitter folder
        shutil.copyfile(tier_2_output_dir + '/' + csvs, tier_3_working_dir+ '/' +  csvs)
        allesfitter_priors.csv_cleaner(tier_3_working_dir + '/' + csvs)
        # Run allesfitter
        allesfitter_priors.flare_params(tier_3_working_dir + '/' +  csvs, params_template_dir, tier_3_working_dir + '/params.csv')
        allesfitter_priors.flare_settings(tier_3_working_dir + '/' +  csvs, settings_template_dir, tier_3_working_dir + '/settings.csv', multi_process=True, cores=MCMC_CPUS)
        allesfitter.mcmc_fit(tier_3_working_dir)
        allesfitter.mcmc_output(tier_3_working_dir)

        #Extract relevant parameters and append them to relevant lists
        # Full Data csv
        parameter_list = allesfitter_priors.return_parameters(tier_3_working_dir + '/results/mcmc_table.csv')
        if len(parameter_list) == 0:
            continue
        energy = allesfitter_priors.flare_energy(parameter_list[3], parameter_list[6], T,  host_radius)
        parameter_list_list[0].append(host_name)
        parameter_list_list[1].append(csvs[:-4])
        for index in range(9):
            parameter_list_list[index+2].append(parameter_list[index])
        parameter_list_list[11].append(energy)
        parameter_list_list[12].append(T)
        parameter_list_list[13].append(host_radius)


        #Per target csv
        #Remove current csv
        os.remove(tier_3_working_dir + '/' +  csvs)
        #Copy relevant graphs and log files
        shutil.copyfile(tier_3_working_dir + '/results/mcmc_corner.pdf', tier_3_output_dir + '/mcmc_corner_' + csvs[:-4] + '.pdf')
        shutil.copyfile(tier_3_working_dir + '/results/mcmc_fit_b.pdf', tier_3_output_dir + '/mcmc_fit_' + csvs[:-4] + '.pdf')
        result_dir = os.listdir(tier_3_working_dir + '/results')
        for dirs in result_dir:
            if dirs[-3] == 'l':
                shutil.copyfile(tier_3_working_dir + '/results/' + dirs, tier_3_output_dir + '/mcmc_' + csvs[:-4] + '.log')
        #Delete results folder
        for files in os.listdir(tier_3_working_dir + '/results'):
            os.remove(tier_3_working_dir + '/results/' + files)
        # except:
        #     print(1)
        #     parameter_list_list[0].append(host_name)
        #     parameter_list_list[1].append(csvs[:-4])
        #     for index in range(9):
        #         parameter_list_list[index+2].append(np.nan)
        #     parameter_list_list[11].append(np.nan)
        #     parameter_list_list[12].append(T)
        #     parameter_list_list[13].append(host_radius)
        ZZ = np.column_stack((parameter_list_list[0], parameter_list_list[1], parameter_list_list[2], parameter_list_list[3], parameter_list_list[4], parameter_list_list[5], parameter_list_list[6], parameter_list_list[7], parameter_list_list[8], parameter_list_list[9], parameter_list_list[10], parameter_list_list[11], parameter_list_list[12], parameter_list_list[13]))
        with open(tier_3_output_dir + '/All_TOI_MCMC_Flares.csv', "a") as f:
            np.savetxt(f, ZZ, delimiter=",", fmt='%s')
            f.close()
TESS_Folder_ID = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/')
TOI_Catalog = pd.read_csv('/data/whitsett.n/Reference_Files/All_TOI_12_17_23.csv')
total_flares = 0
total_possible_flares = 0
total_observation_time = 0
TOI_ID_list = []
flare_number = []
peak_time = []
amplitude = []
time_scale = []
Teff = []
radius = []
flare_phase = []
TOI_period = []
total_flare_energies = []
total_flare_phases = []
item_count = 0
total_periastron_list = []
total_periastron_epoch_list = []
total_epoch_list = []
list_index = 0

tier0_tau = []
tier1_tau = []
tier2_tau = []
for M_dwarves in TESS_Folder_ID:
    ##Iteration Scheme
    TOI_ID = M_dwarves
    print(TOI_ID)

    a = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves)
    os.makedirs('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/', exist_ok=True)
    print(item_count, M_dwarves)
    ##Relevant parameters from TOI catalog
    period = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_orbper'])[0]
    epoch = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_tranmid'])[0]
    stellar_radius = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_rad'])[0]
    T = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_teff'])[0]
    # periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbper'])[0]
    # epoch_periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbtper'])[0]
    if T == '':
        T = np.NAN
    if stellar_radius == '':
        stellar_radius = np.NAN
    # if epoch_periastron == '':
    #     epoch_periastron = np.NAN
    # if periastron == '':
    #     periastron = np.NAN
    if epoch == '':
        epoch = np.NAN
    flare_count = 1
    
    ##Trackable values per star
    phase_folded_flare_list = []
    flare_amplitude = []
    flare_time_scale = []
    flare_energy = []
    accepted_flare_index = []
    flare_phase = []
    accepted_flare_number = []
    observation_time = 0
    possible_flares = 0
    # periastron_list = []
    # periastron_epoch_list = []
    epoch_list = []
    
    ##Total Trackable Lists for all Data
    TOI_ID_list = []
    flare_number = []
    peak_time = []
    amplitude = []
    time_scale = []
    Teff = []
    radius = []
    flare_phase = []
    TOI_period = []
    total_flare_energies = []
    total_flare_phases = []
    item_count = 0
    # total_periastron_list = []
    # total_periastron_epoch_list = []
    total_epoch_list = []
    for folders in a:
        t1 = timer.time()
        try:
            b, pdcsap_flux, pdcsap_error = TESS_data_extract('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves + '/' + folders, PDCSAP_ERR=True)
        except:
            continue
        if folders.endswith('a_fast-lc.fits') == True:
            observation_time += len(pdcsap_flux)*(0.33333333)
            total_observation_time += len(pdcsap_flux)*(0.33333333)
            fast = True
        elif folders.endswith('a_fast-lc.fits') == False:  
            observation_time += len(pdcsap_flux)*2
            total_observation_time += len(pdcsap_flux)*2
            fast = False
        time, flux, pdcsap_error = delete_nans(b, pdcsap_flux, pdcsap_error)
        detrend_flux = SMA_detrend(time, flux, 80)
        t2 = timer.time()
        tier0_tau.append(t2-t1)
        flares, lengths = flare_ID(detrend_flux, 3, fast=fast)
        index = 0
        possible_flares += len(flares)
        total_possible_flares += len(flares)

        for flare_events in flares:
            if flare_events >= 50 and len(flux) - flare_events > 50:
                new_time = time[flare_events-50:flare_events+50]
                new_data = flux[flare_events-50:flare_events+50]
                new_error = pdcsap_error[flare_events-50:flare_events+50]
                recenter = np.max(new_data[45:50+lengths[index]])
            elif flare_events < 50:
                new_time = time[0:flare_events+50]
                new_data = flux[0:flare_events+50]
                new_error = pdcsap_error[flare_events:flare_events+50]
                try:
                    recenter = np.max(new_data[flare_events-5:flare_events+lengths[index]])
                except:
                    continue
            elif len(flux) - flare_events < 50:
                new_time = time[flare_events:]
                new_data = flux[flare_events:]
                new_error = pdcsap_error[flare_events:]
                try:
                    recenter = np.max(new_data[flare_events:])
                except:
                    continue
            norm_time = time[flare_events]
            events = np.where(new_data == recenter)[0][0]
            if period != '':
                phase = phase_folder(new_time, period, epoch)
                flare_phase_value = phase[events]
            elif period == '':
                phase = []
                flare_phase_value = np.NAN
                period = np.NAN
            t3 = timer.time()
            tier1_tau.append(t3-t2)
            new_time = (new_time - new_time[events])*24*60
            criteria1 = False
            if recenter > np.mean(new_data)+3*(np.std(new_data)):
                criteria1 = True
            try:
                if criteria1 == True and new_data[events+1] > np.mean(new_data)+(np.std(new_data)):
                    if lengths[index] >= 25:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
                        squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
                        chi2_cutoff = 20.514
                    elif lengths[index] >= 15 and lengths[index] < 25:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
                        squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
                        chi2_cutoff = 12.554
                    elif lengths[index] > 5 and lengths[index] < 15:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
                        squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
                        chi2_cutoff = 3.55
                    elif lengths[index] <= 5:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
                        squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
                        chi2_cutoff = 1.455
                    chi_squared = np.sum(squares)
                    if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
                        time_scale.append(popt[1])
                        flare_time_scale.append(popt[1])
                        amplitude.append(popt[0])
                        flare_amplitude.append(popt[0])
                        peak_time.append(norm_time)
                        
                        TOI_ID_list.append(TOI_ID)
                        flare_number.append(flare_count)
                        try:
                            energy = bolo_flare_energy(popt, T, stellar_radius, pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000), t_unit='days')
                        except:
                            energy = np.NaN
                        flare_energy.append(energy)
                        total_flare_energies.append(energy)
                        Teff.append(T)
                        radius.append(stellar_radius)
                        try:
                            X = np.column_stack((new_time, alles_data, error))
                        except:
                            X = np.column_stack(([0],[0],[0]))
                        baseline = st.median(new_data)*(lengths[index])*2
                        median = st.median(new_data)
                        flare_phase.append(flare_phase_value)
                        accepted_flare_index.append(flares[index])
                        accepted_flare_number.append(flare_count)
                        total_flare_phases.append(flare_phase_value)
                        TOI_period.append(period)
                        # total_periastron_list.append(periastron)
                        # total_periastron_epoch_list.append(epoch_periastron)
                        # periastron_list.append(periastron)
                        # periastron_epoch_list.append(epoch_periastron)
                        total_epoch_list.append(epoch)
                        epoch_list.append(epoch)
                        print(flare_count)
                        if lengths[index] > 5:
                            print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
                        np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
                        flare_count += 1
                        total_flares += 1
                        t4 = timer.time()
                        tier2_tau.append(t4-t3)
            except:
                continue
                        
                            
                index += 1
        ZZ = np.column_stack((np.array(np.array(tier0_tau).mean()), np.array(np.array(tier1_tau).mean()/len(flares)), np.array(np.array(tier2_tau).mean())))
        with open("/data/whitsett.n/Pipeline_Validation/Time_Stats.csv", "a") as f:
            np.savetxt(f, ZZ, delimiter=",", fmt='%s')
            f.close()
        tier0_tau = []
        tier1_tau = []
        tier2_tau = []
    if observation_time == 0:
        observation_time = -1
    if possible_flares == 0 :
        possible_flares = -1
    Y = np.column_stack((flare_amplitude, flare_time_scale, flare_energy))
    np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/All_Flares.csv', Y, delimiter=',')
    f = open('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Host_Statistics.txt', 'w')
    f.write('Flares/Day: ' + str(flare_count*60*24/observation_time) + '\n' + 'Possible Flares: ' + str(possible_flares) + '\n' +
            'Accepted Flares: ' + str(flare_count) + '\n' + 'Accepted/Possible Flares: ' + str(flare_count/possible_flares) + '\n' +
            'Observation Time (min): ' + str(observation_time))
    f.close()
    item_count += 1
    ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, total_flare_energies, Teff, radius, TOI_period))
    with open("/data/whitsett.n/Tier_3/All_TOI/All_TOI_Flares2.csv", "a") as f:
        np.savetxt(f, ZZ, delimiter=",", fmt='%s')
        f.close()
    print(len(TOI_ID_list), len(flare_number), len(peak_time), len(amplitude), len(time_scale), len(total_flare_energies), len(Teff), len(radius), len(TOI_period))
    print(list_index)
g = open('/data/whitsett.n/Tier_3/All_TOI/All_statistic_TOI.txt', 'w')
g.write('Total Flares: ' + str(total_flares) + '\n' + 'Net Observation Time (Days): ' + str(total_observation_time/(60*24)))
g.close()
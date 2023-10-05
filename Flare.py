# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:20:10 2023

@author: Nate Whitsett
"""

from astropy.io import fits
from astropy.timeseries import LombScargle as LS
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import statistics as st

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
    pdcsap_flux_method = hdul[1].header['PDCMETHD']
    if SAP_ERR == False and PDCSAP_ERR == False:
        return time, sap_flux, pdcsap_flux, pdcsap_flux_method
    if SAP_ERR == False and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux, pdcsap_flux_error, pdcsap_flux_method
    if SAP_ERR == True and PDCSAP_ERR == False:
        return time, sap_flux, sap_flux_error, pdcsap_flux, pdcsap_flux_method
    if SAP_ERR == True and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux, sap_flux_error, pdcsap_flux, pdcsap_flux_method


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
    min_time : int, optional
        The cadence of expected flares, in units of the time interval of the 
        provided data. This is to avoid counting multiple high sigma data 
        points from the same flare as multiple flares. The default is 3.

    Returns
    -------
    flare_indices : numpy array 
        Outputs a list of potential flares as lists of indices in the provided
        data. The index begins at the triggering data point.

    '''
    std_dev = np.std(data)
    sigma = sigma*std_dev
    mu = np.mean(data)
    count = 0
    flare_indices = []
    flare_length = 0
    flare_length_list = []
    flare = False
    for flux in data:
        if flux > (mu + sigma) and flare == False:
            flare = True
        if flare == True and flux > (mu+ 4*sigma):
            flare = False
        if flare == True and data[count+1] > (mu +sigma)/3:
            flare_length += 1
        if flare == True and data[count+1] < (mu +sigma)/3:
            flare = False
            flare_indices.append(count-flare_length)
            flare_length_list.append(flare_length)
            flare_length = 0
        count += 1
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

def SMA_detrend(data, time_scale):
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
    SMA = np.array(mov_average)
    return (data - SMA), SMA
        
def flare_phase_folded_ID(phase, flare_array, period, epoch):
    new_ID_list = []
    for indices in flare_array:
        new_ID_list.append(((phase[indices] - (epoch+period/2)) % period)-period/2)
    return np.array(new_ID_list)



file_path = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/M Dwarf Hosts/TOI 122.01/tess2018206045859-s0001-0000000231702397-0120-s/tess2018206045859-s0001-0000000231702397-0120-s_lc.fits'
a, sap_flux, pdcsap_flux, method = TESS_data_extract(file_path)
ID = file_path[::-1]
ID2 = ID[15:24]
TESS_ID = ID2[::-1]
time, flux = delete_nans(a, pdcsap_flux)
detrend_flux, detrend_array = SMA_detrend(flux, 100)
median = st.median(detrend_flux.tolist())
detrend_flux = detrend_flux / median
# phase = phase_folder(time, 2.7299, 1413.2)
flares, lengths = flare_ID(detrend_flux, 3)
# new_flares = flare_phase_folded_ID(phase, flares, 2.7299, 1413.2)

parameter_list = []
chi_square = []
# plt.plot(time, detrend_flux)
# for events in flares:
#     plt.axvline(time[events], c='green', linewidth=0.5)
count = 0
acceptable_flare = 0
possible_flares = len(flares)
dataframe_export = pd.DataFrame()
for events in flares:
    detrend_correction = detrend_array
    normalized_time = time - time[events]
    new_time = normalized_time.tolist()
    new_data = detrend_flux.tolist()
    if lengths[count] >= 25 and lengths[count] <= 35:
        new_time = np.array(new_time[events-5:events+30])*24*60
        new_data = np.array(new_data[events-5:events+30])
        detrend_correction = detrend_array[events-5:events+30]
        chi2_cutoff = 20.599
    if lengths[count] >= 15 and lengths[count] < 25:
        new_time = np.array(new_time[events-5:events+20])*24*60
        new_data = np.array(new_data[events-5:events+20])
        detrend_array = detrend_array[events-5:events+20]
        chi2_cutoff = 10.085
    if lengths[count] > 5 and lengths[count] < 15:
        new_time = np.array(new_time[events-5:events+10])*24*60
        new_data = np.array(new_data[events-5:events+10])
        detrend_correction = detrend_array[events-5:events+10]
        chi2_cutoff = 2.833
    else:
        new_time = np.array(new_time[events-5:events+7])*24*60
        new_data = np.array(new_data[events-5:events+7])
        detrend_correction = detrend_array[events-5:events+7]
        chi2_cutoff = 1.064
    new_data = new_data - np.mean(new_data)
    popt, pcov = curve_fit(exp_decay, new_time[5:15], new_data[5:15], maxfev=5000)
    half_max = new_data[5]/2
    amp = new_data.max()
    FWHMt = np.log((half_max - popt[2])/popt[0])/(-popt[1]) 
    squares = (new_data[5:15] - exp_decay(new_time[5:15], *popt))**2/(np.var(new_data[5:15]))
    chi_squared = np.sum(squares)
    if chi_squared < chi2_cutoff:
        dataframe_export['Flare ' + str(events) + ' Time'] = new_time
        dataframe_export['Flare ' + str(events) + ' Data'] = new_data
        dataframe_export['Flare ' + str(events) + ' Detrend Correction'] = detrend_correction
        dataframe_export.reset_index()
        dataframe_export.to_csv('C:/Users/Nate Whitsett/Desktop/Flares/TOI 122.01/Flares.csv', index=False)
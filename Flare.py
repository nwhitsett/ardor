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
# from aflare.py import aflare1 as af

def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def flare_fit(time, data):
    popt, pcov = curve_fit(exp_decay, time, data)
    plt.plot(time, exp_decay(time, *popt))
    plt.plot(time, data, 'b-', label='data')
    return popt


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

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
    if SAP_ERR == False and PDCSAP_ERR == False:
        return time, sap_flux, pdcsap_flux
    if SAP_ERR == False and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux, pdcsap_flux_error
    if SAP_ERR == True and PDCSAP_ERR == False:
        return time, sap_flux, sap_flux_error, pdcsap_flux
    if SAP_ERR == True and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux, sap_flux_error, pdcsap_flux


def phase_folder(time, period, epoch):
    '''

    Parameters
    ----------
    time : np.array
        The time array to be phase folded
    period : float
        The period of the planet, in the same units as the time array
    epoch : float
        The time of the first transit (i.e. transit epoch), same unit as time

    Returns
    -------
    phase : np.array
        Numpy array of the same dimension as the input time array, but phase
        folded to the period of the planet. Transit centered at 0.


    '''
    phase = (time - (epoch+period/2)) % period
    phase = phase - period/2
    return phase

def flare_ID(data, sigma, min_time = 3):
    std_dev = np.std(data)
    sigma = sigma*std_dev
    mu = np.mean(data)
    count = 0
    flare_indices = []
    flare_length = 0
    flare = False
    for flux in data:
        if flux > (mu + sigma) and flare == False:
            flare = True
        if flare == True and flux > (mu+ 4*sigma):
            flare = False
        if flare == True and data[count: (count+ min_time)].sum()/min_time > (mu +sigma):
            flare_length += 1
        if flare == True and data[count: (count+min_time)].sum()/min_time < (mu +sigma):
            flare = False
            flare_indices.append(count-flare_length)
            flare_length = 0
        count += 1
    return flare_indices

def delete_nans(time, data):
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
    return data - SMA
        
def flare_phase_folded_ID(phase, flare_array, period, epoch):
    new_ID_list = []
    for indices in flare_array:
        new_ID_list.append(((time[indices] - (epoch+period/2)) % period)-period/2)
    return np.array(new_ID_list)


        
file_path = 'C:/Users/Nathan/OneDrive - Washington University in St. Louis/Grad School/Spring 2023/Research/HST Cycle 31/MAST Data/K2-25 b/TESS/tess2021284114741-s0044-0000000434226736-0215-s/tess2021284114741-s0044-0000000434226736-0215-s_lc.fits'
a, sap_flux, pdcsap_flux = TESS_data_extract('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Grad School/Spring 2023/Research/HST Cycle 31/MAST Data/K2-25 b/TESS/tess2021284114741-s0044-0000000434226736-0215-s/tess2021284114741-s0044-0000000434226736-0215-s_lc.fits')
ID = file_path[::-1]
ID2 = ID[15:24]
TESS_ID = ID2[::-1]
time, flux = delete_nans(a, pdcsap_flux)
detrend_flux = SMA_detrend(flux, 80)
median = st.median(detrend_flux.tolist())
detrend_flux = detrend_flux / median
phase = phase_folder(time, 2.7299, 1413.2)
flares = flare_ID(detrend_flux, 5, min_time=4)
new_flares = flare_phase_folded_ID(phase, flares, 2.7299, 1413.2)

parameter_list = []
RSS_list = []
# for events in flares:
#     plt.axvline(time[events], c='green', linewidth=0.5)

for events in flares:
    count = 1
    normalized_time = time - time[events]
    new_time = normalized_time.tolist()
    new_data = detrend_flux.tolist()
    new_time = np.array(new_time[events-5:events+10])
    new_data = np.array(new_data[events-5:events+10])
    popt, pcov, infodict, mesg, ier  = curve_fit(exp_decay, new_time[6::], new_data[6::], full_output=True, maxfev=5000)
    plt.plot(new_time[4::]*2*24*60+4, exp_decay(new_time[4::], *popt), linestyle='dotted', label = 'Flare candidate ' + str(count) + ' model')
    plt.scatter(new_time*2*24*60+4, new_data, s=7, label='Flare candidate ' + str(count) + ' data', color='red')
    plt.legend()
    plt.xlabel('Time Since Flare (min)')
    plt.ylabel('Standardized Relative Flux')
    plt.title('TOI ' + TESS_ID + ' Flare Index ' + str(events))

    plt.show()
    plt.clf()
    check = input('Continue?')
    parameter_list.append(str("a = " + str(popt[0]) + " +/- " + str(pcov[0,0]**0.5)))
    parameter_list.append(str("b = " + str(popt[1]) + " +/- " + str(pcov[1,1]**0.5)))
    parameter_list.append(str("c = " + str(popt[2]) + " +/- " + str(pcov[2,2]**0.5)))
    for values in infodict['fvec']:
        RSS = 0
        RSS += values**2
    RSS_list.append(list((events, RSS)))
    count += 1
print(RSS_list)

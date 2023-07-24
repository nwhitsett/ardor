# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:20:10 2023

@author: Nate Whitsett
"""

from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def TESS_FITS_csv(input_file, csv_directory, csv_name=None):
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


def phase_folder(time, data, period, epoch):
    phase = (time - epoch) % period
    return phase, data

time, sap_flux, pdcsap_flux = TESS_data_extract('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Spring 2023/Research/HST Cycle 31/MAST Data/K2-25 b/TESS/tess2021284114741-s0044-0000000434226736-0215-s/tess2021284114741-s0044-0000000434226736-0215-s_lc.fits')
phase, flux = phase_folder(time, sap_flux, 3.48475, 2501.7741)
reduced = flux - 767
plt.scatter(time, sap_flux, s=0.5)
plt.title('SAP_FLUX')
plt.plot()

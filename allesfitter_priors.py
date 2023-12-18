# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:48:26 2023

@author: Nathan
"""

import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import os
import aflare
import math

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx
def flare_settings(data_file_dir, settings_template_dir, output_dir, flares=1, multi_process = False, cores = 1):
    settings = pd.read_csv(settings_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    settings.at[5, 'value'] = name
    settings.at[10, 'value'] = multi_process
    if multi_process == True:
        settings.at[11, 'value'] = cores
    settings.at[33, '#name'] = 'baseline_flux_' + str(name)
    settings.at[37, '#name'] = 'error_flux_' + str(name)
    settings.at[41, 'value'] = flares
    settings.to_csv(output_dir, index=False)
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
def allesfitter_priors(flare_csv, time_unit='min'):
    if time_unit == 'min':
        time_unit = 1140
    if time_unit == 'days':
        time_unit = 1
    flare = pd.read_csv(flare_csv, header=None)
    time = flare[0]
    flux_median = np.median(flare[1][:25])
    flux = flare[1] - flux_median
    params, pcov = curve_fit(exp_decay, time[30:], flux[30:])
    param_error = np.sqrt(np.diag(pcov))
    x = np.linspace(0, 10, num = 2000)/(time_unit)
    y = exp_decay(x, params[0], params[1]*time_unit, params[2])
    amp, idx = find_nearest(y, y.max()/2)
    amp = 2*amp
    tau = x[idx]
    plt.scatter(time/time_unit, flux)
    error = param_error[1]/time_unit
    x_aflare = np.linspace(-0.02, 0.04, num=1000)
    y_aflare = aflare.aflare1(x_aflare, 0, tau, amp)
    plt.plot(x,y)
    plt.plot(x_aflare, y_aflare)
    return amp, tau, error*10

def flare_params(data_file_dir, params_template_dir, output_dir, flares=1):
    amp, tau, tau_error = allesfitter_priors(data_file_dir)
    data = pd.read_csv(data_file_dir, index_col = False, header=None)
    error = data[2]
    for index in range(len(error)):
        if math.isnan(error[index]) == True:
            error[index] = 1e-3
    flux_err = np.log10(np.median(error))
    params = pd.read_csv(params_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    params.at[1, 'value'] = 0
    params.at[2, 'value'] = tau
    params.at[3, 'value'] = amp
    params.at[5, 'value'] = np.log10(np.median(error))
    
    params.at[1, 'bounds'] = 'uniform ' + str(-0.0001) + ' ' + str(0.0001)
    params.at[2, 'bounds'] = 'uniform ' + str(amp-0.5) + ' ' + str(amp+0.05)
    params.at[3, 'bounds'] = 'uniform ' + str(tau-tau_error) + ' ' + str(tau + tau_error)
    params.at[5, 'bounds'] = 'uniform ' + str(flux_err - 1) + ' ' + str(flux_err + 1)
    # if flares > 1:
    #     for number in range(2, flares+1):
    #         params.at[1+3*(number-1), '#name'] = 'flare_tpeak_' + str(number)
    #         params.at[2+3*(number-1), '#name'] = 'flare_fwhm_' + str(number)
    #         params.at[3+3*(number-1), '#name'] = 'flare_ampl_' + str(number)
            
    #     params.at[5+3*(flares-1), '#name'] = '#errors (overall scaling) per instrument'
    #     params.at[6+3*(flares-1), '#name'] = 'log_err_flux_' + str(name)
    #     params.at[6+3*(flares-1), 'label'] = '$\log{\sigma_\mathrm{' + str(name) + '}}$'
    # else:
    params.at[5, '#name'] = '#errors (overall scaling) per instrument'
    params.at[5, '#name'] = 'log_err_flux_' + str(name)
    params.at[5, 'label'] = '$\log{\sigma_\mathrm{' + str(name) + '}}$'
    params.to_csv(output_dir, index=False)


allesfitter_priors('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Flares_New/TOI 136.01/Flare5.csv')

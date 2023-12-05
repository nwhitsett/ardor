# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:48:26 2023

@author: Nathan
"""

import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simpson
import numpy as np
import pandas as pd
import os
import aflare
import math
import planck_law

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
    flare = pd.read_csv(flare_csv, header=None)
    count = 0
    for vals in flare[0]:
        if np.log10(np.abs(vals)) > 2:
            count += 1
    flare.drop(index=flare.index[:count], axis=0, inplace=True)
    flare.reset_index(inplace=True, drop=True)
    time = flare[0]
    flux = flare[1] - 1.0
    print(flare)
    try:
        if np.abs(time[0]) > 1:
            params, pcov = curve_fit(exp_decay, time[(30-count):], flux[(30-count):], maxfev=5000)
            time_unit = 1140
        elif np.abs(time[0]) < 1:
            params, pcov = curve_fit(exp_decay, time[(30-count):]*1140, flux[(30-count):], maxfev=5000)
            time_unit = 1
    except:
        params = (0, 5.0, 0.01)
    x = np.linspace(0, 10, num = 1500)/(1140)
    y = exp_decay(x, params[0], params[1]*time_unit, params[2])
    amp, idx = find_nearest(y, y.max()/2)
    amp = flux.max()
    tau = x[idx]
    for index in range(len(flare[2])):
        if math.isnan(flare[2][index]) == True:
            flare[2][index] = 1e-3
    if np.abs(flare[0][0]) > 1:
        flare[0] = flare[0]/time_unit
    flare.to_csv(flare_csv, index=False, header=False)
    return amp, tau

def flare_params(data_file_dir, params_template_dir, output_dir, flares=1):
    amp, tau = allesfitter_priors(data_file_dir)
    data = pd.read_csv(data_file_dir, index_col = False, header=None)
    error = data[2]
    for index in range(len(error)):
        if math.isnan(error[index]) == True:
            error[index] = 1e-3
    flux_err = np.log(np.median(error))
    params = pd.read_csv(params_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    params.at[1, 'value'] = 0
    params.at[2, 'value'] = tau
    params.at[3, 'value'] = amp
    params.at[5, 'value'] = flux_err
    params.at[1, 'bounds'] = 'uniform ' + str(-0.001) + ' ' + str(0.001)
    params.at[2, 'bounds'] = 'uniform ' + str(0) + ' ' + str(2.5*tau)
    params.at[3, 'bounds'] = 'uniform ' + str(amp*0.95) + ' ' + str(3.5*amp)
    params.at[5, 'bounds'] = 'uniform ' + str(-1 + flux_err) + ' ' + str(1 + flux_err)
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
    params.at[5, '#name'] = 'ln_err_flux_' + str(name)
    params.at[5, 'label'] = '$\ln{\sigma_\mathrm{' + str(name) + '}}$'
    params.to_csv(output_dir, index=False)

def return_parameters(mcmc_table_dir):
    data = pd.read_csv(mcmc_table_dir)
    t_peak = data['median'][1]
    t_peak_m = data['lower_error'][1]
    t_peak_p = data['upper_error'][1]
    fwhm = data['median'][2]
    fwhm_m = data['lower_error'][2]
    fwhm_p = data['upper_error'][2]
    amp = data['median'][3]
    amp_m = data['lower_error'][3]
    amp_p = data['upper_error'][3]

    return [t_peak, t_peak_m, t_peak_p, fwhm, fwhm_m, fwhm_p, amp, amp_m, amp_p]

def flare_energy(fwhm, ampl, Teff, R_stellar):
    x = np.linspace(0, 0.02, num = 2000)
    y = aflare.aflare1(x, 0.01, fwhm, ampl)
    flare_area = simpson(y, x)
    print(fwhm, ampl, Teff, R_stellar)
    color_factor = planck_law.planck_integrator(600e-6, 1000e-6, Teff)/planck_law.planck_integrator(600e-6, 1000e-6, 9000)
    energy = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
    return energy
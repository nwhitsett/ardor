# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:48:26 2023

@author: Nathan
"""

import pandas as pd


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

def flare_params(data_file_dir, params_template_dir, output_dir, param_list=None, error_list=None, flares=1):
    params = pd.read_csv(params_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    if flares > 1:
        for number in range(2, flares+1):
            params.at[1+3*(number-1), '#name'] = 'flare_tpeak_' + str(number)
            params.at[2+3*(number-1), '#name'] = 'flare_fwhm_' + str(number)
            params.at[3+3*(number-1), '#name'] = 'flare_ampl_' + str(number)
            
        params.at[5+3*(flares-1), '#name'] = '#errors (overall scaling) per instrument'
        params.at[6+3*(flares-1), '#name'] = 'log_err_flux_' + str(name)
        params.at[6+3*(flares-1), 'label'] = '$\log{\sigma_\mathrm{' + str(name) + '}}$'
    else:
        params.at[5, '#name'] = '#errors (overall scaling) per instrument'
        params.at[6, '#name'] = 'log_err_flux_' + str(name)
        params.at[6, 'label'] = '$\log{\sigma_\mathrm{' + str(name) + '}}$'
    params.to_csv(output_dir, index=False)
    
flare_params('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/settings.csv', 'C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/params.csv', 'C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/params1.csv', flares = 3)
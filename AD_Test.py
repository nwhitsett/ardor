# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
from skgof import ks_test, cvm_test, ad_test
from matplotlib import pyplot as plt
from scipy.stats import uniform
from scipy.stats import ks_1samp
from scipy.stats import ksone
import csv
from itertools import zip_longest

def ks_critical_value(n_trials, alpha):
    return ksone.ppf(1-alpha/2, n_trials)
TOI_flares = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv")
Exoplanet_flares = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv")
count = 0
KS_p_values = []
KS_center_statistic = []
AD_p_values = []
AD_center_statistic = []
test = []
count1 = 0
# fig, ax = plt.subplots(2,2)
N = 100

KS_p_values_per_host = []
KS_center_statistic_per_host = []
AD_p_values_per_host = []
AD_center_statistic_per_host = []

interest_hosts = []
host_num = 0
All_CDF = []
for hosts in np.array(TOI_flares['Host_Name']):
    name = hosts
    planet_period = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0]
    phases = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period_Phase'])/planet_period
    phases1 = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period_Phase'])/planet_period
    obs_time = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Obs_Time'])[0]
    # peri_phases = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period_Phase'])/planet_period
    r = uniform.rvs(size=50)
    count += 1
    if count == len(phases) and np.isnan(planet_period) == False and hosts == 'WASP-121':
        
        a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
        b = ad_test(phases, uniform(0,1))
        KS_orig_p = a.pvalue
        KS_orig_D = a.statistic_location
        AD_orig_p = b.pvalue
        AD_orig_D = b.statistic
        
        # if planet_period < 20 and (KS_orig_p > 0.95 or AD_orig_p > 0.95):
        #     interest_hosts.append((hosts, KS_orig_p, AD_orig_p, len(phases), np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0]))
        #     print(hosts, KS_orig_p, AD_orig_p, len(phases), np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0], obs_time)
        #     x = np.sort(phases1)
        #     y = np.arange(len(x))/float(len(x))
        #     plt.plot(x, y)
        #     plt.plot(np.linspace(0,1), np.linspace(0,1))
        #     plt.title(str(hosts) + ' CDF')
        #     plt.show()
        #     input('e')
        #     plt.clf()
        # KS_p_values.append(KS_orig_p)
        # KS_center_statistic.append(KS_orig_D)
        # AD_p_values.append(AD_orig_p)
        # AD_center_statistic.append(AD_orig_D)
        # phases -= 0.1
        for phase_sample in range(N):
            phases = phases + (1/N)
            for index in range(len(phases)):
                if phases[index] > 1:
                    phases[index] = phases[index] - 1  
            phases = np.sort(phases)
            a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
            b = ad_test(phases, uniform(0,1), assume_sorted=True)
            
            KS_p_values.append(a.pvalue)
            KS_center_statistic.append(a.statistic_location)
            
            AD_p_values.append(b.pvalue)
            AD_center_statistic.append(b.statistic)
            count1 += 1
        KS_center_statistic_per_host.append(np.median(KS_center_statistic))
        KS_p_values_per_host.append(np.median(KS_p_values))
        AD_center_statistic_per_host.append(np.median(AD_center_statistic))
        AD_p_values_per_host.append(np.median(AD_p_values))
        if planet_period < 20 and (np.median(KS_p_values) < 0.05 or np.median(AD_p_values) < 0.05) and len(phases) >= 3:
            interest_hosts.append((hosts, np.median(KS_p_values), np.median(AD_p_values), len(phases), np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0]))
            print(hosts, np.median(KS_p_values), np.median(AD_p_values), len(phases), np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0])
            All_CDF.append(phases1)
            # x = np.sort(phases1)
            # y = np.arange(len(x))/float(len(x))
            # plt.plot(x, y)
            # plt.plot(np.linspace(0,1), np.linspace(0,1))
            # plt.show()
            # print(planet_period)
            # input('e')
            # plt.clf()
        KS_p_values = []
        KS_center_statistic = []
        AD_p_values = []
        AD_center_statistic = []
        test = []
        count = 0
        print(host_num)
        host_num += 1
        phases = phases1

with open("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_CDFs.csv","w+") as f:
    writer = csv.writer(f)
    for values in zip_longest(*All_CDF):
        writer.writerow(values)

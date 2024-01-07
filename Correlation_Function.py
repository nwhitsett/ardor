# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 20:41:03 2023

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import uniform
from scipy.stats import ks_1samp
flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare Output Files/All Exoplanets/All_Exoplanet_Flares.csv")
count = 0
p_values = []
center_statistic = []
test = []

for hosts in np.array(flares['# ID']):
    name = hosts
    planet_period = np.array(flares.loc[flares['# ID'] == hosts, 'Period (Days)'])[0]
    phases = np.array(flares.loc[flares['# ID'] == hosts, 'Phase'])/planet_period
    count += 1
    r = uniform.rvs(size=50)
    if count == len(phases):
        c, d = ks_1samp(r, uniform.cdf, args=(0, 1), alternative='two-sided')
        a, b = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
        p_values.append(b)
        center_statistic.append(a)
        # if count > 5 and b < 0.01 and planet_period < 100:
        #     print(1)
        #     p_values.append(b)
        #     center_statistic.append(a)
            # x = np.sort(phases)
            # y = np.arange(len(phases))/float(len(phases))
            # plt.plot(x,y, label= hosts)
            # print(planet_period)
            
        count = 0
# plt.hist(center_statistic, bins=np.linspace(0,1, num=15))
plt.hist(p_values, bins=np.linspace(0,1, num=15))
# r = uniform.cdf(np.linspace(0, 1))
# x = np.sort(r)
# y = np.arange(len(x))/float(len(x))
# plt.plot(x, y, label='Uniform Distribution')
# plt.axvline(0.5, linewidth=1, linestyle='--')
# plt.title('CPF of Non-Uniform Flaring Exoplanet Hosts')
# plt.xlabel('Phase (Periastron Centered)')
# plt.ylabel('Cumulative Probability')

# plt.figure(dpi=1000)
# plt.savefig('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/CDF.png', dpi=1000)
# plt.show()
# periastron_phase = np.array(flares.loc[flares['# ID'] == "HD 41004 B", 'Flare Phase (Periastron Centered)'])

# planet_period = np.array(flares.loc[flares['# ID'] == "HD 41004 B", 'Period (Days)'])[0]
# periastron_phase = periastron_phase
# ks_1samp(periastron_phase, uniform.cdf, args=(-planet_period/2, planet_period))
# r = uniform.cdf(np.linspace(-planet_period/2, planet_period/2), loc = -planet_period/2, scale = planet_period)
# x = np.sort(periastron_phase)
# y = np.arange(len(x))/float(len(x))
# plt.plot(x, y)
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
flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare Output Files/Periastron_Flares.csv")
count = 0
for hosts in np.array(flares['# ID']):
    name = hosts
    planet_period = np.array(flares.loc[flares['# ID'] == hosts, 'Period'])[0]
    phases = np.array(flares.loc[flares['# ID'] == hosts, 'Flare Phase (Periastron Centered)'])
    count += 1
    if count == len(phases):
        a, b = ks_1samp(phases, uniform.cdf, args=(-planet_period/2, planet_period))
        if b < 0.05 and len(phases) > 5 and planet_period < 50:
            print(b, hosts, planet_period, len(phases))
        count = 0

    
periastron_phase = np.array(flares.loc[flares['# ID'] == "HD 41004 B", 'Flare Phase (Periastron Centered)'])

planet_period = np.array(flares.loc[flares['# ID'] == "HD 41004 B", 'Period'])[0]
periastron_phase = periastron_phase
ks_1samp(periastron_phase, uniform.cdf, args=(-planet_period/2, planet_period))
r = uniform.cdf(np.linspace(-planet_period/2, planet_period/2), loc = -planet_period/2, scale = planet_period)

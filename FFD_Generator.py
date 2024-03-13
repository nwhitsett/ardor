# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 13:26:00 2024

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def linear(x, alpha, beta):
    return alpha*x + beta
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
data = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv')

hosts = np.array(data['Host_Name'])

fig, ax = plt.subplots(2, figsize=(8,8))
host_ = 0
for host in hosts:
    if host != host_:
        energy_count_list = []
        energy_count_list2 = []
        obs_time = np.array(data.loc[data['Host_Name'] == host, 'Obs_Time'])[0]
        # energies = np.array(data.loc[(data['Host_Name'] == host) & (d), 'Flare_Energy'])
        energies = np.array(data.loc[(data['Host_Name'] == host) & (data['Vet'] != 3) & (data['Norm_Phase'] < 0.6) & (data['Norm_Phase'] > 0.4), 'Flare_Energy'])
        energies2 = np.array(data.loc[(data['Host_Name'] == host) & (data['Vet'] != 3), 'Flare_Energy'])

        energies[::-1].sort()
        energies = energies[~np.isnan(energies)]
        energies = np.log10(energies)
        energies2[::-1].sort()
        energies2 = energies2[~np.isnan(energies2)]
        energies2 = np.log10(energies2)
        if len(energies) > 10:
            a = np.linspace(energies[0], energies[-1], num = 10)
            c = np.linspace(30, 37, num = 10)
            for bins in a:
                b = sum(x >= bins for x in energies)
                energy_count_list.append(np.log10(b/obs_time))
            coeff, pcov = curve_fit(linear, a, energy_count_list)
            print(coeff)
            ax[0].scatter(a, energy_count_list, marker = 'o', label=host)
            ax[0].plot(c, linear(c, coeff[0], coeff[1]), ':')
        if len(energies2) > 40:
            a = np.linspace(energies2[0], energies2[-1], num = 10)
            c = np.linspace(30, 37, num = 10)
            for bins in a:
                b = sum(x >= bins for x in energies2)
                energy_count_list2.append(np.log10(b/obs_time))
            coeff, pcov = curve_fit(linear, a, energy_count_list2)
            print(coeff)
            ax[1].scatter(a, energy_count_list2, marker = 'o', label=host)
            ax[1].plot(c, linear(c, coeff[0], coeff[1]), ':')
    host_ = host
ax[0].set_ylim(bottom = -4, top = 1)
ax[0].set_title('FFD Plot: Constrained Periastron Hosts (N $\geq$ 10, Flare Phase $\pm$ 0.1)')
ax[0].set_ylabel('Flare Rate $\geq E_{bol} (d^{-1}) $', fontsize='large')
ax[0].set_xlabel('$Log_{10}(E_{bol})$', fontsize='large')
ax[1].set_ylim(bottom = -4, top = 1)
ax[1].set_title('FFD Plot: All Exoplanet Hosts (N$\geq$40)')
ax[1].set_ylabel('Flare Rate $\geq E_{bol} (d^{-1}) $', fontsize='large')
ax[1].set_xlabel('$Log_{10}(E_{bol})$', fontsize='large')
plt.subplots_adjust(wspace = 0.15, hspace = 0.3)
fig.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Publication Documents/FFD.png', dpi=fig.dpi)
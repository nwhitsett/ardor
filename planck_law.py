# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:59:11 2023

@author: Nathan
"""
import numpy as np
from scipy.integrate import simpson
from matplotlib import pyplot as plt
import pandas as pd
def planck_law(lambda0, T):
    h = 6.626e-34
    k = 1.38e-23
    c = 2.997e8
    return ((8*np.pi*h*c)/(lambda0**5))/(np.exp((h*c)/(lambda0*k*T))-1)

def planck_integrator(lambda_min, lambda_max, Teff):
    lambda_min = lambda_min
    lambda_max = lambda_max
    lambda_range = np.linspace(lambda_min, lambda_max, num=50)
    # R_lambda = pd.read_csv('C:/Users/Nate Whitsett/Desktop/tess-response-function-v2.0.csv', index_col=False)
    # R_B = []
    R_A = []
    for wavelengths in lambda_range:
    #     R_B.append(R_lambda.loc[R_lambda['Wavelength'] == wavelengths, 'lambda'].iloc[0]*planck_law(wavelengths*10**-9, Teff))
       R_A.append(planck_law(wavelengths, Teff))
    # integral = simpson(np.array(R_B), lambda_range)
    integral2 = simpson(np.array(R_A), lambda_range)
    return integral2

TOI_List = pd.read_csv('C:/Users/natha/OneDrive/Desktop/Kepler_Field.csv')
TOI_Teff = TOI_List['st_teff']
TOI_TESS_Mag = TOI_List['sy_tmag']
TOI_AB_Mag = []
T_eff_range = np.linspace(3000, 8000)
ratio = []
for index in T_eff_range:
    ratio.append(planck_integrator(0.23e-6, 0.29e-6, index)/planck_integrator(0.6e-6, 1e-6, index))
# for index in range(len(TOI_Teff)):
#     TOI_AB_Mag.append(-2.5*np.log10((2461*10**(-.4*TOI_TESS_Mag[index]))*(planck_integrator(0.23e-6, 0.29e-6, TOI_Teff[index])/planck_integrator(0.6e-6, 1e-6, TOI_Teff[index]))) + 8.9)
# TOI_List['AB Mag'] = np.array(TOI_AB_Mag)
# TOI_List.to_csv('C:/Users/natha/OneDrive/Desktop/All_Exoplanets.csv', index=False)
TOI_AB_Mag = TOI_List['AB Mag']
TOI_Dist = TOI_List['sy_dist']
bins = np.linspace(5, 25, num=15)
plt.title('Kepler Field TESS Mag vs ULTRASAT AB Mag')
plt.hist(TOI_TESS_Mag, bins, alpha=0.5, label='TESS Magnitude')
plt.hist(TOI_AB_Mag, bins, alpha=0.5, label='AB Magnitude')
# plt.hist(ratio, bins)
plt.legend()
plt.xlabel('Magnitude')
plt.ylabel('Count')
plt.axvline(21.8, linestyle='dashed', linewidth=1)
plt.show()


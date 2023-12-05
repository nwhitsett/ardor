# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 18:52:27 2023

@author: Nate Whitsett
"""
import numpy as np
import pandas as pd
import csv

data = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Injection_Test (Late Type M).csv')
amp = list(data['Amplitude'])
FWHM = list(data['FWHM'])
error = list(data['Error'])
Bool = list(data['Accepted?'])
amp_bins = np.linspace(0, 1, num=17)
FWHM_bins = np.linspace(0.001388888, 0.041, num=17)
error_bins = np.linspace(0.0001, 0.25, num=16)

x = []
y = []

for i in range(len(amp_bins)):
    tmp = []
    for index in range(len(amp)):
        if amp[index] > amp_bins[i] and amp[index] < amp_bins[i+1]:
            tmp.append([amp[index], FWHM[index], Bool[index]])
    x.append(tmp)
total = 0
for i in range(len(FWHM_bins)):
    tmp = []
    for cells in x:
        count = 0
        pos = 0
        for flares in cells:
            if flares[1] > FWHM_bins[i] and flares[1] < FWHM_bins[i+1]:
                count += 1
                total += 1
                if flares[2] == 1:
                    pos += 1
        if i == 16:
            continue
        if count == 0:
            tmp.append(np.nan)
        elif count != 0:
            tmp.append(pos/count * 100)
    y.append(tmp)

with open('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Injection_Test_Grid (Late M Dwarf).csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(y)
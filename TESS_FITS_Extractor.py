# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:20:10 2023

@author: Nate Whitsett
"""

from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
hdul = fits.open('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Spring 2023/Research/MAST_2023-05-31T1318/TESS/tess2018206190142-s0001-s0003-0000000100100827/tess2018206190142-s0001-s0003-0000000100100827-00129_dvt.fits')
#hdul.info()
raw_data = hdul[1].data

time = []
phase = []
data = []
error = []

for tuples in raw_data:
    time.append(tuples[0])
    data.append(tuples[4])
    error.append(tuples[5])
for values in time:
    phase.append(values % 0.94145455)
    

plt.plot(range(len(time)), time)
grand_list = pd.DataFrame({'time': time, 'phase': phase, 'data': data, 'error': error})
print(grand_list)
## Test
grand_list.to_csv('C:/Users/Nate Whitsett/Desktop/WASP18b_Phase_Curve')

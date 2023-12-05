# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:06:28 2023

@author: Nathan
"""

from astroquery.mast import Observations
import pandas as pd

target_list = pd.read_csv('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Grad School/Spring 2023/Research/HST Cycle 31/MAST Data/csv-file-toi-catalog.csv')
TOI_list = target_list['Full TOI ID']
RA = target_list['TIC Right Ascension']
DEC = target_list['TIC Declination']
M_dwarf = 0

for temps in target_list['Effective Temperature Value']:
    if temps < 4000:
        M_dwarf += 1
for index in range(M_dwarf):
    try:
        TOI_Name = 'TOI ' + str(TOI_list[index])
        obs_table = Observations.query_criteria(s_ra=RA[index], s_dec=DEC[index], project="TESS", dataproduct_type = "TIMESERIES", radius="0.02 deg", calib_level = 3, target_name= str(target_list['TIC'][index]))
        data_products = Observations.get_product_list(obs_table)
        for data in range(len(data_products)):
            ext_check = data_products[data][6][-6]
            if ext_check == 'c':
                Observations.download_products(data_products[data], download_dir = "C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/TESS Data/" + TOI_Name)
    except:
        continue
            
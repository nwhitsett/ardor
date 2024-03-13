# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:06:28 2023

@author: Nathan
"""

from astroquery.mast import Observations
import pandas as pd
import numpy as np
target_list = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare Output Files/All TOIs/csv-file-toi-catalog.csv")
TOI_list = target_list['Full TOI ID']
RA = target_list['TIC Right Ascension']
DEC = target_list['TIC Declination']
target = len(target_list)
TIC_ID = target_list['TIC']
# for temps in target_list['Effective Temperature Value']:
#     if temps < 4000:
#         M_dwarf += 1
for index in range(100):
    try:
        if str(TOI_list[index]).endswith('.01') == True:
            TOI_Name = str(TOI_list[index])
            ID = TIC_ID[index]
            # print(ID)
            # obs_table = Observations.query_criteria(target_name = ID, dataproduct_type = "TIMESERIES",  calib_level = 3)
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.02, float(RA[index])+0.02], s_dec = [float(DEC[index]) - 0.02, float(DEC[index])+0.02],dataproduct_type = ["TIMESERIES"], calib_level = 3, obs_collection='TESS')
            data_products = Observations.get_product_list(obs_table)
            print(np.log(np.sum(data_products['size'])/(1e6)))
            # for data in range(len(data_products)):
            #     ext_check = data_products[data][6][-6]
            #     if ext_check == 'c':
            #         Observations.download_products(data_products[data], download_dir = "/data/whitsett.n/TESS_Light_Curves/All_TOI/" + TOI_Name)
        else:
            continue
    except:
        print('No associated products')
        continue


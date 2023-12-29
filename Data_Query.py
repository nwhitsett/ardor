# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:06:28 2023

@author: Nathan
"""

from astroquery.mast import Observations
import pandas as pd

target_list = pd.read_csv("/data/whitsett.n/Reference_Files/All_TOI_12_17_23.csv")
TOI_list = target_list['toi']
RA = target_list['ra']
DEC = target_list['dec']
target = len(target_list)
TIC_ID = target_list['tid']
# for temps in target_list['Effective Temperature Value']:
#     if temps < 4000:
#         M_dwarf += 1
for index in range(target):
    try:
        TOI_Name = str(TOI_list[index])
        ID = TIC_ID[index]
        print(ID)
        obs_table = Observations.query_criteria(target_name = ID, dataproduct_type = "TIMESERIES", radius="0.02 deg", calib_level = 3)
        data_products = Observations.get_product_list(obs_table)
        for data in range(len(data_products)):
            ext_check = data_products[data][6][-6]
            if ext_check == 'c':
                Observations.download_products(data_products[data], download_dir = "/data/whitsett.n/TESS_Light_Curves/All_TOI/" + TOI_Name)
    except:
        print('No associated products')
        continue

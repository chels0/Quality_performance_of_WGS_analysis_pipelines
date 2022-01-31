#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:05:08 2022

@author: chelsea
"""

import pandas as pd
df_files= pd.read_csv('/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/FastQC_Raw_reads/PT28-1-18/PT28-1-18_R1_fastqc/xx03_updated.csv', sep=' ')

small = df_files[df_files['Mean'] < 28]

to_cut = small['#Base'].tolist()

for value in to_cut:
    to_cut = value.split('-')

kuk = []

for i in range(len(to_cut)):
    kuk.append(int(to_cut[i]))

cut = min(kuk)    

print(cut)    

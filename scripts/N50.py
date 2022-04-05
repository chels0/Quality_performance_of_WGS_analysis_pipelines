#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 14:45:05 2022

@author: chelsea
"""

import pandas as pd
import sys
import os
import pathlib

directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca_quast_tables/'

list_N50 = []

for filename in os.listdir(directory):
    df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
    df = df[['N50']]
    df = df.iloc[1: , :]    
    list_N50.append(df)
    

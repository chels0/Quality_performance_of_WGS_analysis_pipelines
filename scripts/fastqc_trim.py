#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:05:08 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
df_files= pd.read_csv('xx03.csv', sep=' ')

splt = np.array_split(df_files, 2)

small_end = splt[1]
small_end = small_end[small_end['Mean'] < 28]

to_cut_end = small_end['#Base'].tolist()

list_of_values= []

for value in to_cut_end:
    to_cut_end = value.split('-')
    list_of_values.append(to_cut_end[0])
    list_of_values.append(to_cut_end[1])    

list_end = []
    
for i in range(len(list_of_values)):
    list_end.append(int(list_of_values[i]))

if not list_end:
    trail = 1
else:
    min_end = min(list_end)
    max_end = max(list_end)
    trail = max_end - min_end

if trail > 20:
    trail = 20
else:
    trail = 28
    
print(trail)

    

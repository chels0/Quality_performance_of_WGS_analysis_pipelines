#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:22:38 2022

@author: chelsea
"""
import pandas as pd
import numpy as np
df_files= pd.read_csv('xx03.csv', sep=' ')

splt = np.array_split(df_files, 2)

small_beginning = splt[0]
small_beginning = small_beginning[small_beginning['Mean'] < 28]

to_cut_beginning = small_beginning['#Base'].tolist()

list_of_values = []

for value in to_cut_beginning:
    to_cut_beginning = value.split('-')
    list_of_values.append(to_cut_beginning[0])
    list_of_values.append(to_cut_beginning[1])  
    

list_beginning = []


for i in range(len(list_of_values)):
    list_beginning.append(int(list_of_values[i]))

if not list_beginning:
    lead = 1
else:
    min_beg = min(list_beginning)
    max_beg = max(list_beginning)
    lead = max_beg - min_beg

if lead > 20:
    lead = 20
else:
    lead = 28
    
print(lead)
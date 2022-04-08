#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:24:15 2022

@author: chelsea
"""
import pandas as pd
import numpy as np
import sys
import os

directory = '/mnt/bigdisk/'

df_T = pd.read_csv('/mnt/bigdisk/T.tsv' , header = [0,1],index_col = 0, sep='\t')
df_F = pd.read_csv('/mnt/bigdisk/F.tsv' , header = [0,1], index_col = 0, sep='\t')

df_T_index = []
df_F_index = []

new_index_F = []
new_index_T = []

new_index_tot = []

for index in df_T.index:
    to_use = index.split(' ')
    df_T_index.append(to_use[len(to_use)-1])

for index in df_F.index:
    to_use = index.split(' ')
    df_F_index.append(to_use[len(to_use)-1])
    
for i in range(len(df_T_index)):
    first_char = df_T_index[i][0]
    word = df_T_index[i].replace(first_char, '')
    index_in_F = df_F_index.index(df_F_index[0][0]+word)
    new_index_T.append((df_T_index[i], i))
    new_index_F.append((df_F_index[0][0]+word, index_in_F))
    #new_index_tot.append(df_T_index[i] + ' vs ' + df_F_index[0][0]+word)


# sorted_ = [None]*len(new_index_F)

# for i in range(len(new_index_F)):
#     index = new_index_F[i][1]
#     element = new_index_F[i][0]
#     sorted_.insert(index, element)

# sorted_ = list(filter(None, sorted_))

temp = df_F.iloc[[1]]
temp = temp.iloc[1: ,:]                

#change order of dataframe F to math T
for value in new_index_F:
    f_row = df_F.iloc[[value[1]]]
    bajs = f_row
    temp = pd.concat([temp, bajs])

temp2 = temp.iloc[[value[1]]]
temp2 = temp2.iloc[1:, :]
for i in range(len(df_T)):
    t_row = df_T.iloc[[i]]      #row to be compared to reference along with reference
    f_row = temp.iloc[[i]]

    bajs = pd.concat([t_row, f_row])
    bajs = bajs.diff()
    test = bajs.iloc[[1]]
    temp2 = pd.concat([temp2, test])
    
temp2.index = temp2.index.str.replace('N', 'T')


    

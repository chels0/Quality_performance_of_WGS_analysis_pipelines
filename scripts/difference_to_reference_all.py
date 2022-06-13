#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:25:39 2022

@author: chelsea
"""

import pandas as pd
import os

directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results (kopia)/chewbbaca_quast_tables'
list_of_files = []
for filename in os.listdir(directory):
    if 'samples.txt' in filename:
        list_of_files.append(filename)

list_ = []
list_2 = []
coverage = ['20x', '50x', '100x']
for file in list_of_files:
    with open(directory+'/'+file, "r") as f:
        seq = f.readlines()
    list_ = []
    for sentence in seq:
        seq2 = sentence.split('\t')
        list_.append(seq2)
    for cov in coverage:
        df = pd.DataFrame(list_)
        df.columns = ['Sample', 'Differences', 'Differing groups']
        df = df[df.Sample.str.contains(cov)] #keep only samples of certain coverage
        df = df[df.Sample.str.contains("PT28-1-20")==False]
        df = df[df.Sample.str.contains("PT28-3-20")==False]
        df['Differences'] = df['Differences'].astype(int)
        sum_ = df['Differences'].sum()
        print(sum_)
        list_2.append((file.split('_')[0], cov, sum_))

df2 = pd.DataFrame(list_2)
df2.columns = ['Sample', 'Coverage', 'Differences']
df3 = df2.pivot(index='Sample', columns='Coverage', values='Differences')
df3 = df3.sort_values(by='Sample')

df3.to_csv('Results/chewbbaca_quast_tables/referene_differences.csv', sep=',', encoding='utf-8') #save to csv with name of run



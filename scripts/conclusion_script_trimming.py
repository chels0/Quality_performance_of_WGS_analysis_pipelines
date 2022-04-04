#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 08:56:13 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import itertools
import pathlib
import copy

directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/'

list_of_files = []
list_of_files1 = []
list_of_files2 = []
list_ = []
list2 = []

diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'


for filename in os.listdir('/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/Comparisons/Chewbbaca_comparisons'):
    split_files = filename.split('_')
    tupls = (split_files, filename)
    list_of_files.append(tupls)
    
    list_of_files1.append(split_files[0]) #split on . to
    list_of_files2.append(split_files[2])

for element in list_of_files:
    if len(element[0][0]) == len(element[0][2]):
        list_.append((element[0][0], element[0][2]))

n_char = 'N'

characters = ['F', 'T']
# combination = list(itertools.product(n_char, characters))
# lists = [[] for _ in range(len(combination))]

# for comb in combination:
#     list2 = []
#     for element in list_:
#         if (element[0][0] == comb[0] and element[1][0] == comb[1]) or (element[0][0] == comb[1] and element[1][0] == comb[0]):
#             tuples = (element[0], element[1])
#             list2.append(tuples)
#     for lista in lists:
#         lista.append(list2)
    


template = {diff: [], corr: [], wrong: [], change: []}
list_of_things = [diff, corr, wrong, change]

characters = ['F', 'T']
for char in characters:
    list2 = []
    for element in list_:
        if (char == element[0][0] and n_char == element[1][0]) or (n_char == element[0][0] and char == element[1][0]):
            tuples = (element[0], element[1])
            list2.append(tuples)
    
    
    x20x = pd.MultiIndex.from_product([['20x'], list_of_things], names=["Coverage", ""])
    x50x = pd.MultiIndex.from_product([['50x'], list_of_things], names=["Coverage", ""])
    x100x = pd.MultiIndex.from_product([['100x'], list_of_things], names=["Coverage", ""])
    
    df_20x_template = pd.DataFrame(template, columns = x20x)
    df_50x_template = pd.DataFrame(template, columns = x50x)
    df_100x_template = pd.DataFrame(template, columns = x100x)
    
    test = []
    test1 = []
    test2 = []
    test3 = []
    
    for tuple_ in list2:
        for value in list_of_files:
            if ((tuple_[0] == value[0][0] and tuple_[1] == value[0][2]) or (tuple_[0] == value[0][2] and tuple_[1] == value[0][0])):
                df = pd.read_csv('/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/Comparisons/Chewbbaca_comparisons/'+value[1], index_col = 0, sep='\t')
                df.index = df.index.fillna('No label')
                df_20x = df[df.index.str.contains('20x|No label', regex = True)]
                df_50x = df[df.index.str.contains('50x|No label', regex = True)]
                df_100x = df[df.index.str.contains('100x|No label', regex = True)]
                dfs = []
                dfs.append((df_20x, '20x'))
                dfs.append((df_50x, '50x' ))
                dfs.append((df_100x, '100x'))
                
                for element in dfs:
                    dataframe = element[0]
                    coverage  = element[1]
                    test = []
                    test1 = []
                    test2 = []
                    for col in dataframe.columns:
                        test.append(col)
                    for i in range(len(test)):
                        if i % 2 == 0:
                            test1.append(test[i])
                        else:
                            test2.append(test[i])
                    count_right = 0
                    count_wrong = 0
                    count_change = 0
                    count_diffs = 0
                    test3 = []
                                        
                    
                    if n_char in dataframe.iloc[0,0]:
                        for i in range(len(test1)):
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test1[i]] == '0') & (dataframe[test2[i]] != '0'), dataframe[test1[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
                            df_copy = dataframe.copy()
                            df_copy[wrong]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] == '0'), dataframe[test1[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy[wrong].dropna())
                            df_copy = dataframe.copy()
                            df_copy[change]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] != '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_change = count_change + len(df_copy[change].dropna())
                            df_copy = dataframe.copy()
                            df_copy[diff]  = np.where((dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy[diff].dropna())
                        
                    elif n_char in dataframe.iloc[0,1]:
                        for i in range(len(test1)):
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test2[i]] == '0') & (dataframe[test1[i]] != '0'), dataframe[test2[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
                            df_copy = dataframe.copy()
                            df_copy[wrong]  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0'), dataframe[test2[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy[wrong].dropna())
                            df_copy = dataframe.copy()
                            df_copy[change]  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test2[i]], np.nan)            
                            count_change = count_change + len(df_copy[change].dropna())
                            df_copy = dataframe.copy()
                            df_copy[diff]  = np.where((dataframe[test2[i]] != dataframe[test1[i]]), dataframe[test2[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy[diff].dropna())   
                
                    if n_char in tuple_[0]:
                        n_char_name = tuple_[0]
                        not_n_char_name = tuple_[1]
                    else:
                        n_char_name = tuple_[1]
                        not_n_char_name = tuple_[0]

                    right_result = {(coverage,diff) : count_diffs, (coverage,corr) : count_right, (coverage, wrong) : count_wrong, (coverage, change) : count_change}                        
                    df_test = pd.DataFrame(right_result, index = [n_char_name+ ' vs ' + not_n_char_name])

                        #right_result = {diff : count_diffs, corr : count_right, wrong : count_wrong, change : count_change}
                    if coverage == '20x':
                        df_20x_template= pd.concat([df_20x_template, df_test])
                    elif coverage == '50x':
                        df_50x_template = pd.concat([df_50x_template, df_test])
                    else:
                        df_100x_template = pd.concat([df_100x_template, df_test])
    
    final_result = df_20x_template.merge(df_50x_template, left_index = True, right_index = True)
    final_result = final_result.merge(df_100x_template, left_index = True, right_index = True)
    final_result = final_result.astype(int)
    final_result.to_csv('/mnt/bigdisk/'+char+'.tsv', sep='\t', encoding='utf-8')

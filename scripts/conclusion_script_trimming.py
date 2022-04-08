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


directory = '/Results/Comparisons/Chewbbaca_comparisons


list_of_files = []
list_of_files1 = []
list_of_files2 = []
list_ = []
list2 = []
word_count = []
word_count2 = []

diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'


for filename in os.listdir(directory):
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
    word_list = []
    word_list2 = []
    
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
                    word_list = []
                    word_list2 = []
                    word_count = []
                    word_count2 = []
                                        
                    print(coverage)
                    if n_char in dataframe.iloc[0,1]:
                        for i in range(len(test1)):
                            word_list = []
                            df_copy = dataframe.copy()
                            df_copy[wrong]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] == '0'), dataframe[test1[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy[wrong].dropna())
                            no_top_row = df_copy.iloc[1:, :]
                            
                            mask = (no_top_row[wrong].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False))
                            new_mask = mask.dropna(axis=0)
                            if len(new_mask) != 0:
                                word_list.append(new_mask[0])
                            
                            uniques = list(set(word_list))
                            for words in uniques:
                                counter = word_list.count(words)
                                word_count.append((words, counter))
                            
                            
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test1[i]] == '0') & (dataframe[test2[i]] != '0'), dataframe[test1[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
                            df_copy = dataframe.copy()
                            df_copy[change]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] != '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_change = count_change + len(df_copy[change].dropna())
                            df_copy = dataframe.copy()
                            df_copy[diff]  = np.where((dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy[diff].dropna())
                        
                    elif n_char in dataframe.iloc[0,0]:
                        for i in range(len(test1)):
                            
                            word_list2 = []
                            df_copy = dataframe.copy()
                            df_copy[wrong]  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0'), dataframe[test2[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy[wrong].dropna())
                            no_top_row = df_copy.iloc[1:, :]
                            
                            mask = (no_top_row[wrong].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False))
                            new_mask = mask.dropna(axis=0)
                            
                            if len(new_mask) != 0:
                                word_list2.append(new_mask[0])
                            
                            uniques2 = list(set(word_list2))
                            for words in uniques2:
                                counter = word_list2.count(words)
                                word_count2.append((words, counter))                          
                            
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test2[i]] == '0') & (dataframe[test1[i]] != '0'), dataframe[test2[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
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
                    
                    uniques = list(set(word_count))
                    uniques2 = list(set(word_count2))
                    uniques3 = list(set(uniques + uniques2))
                    
                    if word_count:
                        for unique_value in uniques:
                            df_test[coverage, '% ' + unique_value[0]] = 0
                            
                        for words in word_count:
                            df_test[coverage, '% ' + words[0]] = words[1] + (df_test[coverage, '% ' + words[0]][0])
                        

                    
                    if word_count2 and not word_count:
                        for unique_value in uniques2:
                            df_test[coverage, '% ' + unique_value[0]] = 0
                        for words2 in word_count2:
                            df_test[coverage, '% ' + words2[0]] = words2[1] + (df_test[coverage, '% ' + words2[0]][0])
          
                            
                    if word_count2 and word_count:
                        for words in word_count2:
                            df_test[coverage, '% ' + words[0]] = (words[1] + df_test[coverage, words[0]][0])

                    
                    for i in range(len(df_test)):
                        if df_test[coverage, wrong][i] != 0:
                            if not word_count2 and not word_count:
                                df_test[coverage, '% Wrong-called'] = count_wrong
                    
                    column = int(df_test[coverage, wrong][0])
                    to_subtract_list = []
                    for value in uniques3:

                        to_subtract = int(df_test[coverage, '% ' + value[0]][0])
                        to_subtract_list.append(to_subtract)
                        sum_ = sum(to_subtract_list)
                        wrong_called = column - sum_
                        
                        df_test[coverage, '% Wrong-called'] = wrong_called
                        

                    if coverage == '20x':
                        df_20x_template= pd.concat([df_20x_template, df_test])
                    elif coverage == '50x':
                        df_50x_template = pd.concat([df_50x_template, df_test])
                    else:
                        df_100x_template = pd.concat([df_100x_template, df_test])
    
    final_result = df_20x_template.merge(df_50x_template, left_index = True, right_index = True)
    final_result = final_result.merge(df_100x_template, left_index = True, right_index = True)
    final_result = final_result.replace('nan', 0)
    final_result = final_result.fillna(0)
    final_result = final_result.astype(int)
    final_result.to_csv('Results/Conclusions/'+n_char+'+'char+'_results.tsv', sep='\t', encoding='utf-8')
    
    

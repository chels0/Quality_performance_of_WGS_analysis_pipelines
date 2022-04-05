#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:41:35 2022

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
#append filenames into list_of_files

for filename in os.listdir('/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/Comparisons/Chewbbaca_comparisons'):
    split_files = filename.split('_')
    tupls = (split_files, filename)
    list_of_files.append(tupls)
    list_of_files1.append(split_files[0]) #split on . to
    list_of_files2.append(split_files[2])
    
all_lists = list_of_files1 + list_of_files2

files = []
trimmomatic_files = []
no_trim_files = []
fastp_files = []
filter_files = []
spades_files = []
skesa_files = []

#characters = ['f', 'T', 'N', 'P', 'F']
df_coverages = []

characters = ['P', 'f']
for char in characters:
    for filename in all_lists:
        if char in filename:
            files.append(filename)

    no_dup = list(set(files))
    
    not_char = []
    
    for value in no_dup:
        hej = value.replace(char, '')
        tupls = (value, hej)
        not_char.append(tupls)
    
    template = {'Differences': [], 'Correct': [], 'Wrong': [], 'Change': []}
    list_of_things = ['Differences', 'Correct', 'Wrong', 'Change']
    # upper_level_20 = ['20x', '50x', '100x']

    x20x = pd.MultiIndex.from_product([['20x'], list_of_things], names=["Coverage", ""])
    x50x = pd.MultiIndex.from_product([['50x'], list_of_things], names=["Coverage", ""])
    x100x = pd.MultiIndex.from_product([['100x'], list_of_things], names=["Coverage", ""])

    df_20x_template = pd.DataFrame(template, columns = x20x)
    df_50x_template = pd.DataFrame(template, columns = x50x)
    df_100x_template = pd.DataFrame(template, columns = x100x)
    
    
    # upper_level = ['20x', '50x', '100x']
    # hej = pd.MultiIndex.from_product([upper_level, list_of_things], names=["Coverage", ""])
    # df_template = pd.DataFrame(template, index = ['Run'], columns = hej)
    # test_dic = {'Diff':[]}
    # test_upper = {'20x', '50x'}
    # test_list = ['Differences']
    # test_hej = pd.MultiIndex.from_product([test_upper, test_list], names=["first", "second"])
    # test_df = pd.DataFrame(test_dic, columns = test_hej)
    # test_result = {('20x','Differences'):5}
    # test_df = test_df.append(test_result, ignore_index=True)
    
    
    test = []
    test1 = []
    test2 = []
    test3 = []
    for tuple_ in not_char:
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
                    if char in dataframe.iloc[0,0]:
                        for i in range(len(test1)):
                            df_copy = dataframe.copy()
                            df_copy['correct']  = np.where((dataframe[test1[i]] == '0') & (dataframe[test2[i]] != '0'), dataframe[test1[i]], np.nan)            
                            count_right = count_right + len(df_copy['correct'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['wrong']  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] == '0'), dataframe[test1[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy['wrong'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['change']  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] != '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_change = count_change + len(df_copy['change'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['diffs']  = np.where((dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy['diffs'].dropna())
                    
                    if char in dataframe.iloc[0,1]:
                        for i in range(len(test1)):
                            df_copy = dataframe.copy()
                            df_copy['correct']  = np.where((dataframe[test2[i]] == '0') & (dataframe[test1[i]] != '0'), dataframe[test2[i]], np.nan)            
                            count_right = count_right + len(df_copy['correct'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['wrong']  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0'), dataframe[test2[i]], np.nan)            
                            count_wrong = count_wrong + len(df_copy['wrong'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['change']  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test2[i]], np.nan)            
                            count_change = count_change + len(df_copy['change'].dropna())
                            df_copy = dataframe.copy()
                            df_copy['diffs']  = np.where((dataframe[test2[i]] != dataframe[test1[i]]), dataframe[test2[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy['diffs'].dropna())   
                        
                    right_result = {(coverage,'Differences') : count_diffs, (coverage,'Correct') : count_right, (coverage, 'Wrong') : count_wrong, (coverage, 'Change') : count_change}                        
                    df_test = pd.DataFrame(right_result, index = [tuple_[0]+ ' vs ' + tuple_[1]])

                        #right_result = {'Differences' : count_diffs, 'Correct' : count_right, 'Wrong' : count_wrong, 'Change' : count_change}
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


characters = ['T', 'F', 'N']
for char in characters:
    for filename in all_lists:
        if char in filename:
            files.append(filename)

    no_dup = list(set(files))
    
    not_char = []
    
    for value in no_dup:
        hej = value.replace(char, '')
        tupls = (value, hej)
        not_char.append(tupls)

    
    #([df_20x_template, df_50x_template, df_100x_template], inde)            
    
                #df_new = df_template.copy().dropna()
                    #df_new.set_index('Run', inplace = True)
                    
                    #df = df_new.set_index('Coverage')
                    
                #df_new.to_csv('/mnt/bigdisk/'+char+'.tsv', sep='\t', encoding='utf-8')

    # if char == 'P':
    #     df_pilon = df_template.copy()
    #     df_pilon.set_index('Run', inplace =True)
        
        
        #elif (tuple_[0] == value[0][2] and tuple_[1] == value[0][0]):
        #    df = pd.read_csv('/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/Comparisons/Chewbbaca_comparisons/'+value[1], index_col = 0, sep='\t')

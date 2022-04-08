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

directory = '/Results/Comparisons/Chewbbaca_comparisons
pathlib.Path('Results/Conclusions/').mkdir(parents=True, exist_ok=True)
diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'


list_of_files = []
list_of_files1 = []
list_of_files2 = []
#append filenames into list_of_files

for filename in os.listdir(directory):
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
word_count = []
word_count2 = []

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
    
    template = {diff: [], corr: [], wrong: [], change: []}
    list_of_things = [diff, corr, wrong, change]
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
    # test_list = [diff]
    # test_hej = pd.MultiIndex.from_product([test_upper, test_list], names=["first", "second"])
    # test_df = pd.DataFrame(test_dic, columns = test_hej)
    # test_result = {('20x',diff):5}
    # test_df = test_df.append(test_result, ignore_index=True)
    
    
    test = []
    test1 = []
    test2 = []
    test3 = []
    word_list = []
    word_list2 = []
    
    for tuple_ in not_char:
        for value in list_of_files:
            if ((tuple_[0] == value[0][0] and tuple_[1] == value[0][2]) or (tuple_[0] == value[0][2] and tuple_[1] == value[0][0])):
                df = pd.read_csv(directory+'/'+value[1], index_col = 0, sep='\t')
                df.index = df.index.fillna('No label')
                df_20x = df[df.index.str.contains('20x|No label', regex = True)]
                df_50x = df[df.index.str.contains('50x|No label', regex = True)]
                df_100x = df[df.index.str.contains('100x|No label', regex = True)]
                dfs = []
                dfs.append((df_20x, '20x'))
                dfs.append((df_50x, '50x' ))
                dfs.append((df_100x, '100x'))
                errors = []
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
                    errors = []
                    word_list = []
                    word_list2 = []
                    word_count = []
                    word_count2 = []
                    
                    if char in dataframe.iloc[0,0]:
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
                            df_copy[diff]  = np.where((dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_diffs = count_diffs + len(df_copy[diff].dropna())
                            
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test1[i]] == '0') & (dataframe[test2[i]] != '0'), dataframe[test1[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
                    #         df_copy = dataframe.copy()
                    #         df_copy[wrong]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] == '0'), dataframe[test1[i]], np.nan)
                            
                    #         no_top_row = df_copy.iloc[1: , :]
                    #         mask = ([no_top_row[col].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False) for col in test1])
                    #         for value2 in mask:
                    #             mask2 = value2.dropna().drop_duplicates()
                    #             for col in mask2:
                    #                 errors.append(col) #append letters from each column
                    #         no_dup_errors = list(set(errors))
                    #         word_count = []
    
                    #         if no_dup_errors:
                    #             for value in no_dup_errors:
                    #                 nr_of_words = df_copy.isin([value]).sum(axis=0).sum(axis=0)
                    #                 word_count.append((value, nr_of_words))
                                
                    #         count_wrong = count_wrong + len(df_copy[wrong].dropna())
                            df_copy = dataframe.copy()
                            df_copy[change]  = np.where((dataframe[test1[i]] != '0') & (dataframe[test2[i]] != '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test1[i]], np.nan)            
                            count_change = count_change + len(df_copy[change].dropna())
                            
                    
                    if char in dataframe.iloc[0,1]:
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
                            df_copy[diff]  = np.where((dataframe[test2[i]] != dataframe[test1[i]]), dataframe[test2[i]], np.nan)            
                            skesa_files.append(df_copy[diff])
                            count_diffs = count_diffs + len(df_copy[diff].dropna())
                            
                            df_copy = dataframe.copy()
                            df_copy[corr]  = np.where((dataframe[test2[i]] == '0') & (dataframe[test1[i]] != '0'), dataframe[test2[i]], np.nan)            
                            count_right = count_right + len(df_copy[corr].dropna())
                            
                            df_copy = dataframe.copy()
                            df_copy[change]  = np.where((dataframe[test2[i]] != '0') & (dataframe[test1[i]] == '0') & (dataframe[test1[i]] != dataframe[test2[i]]), dataframe[test2[i]], np.nan)            
                            count_change = count_change + len(df_copy[change].dropna())
                               
                            
                    pilon = tuple_[0]
                    not_pilon = tuple_[1]
                    right_result = {(coverage,diff) : count_diffs, (coverage,corr) : count_right, (coverage, wrong) : count_wrong, (coverage, change) : count_change}                        
                    df_test = pd.DataFrame(right_result, index = [pilon+ ' vs ' + not_pilon])
                    
                    uniques = list(set(word_count))
                    uniques2 = list(set(word_count2))
                    uniques3 = list(set(uniques + uniques2))
                    print(word_count2)
                    if word_count:
                        for unique_value in uniques:
                            df_test[coverage, '% ' + unique_value[0]] = 0
                            
                        for words in word_count:
                            df_test[coverage, '% ' + words[0]] = words[1] + (df_test[coverage, '% ' + words[0]][0])
                        
                        #for unique_value in uniques:
                         #   df_test[coverage, '% Wrong-called'] = (df_test[coverage, wrong][0] - df_test[coverage, '% ' + unique_value])
                    
                    if word_count2 and not word_count:
                        for unique_value in uniques2:
                            df_test[coverage, '% ' + unique_value[0]] = 0
                        for words2 in word_count2:
                            df_test[coverage, '% ' + words2[0]] = words2[1] + (df_test[coverage, '% ' + words2[0]][0])
                            print(words2)
                        #df_test[coverage, '% Wrong-called'] = (count_wrong - (words2[1] + df_test[coverage, '% ' + words2[0]][0]))
                        #for unique_value in uniques2:
                         #   df_test[coverage, '% Wrong-called'] = (df_test[coverage, wrong][0] - df_test[coverage, '% ' + unique_value[0]])                
                            
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
                        df_20x_template= pd.concat([df_20x_template, df_test], sort = False)
                    elif coverage == '50x':
                        df_50x_template = pd.concat([df_50x_template, df_test], sort = False)
                    else:
                        df_100x_template = pd.concat([df_100x_template, df_test], sort = False)
                        
    
    final_result = df_20x_template.merge(df_50x_template, left_index = True, right_index = True)
    final_result = final_result.merge(df_100x_template, left_index = True, right_index = True)
    final_result = final_result.replace('nan', 0)
    final_result = final_result.fillna(0)
    final_result = final_result.astype(int)
    final_result.to_csv('Results/Conclusions/'+char+'results.tsv', sep='\t', encoding='utf-8')

#append letters from each column
                            # print(mask2)
                            # errors2 = errors
                            # no_dup_errors = list(set(errors))
                            # word_count2 = []
                            
                            
                            
                            
                            
                            # no_top_row = df_copy.iloc[1: , :-1]
                            # mask = ([no_top_row[col].str.extract(('((?!^\d+$)+(?!^\W)^.+$)'), expand=False) for col in test2])
                            # for value2 in mask:
                            #     mask2 = value2.dropna().drop_duplicates()
                            #     for col in mask2:
                            #         test = col
                            #         errors.append(test) #append letters from each column
                            # errors2 = errors
                            # no_dup_errors = list(set(errors))
                            # word_count2 = []
                            
                            # if no_dup_errors:
                            #     for value in no_dup_errors:
                            #         nr_of_words = df_copy.isin([value]).sum(axis=0).sum(axis=0)
                            #         word_count2.append((value, nr_of_words))

    
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

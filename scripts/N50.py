#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 14:45:05 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns


#directory = 'Results/chewbbaca_quast_tables/placeholder/'
directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca_quast_tables/placeholder_test/'

list_N50 = []

columns = ['# contigs', 'Largest contig', 'Total length', 
               'Genome fraction (%)', 'GC (%)', 'N50', 'NG50', '# misassemblies']

#columns = ['GC (%)']

for col in columns:
    coverages = ['20x', '50x', '100x']
    empty_df = pd.DataFrame({})
    for cov in coverages:
        empty_df = pd.DataFrame({})
        for filename in os.listdir(directory):
           df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
           
           df = df[[col]]
           df = df.iloc[1: , :]
           mean = df[df.index.str.contains(cov)].mean()
           index = filename.split('_')[0]
           mean.index = [index]
           mean = mean.rename(col)
           df_mean = pd.DataFrame(mean)
           empty_df = pd.concat([df_mean, empty_df]).sort_index()
          

        test = empty_df[col].values - empty_df[col].values[:, None]
        test_df = pd.DataFrame(test, columns = empty_df.index, index = empty_df.index)
        test_df.index.names = [cov]
        

        if col != 'GC (%)' and col != 'Genome fraction (%)':
            test_df = test_df.applymap(int)
        else:
            test_df = test_df.round(2)
        #test_df.to_csv('Results/Conclusions/'+ col+'_'+cov+'_results.tsv', sep='\t', encoding='utf-8') #save to csv with name of run
        
        sns.pairplot(test_df, kind="scatter")


        test_df.plot(kind='scatter', title='20x')
        plt.legend(title='Labels', loc='upper right')
        plt.savefig('Results/Conclusions/'+n_char+char+'_plot_all.png', bbox_inches='tight')

# test = []
# dict_ = { }
# df2 = pd.DataFrame(data=dict_)

# for cov in coverages:
#     for value in list_N50:
#         if value[2] == cov:
#             filename = value[0]
#             mean = value[1]
#             dict2 = {filename : mean}
#             df3 = pd.DataFrame(data=dict2, index = [filename] )
#             df2 = pd.concat([df2, df3])
#     test.append(df2)
#     df2 = pd.DataFrame(data={})
        
        




    
    
        
        
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:28:44 2022

@author: chelsea
"""

import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns


#Script for plotting QUAST metrics, one boxplot for each metric at different coverages
#_____________________________________________________________________________________

directory = '/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca_quast_tables/placeholder/'

#Define which assemblers have been used
assembler_1 = 'SpC'
assembler_2 = 'SpI'
#Define what name should be plotted
assembler_1_full = 'SPAdes --careful'
assembler_2_full = 'SPAdes --isolate'


only_trimming = []
list_of_trimming_column_vales = []
list_of_software_combinations = []
list_of_files = []
list_of_column_values_lists = []
cov_list = []

#columns to plot
columns = ['# contigs', 'Largest contig', 'Total length', 
               'Genome fraction (%)', 'GC (%)', 'N50', 'NG50', '# misassemblies']


columns = ['# contigs', 'Genome fraction (%)', 'N50', '# misassemblies']

#coverages
coverages = ['20x', '50x', '100x']
empty_df = pd.DataFrame({}) #create empty dataframe

#add filenames to list
for filename in os.listdir(directory):
    list_of_files.append(filename)
list_of_files.sort() #sort list

count=0 #counter
#plot all quast metrics (columns) for each coverage
list_100 = []
list_50 = []
list_20 = []



for col in columns:
    median = []
    count = count + 1 #increase counter

    for cov in coverages:
        #empty lists
        only_trimming = []
        list_of_trimming_column_vales = []
        list_of_software_combinations = []
        list_of_column_values_lists = []
        #extract data from files
        for filename in list_of_files:
            df = pd.read_csv(directory+filename, sep='\t', index_col=0) #create tab separated dataframe
            df = df[columns] #keep only relevant columns
            df = df.iloc[1: , :] #remove reference row
            df = df[df.index.str.contains(cov)] #keep only samples of certain coverage
            median.append((filename.split('_')[0], cov, col, df[col].median()))

    df2 = pd.DataFrame(median)
    df2.columns = ['Combination', 'Coverage', 'QUAST metric', 'median']
    
    df2['Combination'] = df2['Combination'].str.replace('F','')
    df2['Combination'] = df2['Combination'].str.replace('N','')
    df2['Combination'] = df2['Combination'].str.replace('T','')
    df2['Combination'] = df2['Combination'].str.replace('P','')
    df2['Combination'] = df2['Combination'].str.replace('200f','')
    df2['Combination'] = df2['Combination'].str.replace('500f','')
    df2['Combination'] = df2['Combination'].str.replace('Ske','SKESA')
    df2['Combination'] = df2['Combination'].str.replace('SpI','SPAdes --isolate')
    df2['Combination'] = df2['Combination'].str.replace('SpC','SPAdes --careful')


    #df2['Coverage'] = df2['Coverage'].str.replace('x','')
    #df2['Coverage'] = df2['Coverage'].astype(int)
    if col == '# contigs':
        fig_nr = 'A'
    if col == 'N50':
        fig_nr = 'B'
    if col == 'Genome fraction (%)':
        fig_nr = 'C'
    if col == '# misassemblies':
        fig_nr = 'D'
    
    print(count)
    
    if count == 4:
        plt.figure()
        sns.set_theme()
        sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
        ax = sns.lineplot(data=df2, x='Coverage', y='median', hue='Combination')
        ax.set_title(fig_nr+': '+col, fontsize = 15)
        plt.ylabel("Median "+col, labelpad=20, fontsize=15)
        plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
    
    if count < 4:
        plt.figure()
        sns.set_theme()
        sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
        ax = sns.lineplot(data=df2, x='Coverage', y='median', hue='Combination', legend=False )
        ax.set_title(fig_nr+': '+col, fontsize = 15)
        plt.ylabel("Median "+col, labelpad=20, fontsize=15)
        plt.tight_layout()
    
    plt.savefig('Results/Conclusions/coverage_plots/'+col+'_'+cov+'_.png', bbox_inches="tight")

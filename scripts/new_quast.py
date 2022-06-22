#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:46:23 2022

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

directory = 'Results/chewbbaca_quast_tables/placeholder/'

#Define which assemblers have been used
assembler_1 = sys.argv[1]
assembler_2 = 'SpI'
#Define what name should be plotted

if assembler_1 == 'Ske' or assemblier_1 == 'Skesa' or assembler_1 == 'skesa':
    assembler_1_full = 'Skesa'
if assembler_1 == 'SpC' or assembler_1 == 'SPAdes --careful' or assembler_1 == '--careful':
    assembler_1_full = 'SPAdes --careful'
assembler_2_full = 'SPAdes --isolate'


only_trimming = []
list_of_trimming_column_vales = []
list_of_software_combinations = []
list_of_files = []
list_of_column_values_lists = []

#columns to plot
columns = ['# contigs', 'Largest contig', 'Total length', 
               'Genome fraction (%)', 'GC (%)', 'N50', 'NG50', '# misassemblies']


columns = ['# contigs', 'Largest contig', 
               'Genome fraction (%)', 'N50', '# misassemblies']

#coverages
coverages = ['20x', '50x', '100x']
empty_df = pd.DataFrame({}) #create empty dataframe

#add filenames to list
for filename in os.listdir(directory):
    list_of_files.append(filename)
list_of_files.sort() #sort list

count=0 #counter
#plot all quast metrics (columns) for each coverage
for col in columns:
    for cov in coverages:
        count = count + 1 #increase counter
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
            df_assembler_1olumn_values_list = df[col].to_list() #turn dataframe column into list
            list_of_software_combinations.append(filename.split('_')[0]) #append software combination in list
            list_of_column_values_lists.append(df_assembler_1olumn_values_list) #append dataframe in list form to list
            
            if 'P' not in filename and 'f' not in filename:
                only_trimming.append(filename.split('_')[0]) #append software combinations without post assembly improvements
                list_of_trimming_column_vales.append(df_assembler_1olumn_values_list) #append column values of files with no post assembly improvements
                    
        #transform dataframe for plotting on the bassis of assembler   
        df_2 = pd.DataFrame(list_of_trimming_column_vales) #create dataframe
        df_2 = df_2.transpose() #transpose
        df_2.columns = only_trimming #set software combinations without psot assembly improvements as columns 
        df_assembler_1 = df_2[['N'+assembler_1, 'F'+assembler_1, 'T'+assembler_1]] #create dataframe with only values from first assembler
        df_assembler_2 = df_2[['N'+assembler_2, 'F'+assembler_2, 'T'+assembler_2 ]] # create dataframe with only values from second assembler
        
        df_assembler_1_melt = df_assembler_1.melt().assign(assembler=assembler_1) #unpivot dataset for plotting
        df_assembler_2_melt = df_assembler_2.melt().assign(assembler=assembler_2) #unpivot dataset for plotting
        
        concated = pd.concat([df_assembler_1_melt, df_assembler_2_melt]) #combine dataframes
        
        #Replace denotations to full names for plotting
        concated['variable'] = concated['variable'].str.replace('SpI','')
        concated['variable'] = concated['variable'].str.replace('SpC','')
        concated['variable'] = concated['variable'].str.replace('Ske','')
        concated['variable'] = concated['variable'].str.replace('N','No trim '+cov)
        concated['variable'] = concated['variable'].str.replace('F','Fastp '+cov)
        concated['variable'] = concated['variable'].str.replace('T','Trimmomatic '+cov)
        concated = concated.rename(columns = {'variable': 'Trimming software', 'value': col })
        
        #set figure labels        
        if col == '# contigs':
            fig_nr = 'A'
        if col == 'Largest contig':
            fig_nr = 'B'
        if col == 'N50':
            fig_nr = 'C'
        if col == 'Genome fraction (%)':
            fig_nr = 'D'
        if col == '# misassemblies':
            fig_nr = 'E'
            
        if count == 3: #if statement to print legend on only the 100x figure
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.violinplot(data=concated, x='Trimming software', y=col, hue='assembler', split = True, cut = 0, palette="vlag", height=8, aspect=15/8, inner='quartiles')
            ax = sns.swarmplot(data=concated, x='Trimming software', y=col, hue='assembler', dodge = True, s=2.5)
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.tight_layout()
            plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
        
        else:
            plt.figure()
            sns.set_theme()
            sns.set(rc={"figure.figsize":(7, 5), "figure.dpi":300, 'savefig.dpi':300})
            ax = sns.violinplot(data=concated, x='Trimming software', y=col, hue='assembler', split = True, cut = 0, palette="vlag", height=8, aspect=15/8, inner='quartiles').set(title=cov)
            plt.legend([],[], frameon=False)
            ax = sns.swarmplot(data=concated, x='Trimming software', y=col, hue='assembler', dodge = True, s=2.5)
            ax.set_title(fig_nr+cov+': '+col+' at coverage '+cov, fontsize = 15)
            plt.legend([],[], frameon=False)
            plt.tight_layout()

        plt.savefig('Results/Conclusions/quast_plots_seaborn/'+col+'_'+cov+'_.png', bbox_inches="tight")

    
        # df_2 = pd.DataFrame(list_of_column_values_lists)
        # df_2 = df_2.transpose()
        # df_2.columns = list_of_software_combinations
        # #df_assembler_1 = df_2[['NSpC', 'FSpC', 'TSpC' ]]
        # skesa = []
        # spades = []
        # for element in list_of_software_combinations:
        #     if 'SpI' in element:
        #         spades.append(element)
        #     if 'Ske' in element:
        #         skesa.append(element)
        # df_assembler_2 = df_2[spades]
        # df_s = df_2[skesa]
        
        # #df_assembler_1_melt = df_assembler_1.melt().assign(assembler='SPAdes --careful')
        # df_assembler_2_melt = df_assembler_2.melt().assign(assembler='SPAdes --isolate')
        # df_s_melt = df_s.melt().assign(assembler='SKESA')
        
        # concated = pd.concat([df_s_melt, df_assembler_2_melt])
        # concated['variable'] = concated['variable'].str.replace('SpI','')
        # concated['variable'] = concated['variable'].str.replace('SpC','')
        # concated['variable'] = concated['variable'].str.replace('Ske','')
        # concated = concated.rename(columns = {'variable': 'Trimming software', 'value': col })
        
        # if count == 3:
        #     plt.figure()
        #     sns.set_theme()
        #     sns.set(rc={"figure.figsize":(10, 5), "figure.dpi":300, 'savefig.dpi':300})
        #     ax = sns.violinplot(data=concated, x='Trimming software', y=col, hue='assembler', split = True, cut = 0, palette="vlag", height=8, aspect=15/8, inner='quartiles').set(title=cov)
        #     ax = sns.swarmplot(data=concated, x='Trimming software', y=col, hue='assembler', dodge = True, s=2.5)
        #     plt.tight_layout()
        #     plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
        
        # else:
        #     plt.figure()
        #     sns.set_theme()
        #     sns.set(rc={"figure.figsize":(10, 5), "figure.dpi":300, 'savefig.dpi':300})
        #     ax = sns.violinplot(data=concated, x='Trimming software', y=col, hue='assembler', split = True, cut = 0, palette="vlag", height=8, aspect=15/8, inner='quartiles').set(title=cov)
        #     plt.legend([],[], frameon=False)
        #     ax = sns.swarmplot(data=concated, x='Trimming software', y=col, hue='assembler', dodge = True, s=2.5)
        #     plt.legend([],[], frameon=False)
        #     plt.tight_layout()

        # plt.savefig('Results/Conclusions/quast_plots_seaborn/'+col+'_'+cov+'_FULL.png', bbox_inches="tight")

        
        # data = [list_of_trimming_column_vales[0], list_of_trimming_column_vales[2], list_of_trimming_column_vales[4], list_of_trimming_column_vales[1], list_of_trimming_column_vales[3], list_of_trimming_column_vales[5]]
        # labels = [only_trimming[0], only_trimming[2], only_trimming[4], only_trimming[1], only_trimming[3], only_trimming[5]]
        # labels = ['Fastp', 'Fastp', 'No trim', 'No trim', 'Trimmomatic', 'Trimmomatic']
        # hu = ['--careful', '--isolate', '--careful', '--isolate', '--careful', '--isolate']
        # #data = list_of_column_values_lists
        # #labels = list_of_software_combinations
        # length = len(labels)
        # values = [*range(1,length+1,1)]
        # #figure = plt.figure(figsize=(10,5))
        # sns.set_theme()
        
        # ax = sns.violinplot(data=df_2, inner=None, color=".8", cut=0, split = True)

        # #ax = sns.boxplot(data=df_2)
        # #ax = sns.swarmplot(data=df_2)
        # ax = sns.stripplot(data=df_2)
        
        # B = plt.violinplot(data)
        # plt.xticks(values, labels = labels)
        # plt.annotate('SPAdes --careful', (0,0), (100, -25), xycoords='axes fraction', textcoords='offset points', va='top')
        # plt.annotate('SPAdes --isolate', (0,0), (380, -25), xycoords='axes fraction', textcoords='offset points', va='top')
        # plt.title('Quast parameter: ' + col + ' at coverage ' +cov)
        # plt.xlabel('Run')
        # plt.ylabel(col)
        # basj = get_box_plot_data(labels, B)
        # basj.to_csv('Results/Conclusions/quast_plots/'+col+'_'+cov+'_.csv', sep='\t', encoding='utf-8')
        #plt.savefig('Results/Conclusions/quast_plots_seaborn/'+col+'_'+cov+'_.png')




       
       

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:11:34 2022

@author: chelsea
"""
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns

directory = 'Results/Conclusions'

#Define column names
diff = '# Differences'
corr = '# Corrections'
wrong = '# Errors'
change = '# Changes'

list_ = ['Trimming', 'Pilon', 'Filtering 200', 'Pilon and Filtering'] #list of plots to be created
coverages = ['20x', '50x', '100x']

list_index = sys.argv[1]
to_plot = list_[int(list_index)] #choose which plot to plot
#to_plot = list_[3]
if to_plot != list_[0]:
    x_axis_label = 'Compared pipelines'
else:
    x_axis_label = 'The trimming software compared to no trimming'  


df_NT = pd.read_csv(directory+'/'+'N+T_results.tsv', index_col = 0, sep='\t', header=[0,1]) #create dataframe of no trimming vs trimming chewbbaca results
df_NF = pd.read_csv(directory+'/'+'N+F_results.tsv', index_col = 0, sep='\t', header=[0,1]) #create datafrmae of no trimming vs Fastp
df_P = pd.read_csv(directory+'/Presults.tsv', index_col= 0, sep='\t', header=[0,1]) #create dataframe of pilon vs no pilon
df_f = pd.read_csv(directory+'/500fresults.tsv', index_col= 0, sep='\t', header=[0,1]) #create dataframe of filtering vs no filterin

to_compare_to_isolate = sys.argv[2] #assembler to compare to --isolate
if to_compare_to_isolate == 'Ske':
    assembler_not_isolate = 'SKESA'
if to_compare_to_isolate == 'SpC':
    assembler_not_isolate = 'SPAdes --careful'

to_remove_from_trimming= ['P', '500f', '200f'] #
to_remove_P = ['500f', '200f']
to_remove_f = ['P']

df_fP = pd.concat([df_P, df_f]) #create a filtering pilon dataframe
#df_fP = df_fP[df_fP.index.str.contains('500fP')]

df_fP = df_fP[df_fP.index.str.contains('500fP')] #extract only Fp
df_fP = df_fP.astype(int)
concated = pd.concat([df_NT, df_NF]) #concatenate the trimming dataframes

#keep only columns without string errors
concated = concated[[('20x', diff), ('20x', corr), ('20x', wrong), ('20x', change), ('50x', diff), ('50x', corr), ('50x', wrong), ('50x', change), ('100x', diff), ('100x', corr), ('100x', wrong), ('100x', change)]] 
pilon_correct = df_P[[('20x', diff), ('20x', corr), ('20x', wrong), ('20x', change), ('50x', diff), ('50x', corr), ('50x', wrong), ('50x', change), ('100x', diff), ('100x', corr), ('100x', wrong), ('100x', change)]]
 
if to_plot != list_[0]: #if dataframe to plot is not trimming
    for cov in coverages:
        if to_plot == list_[1]: #if to plot is pilon
            for val in to_remove_P:
                pilon_correct = pilon_correct[~pilon_correct.index.str.contains(val)] #remove filtering from pilon dataframe
            concated = pilon_correct
        
        if to_plot == list_[2]: #if to plot is filtering
            for val in to_remove_f:
                df_f = df_f[~df_f.index.str.contains(val)] #remove pilon
            concated = df_f
            
        
        if to_plot == list_[3]: #if to plot is both pilon and filtering
            concated = df_fP

        #calculate percenetages
        concated = concated.groupby(concated.index).sum() #group by index
        # concated[(cov, '% Corrections')] = (concated[(cov, corr)])
        # concated[(cov, '% Errors')] = (concated[(cov, wrong)])
        # concated[(cov, '% Changes')] = (concated[(cov, change)])
        concated = concated.round(0) #remove decimals   
        
        if to_plot !=  list_[3]: #if to plot is not combined filtering and pilon
            
            #put each assembler combination in a separate dataframe
            assembler_2 = concated[concated.index.str.contains(to_compare_to_isolate)]
            assembler_2[(cov, 'assembler')] = [assembler_not_isolate, assembler_not_isolate, assembler_not_isolate,] #add assembler as column
            concated_isolate = concated[concated.index.str.contains('I')]
            concated_isolate[(cov,'assembler')] = ['SPAdes --isolate', 'SPAdes --isolate', 'SPAdes --isolate']
            concated = pd.concat([assembler_2, concated_isolate])
        
        else:
            #sama as above
            assembler_2 = concated[concated.index.str.contains(to_compare_to_isolate)]
            assembler_2[(cov, 'assembler')] = [assembler_not_isolate, assembler_not_isolate, assembler_not_isolate, assembler_not_isolate, assembler_not_isolate, assembler_not_isolate]
            concated_isolate = concated[concated.index.str.contains('I')]
            concated_isolate[(cov,'assembler')] = ['SPAdes --isolate', 'SPAdes --isolate', 'SPAdes --isolate', 'SPAdes --isolate', 'SPAdes --isolate', 'SPAdes --isolate']
            concated = pd.concat([assembler_2, concated_isolate])
        
        #Replace denotations to full names for plotting
        concated.index = concated.index.str.replace('SpI','')
        concated.index = concated.index.str.replace(to_compare_to_isolate,'')
        # concated.index = concated.index.str.replace('NP', 'Pilon \nNo trim')
        # concated.index = concated.index.str.replace('TP', 'Pilon \nTrimmomatic')
        # concated.index = concated.index.str.replace('FP', 'Pilon \nFastp')
        # concated.index = concated.index.str.replace('N', '')
        # concated.index = concated.index.str.replace('T', '')
        # concated.index = concated.index.str.replace('F', '')
        # concated.index = concated.index.str.replace('vs', '')
        
        #Plot figures
        plot = concated[[(cov, '# Corrections'), (cov, '# Errors'), (cov, '# Changes'), (cov,'assembler')]]
        diff_values = concated[[(cov, diff)]]
        diff_values.columns = diff_values.columns.droplevel(0) #remove multindex
        diff_list = diff_values[diff].tolist() #set diff values as list
        fig = plt.figure(figsize=(7,5))
        sns.set_theme()
        ax1 = fig.add_subplot(111)

        if to_plot != list_[3]:
            ax1.axvspan(-1, 2.5, facecolor='skyblue', alpha=0.4)
            ax1.axvspan(2.5, 6, facecolor='purple', alpha =0.15)
        else:
            ax1.axvspan(-1, 5.45, facecolor='skyblue', alpha=0.4)
            ax1.axvspan(5.45, 12, facecolor='purple', alpha =0.15)
        
        if cov != '100x':
            plot.plot(ax=ax1, kind='bar', title=cov, legend=False,
                         stacked=True, color=['mediumspringgreen', 'tomato', 'gold'], zorder=10)
            
            if to_plot == list_[3]:
                plt.xticks(rotation=90)
            else:
                plt.xticks(rotation=0)
            plt.xlabel(x_axis_label, labelpad=20, fontsize=15)
            plt.ylabel("Amount", fontsize=15)
    
    
        else:
            hej = plot.plot(ax=ax1, kind='bar', title=cov, stacked=True, 
                               color=['mediumspringgreen', 'tomato', 'gold'], zorder=10)
            if to_plot == list_[3]:
                plt.xticks(rotation=90)
            else:
                plt.xticks(rotation=0)
            plt.xlabel(x_axis_label, labelpad=20, fontsize=15)
            plt.ylabel("Amount", fontsize=15)


            ax1.legend(labels=[assembler_not_isolate, 'SPAdes --isolate','# Corrections', '# Errors',
                               '# Changes'], title='Legend', bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
            #ax1.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
        
        if to_plot == 'Pilon':
            ax1.set_title('B: Proportion of differences  '+cov, fontsize = 15)
        if to_plot == list_[2]:
            ax1.set_title('C: Proportion of differences  '+cov, fontsize = 15)
        if to_plot == 'Pilon and Filtering':
            ax1.set_title('D: Proportion of differences  '+cov, fontsize = 15)
        ax1.bar_label(ax1.containers[2],labels=['n=%g' % e for e in diff_list], label_type = 'center', padding=8, zorder=20) #set number of differences as n in plot

        plt.savefig(directory+'/'+to_plot+'_'+cov+'.png', dpi=300, bbox_inches="tight")


if to_plot == list_[0]:
    for cov in ['20x', '50x', '100x']:
        concated = pd.concat([df_NT, df_NF])
        concated = concated[[('20x', diff), ('20x', corr), ('20x', wrong), ('20x', change), ('50x', diff), ('50x', corr), ('50x', wrong), ('50x', change), ('100x', diff), ('100x', corr), ('100x', wrong), ('100x', change)]]
        for val in to_remove_from_trimming:
            concated = concated[~concated.index.str.contains(val)]
        concated = concated.groupby(concated.index).sum()
        # concated[(cov, '% Corrections')] = (concated[(cov, corr)])
        # concated[(cov, '% Errors')] = (concated[(cov, wrong)])
        # concated[(cov, '% Changes')] = (concated[(cov, change)])
        concated = concated.round(0)    
    
        assembler_2 = concated[concated.index.str.contains(to_compare_to_isolate)]
        assembler_2[(cov, 'assembler')] = [assembler_not_isolate, assembler_not_isolate]
        concated_isolate = concated[concated.index.str.contains('I')]
        concated_isolate[(cov,'assembler')] = ['SPAdes --isolate', 'SPAdes --isolate']
    
        concated = pd.concat([assembler_2, concated_isolate])
        concated.index = concated.index.str.replace('N',' ')
        concated.index = concated.index.str.replace('F','Fastp ')
        concated.index = concated.index.str.replace('T','Trimmomatic ')
        concated.index = concated.index.str.replace('vs',' ')
        concated.index = concated.index.str.replace('SpI','')
        concated.index = concated.index.str.replace(to_compare_to_isolate,'')
            
        plot = concated[[(cov, '# Corrections'), (cov, '# Errors'), (cov, '# Changes'), (cov,'assembler')]]
        
        diff_values = concated[[(cov, diff)]]
        diff_values.columns = diff_values.columns.droplevel(0)
        diff_list = diff_values[diff].tolist()
        fig = plt.figure(figsize=(7,5))
        sns.set_theme()
    
        ax1 = fig.add_subplot(111)
        ax1.axvspan(-1, 1.5, facecolor='skyblue', alpha=0.4)
        ax1.axvspan(1.5, 6, facecolor='purple', alpha =0.15)
        
        if cov != '100x':
            plot.plot(ax=ax1, kind='bar', title=cov, legend=False,
                         stacked=True, color=['mediumspringgreen', 'tomato', 'gold'])
            plt.xticks(rotation=0)
            plt.xlabel(x_axis_label, labelpad=20, fontsize=15)
            plt.ylabel("Amount", fontsize=15)
            ax2 = ax1.twiny()
            ax2.set_xlim(0, 10)
            ax2.set_xticks([2.5, 7.5])
            ax2.set_xticklabels([assembler_not_isolate, 'SPAdes --isolate'])
            ax2.tick_params(pad=17)
            ax2.xaxis.set_ticks_position("bottom")
    
        else:
            hej = plot.plot(ax=ax1, kind='bar', title=cov, stacked=True, 
                               color=['mediumspringgreen', 'tomato', 'gold'])
            plt.xticks(rotation=0)
            plt.xlabel(x_axis_label, labelpad=20, fontsize=15)
            plt.ylabel("Amount", fontsize=15)
            ax2 = ax1.twiny()
            ax2.set_xlim(0, 10)
            ax2.set_xticks([2.5, 7.5])
            ax2.set_xticklabels([assembler_not_isolate, 'SPAdes --isolate'])
            ax2.tick_params(pad=17)
            ax2.xaxis.set_ticks_position("bottom")
            ax1.legend(labels=[assembler_not_isolate, 'SPAdes --isolate','# Corrections', '# Errors',
                               '# Changes'], title='Legend', bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
            #ax1.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
            
        ax1.set_title('A: Proportion of differences  '+cov, fontsize = 15)
        ax1.bar_label(ax1.containers[2],labels=['n=%g' % e for e in diff_list], label_type = 'center', padding=8)
        
        plt.savefig(directory+'/Trimming_'+cov+'.png', dpi=300, bbox_inches="tight")






# to_plot2.plot(kind='bar', title='20x', stacked=True, color=['tomato','lightseagreen', 'purple'])
# to_plot3.plot(kind='bar', title = '50x', legend = False, stacked=True, color=['tomato','lightseagreen', 'purple'])
# to_plot4.plot(kind='bar', title = '100x', stacked = True, color=['tomato','lightseagreen', 'purple'])
# plt.legend(labels = [corr, wrong, change, 'assembler'], title='Labels', loc='upper right')
# plt.xticks(rotation=0)
# labels = [f'{i:.0%}' for i in to_plot2.to_numpy().flatten(order='F')]


# to_plot2 = concated[[('20x', '% Corrections'), ('20x', '% Errors'), ('20x', '% Changes'), ('20x','assembler')]]
# to_plot3 = concated[[('50x', '% Corrections'), ('50x', '% Errors'), ('50x', '% Changes')]]
# to_plot4 = concated[[('100x', '% Corrections'), ('100x', '% Errors'), ('100x', '% Changes')]]

# for i, patch in enumerate(ax.patches):
#     x, y = patch.get_xy()
#     x += patch.get_width() / 2
#     y += patch.get_height() / 2
#     ax.annotate(labels[i], (x, y), ha='center', va='center', c='white')







# fig, axes = plt.subplots(nrows=3, ncols=1, figsize = (8,30))
# plt.setp(axes, xlim=(0,500))
# to_plot2.plot(kind='bar', title='20x', ax = axes[0], legend = False)
# to_plot3.plot(kind='bar', title = '50x', ax=axes[1], legend = False)
# to_plot4.plot(kind='bar', title = '100x', ax=axes[2])

# plt.legend(labels = [diff, corr, wrong, change], title='Labels', loc='upper right')

# to_plot5 = concated[[('20x', diff), ('50x', diff), ('100x', diff)]]
# to_plot6 = concated[[('20x', corr), ('50x', corr), ('100x', corr)]]
# to_plot7 = concated[[('20x', wrong), ('50x', wrong), ('100x', wrong)]]
# to_plot8 = concated[[('20x', change), ('50x', change), ('100x', change)]]


# col_length = len(to_plot5.columns)

# for j in range(col_length):
#     for i in range(len(to_plot5)):
#         diffe = to_plot5.iloc[i,j]
#         to_plot6.iloc[i,j] = to_plot6.iloc[i,j]*(100/diffe)
#         to_plot7.iloc[i,j] = to_plot7.iloc[i,j]*(100/diffe)
#         to_plot8.iloc[i,j] = to_plot8.iloc[i,j]*(100/diffe)

# fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (20,20))
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
# to_plot5.plot(kind='bar', title=' % Differences', ax = axes[0,0])
# to_plot6.plot(kind='bar', title='% Corrections', ax = axes[0,1])
# to_plot7.plot(kind='bar', title='% Errors', ax = axes[1,0])
# to_plot8.plot(kind='bar', title='% Change', ax = axes[1,1])

# concated = concated.transpose()
# concated.columns = concated.columns.droplevel()

# concated = concated.transpose()

# col = concated.columns

# col_i = col[~col.str.contains('C')].to_list()
# col_c = col[~col.str.contains('I')].to_list()
# df_i = concated[col_i]
# df_i = df_i.reset_index()
# df_c = concated[col_c]
# df_c = df_c.reset_index()

# ax = sns.barplot(data=concated, x='variable', y='value', hue='typee')





# for cov in ['20x', '50x', '100x']:

#     concated = concated.groupby(concated.index).sum()
#     concated[(cov, '% Corrections')] = (concated[(cov, corr)]/concated[(cov, diff)]) * 100
#     concated[(cov, '% Errors')] = (concated[(cov, wrong)]/concated[(cov, diff)]) * 100
#     concated[(cov, change)] = (concated[(cov, change)]/concated[(cov, diff)]) * 100

# to_plot2 = concated[[('20x', '% Corrections'), ('20x', '% Errors'), ('20x', change)]]
# to_plot3 = concated[[('50x', '% Corrections'), ('50x', '% Errors'), ('50x', change)]]
# to_plot4 = concated[[('100x', '% Corrections'), ('100x', '% Errors'), ('100x', change)]]

# to_plot2 = to_plot2.reset_index()
# tuplee = [('index', '')]
# tidy = to_plot2.melt(id_vars=tuplee)

# tidy = tidy.rename(columns = {('index',  ''): "sample", "variable_1": "Trimming software"})

# to_plot2.plot(kind='bar', title='20x')

# sns.barplot(data=tidy, x="sample", y= "value", hue = "Trimming software" )

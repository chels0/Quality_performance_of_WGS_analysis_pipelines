#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 08:15:00 2022

@author: chelsea
"""
import itertools

filtering = ['\tfilter_contigs']
improv = ['\tassembly_improvement']
no_trim = ['\tno_trim']
fastp_trim_qc = ['\tfastp_trim_qc']
trimmomatic = ['\ttrimmomatic']
filter_set = ['\tfilter_set']
spades_set = ['\tspades_set']
bools = [' = false\n', ' = true\n']
trues = [' = true\n']
falses = [' = false\n']
spades_variables = [' = --careful\n', ' = --isolate\n']
filter_variables = [' = 100\n', ' = 200\n' , ' = 300\n']

no_trim_set = list(itertools.product(no_trim, trues, filtering, bools, improv, bools, spades_set, spades_variables))
fastp_trim_set = list(itertools.product(fastp_trim_qc, trues, filtering, bools, improv, bools, spades_set, spades_variables))
trimmomatic_set = list(itertools.product(trimmomatic, trues, filtering, bools, improv, bools, spades_set, spades_variables))

filter_list = []

for element in filter_variables:
    filter_list.append(filter_set[0] + element)

no_trim_settings = []
fastp_settings = []
trimmomatic_settings = []

for i in range(len(no_trim_set)):
    no_trim_set[i] = list(no_trim_set[i])
    fastp_trim_set[i] = list(fastp_trim_set[i])
    trimmomatic_set[i] = list(trimmomatic_set[i])
    
    no_fastp_trimmomatic = fastp_trim_qc + falses + trimmomatic + falses
    no_no_trim_fastp = no_trim + falses + fastp_trim_qc + falses
    no_no_trim_trimmomatic = no_trim + falses + trimmomatic + falses
    
    no_trim_set_new = no_trim_set[i]+no_fastp_trimmomatic
    no_trim_settings.append(no_trim_set_new)
    
    trimmomatic_set_new = trimmomatic_set[i] + no_no_trim_fastp
    trimmomatic_settings.append(trimmomatic_set_new)
    
    fastp_trim_set_new = fastp_trim_set[i] + no_no_trim_trimmomatic
    fastp_settings.append(fastp_trim_set_new)

count = 0

length = len(no_trim_settings[0])
for i in range(len(no_trim_settings)):
    no_trim_settings[i][0:length] = [''.join(no_trim_settings[i][0:length])]
    fastp_settings[i][0:length] = [''.join(fastp_settings[i][0:length])]
    trimmomatic_settings[i][0:length] = [''.join(trimmomatic_settings[i][0:length])]

pattern = filtering[0] + trues[0]

no_trim_settings2 = []
no_trim_list = []
for element in no_trim_settings:
    if pattern in element[0]:
        no_trim_complete = list(itertools.product(element, filter_list))
        no_trim_list.append(no_trim_complete)
    else:
        no_trim_settings2.append(element)

no_trim_settings = []
for element in no_trim_list:
    for i in range(len(element)):
        element[i] = list(element[i])
        element[i][0:length] = [''.join(element[i][0:length])]
        no_trim_settings.append(element[i])

fastp_settings2 = []
fastp_list = []
for element in fastp_settings:
     if pattern in element[0]:
         fastp_complete = list(itertools.product(element, filter_list))
         fastp_list.append(fastp_complete)
     else:
         fastp_settings2.append(element)

fastp_settings = []
for element in fastp_list:
    for i in range(len(element)):
        element[i] = list(element[i])
        element[i][0:length] = [''.join(element[i][0:length])]
        fastp_settings.append(element[i])
        
trimmomatic_settings2 = []
trimmomatic_list = []
for element in trimmomatic_list:
     if pattern in element[0]:
         trimmomatic_complete = list(itertools.product(element, filter_list))
         trimmomatic_list.append(fastp_complete)
     else:
         trimmomatic_settings2.append(element)

trimmomatic_settings = []
for element in trimmomatic_list:
    for i in range(len(element)):
        element[i] = list(element[i])
        element[i][0:length] = [''.join(element[i][0:length])]
        trimmomatic_settings.append(element[i])

concat_no_trim = no_trim_settings2 + no_trim_settings
concat_fastp = fastp_settings2 + fastp_settings
concat_trimmomatic = trimmomatic_settings2 + trimmomatic_settings

for element in concat_no_trim:
      count = count + 1
      with open('Run'+ str(count), 'w') as file2:
          for line in element:
              file2.write("%s" % line)
    
for element in concat_fastp:
    count = count + 1
    with open('Run'+ str(count), 'w') as file2:
        for line in element:
            file2.write("%s" % line)

for element in concat_trimmomatic:                
    count = count + 1
    with open('Run'+ str(count), 'w') as file2:
        for line in element:
            file2.write("%s" % line)
    
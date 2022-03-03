#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 09:45:11 2022

@author: chelsea
"""

import itertools

with open("/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/parameters.txt", "r") as f:
    seq = f.readlines()
    
cop = seq.copy()

empty = []
trimming = []
filtering = []
improv = []
no_trim = []
fastp_trim_qc = []
trimmomatic = []

for element in cop:
    empty.append(element.split('='))

list_of_bools = ['true' , 'false']
list_of_parameters = [i[0] for i in empty]




test = ['hej', 'då']

c = list(itertools.product(trimming, list_of_bools, filtering, list_of_bools, improv, list_of_bools,))

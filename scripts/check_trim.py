#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:57:48 2022

@author: chelsea
"""

#script for determining which leading value is needed to cut

import sys

#define variables from nextflow
lead_R1 = sys.argv[1]
lead_R2 = sys.argv[2]

#If lead from reverse read is higher than R1, use R2
if (lead_R2 > lead_R1):
    print(lead_R2)
else:
    print(lead_R1)


    

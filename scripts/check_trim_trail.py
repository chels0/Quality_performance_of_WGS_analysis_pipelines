#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:18:37 2022

@author: chelsea
"""
import sys

trail_R1 = sys.argv[1]
trail_R2 = sys.argv[2]

if (trail_R2 > trail_R1):
    print(trail_R2)
else:
    print(trail_R1)
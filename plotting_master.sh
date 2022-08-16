#!/usr/bin/env bash

#Modify flags if need be

python3 scripts/coverage_plot.py
python3 scripts/boxplot.py Ske no_trim
python3 scripts/boxplot.py Ske Fastp
python3 scripts/boxplot_post_mod.py Ske Pilon
python3 scripts/boxplot_post_mod.py Ske filtering
python3 scripts/boxplot_post_mod.py Ske both
python3 scripts/boxplot.py SpC no_trim
python3 scripts/boxplot.py SpC Fastp
python3 scripts/boxplot_post_mod.py SpC Pilon
python3 scripts/boxplot_post_mod.py SpC filtering
python3 scripts/boxplot_post_mod.py SpC both
python3 scripts/plotting.py 0 Ske
python3 scripts/plotting.py 1 Ske
python3 scripts/plotting.py 2 Ske
python3 scripts/plotting.py 3 Ske
python3 scripts/plotting.py 0 SpC
python3 scripts/plotting.py 1 SpC
python3 scripts/plotting.py 2 SpC
python3 scripts/plotting.py 3 SpC

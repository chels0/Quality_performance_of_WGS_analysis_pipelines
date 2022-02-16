#!/usr/bin/env bash

#Create schema

chewBBACA.py CreateSchema -i /home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Raw_data/cgMLST/jejuni_genomes/ -o jejuni_schema --ptf /home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Raw_data/Campylobacter_jejuni.trn --cpu 6

#Allele calling for wgMLST

chewBBACA.py AlleleCall -i /home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Raw_data/cgMLST/jejuni_genomes/ -g jejuni_schema/ -o jejuni_wgMLST --cpu 6

mv jejuni_wgMLST/result*/* jejuni_wgMLST/.

#Remove paralogs

chewBBACA.py RemoveGenes -i jejuni_wgMLST/results_alleles.tsv -g jejuni_wgMLST/RepeatedLoci.txt -o jejuni_wgMLST/results_alleles_NoParalogs

#Test genome quality

chewBBACA.py TestGenomeQuality -i jejuni_wgMLST/results_alleles_NoParalogs* -n 20 -t 300 -s 5 -o jejuni_wgMLST/genome_quality_32

#Extract cgMLST

chewBBACA.py ExtractCgMLST -i jejuni_wgMLST/results_alleles_NoParalogs* -p 0.95 -o jejuni_cgMLST

#Create that file

cat jejuni_cgMLST/cgMLSTschema.txt |sed -e "s|^|/home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca/jejuni_schema/|g" > jejuni_cgMLST/fullpath_cgMLSTschema.txt

#Allele call cgMLST

chewBBACA.py AlleleCall -i /home/chelsea/Documents/Quality_performance_of_WGS_analysis_pipelines/Results/No_trimming/assemblies/ -g jejuni_cgMLST/fullpath_cgMLSTschema.txt -o results_cgMLST_jejuni --cpu 6

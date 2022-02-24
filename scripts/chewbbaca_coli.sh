#!/usr/bin/env bash

#Create schema

chewBBACA.py CreateSchema -i /mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Raw_data/genome_assemblies_genome_fasta/ncbi-genomes-2022-02-22 -o coli_schema --ptf /mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Raw_data/Campylobacter_coli.trn --cpu 2

#Allele calling for wgMLST

chewBBACA.py AlleleCall -i /mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Raw_data/genome_assemblies_genome_fasta/ncbi-genomes-2022-02-22 -g coli_schema/ -o coli_wgMLST --cpu 2

mv coli_wgMLST/result*/* coli_wgMLST/.

#Remove paralogs

chewBBACA.py RemoveGenes -i coli_wgMLST/results_alleles.tsv -g coli_wgMLST/RepeatedLoci.txt -o coli_wgMLST/results_alleles_NoParalogs

#Test genome quality

chewBBACA.py TestGenomeQuality -i coli_wgMLST/results_alleles_NoParalogs* -n 20 -t 300 -s 5 -o coli_wgMLST/genome_quality_32

#Extract cgMLST

chewBBACA.py ExtractCgMLST -i coli_wgMLST/results_alleles_NoParalogs* -p 0.95 -o coli_cgMLST

#Create that file

cat coli_cgMLST/cgMLSTschema.txt |sed -e "s|^|/mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca/coli_chew/coli_schema/|g" > coli_cgMLST/fullpath_cgMLSTschema.txt

#Allele call cgMLST

chewBBACA.py AlleleCall -i /mnt/bigdisk/Quality_performance_of_WGS_analysis_pipelines/Results/all_assemblies/ -g coli_cgMLST/fullpath_cgMLSTschema.txt -o results_cgMLST_coli --cpu 2

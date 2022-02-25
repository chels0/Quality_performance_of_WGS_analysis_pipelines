#!/usr/bin/env bash

#Create schema

path=/home/chelsea/Documents

#mkdir Results/chewbbaca/Coli

mkdir Results/chewbbaca/Coli

chewBBACA.py CreateSchema -i Raw_data/coli_genomes/ -o Results/chewbbaca/Coli/coli_schema --ptf Raw_data/Campylobacter_coli.trn --cpu 6

#Allele calling for wgMLST

chewBBACA.py AlleleCall -i Raw_data/coli_genomes/ -g Results/chewbbaca/Coli/coli_schema/ -o Results/chewbbaca/Coli/coli_wgMLST --cpu 6

mv Results/chewbbaca/Coli/coli_wgMLST/result*/* Results/chewbbaca/Coli/coli_wgMLST/

#Remove paralogs

chewBBACA.py RemoveGenes -i Results/chewbbaca/Coli/coli_wgMLST/results_alleles.tsv -g Results/chewbbaca/Coli/coli_wgMLST/RepeatedLoci.txt -o Results/chewbbaca/Coli/coli_wgMLST/results_alleles_NoParalogs

#Test genome quality

chewBBACA.py TestGenomeQuality -i Results/chewbbaca/Coli/coli_wgMLST/results_alleles_NoParalogs* -n 20 -t 300 -s 5 -o Results/chewbbaca/Coli/coli_wgMLST/genome_quality

#Extract cgMLST

chewBBACA.py ExtractCgMLST -i Results/chewbbaca/Coli/coli_wgMLST/results_alleles_NoParalogs* -p 0.95 -o Results/chewbbaca/Coli/coli_cgMLST

#Create that file

cat Results/chewbbaca/Coli/coli_cgMLST/cgMLSTschema.txt |sed -e "s|^|/${path}/Quality_performance_of_WGS_analysis_pipelines/Results/chewbbaca/Coli/coli_schema/|g" > Results/chewbbaca/Coli/coli_cgMLST/fullpath_cgMLSTschema.txt

#Allele call cgMLST

chewBBACA.py AlleleCall -i Results/all_assemblies/ -g Results/chewbbaca/Coli/coli_cgMLST/fullpath_cgMLSTschema.txt -o Results/chewbbaca/Coli/Results_cgMLST_coli --cpu 6

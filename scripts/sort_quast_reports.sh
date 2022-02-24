#!/usr/bin/env bash

#script for sorting all .txt quast reports into one folder

ls Results/No_trimming/SPAdes/ > directory_list.txt

mkdir Results/All_quast_reports

path_no_trim='Results/No_trimming'
path_fastp='Results/Trimmed_w_fastp'
path_trimmomatic='Results/Trimmed_w_Trimmomatic'

for dir in $(cat directory_list.txt)
do
        cp ${path_no_trim}/SPAdes/${dir}/*_polished.fasta Results/all_assemblies/.
        cp ${path_no_trim}/SPAdes/${dir}/*_scaffolds.fasta Results/all_assemblies/.
        cp ${path_fastp}/SPAdes/${dir}/*_polished.fasta Results/all_assemblies/.
        cp ${path_fastp}/SPAdes/${dir}/*_scaffolds.fasta Results/all_assemblies/.
        cp ${path_trimmomatic}/SPAdes/${dir}/*_polished.fasta Results/all_assemblies/.
        cp ${path_trimmomatic}/SPAdes/${dir}/*_scaffolds.fasta Results/all_assemblies/.

        cp ${path_no_trim}/SKESA/${dir}/*.fasta Results/all_assemblies/.
        cp ${path_no_trim}/SKESA/${dir}/*.fasta Results/all_assemblies/.
        cp ${path_fastp}/SKESA/${dir}/*.fasta Results/all_assemblies/.
done

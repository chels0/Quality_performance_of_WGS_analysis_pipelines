#!/usr/bin/env bash

#script for placing all assemblies in one folder for chewbbaca



nr_of_files=$(ls Results/No_trimming/SPAdes/ | wc -l)

ls Results/No_trimming/SPAdes/ > directory_list.txt

mkdir Results/all_assemblies

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
	cp ${path_no_trim}/Pilon/SPAdes/${dir}/* Results/all_assemblies/.
	
	cp ${path_no_trim}/SKESA/${dir}/*.fasta Results/all_assemblies/.
	cp ${path_no_trim}/SKESA/${dir}/*.fasta Results/all_assemblies/.
	cp ${path_fastp}/SKESA/${dir}/*.fasta Results/all_assemblies/.
	cp ${path_no_trim}/Pilon/SKESA/${dir}/* Results/all_assemblies

done
#for  i in $(seq1 $nr_of_files)
#do
#	ln -s $dir/*_scaffolds.fasta Results/all_assemblies/.
#	ln -s $dir/*_polished.fasta Results/all

#!/usr/bin/env bash

#script for placing all assemblies in one folder for chewbbaca
#this version is for the following directory order sampleID --> No_trimming osv osv

ls Results/PT* > directory_list.txt

mkdir Results/all_assemblies

path_no_trim='Results/No_trimming'
path_fastp='Results/Trimmed_w_fastp'
path_trimmomatic='Results/Trimmed_w_Trimmomatic'

for dir in $(cat directory_list.txt)
do
	cp Results/${dir}/No_trimming/SPAdes/${dir}* Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_fastp/SPAdes/${dir}* Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_trimmomatic/SPAdes/${dir}/*_polished.fasta Results/all_assemblies/.
	cp Results/${dir}/No_trimming/Pilon/SPAdes/* Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_fastp/Pilon/SPAdes/* Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_trimmomatic/Pilon/SPAdes/* Results/all_assemblies/.
	
	cp Results/${dir}/No_trimming/SKESA/*.fasta Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_fastp/SKESA/*.fasta Results/all_assemblies/.
	cp Results/${dir}/Trimmed_w_trimmomatic/*.fasta Results/all_assemblies/.
	cp Results/${dir}/No_trimming/Pilon/SKESA/* Results/all_assemblies
	cp Results/${dir}/Trimmed_w_fastp/Pilon/SKESA/* Results/all_assemblies
	cp Results/${dir}/Trimmed_w_Trimmomatic/Pilon/SKESA/* Results/all_assemblies
done

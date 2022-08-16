#!/usr/bin/env bash

#Define output directory. If no output directory is chosen the final results will be in the pipeline Results folder. 
while getopts o: flag #define flag
do
    case "${flag}" in
        o) outdir=${OPTARG};;
    esac
done

if [ -z "${outdir}" ]; #if flag is empty, use Results folder.
then
        outdir=Results
fi

if [ ! -d "${outdir}" ] #if the directory does not exist kill the pipeline
then
    echo "Directory: ${outdir} DOES NOT exist. Exiting. Please create directory" 
    exit 9999 # die with error code 9999
fi

echo "your output directory is:" ${outdir} "cancel script if this is wrong"

echo "your output directory is:" ${outdir} "cancel script if this is wrong"

#Generate config files to be used in config_files folder
python3 scripts/automatisation_V4.py

#Add current path to txt file and add backslashes so it can be used with sed to automatically change path of the parameter_settings.txt file 
pwd > path1.txt
sed 's./.\\/.g' path1.txt > path.txt #adding backslashes
path=$(cat path.txt) 
sed "s/PATHS/$path/g" parameter_settings.txt > parameter_path.txt #exchange PATHS in parameter settings to the actual working path of the pipeline

> log.txt #create empty log

source ../../miniconda3/bin/activate env_version1 #activate conda environment

counter=0

#Iterate over each config file generatated and run the pipeline with those settings specified in the config file
for file in config_files/*;
do
	counter=$[counter + 1]

	head -16 parameter_path.txt > nextflow.config #Keep first 15 rows of parameter file which contains paremeters not to be changed and append to nextflow.config
	cat $file >> nextflow.config #append parameters in parameter file which will be changed with each iteration
	grep fastqc_set parameter_path.txt -A 10 >> nextflow.config #grep the last set of parameters not to be changed by the system and append to config file

	cp nextflow.config ./nextflow_config_copy

	date #print todays date
	echo "Starting pipeline with the following parameters"
	echo $file #print which run file is being processed
	cat $file #print which parameters are being used in pipeline this iteration

	#save the above prints to log file
	date >> log.txt
	echo "Starting pipeline with the following parameters" >> log.txt 
	echo $file >> log.txt
	cat $file >> log.txt

	#Run pipeline	
	nextflow pipeline7.nf -profile conda -resume

	#If outdir is not the Results folder in the pipeline directory, move the results to the wanted outdirectory


	if [ "${outdir}" != Results ];
	then
        	mv Results/* $outdir/.
	fi

	if [ "${counter}" == 8 ];
	then
		counter=0
		rm -r work #remove nextflow work directory
	fi

	>nextflow.config #empty the config file
done

bash scripts/sort_assemblies_v3.sh $outdir #sort the assemblies into correct folders for chewbbaca

path2=$(pwd) #define path
#put the path to the cgMLST genes into file
cat Raw_data/Campylobacter_jejuni/Schema/Cjejuni_cgMLST_678_listGenes.txt | sed -e "s|^|${path2}/Raw_data/Campylobacter_jejuni/Schema/Cjejuni_wgMLST_2795_schema/schema_seed_campy_roary_V5/|g" > fullpath_cgMLSTschema.txt

#run chewbbaca on each iteration and multiqc on quast reports
for dir in $outdir/*
do
	cp Raw_data/jejuniref2.fasta ${dir}/Assemblies/.
	chewBBACA.py AlleleCall -i ${dir}/Assemblies/ -g fullpath_cgMLSTschema.txt --cpu 8 -o ${dir}/chewBBACA/cgMLST_results_jejuni
	mv ${dir}/chewBBACA/cgMLST_results_jejuni/results*/* ${dir}/chewBBACA/cgMLST_results_jejuni/.
	multiqc ${dir}/PT*/Quast -o ${dir}/MultiQC
	cp ${dir}/MultiQC/multiqc_data/multiqc_quast.txt ${dir}/MultiQC/multiqc_quast.tsv
done        


source ../../miniconda3/bin/activate python_env_version1 #activate conda environment

rm -r Results/Comparisons
rm -r Results/chewbbaca_quast_tables
rm -r ${outdir}/Comparisons
rm -r ${outdir}/chewbbaca_quast_tables
rm -r Results/Conclusions
rm -r ${outdir}/Conclusions

python3 scripts/chewbbaca_result.py ${outdir}
mkdir Results/chewbbaca_quast_tables/placeholder 
mv Results/chewbbaca_quast_tables/*_results.tsv Results/chewbbaca_quast_tables/placeholder/.
python3 scripts/reduce_chewbbaca_2.py SpI Ske
python3 scripts/reduce_chewbbaca_2.py SpI SpC
python3 scripts/conclusion_script.py
python3 scripts/conclusion_script_trimming.py 
python3 scripts/difference_to_reference_all.py
mv Results/chewbbaca_quast_tables/placeholder/*
rm -r Results/chewbbaca_quast_tables/placeholder
mv Results/chewbbaca_quast_tables ${outdir}/.
mv Results/Comparisons ${outdir}/.
mv Results/Conclusions ${outdir}/.
rm -r Results/chewbbaca_quast_tables
rm -r Results/Comparisons
rm -r Results/Conclusions





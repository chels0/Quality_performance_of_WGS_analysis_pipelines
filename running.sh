#!/usr/bin/env bash

#för att kunna skriva config till det man vill

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
#python3 scripts/reduce_chewbbaca_2.py SpI SpC
python3 scripts/conclusion_script.py
python3 scripts/conclusion_script_trimming.py 
python3 scripts/difference_to_reference_all.py
python3 scripts/boxplot.py Ske no_trim
python3 scripts/boxplot.py Ske Fastp
python3 scripts/boxplot_post_mod.py Ske Pilon
python3 scripts/boxplot_post_mod.py Ske filtering
python3 scripts/boxplot_post_mod.py Ske both
#python3 scripts/boxplot.py SpC no_trim
#python3 scripts/boxplot.py SpC Fastp
#python3 scripts/boxplot_post_mod.py SpC Pilon
#python3 scripts/boxplot_post_mod.py SpC filtering
#python3 scripts/boxplot_post_mod.py SpC both
python3 scripts/plotting.py 0 Ske
python3 scripts/plotting.py 1 Ske
python3 scripts/plotting.py 2 Ske
python3 scripts/plotting.py 3 Ske
#python3 scripts/plotting.py 0 SpC
#python3 scripts/plotting.py 1 SpC
#python3 scripts/plotting.py 2 SpC
#python3 scripts/plotting.py 3 SpC
#python3 scripts/conclusion_TvsF.py
#python3 scripts/N50.py
mv Results/chewbbaca_quast_tables/placeholder/* ..
rm -r Results/chewbbaca_quast_tables/placeholder
mv Results/chewbbaca_quast_tables ${outdir}/.
mv Results/Comparisons ${outdir}/.
mv Results/Conclusions ${outdir}/.
#rm -r Results/chewbbaca_quast_tables
#rm -r Results/Comparisons
#rm -r Results/Conclusions





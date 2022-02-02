#!/usr/bin/env
 
source /home/chelsea/miniconda3/bin/activate env_version1 

mv nextflow.config ./nextflow_change.config

pwd > path1.txt
sed 's./.\\/.g' path1.txt > path.txt

path= "$(cat path.txt)" 
sed "s/PATHS/$path/g" nextflow_change.config > nextflow_config


nextflow pipeline2.nf -profile conda -resume

multiqc Results/. -o Results/MultiQC

rm -r work 


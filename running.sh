#!/usr/bin/env
 
source /home/chelsea/miniconda3/bin/activate env_version1 

nextflow pipeline.nf -profile conda -resume

multiqc Results/. -o Results/MultiQC

rm -r work 


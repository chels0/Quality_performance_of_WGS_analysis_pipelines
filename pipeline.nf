#!/usr/bin/env nextflow

params.str = 'Hello world!'

genomes = Channel.fromFilePairs('/mnt/bigdisk/Dokument/X5/Examensarbete/Raw_data/*{1,2}.fastq.gz')
process spades{
	
	publishDir '/home/chelsea/miniconda3/envs_test', mode: 'copy'

	
	input:
	tuple sampleID, file(reads) from genomes
	
	output:
	file '${sampleID}_scaffolds.fasta' into spades_output 
 
	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} --only-assembler -o $sampleID 	
	"""

}

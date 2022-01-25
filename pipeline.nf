#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
genomes.into { raw_data_for_Fastp; raw_data_for_trimmomatic }

process fastp{
	
	publishDir './Results/Fastp/', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp
	
	output:
	file "${sampleID}/*.gz" into fastp_trimmed
	file "${sampleID}" into fastp_reports
	val "${sampleID}" into fastp_id_name
	stdout result
	
	script:
	
	if ( params.fastp == true && params.fastp_trim == true )

		"""
		mkdir ${sampleID}
		mkdir ${sampleID}/Reports
		cd ${sampleID}
		fastp -i ../${reads[0]} -I ../${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
		mv fastp* ./Reports
		"""
		
	else if ( params.trimmomatic == false && params.fastp == false )

		"""
		echo "Trimming turned off"
		mkdir ${sampleID}
		mv ${reads[0]} ${reads[1]} ${sampleID}/.
		"""		
}

process trimmomatic{

	publishDir './Results/Trimmomatic/', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_trimmomatic 
	
	output:
	file "${sampleID}/*_paired.fq.gz" into trimmomatic_output 
	
	when:
	params.trimmomatic == true
	
	"""
	mkdir ${sampleID}
	cd ${sampleID}
	trimmomatic PE ${reads[0]} ${reads[0]} ${reads[0]}_paired.fq.gz ${reads[0]}_unpaired.fq.gz ${reads[1]}_paired.fq.gz ${reads[1]}_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	"""

}

process spades{
	
	publishDir './Results/SPAdes', mode: 'copy'

	input:
	file sampleID from fastp_trimmed
	
	val sample from fastp_id_name

	
	output:
	file "${sample}/*_scaffolds.fasta" into spades_output
	file "${sample}/*" into spades_all
	val "$sample" into spades_id_name 
 	
 	script:
 	if ( params.trimmomatic == false && params.fastp == false )
		"""
		spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample 
		mv ${sample}/scaffolds.fasta ${sample}/${sample}_scaffolds.fasta
		"""
	else if ( params.fastp_trim == true ) 
	
		"""
		spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample
		mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_fastp_scaffolds.fasta
		"""

}

process quast{
	
	publishDir './Results/Quast', mode: 'copy'

	input:
	file scaffold from spades_output
	val sampleID from spades_id_name
	
	output:
	file "$sampleID" into quast_output
	val "$sampleID" into quast_id_name
	
	"""
	quast.py $scaffold -o $sampleID
	"""
}

result.view { it }




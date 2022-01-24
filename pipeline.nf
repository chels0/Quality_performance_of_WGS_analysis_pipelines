#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
gen = Channel.fromFilePairs(params.path_to_reads)

process fastp{
	
	publishDir './Results/', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from gen
	
	output:
	file "Fastp/$sampleID" optional true into fastp_output
	val "$sampleID" into fastp_id_name
	stdout result
	
	script:
	
	if ( params.fastp == true && params.fastp_trim == true )
		
		"""
		mkdir Fastp
		mkdir Fastp/${sampleID}
		cd Fastp/${sampleID}
		fastp -i ../../${reads[0]} -I ../../${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
		cd ../..
		"""
		
	else if ( params.trimmomatic == false && params.fastp == false )

		"""
		echo "Trimming turned off"
		"""
	
		
}

process spades{
	
	publishDir './Results/SPAdes', mode: 'copy'

	input:
	tuple sampleID, file(reads) from genomes
	file sampleID_trim from fastp_output
	
	output:
	file "${sampleID}/*_scaffolds.fasta" into spades_output
	file "${sampleID}/*" into all_spades
	val "$sampleID" into spades_id_name 
 	
 	
 	script:
 	if ( params.trimmomatic == false && params.fastp == false )
		"""
		spades.py -1 ${reads[0]} -2 ${reads[1]} --only-assembler -o $sampleID
		mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_scaffolds.fasta 
		"""
	else if ( params.fastp_trim == true ) 
	
		"""
		spades.py -1 ${sampleID_trim}_R1_trimmed_fastp.fq.gz -2 ${sampleID_trim}_R2_trimmed_fastp.fq.gz --only-assembler -o $sampleID
		mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_trimmed__fastp_scaffolds.fasta
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




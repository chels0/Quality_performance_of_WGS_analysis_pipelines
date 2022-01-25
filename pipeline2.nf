#!/usr/bin/env nextflow


genomes = Channel.fromFilePairs(params.path_to_reads)
genomes.into { raw_data_for_Fastp; raw_data_for_trimmomatic; raw_data_for_spades }

process fastp{
	
	publishDir './Results/Trimmed_w_Fastp/Fastp/', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_Fastp
	
	output:
	file "${sampleID}/*.gz" into fastp_trimmed
	file "${sampleID}" into fastp_reports
	val "${sampleID}" into fastp_id_name
	stdout result
	
	when:
	params.fastp == true

	"""
	mkdir ${sampleID}
	mkdir ${sampleID}/Reports
	cd ${sampleID}
	fastp -i ../${reads[0]} -I ../${reads[1]} -o ${sampleID}_R1_trimmed_fastp.fq.gz -O ${sampleID}_R2_trimmed_fastp.fq.gz
	mv fastp* ./Reports
	"""
}

process trimmomatic{

	publishDir './Results/Trimmed_w_Trimmomatic/Trimmomatic/', mode: 'copy'
	
	input:
	tuple sampleID, file(reads) from raw_data_for_trimmomatic 
	
	output:
	file "${sampleID}/*_paired.fq.gz" into trimmomatic_output 
	val "${sampleID}" into id_name_trimmomatic
	
	when:
	params.trimmomatic == true
	
	"""
	mkdir ${sampleID}
	cd ${sampleID}
	trimmomatic PE ../${reads[0]} ../${reads[1]} ${sampleID}_R1_paired.fq.gz ${sampleID}_R1_unpaired.fq.gz ${sampleID}_R2_paired.fq.gz ${sampleID}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	"""

}

process spades_after_fastp{
	
	publishDir './Results/Trimmed_w_Fastp/SPAdes', mode: 'copy'

	input:
	file sampleID from fastp_trimmed
	val sample from fastp_id_name

	output:
	file "${sample}/*_scaffolds.fasta" into spades_output_fastp
	file "${sample}/*" into spades_all_fastp
	val "$sample" into spades_id_name_fastp 
 	
 	when:
 	params.fastp == true
 	
 	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_fastp_scaffolds.fasta
	"""

}

process spades_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/SPAdes', mode: 'copy'

	input:
	file sampleID from trimmomatic_output
	val sample from id_name_trimmomatic

	output:
	file "${sample}/*_scaffolds.fasta" into spades_output_trimmomatic
	file "${sample}/*" into spades_all_trimmomatic
	val "$sample" into spades_id_name_trimmomatic
 	
 	when:
 	params.trimmomatic==true

	"""
	spades.py -1 ${sampleID[0]} -2 ${sampleID[1]} --only-assembler -o $sample 
	mv ${sample}/scaffolds.fasta ${sample}/${sample}_trimmed_trimmomatic_scaffolds.fasta
	"""


}

process spades_no_trimming{
	
	publishDir './Results/No_trimming/SPAdes', mode: 'copy'

	input:
	tuple sampleID, file(reads) from raw_data_for_spades

	output:
	file "${sampleID}/*_scaffolds.fasta" into spades_output
	file "${sampleID}/*" into spades_all
	val "$sampleID" into spades_id_name 
 	
 	when:
 	params.no_trim == true

	"""
	spades.py -1 ${reads[0]} -2 ${reads[1]} --only-assembler -o $sampleID 
	mv ${sampleID}/scaffolds.fasta ${sampleID}/${sampleID}_scaffolds.fasta
	"""

}

process quast_no_trimming{
	
	publishDir './Results/No_trimming/Quast', mode: 'copy'

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

process quast_after_fastp{
	
	publishDir './Results/Trimmed_w_Fastp/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_fastp
	val sampleID from spades_id_name_fastp
	
	output:
	file "$sampleID" into quast_output_fastp
	val "$sampleID" into quast_id_name_fastp
	
	"""
	quast.py $scaffold -o $sampleID
	"""
}

process quast_after_trimmomatic{
	
	publishDir './Results/Trimmed_w_Trimmomatic/Quast', mode: 'copy'

	input:
	file scaffold from spades_output_trimmomatic
	val sampleID from spades_id_name_trimmomatic
	
	output:
	file "$sampleID" into quast_output_trimmomatic
	val "$sampleID" into quast_id_name_trimmomatic
	
	"""
	quast.py $scaffold -o $sampleID
	"""
}


result.view { it }

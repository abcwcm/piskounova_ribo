#!/bin/bash

genome_folder=/local/storage/maa7095/references/mouse/GRCm38.p6_M25/transcriptome/
genome_fasta=gencode.vM25.transcripts.fa
kmer_size=31
new_file_name=$(echo "$genome_fasta" | awk -F ".fna" '{print $1".fa"}')
echo $new_file_name
kallisto index \
	-i $genome_folder$new_file_name \
	-k $kmer_size \
	$genome_folder$genome_fasta


#done

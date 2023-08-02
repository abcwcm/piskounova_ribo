#!/bin/bash
# Loop to run kallisto
work_dir=/SSD/maa7095/RIBOseq/Process_xenograft/RIBO_data/trimmed_out/human_reads/non_conding/coding_reads/
genome_dir=/local/storage/maa7095/references/human/MANE/MANE.GRCh38.v1.1/kallisto_index_kmer21/
genome_name=MANE.GRCh38.v1.1.refseq_rna
#gtf_dir=/local/storage/maa7095/references/human/refseq/refseq_v40/
genome_gtf=MANE.GRCh38.v1.1.refseq_genomic.gtf
length=37 #important when its RIBO data the length has to change
SD=0.8

for f in $work_dir*gz
do
common_out="kallisto_"
#common="${work_dir}${common_out}"
sample_name=$(echo "$f" | awk -F "_trimmedUnmapped_outx_R1_001.fastq.gz" '{print $1}')
last_part=$(echo "$sample_name" | rev | cut -d "/" -f 1 | rev)
last_part="${common_out}${last_part}"


        cmd="kallisto quant \
                -i $genome_dir$genome_name \
                -o $work_dir$last_part \
                -t 4 \
		-l $length \
		-s $SD \
                --pseudobam \
		--single \
                --gtf $genome_dir$genome_gtf \
		$f"
        echo $cmd
done

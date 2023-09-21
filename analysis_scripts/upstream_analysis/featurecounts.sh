#!/bin/bash
# Run bbsplit
work_dir=/SSD/maa7095/RIBOseq/Process_xenograft/RNA_data/trimmed_out/
genome_human=/local/storage/maa7095/references/human/refseq/refseq_v40/GCF_000001405.40_GRCh38.p14_genomic.fna
genome_mouse=/local/storage/maa7095/references/mouse/refseq_GRCm39/GCF_000001635.27_GRCm39_genomic.fna
#gtf_dir=/local/storage/maa7095/references/human/refseq/refseq_v40/

for f in $work_dir*.fq.gz
do
#common="${work_dir}${common_out}"
sample_name=$(echo "$f" | awk -F ".fq.gz" '{print $1}')


        cmd="bbsplit.sh \
                build=1 \
                in=${sample_name}.fq.gz \
                ref_x=$genome_human \
                ref_y=$genome_mouse \
                basename=${sample_name}_out%.fq.gz\
	        scafstats=${sample_name}_scaf.txt \
		refstats=${sample_name}_ref.txt	"
               
        echo $cmd
done

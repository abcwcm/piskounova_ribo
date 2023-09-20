#!/bin/bash
# Run bbsplit
work_dir=/SSD/maa7095/RIBOseq/Shino_new_analysis/Xenograft/trimmed_out/
genome_human=/local/storage/maa7095/references/human/gencode_primary/GRCh38.primary_assembly.genome.fa
genome_mouse=/local/storage/maa7095/references/mouse/GRCm39/GRCm39.primary_assembly.genome.fa
#gtf_dir=/local/storage/maa7095/references/human/refseq/refseq_v40/
ambiguous2=toss

for f in $work_dir*_merged.c.tag.fq.gz
do
#common="${work_dir}${common_out}"
sample_name=$(echo "$f" | awk -F "_merged.c.tag.fq.gz" '{print $1}')


        cmd="bbsplit.sh \
                build=1 \
                in=${sample_name}_merged.c.tag.fq.gz \
                ref_human=$genome_human \
                ref_mouse=$genome_mouse \
                basename=${sample_name}_out%.fq.gz\
	        scafstats=${sample_name}_scaf.txt \
		refstats=${sample_name}_ref.txt	\
		path=${work_dir} \
		ambiguous2=${ambiguous2} "
               
        echo $cmd
done

#!/bin/bash

# Default mode is single-end
mode="single"

# Run bbsplit
work_dir=/SSD/maa7095/RIBOseq/Process_xenograft/RIBO_data/trimmed_out/
genome_human=/local/storage/maa7095/references/human/refseq/refseq_v40/GCF_000001405.40_GRCh38.p14_genomic.fna
genome_mouse=/local/storage/maa7095/references/mouse/refseq_GRCm39/GCF_000001635.27_GRCm39_genomic.fna

while getopts "m:" opt; do
  case $opt in
    m)
      mode="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

for f in "$work_dir"*_val_1.fq.gz
do
    # Common options
    sample_name=$(echo "$f" | awk -F "_val_1.fq.gz" '{print $1}')
    common_options="build=1 in=${sample_name}_val_1.fq.gz ref_x=$genome_human ref_y=$genome_mouse basename=${sample_name}_out%.fq.gz scafstats=${sample_name}_scaf.txt refstats=${sample_name}_ref.txt"

    if [ "$mode" = "single" ]; then
        # Single-end mode
        single_cmd="~/miniconda3/envs/bulkrnaseq/bin/bbsplit.sh $common_options"
        echo $single_cmd
    elif [ "$mode" = "paired" ]; then
        # Paired-end mode
        paired_cmd="~/miniconda3/envs/bulkrnaseq/bin/bbsplit.sh $common_options in2=${sample_name}_val_2.fq.gz "
        echo $paired_cmd
    else
        echo "Invalid mode: $mode"
        exit 1
    fi
done

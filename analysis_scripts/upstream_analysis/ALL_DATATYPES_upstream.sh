#!/bin/bash

# Set up raw data
raw_data_directory=""
genome_human=/local/storage/maa7095/references/human/gencode_primary/GRCh38.primary_assembly.genome.fa
genome_mouse=/local/storage/maa7095/references/mouse/GRCm39/GRCm39.primary_assembly.genome.fa
non_coding_folder=/local/storage/maa7095/references/human/non_coding
ribo_mouse_ref=/SSD/maa7095/RIBOseq/Shino_new_analysis/Xenograft/trimmed_out/ambioguos_mouse/ribo_xeno_ref/trimmed_out
MANE_annot=/local/storage/maa7095/references/human/gencode_primary/MANE_ref/MANE.GRCh38.v1.2.ensembl_genomic.gtf
Gencode_Mane=/local/storage/maa7095/references/human/gencode_primary/MANE_ref

# Define allowed AIM values
ALLOWED_AIM=("RIBO_CELL" "RIBO_PDX" "RNA_CELL" "RNA_PDX")

# Function to display usage instructions
show_usage() {
    echo "Usage: $0 --aim <AIM_OPTION> [OPTIONS]"
    echo "Options:"
    echo "  --raw-data-dir <RAW_DATA_DIR>: Specify the raw data directory."
    echo "Available AIM options: ${ALLOWED_AIM[*]}"
    exit 1
}

# Check if the provided --aim option is valid
validate_aim() {
    local aim="$1"
    for valid_aim in "${ALLOWED_AIM[@]}"; do
        if [ "$aim" == "$valid_aim" ]; then
            return 0  # Valid AIM
        fi
    done
    return 1  # Invalid AIM
}

# Process command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --aim)
            AIM_OPTION="$2"
            shift
            ;;
        --raw-data-dir)
            raw_data_directory="$2"
            shift
            ;;
        *)
            show_usage
            ;;
    esac
    shift
done

# Validate AIM option
validate_aim "$AIM_OPTION"
if [ $? -ne 0 ]; then
    echo "Error: Invalid --aim option. Please choose one of the following options: ${ALLOWED_AIM[*]}"
    exit 1
fi

# Check for required options
if [ -z "$raw_data_directory" ] ; then
    echo "Error: Missing required options."
    show_usage
fi

# Rest of your script here...
# Command 1: Trim Galore
function ribo_data_pipeline {
    local file="$1"
    local sample_name
    local f

    sample_name=$(basename "$file" .fastq.gz)
    f=$(basename "$sample_name")

    #Step 1 trimming

    trim_folder="${raw_data_directory}trimmed/"

    if [ ! -d "$trim_folder" ]; then
        mkdir -p "$trim_folder"
        echo "Directory created: $trim_folder"
    fi

    eval "$(conda shell.bash hook)"
    source activate bulkrnaseq

    trim_galore --fastqc --illumina -o  "$trim_folder" --basename "$f" "$file"

    wait

    conda deactivate
    

    collapsed_folder="${trim_folder}collapsed/"

    if [ ! -d "$collapsed_folder" ]; then
        mkdir -p "$collapsed_folder"
        echo "Directory created: $collapsed_folder"
    fi

    #Step 2 remove duplicated

    if [ "$aim" == "RIBO_CELL" ] || [ "$aim" == "RIBO_PDX" ]; then 
        eval "$(conda shell.bash hook)"
        source activate ctk_duplicates
        # Just for RIBO_CELL and RIBO_PDX

        fastq2collapse.pl "${trim_folder}${f}_trimmed.fq.gz" - | gzip -c > "${collapsed_folder}${f}_trimmed.c.fq.gz"
        stripBarcode.pl -format fastq -len 8 "${collapsed_folder}${f}_trimmed.c.fq.gz" - | gzip -c > "${collapsed_folder}${f}_trimmed.c.tag.fq.gz"

        wait

        conda deactivate 
        source activate bulkrnaseq

        if [ "$AIM_OPTION" == "RIBO_CELL" ]; then
            mv "${collapsed_folder}${f}_trimmed.c.tag.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
        else
            bbsplit.sh build=1 in="${collapsed_folder}${f}_trimmed.c.tag.fq.gz" ref_human="$genome_human" ref_mouse="$genome_mouse" basename="${collapsed_folder}${f}_out%.fq.gz" scafstats="${collapsed_folder}${f}_scaf.txt" refstats="${collapsed_folder}${f}_ref.txt" path="$collapsed_folder" ambiguous2=toss
            mv "${collapsed_folder}${f}_outhuman.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
        fi

        wait
    fi

    source activate bulkrnaseq

    if [ "$AIM_OPTION" == "RNA_PDX" ]; then
        bbsplit.sh build=1 in="${trim_folder}${f}_trimmed.fq.gz" ref_human="$genome_human" ref_mouse="$genome_mouse" basename="${collapsed_folder}${f}_out%.fq.gz" scafstats="${collapsed_folder}${f}_scaf.txt" refstats="${collapsed_folder}${f}_ref.txt" path="$collapsed_folder" ambiguous2=toss
        mv "${collapsed_folder}${f}_outhuman.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
    fi

    if [ "$AIM_OPTION" == "RNA_CELL" ]; then
        mv "${trim_folder}${f}_trimmed.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
    fi

    # Run mapping to noncoding reference
    mapping_folder="${collapsed_folder}mapped/"
    if [ ! -d "$mapping_folder" ]; then
        mkdir -p "$mapping_folder"
        echo "Directory created: $mapping_folder"
    fi

    # Remove noncoding RNA using a custom script


    STAR --runThreadN 4 \
        --genomeDir  "${non_coding_folder}" \
        --readFilesIn "${collapsed_folder}${f}_cleaned.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${mapping_folder}non_coding_${f}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode - \
        --outReadsUnmapped Fastx \
        --alignIntronMax 1 \
        --outFilterMultimapNmax 20 \
        --limitBAMsortRAM 3726878714 \
        --alignEndsType Local \
        --winAnchorMultimapNmax 100 \
        --seedSearchStartLmax 20 \
        --sjdbOverhang 43 \
        --outWigType None \
        --sjdbGTFfile -

    gzip "${mapping_folder}non_coding_${f}_"Unmapped.out.mate1

    if [ "$AIM_OPTION" == "RIBO_CELL" ]; then
        mv "${mapping_folder}non_coding_${f}_"Unmapped.out.mate1.gz "${mapping_folder}clean_${f}".fq.gz
       
    else
        # Remove mouse reads based on a custom reference

        STAR --runThreadN 4 \
            --genomeDir  ${ribo_mouse_ref} \
            --readFilesIn "${mapping_folder}non_coding_${f}_"Unmapped.out.mate1.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix "${mapping_folder}non_coding_ribo_mouse_${f}_" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode - \
            --outReadsUnmapped Fastx \
            --alignIntronMax 1 \
            --outFilterMultimapNmax 20 \
            --limitBAMsortRAM 3726878714 \
            --alignEndsType EndToEnd \
            --winAnchorMultimapNmax 100 \
            --seedSearchStartLmax 20 \
            --sjdbOverhang 43 \
            --outWigType None \
            --sjdbGTFfile -

    gzip "${mapping_folder}non_coding_ribo_mouse_${f}_"Unmapped.out.mate1
    mv "${mapping_folder}non_coding_ribo_mouse_${f}_"Unmapped.out.mate1.gz "${mapping_folder}clean_${f}".fq.gz
    fi

    # Third STAR Mapping to mapp to Human Genome

    human_mapping_folder="${mapping_folder}human_mapped/"
    if [ ! -d "$human_mapping_folder" ]; then
        mkdir -p "$human_mapping_folder"
        echo "Directory created: $human_mapping_folder"
    fi

    STAR --runThreadN 4 \
        --genomeDir ${Gencode_Mane} \
        --readFilesIn "${mapping_folder}clean_${f}".fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix  "${human_mapping_folder}human_mapped_${f}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --outReadsUnmapped Fastx \
        --alignIntronMax 1 \
        --outFilterMultimapNmax 20 \
        --limitBAMsortRAM 3726878714 \
        --alignEndsType EndToEnd \
        --winAnchorMultimapNmax 100 \
        --seedSearchStartLmax 20 \
        --sjdbOverhang 43 \
        --outWigType None \
        --sjdbGTFfile ${MANE_annot} \
        --outSAMattributes Standard

    echo "Processed: $file"
}

# Parallelize the ribo_data_pipeline function calls
#export -f ribo_data_pipeline
#find "$raw_data_directory" -type f -name "*.fastq.gz" | parallel ribo_data_pipeline


for file in ${raw_data_directory}*.fastq.gz; do
    ribo_data_pipeline "$file" &
done

wait
eval "$(conda shell.bash hook)"
source activate bulkrnaseq


human_mapping_folder="${raw_data_directory}trimmed/collapsed/mapped/human_mapped"

if [ "$AIM_OPTION" == "RNA_CELL" ] || [ "$AIM_OPTION" == "RNA_PDX" ]; then 
    featureCounts -t CDS -g gene_id -O -s 0 -a "${MANE_annot}" -o "${human_mapping_folder}/CDS_RNA_counts_not_strand.txt" "${human_mapping_folder}"/*Aligned.sortedByCoord.out.bam
fi

echo "All done"
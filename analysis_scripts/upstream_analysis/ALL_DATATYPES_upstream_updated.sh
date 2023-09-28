#!/bin/bash

# Define allowed AIM values
ALLOWED_AIM=("RIBO_CELL" "RIBO_PDX" "RNA_CELL" "RNA_PDX")

# Function to display usage instructions
show_usage() {
    echo "Usage: $0 --config <CONFIG_FILE> --aim <AIM_OPTION> [OPTIONS]"
    echo "Options:"
    echo "  --config <CONFIG_FILE>: Specify a custom configuration file."
    echo "  --aim <AIM_OPTION>: Specify the AIM option."
    echo "  --raw-data-dir <RAW_DATA_DIR>: Specify the raw data directory."
    echo "Available AIM options: ${ALLOWED_AIM[*]}"
    exit 1
}

# Process command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --config)
            config_file="$2"
            shift
            ;;
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

# Check if a configuration file was provided
if [ -z "$config_file" ]; then
    echo "Error: Please specify a custom configuration file using the --config option."
    show_usage
fi

# Include the specified configuration file
if [ -f "$config_file" ]; then
    source "$config_file"
else
    echo "Error: The specified configuration file '$config_file' does not exist."
    exit 1
fi

# Validate AIM option
validate_aim() {
    local aim="$1"
    for valid_aim in "${ALLOWED_AIM[@]}"; do
        if [ "$aim" == "$valid_aim" ]; then
            return 0  # Valid AIM
        fi
    done
    return 1  # Invalid AIM
}

# Check for required options
if [ -z "$raw_data_directory" ] ; then
    echo "Error: Missing required options."
    show_usage
fi


# Access the genome_human variable from the config file
echo "genome_human is set to: $genome_human"
echo "mouse genome is set to: $genome_mouse"
echo "noncoding RNA refernece is set to: $non_coding_folder"
echo "mouse riboseq reference is set to: $ribo_mouse_ref"
echo "annotation file is set to: $MANE_annot"
echo "indexed genome refernece is set to: $Gencode_Mane"




# Start with the pipeline
function ribo_data_pipeline {
    local file="$1"
    local sample_name
    local f

    sample_name=$(basename "$file" .fastq.gz)
    echo " $sample_name"
    sample_name2=$(basename "$sample_name")
    f=$(echo "$sample_name2" | awk -F "_S" '{print $1}')
    echo " *** Start working with sample $f ***"

    #Step 1 trimming
    echo " *** Start Step 1 trimming for $f ***"

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

    

    if [ "$AIM_OPTION" == "RIBO_CELL" ] || [ "$AIM_OPTION" == "RIBO_PDX" ]; then 
        eval "$(conda shell.bash hook)"
        source activate ctk_duplicates
        # Just for RIBO_CELL and RIBO_PDX

        #Step 2 trimming
        echo " *** Start Step 2 PCR duplicates removal for $f ***"

        fastq2collapse.pl "${trim_folder}${f}_trimmed.fq.gz" - | gzip -c > "${collapsed_folder}${f}_trimmed.c.fq.gz"
        wait
        stripBarcode.pl -format fastq -len 8 "${collapsed_folder}${f}_trimmed.c.fq.gz" - | gzip -c > "${collapsed_folder}${f}_trimmed.c.tag.fq.gz"

        conda deactivate 
        source activate bulkrnaseq

        wait

        if [ "$AIM_OPTION" == "RIBO_CELL" ]; then
            mv "${collapsed_folder}${f}_trimmed.c.tag.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
        else

            #Step 3 trimming
            echo " *** Start Step 3 removing mouse reads for $f ***"

            bbsplit.sh build=1 in="${collapsed_folder}${f}_trimmed.c.tag.fq.gz" ref_human="$genome_human" ref_mouse="$genome_mouse" basename="${collapsed_folder}${f}_out%.fq.gz" scafstats="${collapsed_folder}${f}_scaf.txt" refstats="${collapsed_folder}${f}_ref.txt" path="$collapsed_folder" ambiguous2=toss
            wait
            mv "${collapsed_folder}${f}_outhuman.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
        fi
        wait

        
    fi
    wait

    source activate bulkrnaseq

    if [ "$AIM_OPTION" == "RNA_PDX" ]; then

        #Step 3 trimming
        echo " *** Start Step 3 removing mouse reads for $f ***"
        bbsplit.sh build=1 in="${trim_folder}${f}_trimmed.fq.gz" ref_human="$genome_human" ref_mouse="$genome_mouse" basename="${collapsed_folder}${f}_out%.fq.gz" scafstats="${collapsed_folder}${f}_scaf.txt" refstats="${collapsed_folder}${f}_ref.txt" path="$collapsed_folder" ambiguous2=toss
        wait
        mv "${collapsed_folder}${f}_outhuman.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
    fi
    wait

    if [ "$AIM_OPTION" == "RNA_CELL" ]; then
        mv "${trim_folder}${f}_trimmed.fq.gz" "${collapsed_folder}${f}_cleaned.fq.gz"
    fi

    wait

    # Run mapping to noncoding reference
    mapping_folder="${collapsed_folder}mapped/"
    if [ ! -d "$mapping_folder" ]; then
        mkdir -p "$mapping_folder"
        echo "Directory created: $mapping_folder"
    fi

    # Remove noncoding RNA using a custom script
    #Step 4 trimming
    echo " *** Start Step 4 mapping to noncodingRNA for $f ***"

    wait
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

    wait

    if [ "$AIM_OPTION" != "RIBO_PDX" ]; then
        mv "${mapping_folder}non_coding_${f}_"Unmapped.out.mate1.gz "${mapping_folder}clean_${f}".fq.gz
       
    else

    wait
        # Remove mouse reads based on a custom reference
        #Step 5 
        echo " *** Start Step 5 mapping to mouse riboseq reference for $f ***"

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
    wait
    mv "${mapping_folder}non_coding_ribo_mouse_${f}_"Unmapped.out.mate1.gz "${mapping_folder}clean_${f}".fq.gz
    fi

    wait

    # Third STAR Mapping to mapp to Human Genome

    human_mapping_folder="${mapping_folder}human_mapped/"
    if [ ! -d "$human_mapping_folder" ]; then
        mkdir -p "$human_mapping_folder"
        echo "Directory created: $human_mapping_folder"
    fi

    #Step 6 
    echo " *** Start Step 6 mapping to human genome for $f ***"

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
#find "$raw_data_directory" -type f -name "*.fastq.gz" | parallel ribo_data_pipeline "{}"


for file in ${raw_data_directory}*.fastq.gz; do
    ribo_data_pipeline "$file" &
done

wait
eval "$(conda shell.bash hook)"
source activate bulkrnaseq


human_mapping_folder="${raw_data_directory}trimmed/collapsed/mapped/human_mapped"

cd "$human_mapping_folder"

if [ "$AIM_OPTION" == "RNA_CELL" ] || [ "$AIM_OPTION" == "RNA_PDX" ]; then 

    #Step 7 
    echo " *** Start Step 7 featurecounts for $f ***"

    featureCounts -t CDS -g gene_id -O -s 0 -a "${MANE_annot}" -o CDS_RNA_counts_not_strand.txt *Aligned.sortedByCoord.out.bam
fi

echo "All done"
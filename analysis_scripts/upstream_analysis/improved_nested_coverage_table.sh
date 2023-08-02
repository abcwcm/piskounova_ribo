#!/bin/bash
# This script only work with RIBO datasets since read length has to be exactly between 25-45 bp.
# Define two arrays, one for conditions and one for genes. 
bam_dir=/SSD/maa7095/RIBOseq/Process_xenograft/RIBO_data/trimmed_out/human_reads/non_conding/coding_reads/kallisto_output_kmer21/mapped_files/

# Nested loop over conditions and genes.
for file in ${bam_dir}*.bam; do
  for gene in ${bam_dir}*.genome; do
    echo "Processing $gene for $file"
    
    # Create a temporary variable for the file name.
    sample_name=$(echo "$file" | awk -F "[.]bam" '{print $1}')
    echo $sample_name2
    sample_name2=$(echo "$sample_name" | awk -F "/" '{print $NF}')
    gene_name=$(echo "$gene" | awk -F "[.]genome" '{print $1}')
    gene_name2=$(echo "$gene_name" | awk -F "/" '{print $NF}')
    echo $gene_name2
    output_prefix="${gene_name2}_${sample_name2}_"
    echo "Output prefix will be $output_prefix"

    new_folder="${bam_dir}${gene_name2}"

    if test -d "$new_folder"; then
  	echo "Directory already exists"
    else
    # Create the directory and any parent directories that don't exist
    mkdir -p "$new_folder"
   	echo "Directory created"
    fi 

    output_bam="${new_folder}/${output_prefix}filtered.bam"
    output_bed="${new_folder}/${output_prefix}5p.bed"
    output_cov="${new_folder}/${output_prefix}5p.cov"

    # Sort and index the BAM file.
    if test -f "${file}.sorted.bam.bai"; then
    echo "Indexed file  exists"
    else
    samtools sort -o "${file}.sorted.bam" "$file"
    samtools index "${file}.sorted.bam" 
    echo "File sorted and indexed"
    fi
    #Filter for the gene of interest and output to a new BAM file.
    samtools view -h "${file}.sorted.bam" ${gene_name2} | awk 'length($10)>25 && length($10)<45 || $1 ~ /^@/' | samtools view -bS - > "$output_bam"

    # Convert the filtered BAM file to a BED file and calculate genome coverage.
    bedtools bamtobed -i "$output_bam" | bedtools genomecov -5 -dz -g "${gene}" -i - > $output_cov


    #bedtools genomecov -5 -dz -ibam "$output_bam" > $output_cov
    # Print the sum of the coverage table's 3rd column.
    awk '{sum+=$3;} END {print sum;}'  $output_cov
    done
done


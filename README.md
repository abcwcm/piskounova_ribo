
# Selenocysteine tRNA methylation promotes oxidative stress resistance in melanoma metastasis
Scripts related to the analysis of RIBO-seq and bulkRNA-seq data of FTSJ1 KO cell lines and PDX data from primary and metastasis tumors by Nease et al.


![](WCM_MB_LOGO_HZSS1L_CLR_RGB_new.png)

## DATA AVAILABILITY

* raw data (sequencing reads) and processed counts (raw and normalized counts) can be downloaded from [**GEO: X**]

## DATA ANALYSIS

![](RIBO-seq_diagram.png)


## UPSTREAM ANALYSIS

The script to run all upstream steps is found in [ALL_DATATYPES_upstream.sh](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/upstream_analysis/ALL_DATATYPES_upstream.sh). It needs a configuration file where all the necessary reference and annotation files will be pointed. The script works as follows;

    > `./ALL_DATATYPES_upstream.sh --aim RIBO_CELL|RIBO_PDX|RNA_CELL|RNA_PDX --config config.yml  --raw-data-dir data_directory_path`
    
    Example:


    >`./ALL_DATATYPES_upstream_updated.sh --aim RIBO_CELL --config config.yml --raw-data-dir test_Data/ribo_cell/`

- **Step 1**. Reads from bulk RNA-seq and RIBO-seq libraries from cell lines were trimmed using `Trim Galore` v0.6.10 to remove nucleotides with low quality and adaptor contamination.

- **Step 2**. For RIBO-seq data, duplicated reads created by PCR amplification were removed based on UMI and using `fastq2collapse.pl` and `stripBarcode.pl` scripts from `CTK tool kit`.

- **Step 4**. From `Trim Galore` output for RNA-seq and `CTK` output for RIBO-seq, **non-coding RNA** was removed using a custom reference genome composed by miRNA, rRNA, tRNA and lncRNA sequences using `STAR` v2.7.9a with `–alignEndsType Local`.

- **Step 6**. `STAR` was again used with `–alignEndsType EndToEnd` and `–quantMode TranscriptomeSAM`, and using **GRCh38 primary assembly** genome and **MANE v1.2** annotation file to obtain transcriptome and genome mapping coordinates. 

- **Step 7**. Using bam files originated from the mapping of RNA-seq reads to the whole genome, quantification of reads mapping to CDS regions was run using `featureCounts` v 2.0.1.

For Xenograft data, in addition to the steps described above, mouse reads were removed in two different ways:

- **Step 3**. Before non-coding RNA removal `bbsplit.sh` from `BBMap` v38.90 was used with `ambiguous2==”toss”` using gencode **GRCh38 human** and **GRCm39 mouse** references keeping only reads that mapped unambiguously to human reference. 

- **Step 5**. In addition to this, a custom reference using **mouse healthy liver RIBO-seq** data was created and reads that didn’t align to this reference (=Human reads) using `STAR` and `–alignEndsType EndToEnd` were kept for following steps. 

- **Step 8**. Translation efficiency (TE) was calculated based on the CDS gene counts obtained from bulk RNA-seq and counts obtained for RIBO-seq data from `In-frame psite identification and quantification`.([Cell_line_TE.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script5_cell_lines_Translation_Efficiency.Rmd),[PXD_TE.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script5_PDX_Translation_Efficiency.Rmd))





## RIBO P-SITE BASED DOWNSTREAM ANALYSIS

From ribosome profile bam files originated from mapping to the transcriptome `Step 6` in both cell lines and PDX samples in-frame psite coverages were calculated using `riboWaltz` package. This in-frame psite coverages were used then to calculate stalling based on downstream cumulative distribution function (CDF), to identify readthrough by 5’ and 3’ footprints,  to calculate psite codon usage and finally and to quantify CDS in-frame psite counts for all genes to use in translation efficiency (TE) analysis. 


- **Step 7.1: In-frame psite identification and quantification**:

For in-frame psite coverage quantification, reads longer than 25 pb and shorter than 45 pb were kept and P-offsites were calculated for each read length using psite function with default arguments. Based on those P-offsite, CDS in-frame psite coverages were quantified for each gene and save for TE analysis. Once, psite position of each read was identified, only in-frame psites were kept for further stalliing CDF and heatmaps and readthrough analysis.  Psite coverages for each Selenoprotein can be found [here](https://github.com/abcwcm/piskounova_ribo/tree/main/selenoproteins_psite_counts). Codon usage of psite was calculated then using `codon_usage_psite` function from `riboWaltz` (Figure X). ([Cell_line_psite_identification.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script1_cell_lines_inframe_psite_idenitification.Rmd),[PXD_line_psite_identification.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script1_PDX_inframe_psite_identification.Rmd)).



- **Step 7.2: Stalling CDF plots**:

Once the psite coverage through each selenoprotein and nucleotide position was identified, mean values and cumulative fractions across conditions were calculated through all the gene body. ([Cell_line_CDF.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script2_cell_lines_CDF_plots.Rmd),[PXD_CDF.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script2_PDX_CDF_plots.Rmd))

- **Step 7.3: Stalling heatmaps**:

While for CDF plots (Figure X) psite coverage through all the gene body was used for the rest of the analysis only psite coverage falling in CDS regions was used. In this way, using fixed windows, difference values of cumulative fractions between KO and WT (Metastasis and Primary) were calculated (Figure X). Windows sizes were set as follows; 5’ and 3’ lengths from the Sec codon position (UGA) to start (5’) and stop (3’) codons were divide into five equal length segments. Therefore, each fragment would account for the 20 % window of 5’ or 3’ to Sec codon. Then maximum of cumulative fraction in each of those bins was identified for each KO and WT (Met and Primary) conditions and difference was calculated by KO vs WT (Met vs Primary).  ([Cell_line_stalling_heatmap.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script3_cell_lines_stalling_Heatmap_bins.Rmd),[PXD_stalling_heatmap.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script3_PDX_stalling_Heatmap_bins.Rmd))


- **Step 7.4: Ribosome readthrough by footprints plots**:

Ribosome readthrough was determined based on the obtained CDS psite coverage values from the RIBO-seq libraries (Figure X). 5’ and 3’ coverages were calculated on fixed -X and +X bases relative to the Sec codon. Since we were using only CDS psite values and for some of the selenoproteins the distance between Sec and stop codon was limited, these fixed bins were different between selenoproteins. To determine significant differences between KO and WT (Met and Primary) conditions in footprints ratios, unpaired two-sided Student t-test from `rstatix` package was applied with a p < 0.05 cutoff (https://github.com/kassambara/rstatix). 
([Cell_line_footprint.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script4_cell_lines_readthrough_footprints.Rmd),[PXD_footprint.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script4_PDX_readthrough_footprints.Rmd))


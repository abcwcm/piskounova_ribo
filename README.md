
# Selenocysteine tRNA methylation promotes oxidative stress resistance in melanoma metastasis
Scripts related to the analysis of RIBO-seq and bulkRNA-seq data of FTSJ1 KO cell lines and PDX data from primary and metastasis tumors by Nease et al.


![](WCM_MB_LOGO_HZSS1L_CLR_RGB_new.png)

## DATA AVAILABILITY

* raw data (sequencing reads) and processed counts (raw and normalized counts) can be downloaded from [**GEO: X**]


## UPSTREAM ANALYSIS

Reads from bulk RNA-seq and RIBO-seq libraries from cell lines were trimmed using `Trim Galore` v0.6.10 to remove nucleotides with low quality and adaptor contamination. For RIBO-seq data, duplicated reads created by PCR amplification were removed based on UMI and using `fastq2collapse.pl` and `stripBarcode.pl` scripts from `CTK tool kit`. For both type of libraries, non-coding RNA was removed using a custom reference genome composed by miRNA, rRNA, tRNA and lncRNA sequences using `STAR` v2.7.9a. Reads mapping to the custom reference were removed from further analysis. In a second mapping step, `STAR` was again used with `–alignEndsType EndToEnd` and `–quantMode TranscriptomeSAM`, and using **GRCh38 primary assembly** genome and **MANE v1.2** annotation file to obtain transcriptome and genome mapping coordinates. Using bam files originated from the mapping of RNA-seq reads to the whole genome, quantification of reads mapping to CDS regions was calculated using `featureCounts` v 2.0.1. These were used in downstream Translation Efficiency analysis ([corresponding_scrip.sh](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/upstream_analysis/improved_nested_coverage_table.sh). 



For Xenograft data, in addition to the steps described above, mouse reads were removed in two different ways. Before non-coding RNA removal `bbsplit.sh` from `BBMap` v38.90 was used with `ambiguous2==”toss”` using gencode **GRCh38 human** and **GRCm39 mouse** references keeping only reads that mapped unambiguously to human reference. In addition to this, a mouse reference using healthy liver RIBO-seq data was created and reads that didn’t align to this reference (=Human reads) using `STAR` were kept for following steps. 


Scripts used in this steps can be found [here] (https://github.com/abcwcm/piskounova_ribo/tree/main/analysis_scripts/upstream_analysis). 



## PSITE BASED DOWNSTREAM ANALYSIS

From ribosome profile bam files originated from mapping to the transcriptome in both cell lines and PDX samples in-frame psite coverages were calculated using `riboWaltz` package. This in-frame psite coverages were used then to calculate stalling based on downstream cumulative distribution function (CDF), to identify readthrough by 5’ and 3’ footprints,  to calculate psite codon usage and finally, to quantify inframe psite CDS counts for all genes to use in translation efficiency (TE) analysis. 


- **In-frame psite identification and quantification**:

For in-frame psite coverage quantification, reads longer than 25 pb and shorter than 45 pb were kept and P-offsites were calculated for each read length using psite function with default arguments. Once, psite position of each read was identified, only in-frame psites were kept for further analysis. Codon usage of psite was calculated then using `codon_usage_psite` function from `riboWaltz` (Figure X). ([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Cell_lines/Script1_cell_lines_inframe_psite_idenitification.Rmd)[corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/Xenograft/Script1_PDX_inframe_psite_identification.Rmd)).



- **Stalling CDF plots**:

Once the psite coverage through each selenoprotein and nucleotide position was identified, mean values and cumulative fractions across conditions were calculated through all the gene body. 

- **Stalling heatmaps**:

While for CDF plots (Figure X) psite coverage through all the gene body was used for the rest of the analysis only psite coverage falling in CDS regions was used. In this way, using fixed windows, difference values of cumulative fractions between KO and WT were calculated (Figure X). Windows sizes were set as follows; 5’ and 3’ lengths from the Sec codon position (UGA) to start (5’) and stop (3’) codons were divide into five equal length segments. Therefore, each fragment would account for the 20 % window of 5’ or 3’ to Sec codon. Then the maximum of cumulative fraction in each of those bins was identified for each KO and WT conditions and difference was calculated by KO vs WT.  


- **Ribosome readthrough by footprints plots**:
Ribosome readthrough was determined based on the obtained CDS psite coverage values from the RIBO-seq libraries (Figure X). 5’ and 3’ coverages were calculated on fixed -X and +X bases relative to the Sec codon. Since we were using only CDS psite values and for some of the selenoproteins the distance between Sec and stop codon was limited, these fixed bins were different between selenoproteins (refer to GitHub code to identify the number of bases used for each selenoprotein). To determine significant differences between KO and WT conditions in footprints ratios, unpaired two-sided Student t-test from rstatix package was applied with a p < 0.05 cutoff (https://github.com/kassambara/rstatix). 


([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/footprint_No_SE_Suppl.Rmd)). 



- **Delta TE**:
Lastly, gene quantification for TE analysis was based on CDS in-frame psite for every gene. 


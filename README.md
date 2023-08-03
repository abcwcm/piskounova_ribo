
# Selenocysteine tRNA methylation promotes oxidative stress resistance in melanoma metastasis
Scripts related to the analysis of RIBO-seq and bulkRNA-seq data of FTSJ1 KO cell lines and xenograft data from primary and metastasis tumors by Nease et al.


![](WCM_MB_LOGO_HZSS1L_CLR_RGB_new.png)

## DATA AVAILABILITY

* raw data (sequencing reads) and processed counts (raw and normalized counts) can be downloaded from [**GEO: X**]


## UPSTREAM ANALYSIS

Reads from bulk RNA-seq and RIBO-seq libraries from cell lines were trimmed using `Trim Galore` v0.6.10 to remove nucleotides with low quality and adaptor contamination. Non-coding RNA was removed using a custom reference genome composed by miRNA, rRNA, tRNA and lncRNA sequences and using `STAR` v2.7.9a with default parameters. Reads mapping to the custom reference were removed from further analysis. For mapping, **MANE.GRCh38.v1.1** transcriptome reference (DOI: s41586-022-04558-8) was indexed using `Kallisto` pseudoaligner with `--kmer-size=21` option and isoform quantification was performed with `quant` function using **MANE v1.1** annotation file.

For Xenograft data, in addition to the steps described above, mouse reads were removed before non-coding RNA removal using `bbsplit.sh` from `BBMap` v38.90 with default parameters and using **RefSeq GRCh38 v40** human and **GRCm39** mouse references. 

Scripts used in this steps can be found [here] (https://github.com/abcwcm/piskounova_ribo/tree/main/analysis_scripts/upstream_analysis). 

## DOWNSTREAM ANALYSIS

While, otained bam files from ribosome reads were used for downstream cumulative distribution function (CDF), 5’ and 3’ footprints and P and A-site codon usage analysis, quantification counts from `Kallisto` for both RIBO and RNA data were used for DeltaTE analysis. 

- **CDF plots**:
For CDF heatmaps, for each Selenoprotein reads longer than 25 pb and shorter than 45 pb were kept and coverage across each nucleotide was calculated using genomecov function from `bedtools` starting from 5’end of the gene ([corresponding_scrip.sh](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/upstream_analysis/improved_nested_coverage_table.sh)). For each condition, mean values and cumulative fractions across samples were calculated through each Selenoprotein for each nucleotide. Then, using fixed windows, difference of cumulative fractions between KO and WT were calculated. Windows sizes were set as follows; 5’ and 3’ lengths from the Sec codon position (UGA) were divide into five equal length segments. Therefore, each fragment would account for the 20 % window of 5’ or 3’ to Sec codon. Then the maximum of cumulative fraction in each of those bins was identified for each KO and WT conditions and difference was calculated by KO vs WT ([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/CDF_hetamaps_No_SE_Suppl.Rmd)).

- **RIBO occupancy plots**:
The Ribo occupancy were determined based on the obtained coverage values from the RIBO-seq libraries. 5’ and 3’ coverage were calculated on the -200 bases and +200 bases relative to the Sec codon. To determine significant differences between KO and WT conditions in footprints ratios, unpaired two-sided Student t-test from `rstatix` package was applied with a p < 0.05 cutoff ([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/footprint_No_SE_Suppl.Rmd)). 

- **Codon usage**:
Codon usage of P and A-sites was calculated following `riboWaltz` package. Reads shorter than 25 pb were removed and the P and A-site usage identification was performed using psite function with default arguments ([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/codon_usage_No_SE_Suppl.Rmd)).

- **Delta TE**:
Using isoform quantification obtained by `kallisto`, `tximport` package (PMID: 26925227) was used to calculate gene level counts for both ribosome and transcriptome profiles. These values were then used to calculate translation efficiency (TE) measures by `DESeq2` following `DeltaTE` approach ([corresponding_scrip.Rmd](https://github.com/abcwcm/piskounova_ribo/blob/main/analysis_scripts/downstream_analysis/TE_DESEQ_NO_SE.Rmd)). 

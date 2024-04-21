# TCGAPanCan-CRC
Part I of the Athens Comprehensive Cancer Center Project (https://accc.gr/index.html)

#### NOTE: *Core* parts of this REPOSITORY are also part of the broader ACCC-CRC project related to "*Molecular and functional profiling unravels targetable vulnerabilities in colorectal cancer*" upcoming publication.

#### Efstathios-Iason Vlachavas
###### DKFZ-Division of Molecular Genome Analysis (B050)
###### Efstathios-Iason.Vlachavas@dkfz-heidelberg.de
###### svlachavas@eie.gr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10959700.svg)](https://doi.org/10.5281/zenodo.10959700)

Relative public dataset: coadread_tcga_pan_can_atlas_2018

Publication: From The Cancer Genome Atlas (TCGA) consortium. {https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html}

Data were downloaded from the open-source resource cBioPortal web portal (https://www.cbioportal.org/)
and the respective files were utilized (https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018):

- File 1: **data_RNA_Seq_v2_expression_median.txt** which then saved as an .rda file for downstream analysis.

- File 2: **data_mutations_extended.txt** which then was further processed, filtered and saved as maf object and described in the respective Zenodo folder materials.

- File 3: **data_clinical_sample.txt**. Here we modified it to skip the first 4 # header lines and saved it initially as **coadread_tcga_pan_can_atlas_2018_clinical_data.tsv**. Then, after further processing and sample selection for analysis purposes, we ended we 1 phenotype-like object, namely **PanCan.CRC.Merged.Sub.450samples.Clin.Mut.tsv**. This is the main file that will be used in the analyses reproduced in the .qmd files.

###### These initial files, along with modified-processed-additional outputs, are all included in the *Data* & *Mut_Freq* folders for reproducing the complete analysis.

##### **Note**: Any 'save' function in R (as well as other commonly used functions in R to export objects in various formats, such as *write_tsv*) documented along with a hashtag # sign are typically used to demonstrate how an object was created and is not intended to be run as part of regular code execution.


## Description

This public dataset was used as an additional proxy to increase sample size regarding the analysis of the *in-house* ACCC-CRC cohort, to support and generalize molecular findings regarding consistently deregulated pathway/TF activities in MSI and MSS CRC tumors.

## Notes on data acquisition and cohort definition

Metadata utilization
- Clinical data was directly utilized from the cBioPortal web page as explained: Further processing on the selection on patients is described on the main manuscript, when the major criterion was to keep samples with KRAS/BRAF/RAS_RAF_GNAS WT tumors. For the complete criteria and methods please refer to the Materials & Methods section, as well as in the respective .qmd file.

RNASeq data
- The analysis starts essentially with the *semi-processed/batch normalized;RSEM estimated counts*
(Batch normalized from Illumina HiSeq_RNASeqV2), gene annotation and sample phenotype information.

##### For file size, reproducibility and memory considerations, for both datasets we start with .rda files of the raw data versions.

Further details on each of the utilized databases/tools/repositories based on their respective publications:

- progeny: https://doi.org/10.1038/s41467-017-02391-6

- dorothea: https://doi.org/10.1101/gr.240663.118

- decoupleR: https://doi.org/10.1093/bioadv/vbac016

- cBioPortal: https://doi.org/10.1158/2159-8290.CD-12-0095

## Important R packages that need to be installed for reproducing the analysis on the TCGA cohort:

```r

packages = c(
    "edgeR",
    "circlize",
    "ComplexHeatmap",
    "ReactomePA",
    "tidyverse",
    "progeny"
    "decoupleR",
    "dorothea",
    "PreMSIm",
    "limma",
    "here",
    "xml2",
    "downlit",
    "data.table",
    "org.Hs.eg.db",
    "ggplot2",
    "forcats",
    "DOSE",
    "clusterProfiler",
    "ggplot2"
)

if(!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}

library(BiocManager)

for(p in packages) {
    if(!requireNamespace(p)) {
        install(p)
    }
}

```
## Implementation

- The user just needs to download/clone the respective github repository;
- For example `git clone https://github.com/Jasonmbg/TCGAPanCan-CRC.git

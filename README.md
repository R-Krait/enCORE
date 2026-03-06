<p align="center">
<img width="175" src="src/Logo_enCORE.png">
</p>

# enCORE
**Enhancer-Enhancer Network-based Prediction of Clustered Open Regulatory Elements (CORE) using scATAC-seq data**

### Overview
enCORE is a computational framework for identifying highly interactive enhancer clusters from single-cell chromatin accessibility. enCORE uniquely defines such enhancer clusters as CORE (Clustered Open Regulatory Elements). enCORE operates solely on single-cell ATAC-seq data, without requiring single-cell RNA-seq or multimodal measurements (e.g., 10X Multiome RNA/ATAC).

<p align="center">
<img src="src/cover_image.png" width=90%>
</p>

### Installation of enCORE Package
You can install the development version of enCORE from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("R-Krait/enCORE")
```

### Requirements
The packages listed below are required dependencies for enCORE.

- `ArchR (>= 1.0.2)`
- `TxDb.Hsapiens.UCSC.hg38.knownGene`
- `TxDb.Hsapiens.UCSC.hg19.knownGene`
- `TxDb.Mmusculus.UCSC.mm10.knownGene`
- `org.Hs.eg.db`
- `org.Mm.eg.db`
- `GenomicRanges`
- `GenomicFeatures`
- `ChIPseeker`
- `data.table`
- `dplyr`
- `reshape2`
- `AnnotationDbi`
- `progress`
- `igraph`
- `mefa4`
- `parallel`
- `stringr`
- `coop`
- `chromVARmotifs`
- `motifmatchr`
- `kneedle`
- `scales`

The enCORE package also requires command-line tools, STARE & BEDTools.

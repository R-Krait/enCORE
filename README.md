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
- `TxDb.Hsapiens.UCSC.hg38.knownGene (>= 3.16.0)`
- `TxDb.Hsapiens.UCSC.hg19.knownGene (>= 3.2.2)`
- `TxDb.Mmusculus.UCSC.mm10.knownGene (>= 3.10.0)`
- `org.Hs.eg.db (>= 3.16.0)`
- `org.Mm.eg.db (>= 3.16.0)`
- `GenomicRanges (>= 1.50.2)`
- `GenomicFeatures (>= 1.50.4)`
- `ChIPseeker (>= 1.34.1)`
- `data.table (>= 1.16.0)`
- `dplyr (>= 1.1.4)`
- `reshape2 (>= 1.4.4)`
- `AnnotationDbi (>= 1.60.2)`
- `progress (>= 1.2.2)`
- `igraph (>= 1.5.1)`
- `mefa4 (>= 0.3-9)`
- `parallel (>= 4.2.1)`
- `stringr (>= 1.5.1)`
- `coop (>= 0.6-3)`
- `chromVARmotifs (>= 0.2.0)`
- `motifmatchr (>= 1.20.0)`
- `kneedle (>= 1.0.0)`
- `scales (>= 1.3.0)`

The enCORE package also requires command-line tools, STARE & BEDTools.

First, please install mamba as fast alternative to conda for package installation.
```
conda install conda-forge::mamba
```

Then, install STARE-ABC & BEDTools.
```
mamba install -c conda-forge -c bioconda stare-abc bedtools
```

### Example Usage
(...)

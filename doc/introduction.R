## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE, 
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# library(enCORE)
# library(dplyr)
# library(dplyr)
# library(ArchR)
# library(parallel)
# library(motifmatchr)
# library(chromVARmotifs)
# library(BSgenome.Hsapiens.UCSC.hg38)
# setwd("~/PSJ/enCORE_dev")
# h5disableFileLocking()
# addArchRThreads(13)
# addArchRGenome("hg38")
# data("human_pwms_v2")

## -----------------------------------------------------------------------------
# proj4 <- readRDS("~/PSJ/test_enCORE/Save-Proj_r4/Save-ArchR-Project.rds")
# 
# proj5 <- addCoAccessibility(
#   ArchRProj = proj4,
#   reducedDims = "IterativeLSI",
#   maxDist = 500000
# )
# 
# cA <- getCoAccessibility(
#   ArchRProj = proj5,
#   corCutOff = 0.2,
#   resolution = 1,
#   returnLoops = FALSE
# )

## -----------------------------------------------------------------------------
# proj5$Clusters2 <- mapLabels(proj5$Sample, newLabels = remapClust, oldLabels = names(remapClust))

## -----------------------------------------------------------------------------
# data_lump_enCORE <- Extract_initial_enhancer_candidates(proj_atac = proj5, data_coacc = cA, output_dir = "~/PSJ/test_enCORE", organism = "hg38")

## -----------------------------------------------------------------------------
# data_lump_enCORE <- Calculate_gABC_score(data_lump_enCORE = data_lump_enCORE, n_col = 5, STARE_dir = "~/anaconda3/bin", output_dir = "~/PSJ/test_enCORE/results_gABC_cr_f")

## -----------------------------------------------------------------------------
# proj6 <- addMotifAnnotations(ArchRProj = proj5, motifPWMs = human_pwms_v2, name = "Motif")
# proj6 <- addBgdPeaks(proj6)
# proj6 <- addDeviationsMatrix(
#   ArchRProj = proj6,
#   peakAnnotation = "Motif",
#   force = TRUE,
#   threads = 10
# )

## -----------------------------------------------------------------------------
# data_lump_enCORE <- Determine_TF_weight_threshold(proj_atac = proj6, data_lump_enCORE = data_lump_enCORE, list_cluster = c("CRC", "Normal"), use_default_thres = FALSE)

## -----------------------------------------------------------------------------
# data_lump_enCORE <- Distill_CORE_per_cluster(proj_atac = proj6, data_lump_enCORE = data_lump_enCORE, output_dir = "~/PSJ/test_enCORE/test_filtered", option = "active")


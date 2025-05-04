#' Calculate_gABC_score
#'
#' This function will perform gABC scoring using STARE.
#'
#' @param data_lump_enCORE A list containing the data required for running enCORE.
#' @param n_col The number of columns in peak-cluster CPM matrix for gABC scoring.
#' @param STARE_dir Directory containing STARE_ABCpp.
#' @param output_dir Output directory where the gABC scoring results will be saved.
#' @param threshold_score Threshold for the gABC score. Default = 0.01.
#' @export
Calculate_gABC_score <- function(data_lump_enCORE, n_col, STARE_dir, output_dir,
                                 threshold_score = 0.01) {
  data_dir <- data_lump_enCORE[["working_dir"]]
  organism <- data_lump_enCORE[["organism"]]

  if(organism == "hg38") {
    info_org <- "hg38.ncbiRefSeq_filtered_v2.gtf"
  } else if(organism == "hg19") {
    info_org <- "hg19.ncbiRefSeq_filtered_v2.gtf"
  } else if(organism == "mm10") {
    info_org <- "mm10.ncbiRefSeq_filtered_v2.gtf"
  } else {
    stop("enCORE only supports hg38, hg19, and mm10 as reference genome!")
  }

  print("Calculate gABC scores per cluster...")

  system(paste0(STARE_dir, "/STARE_ABCpp -b ", data_dir, "/peak_for_gABC_cr0.2.bed -n 4-",
                n_col, " -a ", data_dir, "/", info_org, " -f false -k 1 -t ",
                threshold_score, " -u ", data_dir, "/genes_for_filtering.txt -o ",
                output_dir, "/results"))

  data_lump_enCORE[["output_dir_gABC"]] <- output_dir

  return(data_lump_enCORE)
}

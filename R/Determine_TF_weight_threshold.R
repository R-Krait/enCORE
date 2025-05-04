#' Determine_TF_weight_threshold
#'
#' This function will determine TF weight threshold.
#'
#' @param proj_atac An ArchRProject object containing chromVAR CISBP v2 motif deviation.
#' @param data_lump_enCORE A list containing the data required for running enCORE.
#' @param list_cluster A list containing cell population information to be used for performing enCORE.
#' @param use_default_thres If TRUE, use predefined_thres as TF weight threshold. If not, use automatically calculated TF weight threshold.
#' @param predefined_thres Pre-defined default threshold. Default = 0.675.
#' @export
Determine_TF_weight_threshold <- function(proj_atac, data_lump_enCORE, list_cluster, use_default_thres = TRUE, predefined_thres = 0.675) {
  organism <- data_lump_enCORE[["organism"]]

  thres_weight <- numeric()
  for(n_clus in list_cluster) {
    tmp_weight_motifs <- ArchR::getMatrixFromProject(
      ArchRProj = proj_atac,
      useMatrix = "MotifMatrix",
      useSeqnames = NULL,
      verbose = TRUE,
      binarize = FALSE,
      threads = getArchRThreads(),
      logFile = createLogFile("getMatrixFromProject")
    )

    if(organism == "hg38") {
      motif_database <- chromVARmotifs::human_pwms_v2
    } else if(organism == "hg19") {
      motif_database <- chromVARmotifs::human_pwms_v2
    } else if(organism == "mm10") {
      motif_database <- chromVARmotifs::mouse_pwms_v2
    } else {
      stop("enCORE only supports hg38, hg19, and mm10 as reference genome!")
    }

    tmp_group_info <- tmp_weight_motifs$Clusters2
    tmp_group_info <- which(tmp_group_info == n_clus)
    tmp_weight_motifs <- tmp_weight_motifs@assays@data$z
    tmp_weight_motifs <- tmp_weight_motifs[, tmp_group_info]

    tmp_weight_motifs <- as.matrix(tmp_weight_motifs)
    tmp_row_id <- rownames(tmp_weight_motifs)
    tmp_weight_motifs <- rowMedians(tmp_weight_motifs)

    tmp_motif_row_name <- c()
    for(i in tmp_row_id) {
      tmp_motif_row_name <- c(tmp_motif_row_name, motif_database@listData[[i]]@name)
    }

    names(tmp_weight_motifs) <- tmp_motif_row_name

    n_sig_motifs <- length(which(tmp_weight_motifs >= 0.675))
    if(n_sig_motifs >= 30) {
      weight_thres <- 0.675
    } else {
      for(thr in seq(0.675, 0.075, -0.1)) {
        n_sig_motifs <- length(which(tmp_weight_motifs >= thr))

        if(n_sig_motifs >= 30) {
          weight_thres <- thr
          break
        } else {
          weight_thres <- -777
        }
      }
    }

    thres_weight <- c(thres_weight, weight_thres)
  }

  if(length(which(thres_weight < 0)) > 0) {
    idx_outlier <- which(thres_weight == -777)
    cpid_outlier <- list_cluster[idx_outlier]

    print(paste0("outlier cell population: [", cpid_outlier, "]!"))
    stop("Please exclude outlier cell population!")
  } else if(isTRUE(use_default_thres)) {
    thres_weight_f <- predefined_thres
    print(paste0("enCORE adopts the default threshold, ", predefined_thres, "!"))
  } else {
    thres_weight_f <- min(thres_weight)
    print(paste0("enCORE adopts the automatically calculated threshold!"))
  }

  data_lump_enCORE[["list_cluster"]] <- list_cluster
  data_lump_enCORE[["threshold_TF_weight"]] <- thres_weight_f

  return(data_lump_enCORE)
}

#' iterative_proximal_filtering
#'
#' This function will execute iterative proximal filtering to extract CORE regions from the active option.
#'
#' @param data_lump_enCORE A list containing the data required for running enCORE.
#' @param n_celltype Cell population information in list_cluster.
#' @param thres_inactive Threshold for the accessibility within promoter to determine inactive genes. Default = 12.5.
#' @param thres_eprox Effective distance of inactive proximal enhancers. Default = 7500 bp.
#' @export
iterative_proximal_filtering <- function(data_lump_enCORE, n_celltype, thres_inactive = 12.5, thres_eprox = 7500) {
  link_info <- data_lump_enCORE[["link_info"]]
  link_info_total_n <- data_lump_enCORE[["link_info_total_n"]]
  data_coacc <- data_lump_enCORE[["info_coacc"]]
  candidate_enh <- data_lump_enCORE[["candidate_enh"]]
  candidate_enh_u <- data_lump_enCORE[["candidate_enh_u"]]
  promoters_info <- data_lump_enCORE[["promoters_info"]]
  tss_info_filt_act <- data_lump_enCORE[["tss_info_filt_act"]]
  promoters_info_filt_act <- data_lump_enCORE[["promoters_info_filt_act"]]
  peak_info <- data_lump_enCORE[["peak_info"]]
  peak_for_addcorr <- data_lump_enCORE[["peak_for_addcorr"]]
  peak_activity <- data_lump_enCORE[["peak_activity"]]
  peak_for_gABC <- data_lump_enCORE[["peak_for_gABC"]]
  tmp_working_dir <- data_lump_enCORE[["working_dir"]]
  dir_scores <- data_lump_enCORE[["output_dir_gABC"]]
  dir_pot <- data_lump_enCORE[["output_dir_CORE"]]
  n_core <- data_lump_enCORE[["n_core"]]
  organism <- data_lump_enCORE[["organism"]]
  tmp_list_cluster <- data_lump_enCORE[["list_cluster"]]
  names_order <- data_lump_enCORE[["names_order"]]
  thres_TF_weight <- data_lump_enCORE[["threshold_TF_weight"]]

  if(organism == "hg38") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    db_org <- org.Hs.eg.db::org.Hs.eg.db
  } else if(organism == "hg19") {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    db_org <- org.Hs.eg.db::org.Hs.eg.db
  } else if(organism == "mm10") {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    db_org <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop("enCORE only supports hg38, hg19, and mm10 as reference genome!")
  }

  list_celltype <- as.list(c(n_celltype))
  names(list_celltype) <- c(n_celltype)

  test_list <- read.table(paste0(dir_pot, "/CORE_potential_", n_celltype, "_f.bed"),
                          sep = "\t", fill = TRUE)
  colnames(test_list) <- c("seqnames", "start", "end")

  test_list <- GenomicRanges::makeGRangesFromDataFrame(test_list)
  test_list <- ChIPseeker::annotatePeak(test_list, TxDb = txdb, level = "gene")

  test_list <- as.data.frame(test_list)

  anno.result <- AnnotationDbi::select(db_org, keys = as.character(test_list$geneId), columns = "SYMBOL", keytype = "ENTREZID")
  test_list$ENTREZID <- as.character(test_list$geneId)
  test_list$geneId <- NULL

  test_list <- dplyr::left_join(test_list, anno.result, by = "ENTREZID")
  test_list$ENTREZID <- NULL
  test_list$gene_id <- test_list$SYMBOL
  test_list$SYMBOL <- NULL

  test_list <- unique(test_list)
  test_list <- na.omit(test_list)
  test_list <- test_list$gene_id

  list_genes <- unique(test_list)

  for(i in 1:length(names(list_celltype))) {
    tmp <- c()
    for(j in 1:length(list_genes)) {
      tmp_peak_ids <- subset(promoters_info_filt_act, subset = gene_id == list_genes[j])

      if(nrow(tmp_peak_ids) == 0) {
        tmp_value <- 0

        tmp <- c(tmp, tmp_value)
      } else {
        tmp_peak_ids <- as.numeric(tmp_peak_ids$peak_idx)
        tmp_peak_ids <- subset(peak_info, subset = idx %in% tmp_peak_ids)

        if(nrow(tmp_peak_ids) == 1) {
          tmp_peak_ids <- rownames(tmp_peak_ids)

          tmp_peak_ids <- peak_activity[tmp_peak_ids, ]

          tmp_value <- sum(as.numeric(tmp_peak_ids[names(list_celltype)[i]]))
        } else {
          tmp_peak_ids <- rownames(tmp_peak_ids)

          tmp_peak_ids <- peak_activity[tmp_peak_ids, ]

          tmp_value <- sum(as.numeric(tmp_peak_ids[, names(list_celltype)[i]]))
        }
        tmp <- c(tmp, tmp_value)
      }
    }

    names(tmp) <- list_genes
    list_celltype[[names(list_celltype)[i]]] <- tmp

    print(i)
  }

  list_celltype <- as.data.frame(list_celltype[[1]])
  colnames(list_celltype) <- c("values")

  genes_inactive <- rownames(list_celltype)[which(list_celltype$values < thres_inactive)]
  tss_info_filt_act_t <- tss_info_filt_act[genes_inactive]

  for(i in 1:length(tss_info_filt_act_t)) {
    tmp_df <- tss_info_filt_act_t[[i]]
    tmp_df <- tmp_df[, c("seqnames", "start", "end")]
    tmp_df$start <- tmp_df$start - thres_eprox
    tmp_df$end <- tmp_df$end + thres_eprox

    tss_info_filt_act_t[[i]] <- tmp_df

    print(i)
  }

  tss_info_filt_act_t <- do.call(rbind, tss_info_filt_act_t)
  write.table(tss_info_filt_act_t, paste0(dir_pot, "/enhc_inactive_", n_celltype, ".bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  ### iterative filtering test
  test_list <- read.table(paste0(dir_pot, "/CORE_potential_", n_celltype, "_f.bed"), sep = "\t", fill = TRUE)
  colnames(test_list) <- c("seqnames", "start", "end")

  ### subtract
  system(paste0("bedtools subtract -a ", paste0(dir_pot, "/CORE_potential_", n_celltype, "_f.bed"),
                " -b ", paste0(dir_pot, "/enhc_inactive_", n_celltype, ".bed"),
                " > ", paste0(dir_pot, "/CORE_active_", n_celltype, "_initial.bed")))

  e_list <- read.table(paste0(dir_pot, "/CORE_active_", n_celltype, "_initial.bed"), sep = "\t", fill = TRUE)
  colnames(e_list) <- c("seqnames", "start", "end")

  e_list_ovl <- GenomicRanges::makeGRangesFromDataFrame(e_list)
  test_list_gr <- GenomicRanges::makeGRangesFromDataFrame(test_list)

  e_list_ovl <- GenomicRanges::findOverlaps(e_list_ovl, test_list_gr)
  e_list_ovl <- as.data.frame(e_list_ovl)

  e_list_ovl$diff_len <- abs((e_list$end[e_list_ovl$queryHits] - e_list$start[e_list_ovl$queryHits]) -
                               (test_list$end[e_list_ovl$subjectHits] - test_list$start[e_list_ovl$subjectHits]))

  e_list_iter <- subset(e_list_ovl, subset = diff_len != 0)
  e_list_iter <- e_list[e_list_iter$queryHits, ]

  rownames(e_list_iter) <- paste(e_list_iter$seqnames, e_list_iter$start, e_list_iter$end,
                                 sep = "_")
  rownames(e_list) <- paste(e_list$seqnames, e_list$start, e_list$end, sep = "_")

  e_list_iter <- GenomicRanges::makeGRangesFromDataFrame(e_list_iter, keep.extra.columns = TRUE)

  e_list_iter <- ChIPseeker::annotatePeak(e_list_iter, TxDb = txdb, level = "gene")

  e_list_iter <- as.data.frame(e_list_iter)

  anno.result <- AnnotationDbi::select(db_org, keys = as.character(e_list_iter$geneId), columns = "SYMBOL", keytype = "ENTREZID")
  e_list_iter$ENTREZID <- as.character(e_list_iter$geneId)
  e_list_iter$geneId <- NULL

  e_list_iter <- dplyr::left_join(e_list_iter, anno.result, by = "ENTREZID")
  e_list_iter$ENTREZID <- NULL
  e_list_iter$gene_id <- e_list_iter$SYMBOL
  e_list_iter$SYMBOL <- NULL

  e_list_iter <- unique(e_list_iter)
  e_list_iter_inactive <- subset(e_list_iter, subset = gene_id %in% genes_inactive)

  e_list_iter_inactive <- paste(e_list_iter_inactive$seqnames, e_list_iter_inactive$start,
                                e_list_iter_inactive$end, sep = "_")

  e_list <- e_list[setdiff(rownames(e_list), e_list_iter_inactive), ]
  write.table(e_list, paste0(dir_pot, "/CORE_active_", n_celltype, "_iter.bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

  ### sort .bed file
  system(paste0("bedtools sort -i ", paste0(dir_pot, "/CORE_active_", n_celltype, "_iter.bed"), " > ",
                paste0(dir_pot, "/CORE_active_", n_celltype, "_iter_sorted.bed")))

  ### merge overlapped regions
  system(paste0("bedtools merge -i ", paste0(dir_pot, "/CORE_active_", n_celltype, "_iter_sorted.bed"), " > ",
                paste0(dir_pot, "/CORE_active_", n_celltype, "_f.bed")))

  ### remove initial .bed files
  system(paste0("rm ", paste0(dir_pot, "/CORE_active_", n_celltype, "_iter.bed")))
  system(paste0("rm ", paste0(dir_pot, "/CORE_active_", n_celltype, "_iter_sorted.bed")))

  print("Write CORE from the active option!")
  print(paste0(n_celltype, ": Done!"))

  data_lump_enCORE <- as.list(c("link_info", "link_info_total_n", "info_coacc", "candidate_enh", "candidate_enh_u",
                                "promoters_info", "tss_info_filt_act", "promoters_info_filt_act",
                                "peak_info", "peak_for_addcorr", "peak_activity", "peak_for_gABC",
                                "working_dir", "output_dir_gABC", "output_dir_CORE", "organism",
                                "names_order", "list_cluster", "threshold_TF_weight"))
  names(data_lump_enCORE) <- c("link_info", "link_info_total_n", "info_coacc", "candidate_enh", "candidate_enh_u",
                               "promoters_info", "tss_info_filt_act", "promoters_info_filt_act",
                               "peak_info", "peak_for_addcorr", "peak_activity", "peak_for_gABC",
                               "working_dir", "output_dir_gABC", "output_dir_CORE", "organism",
                               "names_order", "list_cluster", "threshold_TF_weight")

  data_lump_enCORE[["link_info"]] <- link_info
  data_lump_enCORE[["link_info_total_n"]] <- link_info_total_n
  data_lump_enCORE[["info_coacc"]] <- data_coacc
  data_lump_enCORE[["candidate_enh"]] <- candidate_enh
  data_lump_enCORE[["candidate_enh_u"]] <- candidate_enh_u
  data_lump_enCORE[["promoters_info"]] <- promoters_info
  data_lump_enCORE[["tss_info_filt_act"]] <- tss_info_filt_act
  data_lump_enCORE[["promoters_info_filt_act"]] <- promoters_info_filt_act
  data_lump_enCORE[["peak_info"]] <- peak_info
  data_lump_enCORE[["peak_for_addcorr"]] <- peak_for_addcorr
  data_lump_enCORE[["peak_activity"]] <- peak_activity
  data_lump_enCORE[["peak_for_gABC"]] <- peak_for_gABC
  data_lump_enCORE[["working_dir"]] <- tmp_working_dir
  data_lump_enCORE[["output_dir_gABC"]] <- dir_scores
  data_lump_enCORE[["output_dir_CORE"]] <- dir_pot
  data_lump_enCORE[["n_core"]] <- n_core
  data_lump_enCORE[["organism"]] <- organism
  data_lump_enCORE[["names_order"]] <- names_order
  data_lump_enCORE[["list_cluster"]] <- tmp_list_cluster
  data_lump_enCORE[["threshold_TF_weight"]] <- thres_TF_weight

  return(data_lump_enCORE)
}

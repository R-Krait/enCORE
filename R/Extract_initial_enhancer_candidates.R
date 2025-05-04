#' Extract_initial_enhancer_candidates
#'
#' This function will extract initial enhancer candidates for gABC scoring.
#'
#' @param proj_atac An ArchRProject object.
#' @param data_coacc Co-accessibility information from ArchR::getCoAccessibility function.
#' @param output_dir Output directory where the required files for gABC scoring will be saved.
#' @param organism Organism information. Accepted values are "hg38", "hg19", and "mm10". Default = "hg38".
#' @param n_core An integer specifying the number of CPU-cores to perform enCORE. Default = 4.
#' @export
Extract_initial_enhancer_candidates <- function(proj_atac, data_coacc, output_dir, organism = "hg38", n_core = 4) {
  link_info <- as.data.frame(data_coacc@listData)
  link_info_total_n <- link_info

  peak_info <- proj_atac@peakSet
  peak_info <- data.frame(idx = 1:length(peak_info@seqnames),
                          chr = peak_info@seqnames, start = peak_info@ranges@start,
                          end = (peak_info@ranges@start + 500), midpoint = (peak_info@ranges@start + 250))
  rownames(peak_info) <- paste(peak_info$chr, peak_info$start, peak_info$end, sep = "_")

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

  res <- GenomicRanges::makeGRangesFromDataFrame(peak_info[, c("chr", "start", "end")])
  annotation_resLFC <- ChIPseeker::annotatePeak(res, TxDb = txdb, level = "gene")
  annotation_resLFC <- as.data.frame(annotation_resLFC)

  annotation_resLFC$peak_id <- paste(annotation_resLFC$seqnames, annotation_resLFC$start, annotation_resLFC$end, sep = "_")
  annotation_resLFC <- annotation_resLFC[, c("peak_id", "annotation")]

  peak_info$annotation <- annotation_resLFC$annotation

  peak_for_gABC <- ArchR::getMatrixFromProject(proj_atac, useMatrix='PeakMatrix')
  row_data_gABC <- paste(peak_info$chr, peak_info$start, peak_info$end, sep = "_")
  col_data_gABC <- peak_for_gABC@colData@rownames
  col_celltype_gABC <- peak_for_gABC@colData$Clusters2

  peak_for_gABC <- peak_for_gABC@assays
  peak_for_gABC <- peak_for_gABC@data@listData
  peak_for_gABC <- peak_for_gABC$PeakMatrix

  tmp <- data.table::as.data.table(summary(peak_for_gABC))
  tmp$i <- row_data_gABC[tmp$i]
  tmp$j <- col_celltype_gABC[tmp$j]
  tmp <- tmp[, sum(x), by = list(i, j)]

  gc()

  tmp_mat <- reshape2::dcast(tmp, i ~ j)
  tmp_mat$idx <- peak_info[tmp_mat$i, "idx"]
  tmp_mat <- tmp_mat %>% dplyr::arrange(idx)
  tmp_mat$idx <- NULL

  tmp_mat <- as.data.frame(tmp_mat)
  rownames(tmp_mat) <- tmp_mat$i
  tmp_mat$i <- NULL

  gc()

  tmp_mat[is.na(tmp_mat)] <- 0

  peak_for_gABC <- as.matrix(tmp_mat)
  peak_for_gABC <- apply(peak_for_gABC, 2, function(x) (x/sum(x))*1000000)

  peak_for_addcorr <- peak_for_gABC

  names_order <- colnames(peak_for_gABC)

  candidate_enh <- subset(peak_info, subset = !((annotation == "Promoter (<=1kb)") | (annotation == "Promoter (1-2kb)")))
  candidate_enh <- subset(candidate_enh, subset = idx %in% c(link_info$queryHits, link_info$subjectHits))
  candidate_enh <- paste(candidate_enh$chr, candidate_enh$start, candidate_enh$end, sep = "_")

  promoters_info <- GenomicFeatures::genes(txdb, single.strand.genes.only = FALSE)
  promoters_info <- GenomicRanges::trim(promoters_info)
  promoters_info <- GenomicRanges::resize(promoters_info, width=1, fix='start')
  promoters_info <- as.data.frame(promoters_info)

  anno.result <- AnnotationDbi::select(db_org, keys = as.character(promoters_info$group_name), columns = "SYMBOL", keytype = "ENTREZID")
  promoters_info$ENTREZID <- as.character(promoters_info$group_name)
  promoters_info$group_name <- NULL

  promoters_info <- dplyr::left_join(promoters_info, anno.result, by = "ENTREZID")
  promoters_info$ENTREZID <- NULL
  promoters_info$gene_id <- promoters_info$SYMBOL
  promoters_info$SYMBOL <- NULL
  promoters_info$group <- NULL

  promoters_info <- unique(promoters_info)
  promoters_info <- subset(promoters_info, seqnames %in% c(paste0("chr", 1:22), "chrX"))
  promoters_info <- na.omit(promoters_info)

  promoters_info <- split(promoters_info, f = promoters_info$gene_id)
  tss_info_filt_act <- promoters_info
  peak_ranges <- GenomicRanges::makeGRangesFromDataFrame(peak_info[, c("chr", "start", "end")])

  gc()

  print("Extract genes for filtering...")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(promoters_info),
                                   clear = FALSE,
                                   width= 65)

  for(i in 1:length(promoters_info)) {
    pb$tick()
    tmp_pr <- promoters_info[[i]][, c("seqnames", "start", "end")]
    tmp_pr$start <- tmp_pr$start - 2000
    tmp_pr$end <- tmp_pr$end + 2000
    tmp_pr <- GenomicRanges::makeGRangesFromDataFrame(tmp_pr)
    tmp_pr <- GenomicRanges::findOverlaps(tmp_pr, peak_ranges)

    tmp_pr <- tmp_pr@to
    if(length(tmp_pr) == 0) {
      promoters_info[[i]] <- 0
    } else {
      tmp_pr <- as.numeric(tmp_pr)
      tmp_pr <- peak_info$idx[tmp_pr]

      promoters_info[[i]] <- tmp_pr
    }

    Sys.sleep(1 / round(length(promoters_info) / 10))
  }

  pb$terminate()

  genes_for_filtering <- c()
  for(i in names(promoters_info)) {
    if(promoters_info[[i]][1] != 0) {
      genes_for_filtering <- c(genes_for_filtering, i)
    }
  }

  genes_for_filtering <- as.data.frame(genes_for_filtering)

  promoters_info <- promoters_info[genes_for_filtering$genes_for_filtering]

  print("Construct promoter information...")

  tot_pid <- as.character(unique(c(link_info$queryHits, link_info$subjectHits)))

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(promoters_info),
                                   clear = FALSE,
                                   width= 65)

  for(i in 1:length(promoters_info)) {
    pb$tick()
    tmp_pid <- promoters_info[[i]]
    idx_pid <- which(as.character(tmp_pid) %chin% tot_pid)

    if(length(idx_pid) == 0) {
      promoters_info[[i]] <- 0
    } else {
      promoters_info[[i]] <- tmp_pid[idx_pid]
    }

    Sys.sleep(1 / round(length(promoters_info) / 10))
  }

  pb$terminate()

  gc()

  promoters_info <- reshape2::melt(promoters_info)
  promoters_info <- subset(promoters_info, value != 0)
  colnames(promoters_info) <- c("peak_idx", "gene_id")

  promoters_info_filt_act <- promoters_info
  promoters_info_filt_act$peak_idx <- as.character(promoters_info_filt_act$peak_idx)

  network_link <- link_info[, c("queryHits", "subjectHits", "correlation")]
  colnames(network_link) <- c("V1", "V2", "weight")

  print("Remove duplicated edges...")
  network_link <- network_link[!duplicated(t(apply(network_link[, c(1, 2)], 1, sort))), ]
  network_link$V1 <- as.character(network_link$V1)
  network_link$V2 <- as.character(network_link$V2)

  network_link <- igraph::graph.data.frame(network_link, directed = FALSE)

  link_adj <- igraph::as_adjacency_matrix(network_link, attr = "weight")
  row_id <- intersect(rownames(link_adj), as.character(unique(promoters_info$peak_idx)))
  col_id <- intersect(colnames(link_adj), as.character(peak_info$idx[which(rownames(peak_info) %in% candidate_enh)]))
  link_adj <- link_adj[row_id, col_id]

  gc()

  link_adj <- mefa4::Melt(link_adj)
  colnames(link_adj) <- c("peak_idx", "enhc_idx", "weight")
  link_adj <- subset(link_adj, weight != 0)
  link_adj <- na.omit(link_adj)

  gc()

  link_adj$peak_idx <- as.character(link_adj$peak_idx)
  link_adj$enhc_idx <- as.character(link_adj$enhc_idx)
  link_adj$weight <- as.numeric(link_adj$weight)

  promoters_info$peak_idx <- as.character(promoters_info$peak_idx)
  promoters_info <- subset(promoters_info, peak_idx %in% link_adj$peak_idx)

  genes_for_filtering <- subset(genes_for_filtering, genes_for_filtering %in% promoters_info$gene_id)

  write.table(genes_for_filtering, paste0(output_dir, "/genes_for_filtering.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE,
              sep = "\t")

  print("Write genes for filtering...")

  candidate_enh_u <- unique(link_adj$enhc_idx)
  candidate_enh_u <- rownames(peak_info)[which(peak_info$idx %in% candidate_enh_u)]

  peak_activity <- peak_for_gABC

  peak_for_gABC <- peak_for_gABC[rownames(peak_for_gABC)[which(rownames(peak_for_gABC) %in% candidate_enh_u)], ]

  peak_for_gABC <- cbind(peak_info[rownames(peak_for_gABC), c(2:4)], peak_for_gABC)
  colnames(peak_for_gABC) <- c("#chr", colnames(peak_for_gABC)[2:3], names_order)

  write.table(peak_for_gABC, paste0(output_dir, "/peak_for_gABC_cr0.2.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  print("Write peak-cluster count matrix for gABC scoring...")

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
  data_lump_enCORE[["info_coacc"]] <- cA
  data_lump_enCORE[["candidate_enh"]] <- candidate_enh
  data_lump_enCORE[["candidate_enh_u"]] <- candidate_enh_u
  data_lump_enCORE[["promoters_info"]] <- promoters_info
  data_lump_enCORE[["tss_info_filt_act"]] <- tss_info_filt_act
  data_lump_enCORE[["promoters_info_filt_act"]] <- promoters_info_filt_act
  data_lump_enCORE[["peak_info"]] <- peak_info
  data_lump_enCORE[["peak_for_addcorr"]] <- peak_for_addcorr
  data_lump_enCORE[["peak_activity"]] <- peak_activity
  data_lump_enCORE[["peak_for_gABC"]] <- peak_for_gABC
  data_lump_enCORE[["working_dir"]] <- output_dir
  data_lump_enCORE[["output_dir_gABC"]] <- NULL
  data_lump_enCORE[["output_dir_CORE"]] <- NULL
  data_lump_enCORE[["n_core"]] <- n_core
  data_lump_enCORE[["organism"]] <- organism
  data_lump_enCORE[["names_order"]] <- names_order
  data_lump_enCORE[["list_cluster"]] <- NULL
  data_lump_enCORE[["threshold_TF_weight"]] <- NULL

  return(data_lump_enCORE)
}

#' extract_raw_CORE_per_cluster
#'
#' This function will extract raw CORE regions from each cell population.
#'
#' @param proj_atac An ArchRProject object containing chromVAR CISBP v2 motif deviation.
#' @param data_lump_enCORE A list containing the data required for running enCORE.
#' @param n_celltype Cell population information in list_cluster.
#' @param thres_cos_initial Threshold for the cosine similarity to extract interaction-shared decomposed components. Default = 0.4.
#' @param thres_cos_link_loss Threshold for the cosine similarity to perform Link Loss Correction. Default = 0.4.
#' @param thres_cos_re_merge Threshold for the cosine similarity to perform Sparsity Correction. Default = 0.92.
#' @param thres_corr_add Threshold for the co-accessibility value to incorporate additional nodes with high accessibility. Default = 0.49.
#' @param thres_topic_add Threshold for the accessibility to incorporate additional nodes with high accessibility. Default value is top 0.5% accessibility of original enhancer candidates. In other words, Default = 0.5.
#' @param max_dist Maximum effective distance between enhancer candidates. Default = 500000.
#' @param max_nCORE The maximum number of CORE regions. Default = 0.1*26000/1.31.
#' @param max_iter_elbow The maximum number of iterations for elbow point determination. Default = 10.
#' @param min_count The minimum CPM count for CORE constituents. Only constituents with a CPM above this threshold will be retained.
#' @export
extract_raw_CORE_per_cluster <- function(proj_atac, data_lump_enCORE, n_celltype,
                                         thres_cos_initial = 0.4, thres_cos_link_loss = 0.4,
                                         thres_cos_re_merge = 0.92, thres_corr_add = 0.49, thres_topic_add = 0.5, max_dist = 500000,
                                         max_nCORE = 0.1*26000/1.31, max_iter_elbow = 10, min_count = 5.0) {
  link_info <- data_lump_enCORE[["link_info"]]
  link_info_total_n <- data_lump_enCORE[["link_info_total_n"]]
  info_coacc <- data_lump_enCORE[["info_coacc"]]
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
  output_dir <- data_lump_enCORE[["output_dir_CORE"]]
  n_core <- data_lump_enCORE[["n_core"]]
  organism <- data_lump_enCORE[["organism"]]
  tmp_list_cluster <- data_lump_enCORE[["list_cluster"]]
  names_order <- data_lump_enCORE[["names_order"]]
  thres_TF_weight <- data_lump_enCORE[["threshold_TF_weight"]]

  gABC_per_cellType <- as.list(names_order)
  names(gABC_per_cellType) <- c(names_order)

  for(i in 1:length(gABC_per_cellType)) {
    tmp_gABC <- read.table(gzfile(paste0(dir_scores, "/results_ABCpp_scoredInteractions_", gABC_per_cellType[[i]], ".txt.gz")),
                           sep = "\t", fill = TRUE)

    colnames(tmp_gABC) <- c("chr", "start", "end", "ensembl_id", "gene_name", "peak_id", "signal_value",
                            "contact_frequency", "adapted_activity", "scaled_contact", "intergenic_score", "dist_TSS", "ABC_score")

    tmp_gABC$chr <- paste0("chr", tmp_gABC$chr)
    tmp_gABC$peak_id <- paste(tmp_gABC$chr, tmp_gABC$start, tmp_gABC$end, sep = "_")

    gABC_per_cellType[[i]] <- tmp_gABC
  }

  results_gABC <- gABC_per_cellType[[n_celltype]]

  print("Collect filtered enhancer candidates...")

  enhancer_predicted <- c()
  for(i in 1:length(gABC_per_cellType)) {
    tmp_enhp <- gABC_per_cellType[[i]]$peak_id
    enhancer_predicted <- c(enhancer_predicted, tmp_enhp)
  }

  enhancer_predicted <- unique(enhancer_predicted)

  if(organism == "hg38") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    db_org <- org.Hs.eg.db::org.Hs.eg.db
    motif_database <- chromVARmotifs::human_pwms_v2
  } else if(organism == "hg19") {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    db_org <- org.Hs.eg.db::org.Hs.eg.db
    motif_database <- chromVARmotifs::human_pwms_v2
  } else if(organism == "mm10") {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    db_org <- org.Mm.eg.db::org.Mm.eg.db
    motif_database <- chromVARmotifs::mouse_pwms_v2
  } else {
    stop("enCORE only supports hg38, hg19, and mm10 as reference genome!")
  }

  tss_info <- GenomicFeatures::genes(txdb, single.strand.genes.only = FALSE)
  tss_info <- GenomicRanges::trim(tss_info)
  tss_info <- GenomicRanges::resize(tss_info, width=1, fix='start')
  tss_info <- as.data.frame(tss_info)

  anno.result <- AnnotationDbi::select(db_org, keys = as.character(tss_info$group_name), columns = "SYMBOL", keytype = "ENTREZID")
  tss_info$ENTREZID <- as.character(tss_info$group_name)
  tss_info$group_name <- NULL

  tss_info <- dplyr::left_join(tss_info, anno.result, by = "ENTREZID")
  tss_info$ENTREZID <- NULL
  tss_info$gene_id <- tss_info$SYMBOL
  tss_info$SYMBOL <- NULL
  tss_info$group <- NULL

  tss_info <- unique(tss_info)
  tss_info <- subset(tss_info, seqnames %in% c(paste0("chr", 1:22), "chrX"))
  tss_info <- na.omit(tss_info)

  tss_info <- data.frame(gene_name = tss_info$gene_id, TSS = tss_info$start)
  tss_info <- subset(tss_info, gene_name %in% results_gABC$gene_name)

  results_gABC <- dplyr::left_join(results_gABC, tss_info, by = "gene_name")
  results_gABC$TSS <- as.numeric(results_gABC$TSS)

  results_gABC_list <- split(results_gABC, f = results_gABC$chr)
  total_gene_id <- results_gABC_list

  for(i in names(total_gene_id)) {
    total_gene_id[[i]] <- unique(total_gene_id[[i]]$gene_name)
  }

  results_gABC_list <- lapply(results_gABC_list, split_by_peak_id)

  print("Convert gABC scores to valid vector...")
  for(i in names(results_gABC_list)) {
    results_gABC_list[[i]] <- lapply(results_gABC_list[[i]], convert_valid_vector, y = total_gene_id[[i]])
  }

  # summary
  link_info$distance <- NA
  link_info$ABC_similarity <- NA

  print("Remove duplicated edges...")
  link_info <- link_info[!duplicated(t(apply(link_info[, c(1, 2)], 1, sort))), ]
  link_info <- link_info[, c("queryHits", "subjectHits", "seqnames", "distance", "ABC_similarity")]

  within_rnet <- sort(unique(as.integer(c(link_info$queryHits, link_info$subjectHits))))
  in_silico_self_loops <- data.frame(queryHits = within_rnet, subjectHits = within_rnet,
                                     seqnames = peak_info$chr[within_rnet],
                                     distance = 0, ABC_similarity = 1)

  link_info <- rbind(link_info, in_silico_self_loops)
  link_info <- link_info %>% dplyr::arrange(queryHits, subjectHits)

  rownames(link_info) <- 1:length(link_info$queryHits)

  link_info$query_id <- rownames(peak_info)[link_info$queryHits]
  link_info$subject_id <- rownames(peak_info)[link_info$subjectHits]

  link_info <- subset(link_info, subset = ((query_id %in% enhancer_predicted) & (subject_id %in% enhancer_predicted)))

  gc()

  print("Assign cosine similarity to network edge...")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(link_info$seqnames),
                                   clear = FALSE,
                                   width= 65)

  peak_id_incl <- results_gABC$peak_id
  store_x <- match(link_info$queryHits, peak_info$idx)
  store_y <- match(link_info$subjectHits, peak_info$idx)

  store_id_x <- paste(peak_info$chr[store_x], peak_info$start[store_x], peak_info$end[store_x], sep = "_")
  store_id_y <- paste(peak_info$chr[store_y], peak_info$start[store_y], peak_info$end[store_y], sep = "_")

  for(i in 1:length(link_info$seqnames)) {
    pb$tick()
    x <- store_x[i]
    y <- store_y[i]

    id_x <- store_id_x[i]
    id_y <- store_id_y[i]

    if((id_x %in% peak_id_incl) & (id_y %in% peak_id_incl)) {
      if(is.na(link_info$ABC_similarity[i]) == TRUE) {
        chr_info <- as.character(link_info$seqnames[i])
        link_info$distance[i] <- abs(peak_info$midpoint[x] - peak_info$midpoint[y])

        link_info$ABC_similarity[i] <- coop::cosine(results_gABC_list[[chr_info]][[id_x]], results_gABC_list[[chr_info]][[id_y]])
      }
    } else {
      if(is.na(link_info$ABC_similarity[i]) == FALSE) {
        link_info$ABC_similarity[i] <- NA
      }
    }

    Sys.sleep(1 / round(length(link_info$seqnames) / 10))
  }

  pb$terminate()

  link_info <- na.omit(link_info)

  link_info$ABC_similarity[which(link_info$ABC_similarity < thres_cos_initial)] <- 0

  link_info_list <- split(link_info, f = link_info$seqnames)

  print("Execute link loss correction...")
  link_info_list <- lapply(link_info_list, extract_decomposed_components, results_gABC_list = results_gABC_list,
                           thres_cos = thres_cos_link_loss, n_core = n_core, max_dist = max_dist)

  print("Distance threshold optimization...")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(link_info_list),
                                   clear = FALSE,
                                   width= 65)

  dist_intra <- numeric()
  dist_inter <- numeric()
  for(i in 1:length(link_info_list)) {
    pb$tick()
    tmp_list <- link_info_list[[i]]
    tmp_list <- tmp_list[["changed_point"]]
    tmp_list$distance_diff_t <- c(0, diff(tmp_list$midpoint))

    for(j in 2:nrow(tmp_list)) {
      if(tmp_list$component[j] == tmp_list$component[(j-1)]) {
        dist_intra <- c(dist_intra, tmp_list$distance_diff_t[j])
      } else {
        dist_inter <- c(dist_inter, tmp_list$distance_diff_t[j])
      }
    }

    Sys.sleep(1 / round(length(link_info_list) / 10))
  }

  pb$terminate()

  ### plot distance intra-components
  hist(dist_intra)

  info_dist <- seq(500, 100000, by = 500)
  p_int_ita <- numeric()
  p_int_ite <- numeric()

  for(i in 1:length(info_dist)) {
    tmp_intra <- length(which(dist_intra <= info_dist[i])) / length(dist_intra)
    tmp_inter <- length(which(dist_inter > info_dist[i])) / length(dist_inter)

    p_int_ita <- c(p_int_ita, tmp_intra)
    p_int_ite <- c(p_int_ite, tmp_inter)
  }

  lld <- p_int_ita - p_int_ite

  idx_dist_thres <- which(lld >= 0)
  idx_dist_thres <- as.numeric(idx_dist_thres[1])
  idx_dist_thres <- (idx_dist_thres * 500) + 500

  print(paste0("Distance threshold is ", idx_dist_thres/1000, "kb"))
  print("Reduce sharing effect based on distance threshold...")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(link_info_list),
                                   clear = FALSE,
                                   width= 65)

  for(i in 1:length(link_info_list)) {
    pb$tick()
    link_per_cluster <- link_info_list[[i]][["decomposed_components"]]
    chr_info <- link_info_list[[i]][["changed_point"]]$chr[1]

    for(j in 1:length(link_per_cluster)) {
      link_cl <- link_per_cluster[[j]]
      link_cl$dist_diff <- c(0, diff(link_cl$midpoint))
      idx_decomposed_point <- which(link_cl$dist_diff > idx_dist_thres)

      if(length(idx_decomposed_point) == 0) {
        component_updated <- rep(paste(chr_info, link_cl$component[1], 0, sep = "_"), nrow(link_cl))
        link_cl$component_stitching <- component_updated
      } else {
        idx_decomposed_point <- c((idx_decomposed_point[1] - 1), diff(idx_decomposed_point), (nrow(link_cl) - idx_decomposed_point[length(idx_decomposed_point)] + 1))

        component_updated <- c()
        for(k in 1:length(idx_decomposed_point)) {
          tmp_decompose <- rep(paste(chr_info, link_cl$component[1], (k - 1), sep = "_"), idx_decomposed_point[k])
          component_updated <- c(component_updated, tmp_decompose)
        }

        link_cl$component_stitching <- component_updated
      }

      link_per_cluster[[j]] <- link_cl
    }

    link_per_cluster <- do.call(rbind, link_per_cluster)

    link_info_list[[i]] <- link_per_cluster

    Sys.sleep(1 / round(length(link_info_list) / 10))
  }

  pb$terminate()

  print("Re-merging step... (Sparsity correction)")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(link_info_list),
                                   clear = FALSE,
                                   width= 65)

  for(i in 1:length(link_info_list)) {
    pb$tick()
    chr_info <- names(link_info_list)[i]

    test_env <- link_info_list[[chr_info]]
    test_env$midpoint <- as.numeric(test_env$midpoint)
    test_env <- test_env %>% dplyr::arrange(midpoint)
    test_env <- split(test_env, f = test_env$component)

    test_env <- lapply(test_env, remerge_component, chr_info = chr_info, thres_cos = thres_cos_re_merge, results_gABC_list = results_gABC_list)
    test_env <- do.call(rbind, test_env)
    link_info_list[[i]] <- test_env

    Sys.sleep(1 / round(length(link_info_list) / 10))
  }

  pb$terminate()

  link_info_list <- do.call(rbind, link_info_list)
  link_info_list$topic <- peak_for_gABC[link_info_list$peak_id, n_celltype]
  rownames(link_info_list) <- link_info_list$peak_id

  ### add additional node
  link_info_total_corr <- link_info_total_n

  link_info_total_n$query_id <- rownames(peak_info)[link_info_total_n$queryHits]
  link_info_total_n$subject_id <- rownames(peak_info)[link_info_total_n$subjectHits]

  link_info_total_n <- link_info_total_n[which(link_info_total_n$correlation > thres_corr_add), ]
  link_info_total_n <- link_info_total_n[which((link_info_total_n$query_id %in% link_info_list$peak_id) | (link_info_total_n$subject_id %in% link_info_list$peak_id)), ]

  additional_nodes_cand <- setdiff(candidate_enh, candidate_enh_u)
  additional_nodes <- c()

  print("Check additional node candidates...")

  pb <- progress::progress_bar$new(format = " Progress: [:bar] :percent, ETA: :eta",
                                   total = length(additional_nodes_cand),
                                   clear = FALSE,
                                   width= 65)

  for(i in 1:length(additional_nodes_cand)) {
    pb$tick()
    tmp_add_cand <- c(link_info_total_n$subject_id[which(link_info_total_n$query_id == additional_nodes_cand[i])],
                      link_info_total_n$query_id[which(link_info_total_n$subject_id == additional_nodes_cand[i])])
    tmp_add_cand <- setdiff(unique(tmp_add_cand), additional_nodes_cand[i])

    if(length(intersect(tmp_add_cand, link_info_list$peak_id)) > 0) {
      additional_nodes <- c(additional_nodes, additional_nodes_cand[i])
    }

    Sys.sleep(1 / round(length(additional_nodes_cand) / 10))
  }

  pb$terminate()

  thres_topic_add_u <- thres_topic_add / 100

  thres_topic <- sort(link_info_list$topic, decreasing = TRUE)[thres_topic_add_u*nrow(link_info_list)]

  additional_nodes_u <- additional_nodes[which(peak_for_addcorr[additional_nodes, n_celltype] >= thres_topic)]

  additional_info_list <- as.list(additional_nodes_u)
  names(additional_info_list) <- additional_nodes_u

  print("Add nodes with high accessibility...")

  for(i in names(additional_info_list)) {
    tmp_add <- c(link_info_total_n$subject_id[which(link_info_total_n$query_id == additional_info_list[[i]])],
                 link_info_total_n$query_id[which(link_info_total_n$subject_id == additional_info_list[[i]])])
    tmp_add <- setdiff(unique(tmp_add), additional_info_list[[i]])
    tmp_add_id <- tmp_add

    tmp_add <- link_info_list[which(link_info_list$peak_id %in% tmp_add), ]

    tmp_mid_add <- stringr::str_split(additional_info_list[[i]], pattern = "_")
    tmp_mid_add <- do.call(rbind, tmp_mid_add)
    tmp_mid_add <- as.numeric(tmp_mid_add[1, 2]) + 250

    tmp_add$dist_add <- abs(tmp_add$midpoint - tmp_mid_add)

    tmp_test_n <- link_info_total_n[which((link_info_total_n$query_id == additional_info_list[[i]]) | (link_info_total_n$subject_id == additional_info_list[[i]])), ]
    tmp_test_n <- tmp_test_n[which(tmp_test_n$query_id %in% tmp_add_id), ]
    tmp_test_n <- data.frame(peak_id = tmp_test_n$query_id, corr = tmp_test_n$correlation)
    rownames(tmp_test_n) <- tmp_test_n$peak_id

    tmp_add$corr <- tmp_test_n[tmp_add$peak_id, "corr"]

    tmp_add$score_assign <- (tmp_add$corr * tmp_add$topic) / tmp_add$dist_add

    tmp_add_id <- tmp_add$component_stitching[which(tmp_add$score_assign == max(tmp_add$score_assign))]

    if(length(tmp_add_id) == 0) {
      tmp_add_id_add <- NA
    } else {
      tmp_add_id_add <- tmp_add_id
    }

    additional_info_list[[i]] <- tmp_add_id_add
  }

  link_info_list_add <- data.frame(peak_id = names(additional_info_list),
                                   component_stitching = unlist(additional_info_list))

  link_info_list_add <- na.omit(link_info_list_add)

  if(nrow(link_info_list_add) == 0) {
    link_info_list <- link_info_list
  } else {
    link_info_list_add$component <- "additional"
    tmp_peak_id_info_add <- stringr::str_split(link_info_list_add$peak_id, pattern = "_")
    tmp_peak_id_info_add <- do.call(rbind, tmp_peak_id_info_add)
    link_info_list_add$chr <- tmp_peak_id_info_add[, 1]
    link_info_list_add$start <- as.numeric(tmp_peak_id_info_add[, 2])
    link_info_list_add$end <- as.numeric(tmp_peak_id_info_add[, 3])
    link_info_list_add$difference <- "additional"
    link_info_list_add$midpoint <- link_info_list_add$start + 250
    link_info_list_add$dist_diff <- "additional"
    link_info_list_add$topic <- peak_for_addcorr[link_info_list_add$peak_id, n_celltype]

    link_info_list_add <- link_info_list_add[, colnames(link_info_list)]
    link_info_list <- rbind(link_info_list, link_info_list_add)
  }

  ### write total final enhancer candidates
  write.table(link_info_list[, c("chr", "start", "end")], paste0(output_dir, "/total_enhc_", n_celltype, ".bed"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

  print("Write total final enhancer candidates...")

  link_info_total_corr <- link_info_total_corr[!duplicated(t(apply(link_info_total_corr[, c(1, 2)], 1, sort))), ]
  link_info_total_corr <- link_info_total_corr[, c("queryHits", "subjectHits", "seqnames")]

  within_rnet_corr <- sort(unique(as.integer(c(link_info_total_corr$queryHits, link_info_total_corr$subjectHits))))
  in_silico_self_loops_corr <- data.frame(queryHits = within_rnet_corr, subjectHits = within_rnet_corr,
                                          seqnames = peak_info$chr[within_rnet_corr])

  link_info_total_corr <- rbind(link_info_total_corr, in_silico_self_loops_corr)
  link_info_total_corr <- link_info_total_corr %>% dplyr::arrange(queryHits, subjectHits)

  rownames(link_info_total_corr) <- 1:length(link_info_total_corr$queryHits)

  link_info_total_corr$query_id <- rownames(peak_info)[link_info_total_corr$queryHits]
  link_info_total_corr$subject_id <- rownames(peak_info)[link_info_total_corr$subjectHits]

  link_info_total_corr <- link_info_total_corr[which((link_info_total_corr$query_id %in% link_info_list$peak_id) & (link_info_total_corr$subject_id %in% link_info_list$peak_id)), ]

  network_reg <- link_info_total_corr[, c("query_id", "subject_id", "queryHits", "subjectHits")]
  network_reg <- subset(network_reg, subset = ((query_id %in% link_info_list$peak_id) & (subject_id %in% link_info_list$peak_id)))

  # distance
  data_coacc <- as.data.frame(info_coacc@listData)
  network_reg <- dplyr::left_join(network_reg, data_coacc, by = c("queryHits", "subjectHits"))

  network_reg <- network_reg[, c("query_id", "subject_id", "correlation")]
  network_reg$correlation <- as.numeric(network_reg$correlation)

  peak_for_nreg <- peak_for_addcorr[, c(1:(length(colnames(peak_for_addcorr))))]
  idx_cts <- which(colnames(peak_for_nreg) == n_celltype)

  options(digits=22)

  network_reg_1 <- network_reg[!(is.na(network_reg$correlation)), ]
  network_reg_1 <- network_reg_1[, c(2, 1, 3)]
  colnames(network_reg_1) <- c("query_id", "subject_id", "correlation")

  network_reg <- rbind(network_reg, network_reg_1)

  print("Pseudo-directed construction!")
  print("[0] Assign edge weights...")

  topic_sum_comp <- link_info_list %>% dplyr::group_by(component_stitching) %>% dplyr::summarise(topic_sum = sum(topic))
  topic_sum_comp <- as.data.frame(topic_sum_comp)
  rownames(topic_sum_comp) <- topic_sum_comp$component_stitching

  network_reg$component_q <- link_info_list[network_reg$query_id, "component_stitching"]
  network_reg$component_s <- link_info_list[network_reg$subject_id, "component_stitching"]

  network_reg$tmp_weight_0 <- as.numeric(peak_for_nreg[network_reg$query_id, n_celltype])
  network_reg$tmp_weight_1 <- as.numeric(peak_for_nreg[network_reg$subject_id, n_celltype])

  network_reg$tmp_weight_r0 <- topic_sum_comp[network_reg$component_q, "topic_sum"]
  network_reg$tmp_weight_r1 <- topic_sum_comp[network_reg$component_s, "topic_sum"]

  tmp_reg_0 <- stringr::str_split(network_reg$query_id, pattern = "_")
  tmp_reg_1 <- stringr::str_split(network_reg$subject_id, pattern = "_")
  tmp_reg_0 <- do.call(rbind, tmp_reg_0)
  tmp_reg_1 <- do.call(rbind, tmp_reg_1)

  network_reg$distance_bw <- abs(as.numeric(tmp_reg_0[, 2]) - as.numeric(tmp_reg_1[, 2]))

  network_reg$tmp_weight_0 <- (network_reg$tmp_weight_0 * network_reg$tmp_weight_1 * network_reg$tmp_weight_r1) / (network_reg$tmp_weight_r0 + network_reg$tmp_weight_r1)

  network_reg$correlation <- network_reg$correlation * network_reg$tmp_weight_0 / network_reg$distance_bw

  network_reg[is.na(network_reg$correlation), "correlation"] <- 1
  network_reg <- network_reg[, c("query_id", "subject_id", "correlation")]
  colnames(network_reg) <- c("query_id", "subject_id", "weight")
  #####

  network_reg <- igraph::graph.data.frame(network_reg, directed = TRUE)
  network_reg <- igraph::simplify(network_reg)

  weight_motifs <- ArchR::getMatrixFromProject(
    ArchRProj = proj_atac,
    useMatrix = "MotifMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
  )

  group_info <- weight_motifs$Clusters2
  group_info <- which(group_info == n_celltype)
  weight_motifs <- weight_motifs@assays@data$z
  weight_motifs <- weight_motifs[, group_info]

  weight_motifs <- as.matrix(weight_motifs)
  tmp_row_id <- rownames(weight_motifs)
  weight_motifs <- rowMedians(weight_motifs)

  motif_row_name <- c()
  for(i in tmp_row_id) {
    motif_row_name <- c(motif_row_name, motif_database@listData[[i]]@name)
  }

  names(weight_motifs) <- motif_row_name

  if(length(which(weight_motifs >= thres_TF_weight)) == 0) {
    stop("Please check your thres_TF_weight value!")
  }

  weight_motifs[weight_motifs < thres_TF_weight] <- 0
  test_weight <- weight_motifs

  peaks_for_link <- GenomicRanges::makeGRangesFromDataFrame(link_info_list[, c("chr", "start", "end")])
  results_motifmatchr <- motifmatchr::matchMotifs(motif_database, peaks_for_link, genome = organism, bg = "genome", out = "scores")
  results_motifmatchr <- motifmatchr::motifMatches(results_motifmatchr)
  results_motifmatchr <- as.matrix(results_motifmatchr)
  results_motifmatchr <- as.data.frame(results_motifmatchr)

  colnames(results_motifmatchr) <- motif_row_name

  results_motifmatchr[results_motifmatchr == TRUE] <- 1
  results_motifmatchr[results_motifmatchr == FALSE] <- 0

  results_motifmatchr <- results_motifmatchr[, intersect(colnames(results_motifmatchr), names(test_weight))]

  motif_weights <- rep(0, length(rownames(results_motifmatchr)))
  names(motif_weights) <- rownames(results_motifmatchr)

  results_motifmatchr <- as.matrix(results_motifmatchr)
  real_col_id <- colnames(results_motifmatchr)
  test_weight <- test_weight[real_col_id]

  for(i in 1:length(results_motifmatchr[, 1])) {
    motif_weights[i] <- sum(results_motifmatchr[i, ] * test_weight)
  }

  link_info_list$motif_weights <- motif_weights

  print("[1] Assign node weights...")
  link_info_list$topic_updated <- 100 * scales::rescale(link_info_list$topic, to = c(0, 1)) * scales::rescale(link_info_list$motif_weights, to = c(0, 1))

  topic_weight <- link_info_list[names(igraph::V(network_reg)), "topic_updated"]

  pagerank_topic <- igraph::page.rank(network_reg, directed = TRUE, personalized = topic_weight)
  pagerank_topic <- as.data.frame(pagerank_topic$vector)

  colnames(pagerank_topic) <- "topic_sensitive_pagerank"

  test_comp <- link_info_list
  test_comp$topic_sensitive_pagerank <- pagerank_topic[test_comp$peak_id, "topic_sensitive_pagerank"]

  test_comp <- test_comp %>% dplyr::group_by(component_stitching) %>% dplyr::summarise(importance = sum(topic_sensitive_pagerank))
  test_comp$importance <- 100 * (test_comp$importance / max(test_comp$importance))
  test_comp <- test_comp %>% dplyr::arrange(importance)

  plot(1:nrow(test_comp), test_comp$importance, xlab = "rank", ylab = "summarized centrality")

  test <- test_comp
  test$component_stitching <- NULL
  colnames(test) <- "y"
  test$x <- 1:nrow(test_comp)
  test <- test[, c("x", "y")]

  idx_elbow <- kneedle::kneedle(test$x, test$y, sensitivity = 2)
  idx_elbow <- idx_elbow[1]
  idx_elbow <- ceiling(idx_elbow)
  idx_ref <- nrow(test_comp) - floor(max_nCORE) + 1

  if(idx_elbow >= idx_ref) {
    idx_elbow <- idx_elbow
  } else {
    tmp_idx_ref <- numeric()
    for(i in c(3:max_iter_elbow)) {
      tmp_idx_ref_i <- kneedle::kneedle(test$x, test$y, sensitivity = i)

      if(tmp_idx_ref_i[1] >= idx_ref) {
        tmp_idx_ref <- c(tmp_idx_ref, tmp_idx_ref_i[1])
        break
      }
    }

    idx_elbow <- tmp_idx_ref[length(tmp_idx_ref)]
  }

  test_list <- subset(link_info_list, subset = component_stitching %in% test_comp$component_stitching[idx_elbow:nrow(test_comp)])
  test_list <- split(test_list, f = test_list$component_stitching)

  for(i in 1:length(test_list)) {
    test_list[[i]] <- extract_ranges_for_bed(test_list[[i]], min_count = min_count)
  }

  test_list <- do.call(rbind, test_list)
  test_list <- na.omit(test_list)

  write.table(test_list, paste0(output_dir, "/CORE_potential_", n_celltype, "_i.bed"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  ### sort .bed file
  system(paste0("bedtools sort -i ", paste0(output_dir, "/CORE_potential_", n_celltype, "_i.bed"), " > ",
                paste0(output_dir, "/CORE_potential_", n_celltype, "_i_sorted.bed")))

  ### merge overlapped regions
  system(paste0("bedtools merge -i ", paste0(output_dir, "/CORE_potential_", n_celltype, "_i_sorted.bed"), " > ",
                paste0(output_dir, "/CORE_potential_", n_celltype, "_f.bed")))

  ### remove initial .bed files
  system(paste0("rm ", paste0(output_dir, "/CORE_potential_", n_celltype, "_i.bed")))
  system(paste0("rm ", paste0(output_dir, "/CORE_potential_", n_celltype, "_i_sorted.bed")))

  print("Write CORE from the potential option!")
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
  data_lump_enCORE[["output_dir_CORE"]] <- output_dir
  data_lump_enCORE[["n_core"]] <- n_core
  data_lump_enCORE[["organism"]] <- organism
  data_lump_enCORE[["names_order"]] <- names_order
  data_lump_enCORE[["list_cluster"]] <- tmp_list_cluster
  data_lump_enCORE[["threshold_TF_weight"]] <- thres_TF_weight

  return(data_lump_enCORE)
}

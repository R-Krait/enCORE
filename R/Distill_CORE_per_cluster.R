#' Distill_CORE_per_cluster
#'
#' This function will distill CORE regions per cell population in list_cluster.
#'
#' @param proj_atac An ArchRProject object containing chromVAR CISBP v2 motif deviation.
#' @param data_lump_enCORE A list containing the data required for running enCORE.
#' @param output_dir Output directory where extracted CORE regions will be saved.
#' @param option There are two options in enCORE framwork, potential and active. If option == "potential", distill only raw CORE regions. If option == "active", distill both raw and iterative-proximal-filtered CORE regions.
#' @param thres_cos_initial Threshold for the cosine similarity to extract interaction-shared decomposed components. Default = 0.4.
#' @param thres_cos_link_loss Threshold for the cosine similarity to perform Link Loss Correction. Default = 0.4.
#' @param thres_cos_re_merge Threshold for the cosine similarity to perform Sparsity Correction. Default = 0.92.
#' @param thres_corr_add Threshold for the co-accessibility value to incorporate additional nodes with high accessibility. Default = 0.49.
#' @param thres_topic_add Threshold for the accessibility to incorporate additional nodes with high accessibility. Default value is top 0.5% accessibility of original enhancer candidates. In other words, Default = 0.5.
#' @param max_dist Maximum effective distance between enhancer candidates. Default = 500000.
#' @param max_nCORE The maximum number of CORE regions. Default = 0.1*26000/1.31.
#' @param max_iter_elbow The maximum number of iterations for elbow point determination. Default = 10.
#' @param min_count The minimum CPM count for CORE constituents. Only constituents with a CPM above this threshold will be retained.
#' @param thres_inactive Threshold for the accessibility within promoter to determine inactive genes. Default = 12.5.
#' @param thres_eprox Effective distance of inactive proximal enhancers. Default = 7500 bp.
#' @export
Distill_CORE_per_cluster <- function(proj_atac, data_lump_enCORE, output_dir, option = "potential",
                                     thres_cos_initial = 0.4, thres_cos_link_loss = 0.4,
                                     thres_cos_re_merge = 0.92, thres_corr_add = 0.49, thres_topic_add = 0.5, max_dist = 500000,
                                     max_nCORE = 0.1*26000/1.31, max_iter_elbow = 10, min_count = 5.0,
                                     thres_inactive = 12.5, thres_eprox = 7500) {
  data_lump_enCORE[["output_dir_CORE"]] <- output_dir
  list_cluster <- data_lump_enCORE[["list_cluster"]]

  for(n_celltype in list_cluster) {
    data_lump_enCORE_cp <- data_lump_enCORE
    if(option == "potential") {
      data_lump_enCORE_cp <- extract_raw_CORE_per_cluster(proj_atac = proj_atac, data_lump_enCORE = data_lump_enCORE_cp, n_celltype = n_celltype,
                                                          thres_cos_initial = thres_cos_initial, thres_cos_link_loss = thres_cos_link_loss,
                                                          thres_cos_re_merge = thres_cos_re_merge, thres_corr_add = thres_corr_add, thres_topic_add = thres_topic_add,
                                                          max_dist = max_dist, max_nCORE = max_nCORE, max_iter_elbow = max_iter_elbow, min_count = min_count)
    } else if(option == "active") {
      data_lump_enCORE_cp <- extract_raw_CORE_per_cluster(proj_atac = proj_atac, data_lump_enCORE = data_lump_enCORE_cp, n_celltype = n_celltype,
                                                          thres_cos_initial = thres_cos_initial, thres_cos_link_loss = thres_cos_link_loss,
                                                          thres_cos_re_merge = thres_cos_re_merge, thres_corr_add = thres_corr_add, thres_topic_add = thres_topic_add,
                                                          max_dist = max_dist, max_nCORE = max_nCORE, max_iter_elbow = max_iter_elbow, min_count = min_count)

      data_lump_enCORE_cp <- iterative_proximal_filtering(data_lump_enCORE = data_lump_enCORE_cp, n_celltype = n_celltype,
                                                          thres_inactive = thres_inactive, thres_eprox = thres_eprox)
    } else {
      stop("There are only two options, potential or active!")
    }
  }

  return(data_lump_enCORE)
}

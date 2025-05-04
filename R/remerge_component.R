#' remerge_component
#'
#' This function will re-merge decomposed components with high cosine similarity to perform Sparsity Correction.
#'
#' @param env_T An element of the list containing co-accessibility network edges split by decomposed component ID after Link Loss Correction.
#' @param chr_info Chromosome information.
#' @param results_gABC_list results_gABC_list after generating decomposed components.
#' @param thres_cos Threshold for the cosine similarity. Default = 0.92.
#' @export
remerge_component <- function(env_T, chr_info, results_gABC_list, thres_cos = 0.92) {
  list_comp <- env_T$component_stitching[!duplicated(env_T$component_stitching)]

  if(length(list_comp) == 1) {
    env_T <- env_T
  } else {
    test_list <- rep(0, length(list_comp))

    for(i in 1:(length(list_comp) - 1)) {
      test_comp_0 <- list_comp[i]
      test_comp_1 <- list_comp[(i + 1)]

      cent_0 <- calculate_centroid_remerge(env_T, test_comp_0, chr_info, results_gABC_list)
      cent_1 <- calculate_centroid_remerge(env_T, test_comp_1, chr_info, results_gABC_list)

      metric_cosine <- coop::cosine(cent_0, cent_1)

      if(metric_cosine >= thres_cos) {
        next
      }

      test_list[(i + 1)] <- 1
    }

    idx_cliff <- which(test_list == 1)

    if(length(idx_cliff) == 0) {
      if(sum(test_list) == 0) {
        prefix <- stringr::str_split(list_comp[1], pattern = "_")
        prefix <- do.call(rbind, prefix)
        prefix <- prefix[1, c(1:(length(prefix[1, ]) - 1))]
        prefix <- paste(prefix, collapse = "_")

        tmp_test <- rep(paste0(prefix, "_", 0), nrow(env_T))
        env_T$component_stitching <- tmp_test
      } else {
        env_T <- env_T
      }
    } else {
      idx_cliff <- c(idx_cliff[1] - 1, diff(idx_cliff), (length(test_list) - idx_cliff[length(idx_cliff)] + 1))

      list_comp_u <- c()
      prefix <- stringr::str_split(list_comp[1], pattern = "_")
      prefix <- do.call(rbind, prefix)
      prefix <- prefix[1, c(1:(length(prefix[1, ]) - 1))]
      prefix <- paste(prefix, collapse = "_")

      for(j in 1:length(idx_cliff)) {
        tmp_test <- rep(paste0(prefix, "_", (j - 1)), idx_cliff[j])
        list_comp_u <- c(list_comp_u, tmp_test)
      }

      names(list_comp_u) <- list_comp

      for(k in 1:nrow(env_T)) {
        env_T$component_stitching[k] <- as.character(list_comp_u[env_T$component_stitching[k]])
      }
    }
  }

  return(env_T)
}

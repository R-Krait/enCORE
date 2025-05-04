#' calculate_centroid_remerge
#'
#' This function will calculate centroid of gABC scores for decomposed components to perform Sparsity Correction.
#'
#' @param env_T An element of the list containing co-accessibility network edges split by decomposed component ID after Link Loss Correction.
#' @param x Query or Subject decomposed component to calculate centroid of gABC scores.
#' @param chr_info Chromosome information.
#' @param results_gABC_list results_gABC_list after generating decomposed components.
#' @export
calculate_centroid_remerge <- function(env_T, x, chr_info, results_gABC_list) {
  x <- env_T$peak_id[env_T$component_stitching == x]
  x <- results_gABC_list[[chr_info]][x]
  x <- do.call(cbind, x)

  x <- rowSums(x)

  return(x)
}

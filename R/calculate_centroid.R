#' calculate_centroid
#'
#' This function will calculate centroid of gABC scores for enhancer candidates.
#'
#' @param p An element of the list to which the centroid information will be added.
#' @param chr_info Chromosome information.
#' @param results_gABC_list results_gABC_list after valid vector generation.
#' @export
calculate_centroid <- function(p, chr_info, results_gABC_list) {
  p <- names(igraph::V(p))
  p <- results_gABC_list[[chr_info]][p]
  p <- do.call(cbind, p)

  p <- rowSums(p)

  return(p)
}

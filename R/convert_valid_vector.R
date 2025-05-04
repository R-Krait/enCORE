#' convert_valid_vector
#'
#' This function will convert gABC scores to valid vector to perform enCORE.
#'
#' @param x results_gABC_list element
#' @param y total_gene_id element
#' @export
convert_valid_vector <- function(x, y) {
  z <- rep(0, length(y))
  names(z) <- y
  z[x$gene_name] <- x$ABC_score

  return(z)
}

#' split_by_peak_id
#'
#' This function will split each element of the results_gABC_list by peak_id.
#'
#' @param x results_gABC_list element.
#' @export
split_by_peak_id <- function(x) {
  y <- split(x, f = x$peak_id)

  return(y)
}

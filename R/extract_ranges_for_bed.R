#' extract_ranges_for_bed
#'
#' This function will convert the extracted CORE regions to .bed file.
#'
#' @param x Each raw CORE region in test_list variable.
#' @param min_count The minimum CPM count for CORE constituents. Only constituents with a CPM above this threshold will be retained.
#' @export
extract_ranges_for_bed <- function(x, min_count = 5.0) {
  x$start <- as.integer(x$start)
  x <- x %>% dplyr::arrange(start)

  idx_ab <- which(x$topic > min_count)

  if(length(idx_ab) == 0) {
    x <- data.frame(chr = NA, start = NA, end = NA)
  } else {
    x <- x[idx_ab, ]
    x <- data.frame(chr = x$chr[1], start = x$start[1], end = x$end[nrow(x)])
  }

  return(x)
}

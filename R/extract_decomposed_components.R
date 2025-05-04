#' extract_decomposed_components
#'
#' This function will execute Link Loss Correction.
#'
#' @param link_per_chr A list containing co-accessibility network edges split by chromosome information.
#' @param results_gABC_list results_gABC_list after assigning cosine similarity to co-accessibility network edge.
#' @param thres_cos Threshold for the cosine similarity. Default = 0.4.
#' @param max_dist Maximum effective distance between enhancer candidates. Default = 500000.
#' @param n_core An integer specifying the number of CPU-cores to perform enCORE. Default = 4.
#' @export
extract_decomposed_components <- function(link_per_chr, results_gABC_list, thres_cos = 0.4, max_dist = 500000, n_core = 4) {
  chr_info <- as.character(link_per_chr$seqnames[1])
  link_per_chr <- link_per_chr[, c("query_id", "subject_id", "ABC_similarity")]
  colnames(link_per_chr) <- c("V1", "V2", "weight")
  link_per_chr <- link_per_chr[which(!(link_per_chr$weight == 0)), ]

  network_per_chr <- igraph::graph.data.frame(link_per_chr, directed = FALSE)

  # decompose
  decomposed_components <- igraph::decompose(network_per_chr)

  centroid_list <- decomposed_components

  centroid_list <- parallel::mclapply(centroid_list, calculate_centroid, chr_info = chr_info,
                                      results_gABC_list = results_gABC_list, mc.cores = n_core)

  cosine_cent <- data.frame(query = 1:(length(centroid_list) - 1), idx_subject = NA)

  for(i in 1:nrow(cosine_cent)) {
    list_candidate <- names(igraph::V(decomposed_components[[cosine_cent$query[i]]]))
    list_candidate <- stringr::str_split(list_candidate, pattern = "_")
    list_candidate <- do.call(rbind, list_candidate)
    list_candidate <- as.numeric(list_candidate[, 3])
    list_candidate <- max(list_candidate) + max_dist

    if(i == nrow(cosine_cent)) {
      cosine_cent$idx_subject[i] <- cosine_cent$query[i]
    } else {
      check_cand <- numeric()
      pt_start <- i + 1
      for(j in pt_start:nrow(cosine_cent)) {
        subject_start <- names(igraph::V(decomposed_components[[cosine_cent$query[j]]]))
        subject_start <- stringr::str_split(subject_start, pattern = "_")
        subject_start <- do.call(rbind, subject_start)
        subject_start <- as.numeric(subject_start[, 2])
        subject_start <- min(subject_start)

        if(subject_start <= list_candidate) {
          check_cand <- c(check_cand, j)
        } else {
          break
        }
      }

      if(length(check_cand) == 0) {
        cosine_cent$idx_subject[i] <- cosine_cent$query[i]
      } else {
        subject_cand <- do.call(cbind, centroid_list[cosine_cent$query[check_cand]])
        cosine_object <- cbind(centroid_list[[cosine_cent$query[i]]], subject_cand)
        colnames(cosine_object) <- as.character(c(cosine_cent$query[i], cosine_cent$query[check_cand]))

        cosine_object <- coop::cosine(cosine_object)

        cosine_object <- cosine_object[as.character(cosine_cent$query[i]), as.character(cosine_cent$query[check_cand])]
        idx_subject_tmp <- which(cosine_object >= thres_cos)

        if(length(idx_subject_tmp) == 0) {
          cosine_cent$idx_subject[i] <- cosine_cent$query[i]
        } else {
          idx_subject_tmp <- max(idx_subject_tmp)
          if(length(cosine_object) == 1) {
            cosine_cent$idx_subject[i] <- cosine_cent$query[as.numeric(as.character(cosine_cent$query[check_cand])[idx_subject_tmp])]
          } else {
            cosine_cent$idx_subject[i] <- cosine_cent$query[as.numeric(names(cosine_object)[idx_subject_tmp])]
          }
        }
      }
    }
  }

  idx_merged_point <- data.frame(start = as.numeric(cosine_cent$query),
                                 end = as.numeric(cosine_cent$idx_subject),
                                 pt_inf = c(1, rep(0, (nrow(cosine_cent) - 1))))

  save_point <- numeric()
  for(i in 1:nrow(idx_merged_point)) {
    if(length(save_point) != 0) {
      if(i < max(save_point)) {
        next
      }
    }

    if(i == nrow(idx_merged_point)) {
      break
    } else {
      checkpoint <- i + 1
      savepoint <- i
      tmp_save <- idx_merged_point$end[i]
      for(j in checkpoint:nrow(idx_merged_point)) {
        if(tmp_save < idx_merged_point$start[j]) {
          break
        } else {
          savepoint <- j
          tmp_save <- max(tmp_save, idx_merged_point$end[j])
        }
      }

      if((savepoint + 1) > nrow(idx_merged_point)) {
        break
      } else {
        idx_merged_point$pt_inf[(savepoint + 1)] <- 1
        save_point <- c(save_point, (savepoint + 1))
      }
    }
  }

  idx_mp <- which(idx_merged_point$pt_inf == 1)

  idx_merged_point_u <- data.frame(idx_q = idx_mp, idx_s = (idx_mp + diff(c(idx_mp, (nrow(idx_merged_point) + 1))) - 1),
                                   start_u = NA, end_u = NA)

  for(i in 1:nrow(idx_merged_point_u)) {
    idx_list <- seq(idx_merged_point_u$idx_q[i], idx_merged_point_u$idx_s[i], 1)
    idx_list <- as.numeric(c(idx_merged_point$start[idx_list], idx_merged_point$end[idx_list]))

    idx_merged_point_u$start_u[i] <- min(idx_list)
    idx_merged_point_u$end_u[i] <- max(idx_list)
  }

  idx_merged_point_u <- idx_merged_point_u[, c(3, 4)]

  idx_merged_point_u$mc <- paste("m", 1:nrow(idx_merged_point_u), sep = "_")

  cluster_info <- 1:length(decomposed_components)
  names(cluster_info) <- cluster_info

  for(i in 1:nrow(idx_merged_point_u)) {
    id_m <- seq(idx_merged_point_u$start_u[i], idx_merged_point_u$end_u[i])
    cluster_info[id_m] <- idx_merged_point_u$mc[i]
  }

  changed_point <- decomposed_components

  for(i in 1:length(changed_point)) {
    tmp <- data.frame(component = i, peak_id = names(igraph::V(changed_point[[i]])))
    changed_point[[i]] <- tmp
  }

  changed_point <- do.call(rbind, changed_point)
  id_split <- stringr::str_split(changed_point$peak_id, pattern = "_")
  id_split <- do.call(rbind, id_split)

  changed_point$chr <- id_split[, 1]
  changed_point$start <- id_split[, 2]
  changed_point$end <- id_split[, 3]

  changed_point <- changed_point[order(as.numeric(changed_point$start), decreasing = FALSE), ]
  changed_point$difference <- 0
  changed_point$midpoint <- as.numeric(changed_point$start) + 250

  for(i in 1:nrow(changed_point)) {
    changed_point$component[i] <- cluster_info[changed_point$component[i]]
  }

  for(i in 1:nrow(changed_point)) {
    if(i != 1) {
      if(changed_point$component[i] != changed_point$component[i - 1]) {
        changed_point$difference[i] <- 1
      }
    }
  }

  changed_point$start <- as.integer(changed_point$start)
  changed_point$end <- as.integer(changed_point$end)
  decomposed_components_updated <- split(changed_point, f = changed_point$component)

  results_decompose <- as.list(c(1, 2, 3))
  names(results_decompose) <- c("changed_point", "decomposed_components", "cluster_info")

  results_decompose[["changed_point"]] <- changed_point
  results_decompose[["decomposed_components"]] <- decomposed_components_updated
  results_decompose[["cluster_info"]] <- cluster_info

  return(results_decompose)
}

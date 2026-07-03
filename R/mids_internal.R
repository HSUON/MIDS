#' @importFrom magrittr %>%
run_mids_with_alpha <- function(df_group,
                                frame_col,
                                length_col,
                                precision_col,
                                alpha_used,
                                compare_to = c("earliest", "latest")) {

  compare_to <- match.arg(compare_to)

  df_group <- df_group %>%
    dplyr::arrange(.data[[frame_col]], local_id)

  n <- nrow(df_group)

  if (n < 2) {
    return(df_group %>% dplyr::mutate(individual_id = local_id))
  }

  frames <- sort(unique(df_group[[frame_col]]))
  assigned_rows <- list()

  first_frame <- frames[1]

  df_first <- df_group %>%
    dplyr::filter(.data[[frame_col]] == first_frame) %>%
    dplyr::arrange(local_id) %>%
    dplyr::mutate(individual_id = seq_len(dplyr::n()))

  assigned_rows[[as.character(first_frame)]] <- df_first
  next_individual_id <- max(df_first$individual_id) + 1

  if (length(frames) > 1) {
    for (fr in frames[-1]) {

      current <- df_group %>%
        dplyr::filter(.data[[frame_col]] == fr) %>%
        dplyr::arrange(local_id)

      previous_all <- dplyr::bind_rows(assigned_rows)

      # Full measurement history for each known individual
      known_history <- previous_all %>%
        dplyr::select(
          individual_id,
          previous_local_id = local_id,
          previous_Frame = dplyr::all_of(frame_col),
          prev_Length = dplyr::all_of(length_col),
          prev_PE = dplyr::all_of(precision_col)
        )

      # Compare each current fish to all previous observations of each known individual
      pair_tests <- tidyr::expand_grid(
        individual_id = unique(known_history$individual_id),
        local_id = current$local_id
      ) %>%
        dplyr::left_join(known_history, by = "individual_id") %>%
        dplyr::left_join(
          current %>%
            dplyr::select(
              local_id,
              dplyr::all_of(length_col),
              dplyr::all_of(precision_col)
            ) %>%
            dplyr::rename(
              cur_Length = dplyr::all_of(length_col),
              cur_PE = dplyr::all_of(precision_col)
            ),
          by = "local_id"
        ) %>%
        dplyr::mutate(
          z = abs(prev_Length - cur_Length) / sqrt(prev_PE^2 + cur_PE^2),
          p = 2 * (1 - stats::pnorm(z)),
          non_distinct = p >= alpha_used
        )

      # A current fish can only match an existing individual if it overlaps
      # with all previous observations assigned to that individual
      candidate_pairs <- pair_tests %>%
        dplyr::group_by(individual_id, local_id) %>%
        dplyr::summarise(
          all_history_matches = all(non_distinct),
          min_p = min(p),
          .groups = "drop"
        ) %>%
        dplyr::filter(all_history_matches)

      if (nrow(candidate_pairs) > 0) {

        vertices <- tibble::tibble(
          name = c(
            paste0("ind_", unique(previous_all$individual_id)),
            paste0("fish_", current$local_id)
          ),
          type = c(
            rep(FALSE, length(unique(previous_all$individual_id))),
            rep(TRUE, nrow(current))
          )
        )

        edges <- candidate_pairs %>%
          dplyr::transmute(
            from = paste0("ind_", individual_id),
            to = paste0("fish_", local_id)
          )

        g <- igraph::graph_from_data_frame(
          d = edges,
          directed = FALSE,
          vertices = vertices
        )

        mm <- igraph::max_bipartite_match(g)

        matched <- tibble::tibble(
          from = names(mm$matching),
          to = as.character(mm$matching)
        ) %>%
          dplyr::filter(!is.na(to)) %>%
          dplyr::filter(grepl("^ind_", from) & grepl("^fish_", to)) %>%
          dplyr::transmute(
            individual_id = as.integer(sub("^ind_", "", from)),
            local_id = as.integer(sub("^fish_", "", to))
          )

      } else {
        matched <- tibble::tibble(
          individual_id = integer(),
          local_id = integer()
        )
      }

      current_assigned <- current %>%
        dplyr::left_join(matched, by = "local_id")

      unmatched_n <- sum(is.na(current_assigned$individual_id))

      if (unmatched_n > 0) {
        current_assigned$individual_id[is.na(current_assigned$individual_id)] <-
          seq(from = next_individual_id, length.out = unmatched_n)

        next_individual_id <- next_individual_id + unmatched_n
      }

      assigned_rows[[as.character(fr)]] <- current_assigned
    }
  }

  dplyr::bind_rows(assigned_rows) %>%
    dplyr::arrange(.data[[frame_col]], local_id)
}

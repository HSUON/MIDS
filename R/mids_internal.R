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

      known <- previous_all %>%
        dplyr::group_by(individual_id) %>%
        {
          if (compare_to == "earliest") {
            dplyr::arrange(., .data[[frame_col]], local_id)
          } else {
            dplyr::arrange(., dplyr::desc(.data[[frame_col]]), dplyr::desc(local_id))
          }
        } %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(
          individual_id,
          dplyr::all_of(length_col),
          dplyr::all_of(precision_col)
        ) %>%
        dplyr::rename(
          Length = dplyr::all_of(length_col),
          Precision_Error = dplyr::all_of(precision_col)
        )

      candidate_pairs <- tidyr::expand_grid(
        individual_id = known$individual_id,
        local_id = current$local_id
      ) %>%
        dplyr::left_join(
          known %>%
            dplyr::rename(
              prev_Length = Length,
              prev_PE = Precision_Error
            ),
          by = "individual_id"
        ) %>%
        dplyr::left_join(
          current %>%
            dplyr::select(
              local_id,
              dplyr::all_of(length_col),
              dplyr::all_of(precision_col)
            ) %>%
            dplyr::rename(
              Length = dplyr::all_of(length_col),
              Precision_Error = dplyr::all_of(precision_col)
            ) %>%
            dplyr::rename(
              cur_Length = Length,
              cur_PE = Precision_Error
            ),
          by = "local_id"
        ) %>%
        dplyr::mutate(
          z = abs(prev_Length - cur_Length) / sqrt(prev_PE^2 + cur_PE^2),
          p = 2 * (1 - stats::pnorm(z))
        ) %>%
        dplyr::filter(p >= alpha_used)

      if (nrow(candidate_pairs) > 0) {

        vertices <- tibble::tibble(
          name = c(
            paste0("ind_", known$individual_id),
            paste0("fish_", current$local_id)
          ),
          type = c(
            rep(FALSE, nrow(known)),
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


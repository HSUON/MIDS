utils::globalVariables(c(
  "row_i",
  "row_j",
  "frame_i",
  "frame_j",
  "length_i",
  "length_j",
  "precision_i",
  "precision_j",
  "combined_uncertainty",
  "same_frame",
  "z",
  "p",
  "incompatible",
  "local_id",
  "local_id_i",
  "local_id_j",
  "individual_id",
  "individual_i",
  "individual_j",
  "n_in_frame",
  "Before",
  "After"
))
#' Exact minimum graph colouring
#'
#' Finds the minimum number of colours required so that adjacent vertices
#' receive different colour assignments. Used internally by MIDS to assign
#' the minimum number of inferred individuals consistent with all
#' incompatibility constraints.
#'
#' @param graph An igraph graph object.
#' @param time_limit_seconds Maximum search time in seconds for one group.
#'
#' @return A list containing:
#' \describe{
#'   \item{colors}{A named integer vector of colour assignments.}
#'   \item{n_colors}{The exact minimum number of colours.}
#' }
#'
#' @keywords internal
#' @noRd
exact_minimum_coloring <- function(graph,
                                   time_limit_seconds = 600) {

  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.", call. = FALSE)
  }

  if (!is.numeric(time_limit_seconds) ||
      length(time_limit_seconds) != 1L ||
      !is.finite(time_limit_seconds) ||
      time_limit_seconds <= 0) {
    stop(
      "`time_limit_seconds` must be one positive finite number.",
      call. = FALSE
    )
  }

  n_vertices <- igraph::vcount(graph)
  vertex_names <- igraph::V(graph)$name

  if (n_vertices == 0L) {
    return(
      list(
        colors = stats::setNames(integer(), character()),
        n_colors = 0L
      )
    )
  }

  if (n_vertices == 1L) {
    return(
      list(
        colors = stats::setNames(1L, vertex_names),
        n_colors = 1L
      )
    )
  }

  start_time <- proc.time()[["elapsed"]]

  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(
      graph,
      sparse = FALSE
    )
  )

  storage.mode(adjacency) <- "logical"
  diag(adjacency) <- FALSE

  vertex_degree <- rowSums(adjacency)

  # Lower bound: every vertex in the largest clique requires a separate ID.
  lower_bound <- max(
    1L,
    as.integer(igraph::clique_num(graph))
  )

  # Upper bound only. This greedy result is not used as the final assignment.
  greedy_colors <- igraph::greedy_vertex_coloring(
    graph,
    heuristic = "dsatur"
  )

  upper_bound <- max(as.integer(greedy_colors))

  elapsed_time <- function() {
    proc.time()[["elapsed"]] - start_time
  }

  # Test whether the graph can be coloured using at most k colours.
  solve_k_coloring <- function(k) {

    colors <- integer(n_vertices)

    recursive_search <- function() {

      if (elapsed_time() > time_limit_seconds) {
        stop(
          "Exact MIDS colouring exceeded the time limit of ",
          time_limit_seconds,
          " seconds for one group.",
          call. = FALSE
        )
      }

      uncolored <- which(colors == 0L)

      if (length(uncolored) == 0L) {
        return(TRUE)
      }

      # DSATUR-style vertex selection:
      # highest colour saturation, followed by highest graph degree.
      saturation <- vapply(
        uncolored,
        function(vertex) {

          neighbor_colors <- colors[
            adjacency[vertex, ] & colors > 0L
          ]

          length(unique(neighbor_colors))
        },
        integer(1)
      )

      selected_position <- order(
        -saturation,
        -vertex_degree[uncolored],
        uncolored
      )[1L]

      vertex <- uncolored[selected_position]

      forbidden_colors <- unique(
        colors[adjacency[vertex, ] & colors > 0L]
      )

      available_colors <- setdiff(
        seq_len(k),
        forbidden_colors
      )

      if (length(available_colors) == 0L) {
        return(FALSE)
      }

      used_colors <- sort(
        unique(colors[colors > 0L])
      )

      existing_options <- intersect(
        used_colors,
        available_colors
      )

      unused_options <- setdiff(
        available_colors,
        used_colors
      )

      colors_to_try <- existing_options

      # Only one unused colour label needs to be tested because unused colour
      # labels are otherwise interchangeable.
      if (length(unused_options) > 0L) {
        colors_to_try <- c(
          colors_to_try,
          min(unused_options)
        )
      }

      for (candidate_color in colors_to_try) {

        colors[vertex] <<- candidate_color

        if (recursive_search()) {
          return(TRUE)
        }

        colors[vertex] <<- 0L
      }

      FALSE
    }

    feasible <- recursive_search()

    if (feasible) {
      colors
    } else {
      NULL
    }
  }

  # The first feasible number of colours is the exact minimum.
  for (k in seq.int(lower_bound, upper_bound)) {

    solution <- solve_k_coloring(k)

    if (!is.null(solution)) {

      names(solution) <- vertex_names

      return(
        list(
          colors = solution,
          n_colors = as.integer(k)
        )
      )
    }
  }

  stop(
    "Exact MIDS colouring failed to return a valid solution.",
    call. = FALSE
  )
}

#' Run MIDS using one fixed alpha
#'
#' Evaluates every unique pair of observations within one grouped comparison
#' unit and determines the exact minimum number of inferred individuals.
#'
#' @param df_group Data frame for one grouped comparison unit.
#' @param frame_col Name of the frame column.
#' @param length_col Name of the length column.
#' @param precision_col Name of the precision-error column.
#' @param alpha_used Fixed significance threshold.
#' @param time_limit_seconds Maximum exact-colouring search time.
#'
#' @return The original observations with inferred individual IDs.
#'
#' @keywords internal
run_mids_with_alpha <- function(df_group,
                                frame_col,
                                length_col,
                                precision_col,
                                alpha_used,
                                time_limit_seconds = 600) {

  if (!is.numeric(alpha_used) ||
      length(alpha_used) != 1L ||
      !is.finite(alpha_used) ||
      alpha_used <= 0 ||
      alpha_used >= 1) {
    stop(
      "`alpha_used` must be one finite value between 0 and 1.",
      call. = FALSE
    )
  }

  required_columns <- c(
    frame_col,
    length_col,
    precision_col,
    "local_id"
  )

  missing_columns <- setdiff(
    required_columns,
    names(df_group)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(df_group[required_columns])) {
    stop(
      "Missing values were detected in a required MIDS column.",
      call. = FALSE
    )
  }

  if (!is.numeric(df_group[[length_col]]) ||
      !is.numeric(df_group[[precision_col]])) {
    stop(
      "The length and precision-error columns must be numeric.",
      call. = FALSE
    )
  }

  if (any(!is.finite(df_group[[length_col]])) ||
      any(!is.finite(df_group[[precision_col]]))) {
    stop(
      "Length and precision-error values must be finite.",
      call. = FALSE
    )
  }

  if (any(df_group[[precision_col]] <= 0)) {
    stop(
      "All precision-error values must be greater than zero.",
      call. = FALSE
    )
  }

  df_group <- df_group[
    order(
      df_group[[frame_col]],
      df_group$local_id
    ),
    ,
    drop = FALSE
  ]

  n <- nrow(df_group)

  if (n == 0L) {
    df_group$individual_id <- integer()
    df_group$exact_n_individuals <- integer()
    df_group$alpha_used <- numeric()

    return(df_group)
  }

  if (n == 1L) {
    df_group$individual_id <- 1L
    df_group$exact_n_individuals <- 1L
    df_group$alpha_used <- alpha_used

    return(df_group)
  }

  # Every unique pair of observations within this camera/species group.
  pair_index <- t(
    utils::combn(
      seq_len(n),
      2
    )
  )

  pair_tests <- tibble::tibble(
    row_i = pair_index[, 1],
    row_j = pair_index[, 2]
  ) %>%
    dplyr::mutate(
      local_id_i = df_group$local_id[row_i],
      local_id_j = df_group$local_id[row_j],

      frame_i = df_group[[frame_col]][row_i],
      frame_j = df_group[[frame_col]][row_j],

      length_i = df_group[[length_col]][row_i],
      length_j = df_group[[length_col]][row_j],

      precision_i = df_group[[precision_col]][row_i],
      precision_j = df_group[[precision_col]][row_j],

      same_frame = frame_i == frame_j,

      combined_uncertainty = sqrt(
        precision_i^2 + precision_j^2
      ),

      z = abs(length_i - length_j) /
        combined_uncertainty,

      p = 2 * stats::pnorm(-abs(z)),

      # Same-frame observations must be separate individuals.
      #
      # For cross-frame observations:
      # p < alpha_used  = statistically distinct
      # p >= alpha_used = potentially the same individual
      incompatible =
        same_frame |
        (!same_frame & p < alpha_used)
    )

  # An edge means the two observations cannot share an individual ID.
  incompatibility_edges <- pair_tests %>%
    dplyr::filter(incompatible) %>%
    dplyr::transmute(
      from = as.character(local_id_i),
      to = as.character(local_id_j)
    )

  if (nrow(incompatibility_edges) == 0L) {
    incompatibility_edges <- tibble::tibble(
      from = character(),
      to = character()
    )
  }

  vertices <- tibble::tibble(
    name = as.character(df_group$local_id)
  )

  incompatibility_graph <- igraph::graph_from_data_frame(
    d = incompatibility_edges,
    directed = FALSE,
    vertices = vertices
  )

  exact_solution <- exact_minimum_coloring(
    graph = incompatibility_graph,
    time_limit_seconds = time_limit_seconds
  )

  identity_table <- tibble::tibble(
    local_id = as.integer(
      names(exact_solution$colors)
    ),
    individual_id = as.integer(
      unname(exact_solution$colors)
    )
  )

  result <- df_group %>%
    dplyr::left_join(
      identity_table,
      by = "local_id"
    ) %>%
    dplyr::mutate(
      exact_n_individuals = exact_solution$n_colors,
      alpha_used = alpha_used
    )

  result <- result[
    order(
      result[[frame_col]],
      result$local_id
    ),
    ,
    drop = FALSE
  ]

  # -----------------------------------------------------------------------
  # Internal validation
  # -----------------------------------------------------------------------

  if (anyNA(result$individual_id)) {
    stop(
      "Internal MIDS error: at least one observation received no ID.",
      call. = FALSE
    )
  }

  # Same-frame observations must never share an ID.
  same_frame_conflicts <- result %>%
    dplyr::count(
      .data[[frame_col]],
      individual_id,
      name = "n_in_frame"
    ) %>%
    dplyr::filter(n_in_frame > 1L)

  if (nrow(same_frame_conflicts) > 0L) {
    stop(
      "Internal MIDS error: observations in the same frame received ",
      "the same individual ID.",
      call. = FALSE
    )
  }

  # Statistically incompatible observations must never share an ID.
  incompatible_conflicts <- pair_tests %>%
    dplyr::filter(incompatible) %>%
    dplyr::left_join(
      identity_table %>%
        dplyr::rename(
          local_id_i = local_id,
          individual_i = individual_id
        ),
      by = "local_id_i"
    ) %>%
    dplyr::left_join(
      identity_table %>%
        dplyr::rename(
          local_id_j = local_id,
          individual_j = individual_id
        ),
      by = "local_id_j"
    ) %>%
    dplyr::filter(individual_i == individual_j)

  if (nrow(incompatible_conflicts) > 0L) {
    stop(
      "Internal MIDS error: an incompatible observation pair received ",
      "the same individual ID.",
      call. = FALSE
    )
  }

  result
}

#' Run internal two-pass MIDS
#'
#' @param df_group Data frame for one grouped comparison unit.
#' @param frame_col Name of the frame column.
#' @param length_col Name of the length column.
#' @param precision_col Name of the precision-error column.
#' @param alpha Initial significance threshold.
#' @param use_post_bonferroni Whether to apply the second-pass adjustment.
#' @param time_limit_seconds Maximum exact-colouring search time per pass.
#'
#' @return All observations with final individual IDs and threshold metadata.
#'
#' @keywords internal
run_mids_internal <- function(df_group,
                              frame_col,
                              length_col,
                              precision_col,
                              alpha = 0.05,
                              use_post_bonferroni = TRUE,
                              time_limit_seconds = 600) {

  # Pass 1: standard retrospective assignment.
  first_pass <- run_mids_with_alpha(
    df_group = df_group,
    frame_col = frame_col,
    length_col = length_col,
    precision_col = precision_col,
    alpha_used = alpha,
    time_limit_seconds = time_limit_seconds
  )

  first_pass_n_distinct <- dplyr::n_distinct(
    first_pass$individual_id
  )

  # Second-pass threshold based on the first-pass inferred abundance.
  if (use_post_bonferroni &&
      first_pass_n_distinct > 0L) {

    alpha_second_pass <-
      alpha / first_pass_n_distinct

  } else {

    alpha_second_pass <- alpha
  }

  # Pass 2: repeat the full retrospective assignment.
  final_pass <- run_mids_with_alpha(
    df_group = df_group,
    frame_col = frame_col,
    length_col = length_col,
    precision_col = precision_col,
    alpha_used = alpha_second_pass,
    time_limit_seconds = time_limit_seconds
  )

  final_pass %>%
    dplyr::mutate(
      first_pass_n_distinct = first_pass_n_distinct,
      alpha_first_pass = alpha,
      alpha_second_pass = alpha_second_pass,
      final_n_distinct = dplyr::n_distinct(individual_id),
      assignment_method =
        "Exact retrospective minimum colouring"
    )
}

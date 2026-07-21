#' Run MIDS abundance estimation
#'
#' Applies Multi-frame Individual Distinction by Size (MIDS) to stereo-video
#' observations to estimate the minimum number of distinct individuals
#' consistent with all pairwise length comparisons and same-frame
#' co-occurrence constraints.
#'
#' MIDS is performed retrospectively within each user-defined group. Every
#' unique observation pair is evaluated simultaneously, avoiding assumptions
#' that earlier or later measurements are inherently more reliable.
#'
#' When `use_post_bonferroni = TRUE`, MIDS is first run at `alpha`. The number
#' of first-pass inferred individuals is then used to calculate one fixed
#' second-pass threshold, and the complete retrospective analysis is repeated.
#'
#' @param data A data frame of fish observations.
#' @param group_cols Character vector of grouping columns defining independent
#'   comparison units. These should ordinarily include deployment, camera or
#'   video identifier, and species.
#' @param frame_col Character string naming the frame or time column.
#' @param length_col Character string naming the fish-length column.
#' @param precision_col Character string naming the measurement precision-error
#'   column.
#' @param alpha Initial significance threshold used for pairwise length
#'   comparisons. Default is `0.05`.
#' @param use_post_bonferroni Logical. If `TRUE`, MIDS is rerun using a fixed
#'   second-pass threshold calculated as `alpha` divided by the number of
#'   first-pass inferred individuals.
#' @param time_limit_seconds Maximum time allowed for exact graph colouring
#'   within one grouped comparison unit and pass. Default is 300 seconds.
#'
#' @return A data frame containing one representative row for each inferred
#'   individual. The earliest observation is retained only as the output row;
#'   frame order does not determine identity.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
run_mids <- function(data,
                     group_cols,
                     frame_col,
                     length_col,
                     precision_col,
                     alpha = 0.05,
                     use_post_bonferroni = TRUE,
                     time_limit_seconds = 300) {

  # -----------------------------------------------------------------------
  # Validate arguments
  # -----------------------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  character_arguments <- list(
    group_cols = group_cols,
    frame_col = frame_col,
    length_col = length_col,
    precision_col = precision_col
  )

  if (!is.character(group_cols) ||
      length(group_cols) < 1L ||
      anyNA(group_cols) ||
      any(group_cols == "")) {
    stop(
      "`group_cols` must contain at least one valid column name.",
      call. = FALSE
    )
  }

  for (argument_name in c(
    "frame_col",
    "length_col",
    "precision_col"
  )) {

    argument_value <- character_arguments[[argument_name]]

    if (!is.character(argument_value) ||
        length(argument_value) != 1L ||
        is.na(argument_value) ||
        argument_value == "") {
      stop(
        "`",
        argument_name,
        "` must be one valid column name.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop(
      "`alpha` must be one finite value between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.logical(use_post_bonferroni) ||
      length(use_post_bonferroni) != 1L ||
      is.na(use_post_bonferroni)) {
    stop(
      "`use_post_bonferroni` must be TRUE or FALSE.",
      call. = FALSE
    )
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

  required_cols <- unique(
    c(
      group_cols,
      frame_col,
      length_col,
      precision_col
    )
  )

  missing_cols <- setdiff(
    required_cols,
    names(data)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "The following columns are missing from `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(data[required_cols])) {
    stop(
      "Missing values were detected in required MIDS columns.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[length_col]]) ||
      !is.numeric(data[[precision_col]])) {
    stop(
      "The length and precision-error columns must be numeric.",
      call. = FALSE
    )
  }

  if (any(!is.finite(data[[length_col]])) ||
      any(!is.finite(data[[precision_col]]))) {
    stop(
      "Length and precision-error values must be finite.",
      call. = FALSE
    )
  }

  if (any(data[[precision_col]] <= 0)) {
    stop(
      "All precision-error values must be greater than zero.",
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # Prepare data
  # -----------------------------------------------------------------------

  mids_prepared <- data %>%
    dplyr::arrange(
      dplyr::across(
        dplyr::all_of(
          c(
            group_cols,
            frame_col,
            length_col,
            precision_col
          )
        )
      )
    ) %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(group_cols)
      )
    ) %>%
    dplyr::mutate(
      local_id = dplyr::row_number()
    ) %>%
    dplyr::ungroup()

  # -----------------------------------------------------------------------
  # Run two-pass retrospective MIDS independently within each group
  # -----------------------------------------------------------------------

  mids_assignments <- mids_prepared %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(group_cols)
      )
    ) %>%
    dplyr::group_modify(
      ~ run_mids_internal(
        df_group = .x,
        frame_col = frame_col,
        length_col = length_col,
        precision_col = precision_col,
        alpha = alpha,
        use_post_bonferroni = use_post_bonferroni,
        time_limit_seconds = time_limit_seconds
      )
    ) %>%
    dplyr::ungroup()

  # -----------------------------------------------------------------------
  # Retain one output record per inferred individual
  # -----------------------------------------------------------------------
  #
  # Identity has already been determined retrospectively.
  # The earliest row is retained only as the representative output row.

  mids_distinct <- mids_assignments %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          c(
            group_cols,
            "individual_id"
          )
        )
      )
    ) %>%
    dplyr::arrange(
      .data[[frame_col]],
      local_id,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  mids_distinct
}

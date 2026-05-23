#' Run MIDS abundance estimation
#'
#' Applies Multi-frame Individual Distinction by Size (MIDS) to stereo-video
#' observations to estimate distinct individuals within user-defined groups.
#'
#' @param data A data frame of fish observations.
#' @param group_cols Character vector of grouping columns. These should define
#'   the independent comparison unit, usually deployment/camera/video and species.
#' @param frame_col Column containing frame number or time.
#' @param length_col Column containing fish length.
#' @param precision_col Column containing measurement precision error.
#' @param alpha Significance threshold used for length-overlap comparisons.
#' @param compare_to Whether to compare new observations to the earliest or latest
#'   measurement of a previously assigned individual.
#' @param use_post_bonferroni Logical. If TRUE, reruns MIDS using a Bonferroni
#'   correction based on the first-pass number of distinct individuals.
#'
#' @return A data frame of distinct individual records.
#' @export
#' @importFrom magrittr %>%
run_mids <- function(data,
                     group_cols,
                     frame_col,
                     length_col,
                     precision_col,
                     alpha = 0.05,
                     compare_to = c("earliest", "latest"),
                     use_post_bonferroni = TRUE) {

  compare_to <- match.arg(compare_to)

  required_cols <- c(group_cols, frame_col, length_col, precision_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing from `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  mids_prepared <- data %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(required_cols))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(local_id = dplyr::row_number()) %>%
    dplyr::ungroup()

  mids_assignments <- mids_prepared %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::group_modify(
      ~ {
        first_pass <- run_mids_with_alpha(
          df_group = .x,
          frame_col = frame_col,
          length_col = length_col,
          precision_col = precision_col,
          alpha_used = alpha,
          compare_to = compare_to
        )

        if (use_post_bonferroni) {
          n_distinct_individuals <- dplyr::n_distinct(first_pass$individual_id)
          alpha_used <- if (n_distinct_individuals > 0) {
            alpha / n_distinct_individuals
          } else {
            alpha
          }
        } else {
          alpha_used <- alpha
        }

        run_mids_with_alpha(
          df_group = .x,
          frame_col = frame_col,
          length_col = length_col,
          precision_col = precision_col,
          alpha_used = alpha_used,
          compare_to = compare_to
        )
      }
    ) %>%
    dplyr::ungroup()

  mids_distinct <- mids_assignments %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "individual_id")))) %>%
    dplyr::arrange(.data[[frame_col]], local_id, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  mids_distinct
}

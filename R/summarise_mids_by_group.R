#' Summarise MIDS results by group
#'
#' Summarises original observations, retained MIDS individuals, and removed
#' repeat observations within user-defined groups.
#'
#' @param original_data Original input data frame.
#' @param mids_results Output from run_mids().
#' @param group_cols Character vector of grouping columns.
#'
#' @return A tibble with before, after, and removed counts by group.
#' @export
summarise_mids_by_group <- function(original_data, mids_results, group_cols) {
  
  before <- original_data |>
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "Before")
  
  after <- mids_results |>
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "After")
  
  before |>
    dplyr::left_join(after, by = group_cols) |>
    dplyr::mutate(
      After = tidyr::replace_na(After, 0L),
      Removed = Before - After
    )
}
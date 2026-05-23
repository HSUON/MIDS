#' Summarise MIDS results
#'
#' Summarises the number of original observations, retained MIDS individuals,
#' and removed repeat observations.
#'
#' @param original_data Original input data frame.
#' @param mids_results Output from run_mids().
#'
#' @return A tibble with before, after and removed counts.
#' @export
summarise_mids <- function(original_data, mids_results) {
  tibble::tibble(
    Before = nrow(original_data),
    After = nrow(mids_results),
    Removed = Before - After
  )
}


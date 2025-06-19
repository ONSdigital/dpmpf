
#' Create, run, and destroy a particle filter
#'
#' Create a particle filter, run it, write the output to disk, then destroy
#' the particle filter. The particle filter deals with a single
#' cohort (or, more precisely, a single combination of cohort
#' and sex/gender.)
#'
#' \code{create_run_destroy_pfilter} is an internal function,
#' and end-users should never need to call it directly.
#' \code{create_run_destroy_pfilter} does not do any
#' validity-checking on its arguments: error checking
#' is the responsibility of the calling functions.
#'
#' An individual particle filter can hold millions of
#' values, so we avoid creating multiple particle
#' filters at the same time.
#'
#' \code{create_run_destroy_pfilter} writes files to directory
#' \code{work_dir}.
#'
#' WARNING - THIS DESCRIPTION NEEDS TO BE UPDATED
#' The first file is called \code{tmp-population-<cohort>-<sexgender>.csv},
#' where values for \code{<cohort>} and \code{<sexgender>} are
#' taken from \code{df}. If \code{df} includes a region dimension,
#' then the file has columns \code{"cohort"}, \code{"sexgender"},
#' \code{"time"}, \code{"age"}, \code{"region"},
#' \code{"iteration"}, and \code{"count"}. If \code{df} does not
#' include a region dimension, then the file does not include a
#' \code{"region"} column.
#'
#' Next there is a set of files describing events, each of which
#' has a filename of the form \code{tmp-<event>-<cohort>-<sexgender>.csv}.
#' If \code{df} has a region dimension, and if the value for
#' \code{is_dominant} for \code{df} is \code{TRUE},
#' then \code{<event>} consists of \code{"births"},
#' \code{"deaths"}, \code{"internal_in"}, \code{"internal_out"},
#' \code{"immigration1"}, \code{"emigration1"},
#' \code{"immigration2"}, and \code{"emigration2"}. If
#' If \code{df} does not include a region dimension,
#' then \code{<event>} does not include \code{"internal_in"}
#' or \code{"internal_out"}. If the value for
#' \code{is_dominant} for \code{df} is \code{FALSE},
#' then \code{<event>} does not include \code{"births"}.
#'
#' If \code{df} includes a region dimension,
#' then the file has variables \code{"cohort"}, \code{"sexgender"},
#' \code{"time"}, \code{"age"}, \code{"iteration"},
#' and \code{"count"}. If \code{df} does not include a region dimension,
#' then the file does not include variable \code{"region"}.
#'
#' The final file is called \code{tmp-diagnostic-<cohort>-<sexgender>.csv}.
#' The file has columns \code{"cohort"}, \code{"sexgender"},
#' \code{"time"}, \code{"age"}, \code{"triangle"}, \code{"ess"},
#' \code{"resampled"}, and \code{"n_unique"}.
#'
#' @param df A data frame holding inputs for an individual cohort.
#' @param threshold Threshold for triggering resampling.
#'   A number between 0 and 1.
#' @param work_dir The name of the folder where
#'   intermediate calculations are carried out and
#'   results are stored.
#' @param is_forecast Whether the particle filter
#'   is being used for a forecast.
#' @param quiet Whether to show a message saying calculations
#'   have started.
#'
#' @return Return \code{TRUE} invisibly and
#' writes the files described above to \code{work_dir}.
#'
#' @export
create_run_destroy_pfilter <- function(df,
                                       threshold,
                                       work_dir,
                                       is_forecast,
                                       quiet) {
    if (!quiet)
        show_start_message(df)
    pf <- make_pf(
        df = df,
        threshold = threshold,
        is_forecast = is_forecast
    )
    pf$run(is_forecast = is_forecast)
    pf$write_results(
           is_forecast = is_forecast,
           work_dir = work_dir
       )
    invisible(TRUE)
}




#' Do a forecast for a single cohort
#'
#' Do a forecast for a single cohort. We draw straight from the
#' posterior distribution, so we don't need to use particle filtering.
#'
#' \code{do_forecast} is an internal function,
#' and end-users should never need to call it directly.
#'
#' @inheritParams create_run_destroy_pfilter
#'
#' @return Return \code{TRUE} invisibly and
#' writes the files described above to \code{work_dir}.
#'
#' @export
do_forecast <- function(df,
                        work_dir,
                        quiet) {
    if (!quiet)
        show_start_message(df)
    pf <- make_pf(
        df = df,
        threshold = 0,
        is_forecast = TRUE
    )
    pf$forecast_cohort()
    pf$write_results(
           is_forecast = TRUE,
           work_dir = work_dir
       )
    invisible(TRUE)
}

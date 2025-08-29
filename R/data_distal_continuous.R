#' A synthetic data set of an MRT with continuous distal outcome
#'
#' Simulated longitudinal dataset suitable for illustrating the `dcee()` function.
#' Each row corresponds to one decision point for one subject.
#' The distal outcome `Y` is constant within subject (because it is measured
#' at the end of the study, and here we append it to the long format data
#' as an extra column to conform with the `dcee()` function requirement.
#'
#' @format a data frame with 1500 observations and 11 variables
#' \describe{
#'   \item{userid}{Subject identifier}
#'   \item{dp}{Decision point (1..T)}
#'   \item{X}{Endogenous continuous time-varying covariate}
#'   \item{Z}{Endogenous binary time-varying covariate}
#'   \item{avail}{Availability indicator (0/1)}
#'   \item{A}{Treatment (0/1)}
#'   \item{prob_A}{Randomization probability P(A=1|H_t)}
#'   \item{A_lag1}{Lagged treatment}
#'   \item{Y}{Distal continuous outcome (constant per subject)}
#' }
"data_distal_continuous"

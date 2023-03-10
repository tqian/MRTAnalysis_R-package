#' A synthetic data set of an MRT with binary outcome
#'
#' @format a data frame with 3000 observations and 10 variables
#' \describe{ This random sample uses the baseline model:
#'            log E(Y_{t+1} | A_t = 0, I_t = 1) = alpha_0 + alpha_1 * time / total_T
#'            + alpha_2 * 1(time > total_T/2),
#'            the treatment effect model:
#'            log relative risk = beta_0 + beta_1 * time / total_T,
#'            the probability of treatment assignment p_t: 0.3, 0.5, 0.7 with repetition,
#'            and exogenous probability of availability: 0.8 at all time points.
#'    \item{userid}{individual id number}
#'    \item{time}{decision point index}
#'    \item{time_var1}{time-varying covariate 1, the "standardized time in study", defined as the current decision point index divided by the total number of decision points}
#'    \item{time_var2}{time-varying covariate 2, indicator of "the second half of the study", defined as whether the current decision point index is greater than the total number of decision points divided by 2.}
#'    \item{Y}{binary proximal outcome}
#'    \item{A}{treatment assignment, i.e., whether the intervention is randomized to be delivered (=1) or not (=0) at the current decision point}
#'    \item{rand_prob}{the randomization probability P(A=1) for the current decision point}
#'    \item{avail}{whether the individual is available (=1) or not (=0) at the current decision point}
#' }
"data_binary"

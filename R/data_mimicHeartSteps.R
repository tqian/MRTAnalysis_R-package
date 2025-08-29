#' A synthetic data set that mimics the HeartSteps V1 data structure to illustrate
#' the use of [wcls()] function for continuous proximal outcomes
#'
#' @format a data frame with 7770 observations and 9 variables
#' \describe{
#'     \item{userid}{individual id number}
#'     \item{time}{decision point index}
#'     \item{day_in_study}{day in the study}
#'     \item{logstep_30min}{proximal outcome: the step count in the 30 minutes following the current decision point (log-transformed)}
#'     \item{logstep_30min_lag1}{proximal outcome at the previous decision point (lag-1 outcome): the step count in the 30 minutes following the previous decision point (log-transformed)}
#'     \item{logstep_pre30min}{the step count in the 30 minutes prior to the current decision point (log-transformed); used as a control variable}
#'     \item{is_at_home_or_work}{whether the individual is at home or work (=1) or at other locations (=0) at the current decision point}
#'     \item{intervention}{whether the intervention is randomized to be delivered (=1) or not (=0) at the current decision point}
#'    \item{rand_prob}{the randomization probability P(A=1) for the current decision point}
#'     \item{availability}{whether the individual is available (=1) or not (=0) at the current decision point}
#'     }
"data_mimicHeartSteps"

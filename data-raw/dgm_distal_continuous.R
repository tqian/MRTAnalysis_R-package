# code to prepare data_distal_continuous goes here

# -----------------------------------------------------------------------------
# Data-generating mechanism for distal continuous outcomes in MRTs
# -----------------------------------------------------------------------------
# Features
# (1) X_t is endogenous (depends on past A)
# (2) Y depends on X_t in nonlinear ways
# (3) Interactions present: X_t * A_t and A_t * A_{t-1}
# (4) Additional endogenous binary covariate Z_t
#
# Notes
# - Returns a long-format data.frame with one row per (id, decision point t)
# - The distal outcome Y is constant within each subject id (same value across t)
# - No external dependencies (base R only)
# - Randomization probability can be constant or time-varying via covariates and time
#
# Columns returned
#   userid   : subject id (1..sample_size)
#   dp       : decision point (1..total_T)
#   X        : endogenous continuous covariate at time t
#   Z        : endogenous binary covariate at time t
#   avail    : availability indicator at time t (0/1)
#   A        : treatment at time t (0/1), only drawn when avail == 1
#   prob_A   : P(A_t = 1 | H_t)
#   A_lag1   : lagged treatment A_{t-1}
#   Y        : distal continuous outcome (constant within subject)
#
# This is dgm_distal_outcome_endox_nonlinear_interaction_binaryZ.R from
# Tianchen's DCEE paper simulation, with some simplications (removing the
# "force_A" argument.)
# -----------------------------------------------------------------------------

library(tidyverse)

dgm_distal_continuous <- function(sample_size, total_T,
                                  rand_prob_pattern = c("constant", "timevar"),
                                  availability = c("always", "not-always") # always 1 or not always 1
) {
    expit <- function(x) {
        return(1 / (1 + exp(-x)))
    }

    const_rand_prob <- 0.5
    rand_prob_tuning_param <- 2

    alpha_vec <- seq(from = 1, to = 3, length.out = total_T)
    beta_vec <- seq(from = 1, to = 2, length.out = total_T)
    gamma_vec <- seq(from = 1, to = 1.5, length.out = total_T)
    lambda_vec <- seq(from = -1, to = -2, length.out = total_T)
    xi_vec <- seq(from = 1, to = 2, length.out = total_T)
    theta0 <- -0.5
    theta1 <- 0.5
    theta2 <- 0.5
    zeta0 <- -1
    zeta1 <- 1
    zeta2 <- 1

    # Note:
    # X_t = theta0 + theta1 A_t-1 + theta2 X_t-1 + eta_t
    # Z_t = Bern(expit(zeta0 + zeta1 A_t-1 + zeta2 Z_t-1))
    # Y = sum_t xi_t {g_t(X_t) + Z_t} + sum_t A_t(alpha_t + beta_t X_t + gamma_t Z_t + lambda_t A_t-1) + eps

    rand_prob_pattern <- match.arg(rand_prob_pattern)
    availability <- match.arg(availability)

    df_names <- c("userid", "dp", "X", "Z", "avail", "A", "prob_A", "A_lag1")
    # Note: Even though data is in long format and has Y on every row, Y is
    #       a distal outcome and thus takes the same value across all rows within a person

    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names

    dta$userid <- rep(1:sample_size, each = total_T)
    dta$dp <- rep(1:total_T, times = sample_size)

    # vcov_error <- exponential_vcov(error_tuning_param, total_T)
    # dta$eps <- as.vector(rmvnorm(sample_size, sigma = vcov_error))

    ## Generate X_t, Z_t, and A_t
    Z_lag1 <- rep(0, sample_size)
    X_lag1 <- rep(0, sample_size)
    A_lag1 <- rep(0, sample_size)
    prob_A_lag1 <- rep(0, sample_size)

    for (t in 1:total_T) {
        # row index for the rows corresponding to decision point (dp) t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)

        dp <- dta$dp[row_index]

        ## generate X: X_t = theta0 + theta1 A_t-1 + theta2 X_t-1 + eta_t
        eta <- rnorm(sample_size)
        X <- theta0 + theta1 * A_lag1 + theta2 * X_lag1 + eta

        ## generate Z: Z_t = Bern(expit(zeta0 + zeta1 A_t-1 + zeta2 Z_t-1))
        pZ <- expit(zeta0 + zeta1 * A_lag1 + zeta2 * Z_lag1)
        Z <- rbinom(sample_size, 1, pZ)

        ## generate avail
        if (availability == "always") {
            avail <- rep(1, length(row_index))
        } else if (availability == "not-always") {
            avail <- rbinom(length(row_index), 1, 0.8)
        }

        ## generate A
        if (rand_prob_pattern == "constant") {
            prob_A <- rep(const_rand_prob, sample_size)
        } else if (rand_prob_pattern == "timevar") {
            X_transformed <- X / 6 # it seems that X is mostly symmetric and doesn't exceed 6
            dp_transformed <- (dp - total_T / 2) / total_T
            prob_A <- expit(rand_prob_tuning_param * X_transformed +
                rand_prob_tuning_param * (Z - 0.5) +
                rand_prob_tuning_param * dp_transformed)
            prob_A <- pmax(pmin(prob_A, 0.9), 0.1)
        }
        A <- avail * rbinom(sample_size, 1, prob = prob_A)

        ## Put variables into dta
        dta$X[row_index] <- X
        dta$Z[row_index] <- Z
        dta$avail[row_index] <- avail
        dta$A[row_index] <- A
        dta$prob_A[row_index] <- prob_A
        dta$A_lag1[row_index] <- A_lag1

        # Gather lagged variables to be used
        X_lag1 <- X
        Z_lag1 <- Z
        A_lag1 <- A
        prob_A_lag1 <- prob_A
    }

    dta_Y <- dta %>%
        group_by(userid) %>%
        # Y = sum_t xi_t {g_t(X_t) + Z_t} + sum_t A_t(alpha_t + beta_t X_t + gamma_t Z_t + lambda_t A_t-1) + eps
        summarize(
            expectY =
                sum(xi_vec * (dbeta(X / 12 + 0.5, 2, 2) + Z)) +
                    sum(A * (alpha_vec + beta_vec * X + gamma_vec * Z + lambda_vec * A_lag1))
        ) %>%
        ungroup()
    dta_Y$epsY <- rnorm(sample_size)
    dta_Y$Y <- dta_Y$expectY + dta_Y$epsY

    dta <- dta %>% full_join(dta_Y, by = "userid")

    dta <- dta %>% select(-expectY, -epsY)

    return(dta)
}

set.seed(123)
data_distal_continuous <- dgm_distal_continuous(
    sample_size = 50,
    total_T = 30,
    rand_prob_pattern = "timevar",
    availability = "not-always"
)
usethis::use_data(data_distal_continuous, overwrite = TRUE)

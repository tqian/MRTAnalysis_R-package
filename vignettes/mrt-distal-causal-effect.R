## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(MRTAnalysis)
current_options <- options(digits = 3) # save current options for restoring later
head(data_distal_continuous, 10)

## -----------------------------------------------------------------------------
fit_lm <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~1,
    control_formula = ~ X + Z,
    availability = "avail",
    control_reg_method = "lm"
)
summary(fit_lm)

## -----------------------------------------------------------------------------
fit_mod <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~Z,
    control_formula = ~ Z + X,
    availability = "avail",
    control_reg_method = "lm"
)
summary(fit_mod, lincomb = c(1, 1)) # beta0 + beta1

## -----------------------------------------------------------------------------
fit_gam <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~Z,
    control_formula = ~ s(X) + Z,
    availability = "avail",
    control_reg_method = "gam"
)
summary(fit_gam)

## -----------------------------------------------------------------------------
fit_rf <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~1,
    control_formula = ~ X + Z,
    availability = "avail",
    control_reg_method = "rf" # can replace "rf" with "ranger" for faster implementation
)
summary(fit_rf)

## ----eval = requireNamespace("SuperLearner", quietly = TRUE)------------------
library(SuperLearner)
fit_sl <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~1,
    control_formula = ~ X + Z,
    availability = "avail",
    control_reg_method = "sl"
)
summary(fit_sl)

## -----------------------------------------------------------------------------
fit_cf <- dcee(
    data = data_distal_continuous,
    id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
    moderator_formula = ~1,
    control_formula = ~ X + Z,
    availability = "avail",
    control_reg_method = "gam",
    cross_fit = TRUE, cf_fold = 5
)
summary(fit_cf)

## -----------------------------------------------------------------------------
summary(fit_lm, show_control_fit = TRUE)


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(MRTAnalysis)
current_options <- options(digits = 3) # save current options for restoring later

## -----------------------------------------------------------------------------
head(data_mimicHeartSteps, 10)

## -----------------------------------------------------------------------------
fit1 <- wcls(
    data = data_mimicHeartSteps,
    id = userid,
    outcome = logstep_30min,
    treatment = intervention,
    rand_prob = 0.6,
    moderator_formula = ~1,
    control_formula = ~logstep_pre30min,
    availability = avail
)
summary(fit1)

## -----------------------------------------------------------------------------
fit2 <- wcls(
    data = data_mimicHeartSteps,
    id = userid,
    outcome = logstep_30min,
    treatment = intervention,
    rand_prob = 0.6,
    moderator_formula = ~1,
    control_formula = ~ logstep_pre30min + logstep_30min_lag1 + is_at_home_or_work,
    availability = avail
)
summary(fit2)

## -----------------------------------------------------------------------------
summary(fit2, show_control_fit = TRUE)

## -----------------------------------------------------------------------------
fit3 <- wcls(
    data = data_mimicHeartSteps,
    id = userid,
    outcome = logstep_30min,
    treatment = intervention,
    rand_prob = 0.6,
    moderator_formula = ~is_at_home_or_work,
    control_formula = ~ logstep_pre30min + logstep_30min_lag1 + is_at_home_or_work,
    availability = avail
)
summary(fit3)

## -----------------------------------------------------------------------------
summary(fit3, lincomb = c(1, 1))

## -----------------------------------------------------------------------------
head(data_binary, 30)

## -----------------------------------------------------------------------------
fit4 <- emee(
    data = data_binary,
    id = userid,
    outcome = Y,
    treatment = A,
    rand_prob = rand_prob,
    moderator_formula = ~1,
    control_formula = ~ time_var1 + time_var2,
    availability = avail
)
summary(fit4)

## -----------------------------------------------------------------------------
summary(fit4, show_control_fit = TRUE)

## -----------------------------------------------------------------------------
fit5 <- emee(
    data = data_binary,
    id = userid,
    outcome = Y,
    treatment = A,
    rand_prob = rand_prob,
    moderator_formula = ~time_var1,
    control_formula = ~ time_var1 + time_var2,
    availability = avail
)
summary(fit5)

## -----------------------------------------------------------------------------
summary(fit5, lincomb = rbind(c(1, 0.0333), c(1, 0.5), c(1, 1)))

## -----------------------------------------------------------------------------
options(current_options) # restore old options


## code to prepare `data_mimicHeartSteps` dataset goes here

# A generative model that mimics HeartSteps

# used in the MRT paper at Psychology Methods

# read the empirical distribution of jbsteps30pre.log from HeartSteps
# jbsteps30pre.log.empirical <- readRDS("jbsteps30pre.log.empirical.RDS")

dgm_mimicHeartSteps <- function(sample_size, total_T = 210) {
    #                                  Estimate   95% LCL   95% UCL        SE Hotelling p-value
    # (Intercept)                      1.755916  1.479176  2.032656  0.135310   168.402 < 1e-04 ***
    # jbsteps30pre.log                 0.403241  0.347119  0.459363  0.027440   215.949 < 1e-04 ***
    # jbsteps30.log.lag1               0.063676  0.016708  0.110644  0.022965     7.688 0.00961 **
    # study.day.nogap                 -0.009488 -0.018368 -0.000607  0.004342     4.775 0.03711 *
    # location.homework                0.119917 -0.096898  0.336733  0.106010     1.280 0.26725
    # I(send - 0.6)                    0.442889  0.128992  0.756785  0.153477     8.327 0.00730 **
    # I(send - 0.6):study.day.nogap   -0.017908 -0.030018 -0.005798  0.005921     9.147 0.00517 **
    # I(send - 0.6):location.homework  0.101423 -0.239240  0.442086  0.166565     0.371 0.54733

    alpha_0 <- 1.755916 # intercept
    alpha_1 <- 0.403241 # jbsteps30pre.log
    alpha_2 <- 0.063676 # jbsteps30.log.lag1
    alpha_3 <- -0.009488 # study.day.nogap
    alpha_4 <- 0.119917 # location.homework
    beta_0 <- 0.442889 # I(send - 0.6)
    beta_1 <- -0.017908 # I(send - 0.6):study.day.nogap
    beta_2 <- 0.101423 # I(send - 0.6):location.homework

    df_names <- c(
        "userid", "decision.index.nogap", "study.day.nogap",
        "jbsteps30.log", "jbsteps30.log.lag1", "jbsteps30pre.log",
        "location.homework", "send", "rand_prob", "avail"
    )

    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names

    dta$userid <- rep(1:sample_size, each = total_T)
    dta$decision.index.nogap <- rep(1:total_T, times = sample_size)
    dta$study.day.nogap <- floor((dta$decision.index.nogap - 1) / 5)

    # generate jbsteps30pre.log
    tmp <- rnorm(n = nrow(dta), mean = 20, sd = 50)
    tmp <- pmax(tmp, 0)
    dta$jbsteps30pre.log <- log(tmp + 0.5)
    dta$jbsteps30pre.log <- pmin(dta$jbsteps30pre.log, 8.5)
    # exp(8.5) = 4914.769; max jbstep30 = 5527 in original HeartSteps

    # generate location.homework (using marginal mean of HeartSteps data)
    dta$location.homework <- rbinom(nrow(dta), 1, 0.3624668)

    # generate avail (using marginal mean of HeartSteps data)
    dta$avail <- rbinom(nrow(dta), 1, 0.8038462)

    # generate send
    dta$rand_prob <- 0.6
    dta$send <- rbinom(nrow(dta), 1, 0.6)
    dta$send <- ifelse(dta$avail, dta$send, 0)

    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        if (t == 1) {
            dta$jbsteps30.log.lag1[row_index] <- 0
        } else {
            row_index_pre <- seq(from = t - 1, by = total_T, length = sample_size)
            dta$jbsteps30.log.lag1[row_index] <- dta$jbsteps30.log[row_index_pre]
        }
        this_meanY <- alpha_0 + alpha_1 * dta$jbsteps30pre.log[row_index] + alpha_2 * dta$jbsteps30.log.lag1[row_index] +
            alpha_3 * dta$study.day.nogap[row_index] + alpha_4 * dta$location.homework[row_index] +
            (dta$send[row_index] - 0.6) * (beta_0 + beta_1 * dta$study.day.nogap[row_index] +
                beta_2 * dta$location.homework[row_index])

        # Make it so that the generated logged step count cannot be smaller than log(0.5)
        # Of course this would change the true causal parameter values
        dta$jbsteps30.log[row_index] <- pmax(
            rnorm(n = sample_size, mean = this_meanY, 2.716),
            log(0.5)
        )
        dta$jbsteps30.log[row_index] <- pmin(dta$jbsteps30.log[row_index], 8.5)
        # exp(8.5) = 4914.769; max jbstep30 = 5527 in original HeartSteps
    }

    names(dta)[names(dta) == "decision.index.nogap"] <- "decision_point"
    names(dta)[names(dta) == "study.day.nogap"] <- "day_in_study"
    names(dta)[names(dta) == "jbsteps30.log"] <- "logstep_30min"
    names(dta)[names(dta) == "jbsteps30.log.lag1"] <- "logstep_30min_lag1"
    names(dta)[names(dta) == "jbsteps30pre.log"] <- "logstep_pre30min"
    names(dta)[names(dta) == "location.homework"] <- "is_at_home_or_work"
    names(dta)[names(dta) == "send"] <- "intervention"
    names(dta)[names(dta) == "avail"] <- "avail"
    return(dta)
}

set.seed(123)
data_mimicHeartSteps <- dgm_mimicHeartSteps(sample_size = 37)
usethis::use_data(data_mimicHeartSteps, overwrite = TRUE)

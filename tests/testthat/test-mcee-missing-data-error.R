test_that("mcee errors cleanly on missing/NaN/Inf in core, moderator, and control vars", {
    set.seed(42)

    # ---- base DGM (no missing) ----
    n <- 6
    Ti <- c(5, 7, 6, 6, 8, 5)
    id <- rep(seq_len(n), Ti)
    dp <- unlist(lapply(Ti, seq_len))
    I <- rbinom(length(dp), 1, 0.9)
    A <- rbinom(length(dp), 1, 0.6)
    M <- rbinom(length(dp), 1, plogis(-0.2 + 0.3 * A + 0.1 * scale(dp)))
    Ytmp <- 0.5 * A + 0.6 * M + 0.08 * scale(dp) + rnorm(length(dp), 0, 0.2)
    Y <- ave(Ytmp, id, FUN = function(v) rep(mean(v), length(v)))

    d0 <- data.frame(id, dp, I, A, M, Y, check.names = FALSE)

    # nontrivial weights
    w <- ave(0.3 + 0.7 * (d0$dp / ave(d0$dp, d0$id, FUN = max)),
        d0$id,
        FUN = function(v) v / sum(v)
    )

    # 1) Missing in a core column (dp)
    d1 <- d0
    d1$dp[c(3, 10)] <- NA_integer_
    expect_error(
        mcee(
            data = d1, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf.*dp.*rows\\s+3,\\s+10.*does not support handling missing data",
        ignore.case = TRUE
    )

    # 2) NaN/Inf in moderator vars (here dp is used in moderator)
    d2 <- d0
    d2$dp[c(2, 9)] <- Inf
    expect_error(
        mcee(
            data = d2, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~dp,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf.*dp.*rows\\s+2,\\s+9",
        ignore.case = TRUE
    )

    # 3) Missing in a control-only variable (M here)
    d3 <- d0
    d3$M[c(4, 12, 25)] <- NA_real_
    expect_error(
        mcee(
            data = d3, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf.*M.*rows\\s+4,\\s+12,\\s+25",
        ignore.case = TRUE
    )
})

test_that("mcee_general errors cleanly on missing in config formula variables", {
    set.seed(99)

    # base DGM again
    n <- 6
    Ti <- c(5, 7, 6, 6, 8, 5)
    id <- rep(seq_len(n), Ti)
    dp <- unlist(lapply(Ti, seq_len))
    I <- rbinom(length(dp), 1, 0.9)
    A <- rbinom(length(dp), 1, 0.6)
    M <- rbinom(length(dp), 1, plogis(-0.2 + 0.3 * A + 0.1 * scale(dp)))
    Ytmp <- 0.5 * A + 0.6 * M + 0.08 * scale(dp) + rnorm(length(dp), 0, 0.2)
    Y <- ave(Ytmp, id, FUN = function(v) rep(mean(v), length(v)))
    Z <- rnorm(length(dp)) # extra covariate for configs

    d0 <- data.frame(id, dp, I, A, M, Y, Z, check.names = FALSE)
    w <- ave(0.3 + 0.7 * (d0$dp / ave(d0$dp, d0$id, FUN = max)),
        d0$id,
        FUN = function(v) v / sum(v)
    )

    # Make NA only in a variable used by a config (Z used in q)
    d1 <- d0
    d1$Z[c(6, 13)] <- NA_real_

    cfg_p <- list(method = "glm", formula = ~dp) # binomial auto
    cfg_q <- list(method = "glm", formula = ~ dp + Z) # binomial auto; Z has NA
    cfg_eta <- list(method = "glm", formula = ~dp) # gaussian auto
    cfg_mu <- list(method = "glm", formula = ~ dp + M) # gaussian auto
    cfg_nu <- list(method = "glm", formula = ~dp) # gaussian auto

    expect_error(
        mcee_general(
            data = d1,
            id = "id", dp = "dp", outcome = "Y",
            treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            config_p = cfg_p, config_q = cfg_q,
            config_eta = cfg_eta, config_mu = cfg_mu, config_nu = cfg_nu,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf.*Z.*rows\\s+6,\\s+13",
        ignore.case = TRUE
    )
})

test_that("mcee_userfit_nuisance errors cleanly on missing in supplied nuisance vectors", {
    set.seed(123)

    # base DGM
    n <- 6
    Ti <- c(5, 7, 6, 6, 8, 5)
    id <- rep(seq_len(n), Ti)
    dp <- unlist(lapply(Ti, seq_len))
    I <- rbinom(length(dp), 1, 0.9)
    A <- rbinom(length(dp), 1, 0.6)
    M <- rbinom(length(dp), 1, plogis(-0.2 + 0.3 * A + 0.1 * scale(dp)))
    Ytmp <- 0.5 * A + 0.6 * M + 0.08 * scale(dp) + rnorm(length(dp), 0, 0.2)
    Y <- ave(Ytmp, id, FUN = function(v) rep(mean(v), length(v)))

    d0 <- data.frame(id, dp, I, A, M, Y, check.names = FALSE)
    w <- ave(0.3 + 0.7 * (d0$dp / ave(d0$dp, d0$id, FUN = max)),
        d0$id,
        FUN = function(v) v / sum(v)
    )

    # Create sane nuisance predictions first
    p1 <- plogis(-0.1 + 0.02 * dp) # in (0,1)
    q1 <- plogis(-0.2 + 0.02 * dp + 0.3 * M)
    eta1 <- 0.4 + 0.1 * dp
    eta0 <- 0.3 + 0.05 * dp
    mu1 <- 0.4 + 0.15 * dp + 0.2 * M
    mu0 <- 0.3 + 0.10 * dp + 0.1 * M
    nu1 <- 0.35 + 0.12 * dp
    nu0 <- 0.25 + 0.08 * dp

    # Inject NA into one nuisance vector
    p1_bad <- p1
    p1_bad[c(2, 11, 17)] <- NA_real_

    expect_error(
        mcee_userfit_nuisance(
            data = d0,
            id = "id", dp = "dp", outcome = "Y",
            treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~dp,
            p1 = p1_bad, q1 = q1,
            eta1 = eta1, eta0 = eta0,
            mu1 = mu1, mu0 = mu0,
            nu1 = nu1, nu0 = nu0,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf detected in 'p1'.*rows\\s+2,\\s+11,\\s+17.*does not support handling missing data",
        ignore.case = TRUE
    )

    # Inject Inf into another vector (nu0)
    nu0_bad <- nu0
    nu0_bad[c(3, 9)] <- Inf
    expect_error(
        mcee_userfit_nuisance(
            data = d0,
            id = "id", dp = "dp", outcome = "Y",
            treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            p1 = p1, q1 = q1,
            eta1 = eta1, eta0 = eta0,
            mu1 = mu1, mu0 = mu0,
            nu1 = nu1, nu0 = nu0_bad,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "Missing/NaN/Inf detected in 'nu0'.*rows\\s+3,\\s+9",
        ignore.case = TRUE
    )
})

test_that("mcee_* missing-data message aggregates multiple offenders", {
    set.seed(7)

    n <- 5
    Ti <- rep(5, n)
    id <- rep(seq_len(n), Ti)
    dp <- unlist(lapply(Ti, seq_len))
    I <- rep(1, length(dp))
    A <- rbinom(length(dp), 1, 0.5)
    M <- rbinom(length(dp), 1, 0.5)
    Y <- ave(0.2 * A + 0.3 * M + rnorm(length(dp), 0, .1), id, FUN = function(v) rep(mean(v), length(v)))

    d0 <- data.frame(id, dp, I, A, M, Y)

    # Make both Y and M missing in a few rows (should list both variables)
    d1 <- d0
    d1$Y[c(4, 7)] <- NA_real_
    d1$M[c(3, 8)] <- NA_real_

    w <- rep(1, nrow(d1))

    expect_error(
        mcee(
            data = d1, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        regexp = "(Y.*rows\\s+4,\\s+7|M.*rows\\s+3,\\s+8).*(Y.*rows\\s+4,\\s+7|M.*rows\\s+3,\\s+8)",
        ignore.case = TRUE
    )
})

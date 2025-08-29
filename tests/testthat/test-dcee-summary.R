test_that("summary.dcee: single vector lincomb matches manual L %*% beta", {
    skip_on_cran()

    dat <- data_distal_continuous

    # Need at least 2 betas -> include one moderator (Z)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = FALSE,
        verbose = FALSE
    )

    s <- summary(fit, lincomb = c(0, 1), conf_level = 0.95)

    beta <- fit$fit$beta_hat
    V <- fit$fit$beta_varcov
    L <- matrix(c(0, 1), nrow = 1)
    est <- as.numeric(L %*% beta)
    se <- sqrt(as.numeric(L %*% V %*% t(L)))
    df <- fit$df
    tcrit <- stats::qt(0.975, df = df)

    # pull from summary table
    got <- s$lincomb[1, , drop = FALSE]

    expect_equal(got$Estimate, est, tolerance = 1e-8)
    expect_equal(got$`Std. Error`, se, tolerance = 1e-8)
    expect_equal(got$df, df)
    expect_equal(got$`t value`, est / se, tolerance = 1e-8)
    expect_equal(got[[grep("LCL$", names(got), value = TRUE)]], est - tcrit * se, tolerance = 1e-8)
    expect_equal(got[[grep("UCL$", names(got), value = TRUE)]], est + tcrit * se, tolerance = 1e-8)
})



test_that("summary.dcee: matrix lincomb (multiple rows) and rownames preserved", {
    skip_on_cran()

    dat <- data_distal_continuous

    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        verbose = FALSE
    )

    # L1 = Intercept only; L2 = Moderator only
    L <- rbind(
        Intercept = c(1, 0),
        Z_only    = c(0, 1)
    )

    s <- summary(fit, lincomb = L, conf_level = 0.90) # non-default CL to exercise tcrit change
    expect_true(is.data.frame(s$lincomb))
    expect_equal(rownames(s$lincomb), rownames(L))

    beta <- fit$fit$beta_hat
    V <- fit$fit$beta_varcov
    est <- as.numeric(L %*% beta)
    se <- sqrt(diag(L %*% V %*% t(L)))
    df <- fit$df
    tcrit <- stats::qt(0.95, df = df) # 90% CI

    expect_equal(as.numeric(s$lincomb$Estimate), as.numeric(est), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb$`Std. Error`), as.numeric(se), tolerance = 1e-8)
    expect_equal(s$lincomb$df, rep(df, length(est)))
    expect_equal(as.numeric(s$lincomb$`t value`), as.numeric(est / se), tolerance = 1e-8)

    LCL_col <- grep("LCL$", names(s$lincomb), value = TRUE)
    UCL_col <- grep("UCL$", names(s$lincomb), value = TRUE)
    expect_equal(as.numeric(s$lincomb[[LCL_col]]), as.numeric(est - tcrit * se), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb[[UCL_col]]), as.numeric(est + tcrit * se), tolerance = 1e-8)
})




test_that("summary.dcee: transposes lincomb when user supplies p x k", {
    skip_on_cran()

    dat <- data_distal_continuous

    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        verbose = FALSE
    )

    p <- length(fit$fit$beta_hat)
    expect_equal(p, 2L) # sanity

    # Provide p x k instead of k x p. Here, p x 2 → will be transposed internally.
    L_wrong_shape <- cbind(c(1, 0), c(0, 1)) # 2x2
    colnames(L_wrong_shape) <- NULL
    rownames(L_wrong_shape) <- NULL

    s <- summary(fit, lincomb = L_wrong_shape, conf_level = 0.95)
    expect_true(is.data.frame(s$lincomb))
    expect_equal(nrow(s$lincomb), 2L)
})


test_that("summary.dcee: invalid lincomb shapes error clearly", {
    skip_on_cran()

    dat <- data_distal_continuous

    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        verbose = FALSE
    )

    # vector wrong length
    expect_error(summary(fit, lincomb = c(1, 1, 1)), "`lincomb` vector must have length 2.", ignore.case = TRUE)

    # matrix with wrong cols/rows
    expect_error(summary(fit, lincomb = matrix(1, nrow = 3, ncol = 3)),
        "must have 2 columns",
        fixed = FALSE
    )
})



test_that("summary.dcee: lincomb works with GAM nuisance as well", {
    skip_on_cran()
    skip_if_not_installed("mgcv")

    dat <- data_distal_continuous

    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ s(X) + Z,
        availability = "avail",
        control_reg_method = "gam",
        verbose = FALSE
    )

    s <- summary(fit, lincomb = c(0, 1))
    expect_true(is.data.frame(s$lincomb))
    expect_equal(nrow(s$lincomb), 1L)
    expect_true(all(is.finite(unlist(s$lincomb))))
})

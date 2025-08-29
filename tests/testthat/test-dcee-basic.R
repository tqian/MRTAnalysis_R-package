test_that("dcee: basic lm fit runs and returns expected structure", {
    skip_on_cran()

    data <- data_distal_continuous
    expect_true(is.data.frame(data))

    fit <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = FALSE,
        verbose = FALSE
    )

    expect_s3_class(fit, "dcee_fit")
    expect_type(fit$fit, "list")
    expect_named(
        fit$fit,
        c("beta_hat", "beta_se", "beta_varcov", "conf_int", "conf_int_tquantile", "regfit_a0", "regfit_a1"),
        ignore.order = TRUE
    )
    # estimates are finite
    expect_true(all(is.finite(fit$fit$beta_hat)))
    expect_true(all(is.finite(fit$fit$beta_se)))
    expect_true(all(is.finite(diag(fit$fit$beta_varcov))))

    # df looks sensible: n_subjects - p
    n_id <- length(unique(data$userid))
    p <- length(fit$fit$beta_hat)
    expect_equal(fit$df, n_id - p)

    # summary runs
    s <- summary(fit)
    expect_s3_class(s, "summary.dcee_fit")
    expect_true(is.data.frame(s$excursion_effect))
    expect_equal(nrow(s$excursion_effect), p)
})

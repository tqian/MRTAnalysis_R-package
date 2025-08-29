test_that("dcee: GAM nuisance works (mgcv)", {
    skip_on_cran()
    skip_if_not_installed("mgcv")

    data <- data_distal_continuous

    fit_gam <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~Z, # add a moderator
        control_formula = ~ s(X) + Z, # smooth term allowed
        availability = "avail",
        control_reg_method = "gam",
        cross_fit = FALSE,
        verbose = FALSE
    )

    expect_s3_class(fit_gam, "dcee_fit")
    expect_true(all(is.finite(fit_gam$fit$beta_hat)))
    # summary with lincomb (pull the moderator coefficient)
    s <- summary(fit_gam, lincomb = c(0, 1))
    expect_s3_class(s, "summary.dcee_fit")
    expect_true(is.data.frame(s$lincomb))
    expect_equal(nrow(s$lincomb), 1)
})

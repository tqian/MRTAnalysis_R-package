test_that("dcee: SuperLearner smoke test if installed", {
    skip_on_cran()
    skip_if_not_installed("SuperLearner")

    suppressMessages(library(SuperLearner))

    data <- data_distal_continuous

    fit_sl <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "sl",
        cross_fit = FALSE,
        verbose = FALSE
    )
    expect_s3_class(fit_sl, "dcee_fit")
    expect_true(all(is.finite(fit_sl$fit$beta_hat)))

    # summary should print a short advisory rather than dumping the whole object
    expect_no_error(summary(fit_sl, show_control_fit = TRUE))
})

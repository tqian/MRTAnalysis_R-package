test_that("dcee: default availability (NULL) is treated as all-ones", {
    skip_on_cran()

    data <- data_distal_continuous
    data$avail_backup <- data$avail
    data$avail <- NULL

    fit <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~X,
        availability = NULL, # default path
        control_reg_method = "lm",
        cross_fit = FALSE,
        verbose = FALSE
    )
    expect_s3_class(fit, "dcee_fit")
})

test_that("dcee: trivial numeric weighting function accepted", {
    skip_on_cran()

    data <- data_distal_continuous
    Ttot <- max(data$dp) # or time variable if different
    omega_const <- 1 / Ttot

    fit <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~X,
        availability = "avail",
        control_reg_method = "lm",
        weighting_function = omega_const,
        cross_fit = FALSE,
        verbose = FALSE
    )
    expect_s3_class(fit, "dcee_fit")
})

test_that("dcee: nontrivial column weighting works", {
    skip_on_cran()

    data <- data_distal_continuous
    # make an indicator: weight only the first decision point for each subject
    stopifnot(all(c("userid", "dp") %in% names(data)))
    data$omega <- as.numeric(data$dp == 1)

    fit <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~X,
        availability = "avail",
        control_reg_method = "lm",
        weighting_function = "omega",
        cross_fit = FALSE,
        verbose = FALSE
    )
    expect_s3_class(fit, "dcee_fit")
})

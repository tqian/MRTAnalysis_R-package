test_that("dcee: errors for bad treatment coding", {
    skip_on_cran()

    data <- data_distal_continuous
    data$A_bad <- ifelse(data$A == 1, 2, 0)

    expect_error(
        dcee(
            data = data,
            id = "userid",
            outcome = "Y",
            treatment = "A_bad",
            rand_prob = "prob_A",
            moderator_formula = ~1,
            control_formula = ~X,
            availability = "avail",
            control_reg_method = "lm",
            verbose = FALSE
        ),
        regexp = "must be coded 0/1",
        fixed = FALSE
    )
})

test_that("dcee: s() terms rejected unless method = 'gam'", {
    skip_on_cran()

    data <- data_distal_continuous

    expect_error(
        dcee(
            data = data,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator_formula = ~1,
            control_formula = ~ s(X), # not allowed with lm
            availability = "avail",
            control_reg_method = "lm",
            verbose = FALSE
        ),
        "`s( )` terms are only supported for control_reg_method = 'gam'.",
        fixed = TRUE
    )
})

test_that("dcee: control_formula ~ 1 only allowed with set_to_zero", {
    skip_on_cran()

    data <- data_distal_continuous

    expect_error(
        dcee(
            data = data,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator_formula = ~1,
            control_formula = ~1, # only allowed with set_to_zero
            availability = "avail",
            control_reg_method = "rf",
            verbose = FALSE
        ),
        "`control_formula = ~ 1` is only allowed with control_reg_method in {'set_to_zero','lm','gam'}.",
        fixed = TRUE
    )
})

test_that("rand_prob constraint only enforced when availability == 1", {
    data <- data_distal_continuous

    # Allow invalid rand_prob at avail==0
    data$prob_A_bad <- data$prob_A
    data$prob_A_bad[data$avail == 0] <- 0 # invalid but ignored because avail==0
    expect_no_error(
        dcee(
            data,
            id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A_bad",
            moderator_formula = ~1, control_formula = ~X, availability = "avail",
            control_reg_method = "lm",
            verbose = FALSE
        )
    )

    # Still error when invalid and avail==1
    data$prob_A_bad2 <- data$prob_A
    data$prob_A_bad2[data$avail == 1][1] <- 0
    expect_error(
        dcee(
            data,
            id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A_bad2",
            moderator_formula = ~1, control_formula = ~X, availability = "avail",
            control_reg_method = "lm",
            verbose = FALSE
        ),
        "must lie strictly in \\(0,1\\).*availability = 1",
        fixed = FALSE
    )
})

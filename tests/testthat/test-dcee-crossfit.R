test_that("dcee: cross-fitting executes and structure is intact", {
    skip_on_cran()

    data <- data_distal_continuous

    fit_cf <- dcee(
        data = data,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = TRUE,
        cf_fold = 5,
        verbose = FALSE
    )

    expect_s3_class(fit_cf, "dcee_fit")
    expect_true(all(is.finite(fit_cf$fit$beta_hat)))
    # Stage-1 objects exist (last-fold models)
    expect_true(is.list(fit_cf$fit$regfit_a0) || inherits(fit_cf$fit$regfit_a0, c("lm", "gam")) || is.null(fit_cf$fit$regfit_a0))
    expect_true(is.list(fit_cf$fit$regfit_a1) || inherits(fit_cf$fit$regfit_a1, c("lm", "gam")) || is.null(fit_cf$fit$regfit_a1))

    # summary prints with control fits shown (don’t assert content, just no error)
    expect_no_error(summary(fit_cf, show_control_fit = TRUE))
})


test_that("Stage-1 cross-fitting handles uneven folds and predicts all rows (lm)", {
    skip_on_cran()
    set.seed(123)

    # small synthetic distal-outcome data where n_id is NOT divisible by cf_fold
    make_toy_distal <- function(n_id = 23, T = 4) {
        ids <- rep(seq_len(n_id), each = T)
        dp <- rep(seq_len(T), times = n_id)
        X <- rnorm(n_id * T)
        Z <- rbinom(n_id * T, 1, 0.4)
        avail <- rbinom(n_id * T, 1, 0.9)
        pA <- plogis(0.2 * X - 0.5 * Z) # varies in (0,1)
        A <- rbinom(n_id * T, 1, ifelse(avail == 1, pA, 0)) * avail

        # distal Y: one value per id (repeat across rows)
        # depends on summaries of X, Z, and exposure to A
        Y_id <- tapply(X, ids, mean) * 1.0 + tapply(Z, ids, mean) * 0.5 +
            tapply(A, ids, mean) * 1.2 + rnorm(n_id, sd = 0.7)
        Y <- Y_id[as.character(ids)]

        data.frame(
            userid = ids, dp = dp,
            X = X, Z = Z,
            avail = avail, prob_A = pA,
            A = A, Y = as.numeric(Y)
        )
    }

    dta <- make_toy_distal(n_id = 23, T = 5)
    cf_fold <- 5
    expect_false(nrow(unique(data.frame(id = dta$userid))) %% cf_fold == 0)

    # Call the internal Stage-1 helper via ::: to check out-of-fold predictions
    out1 <- MRTAnalysis:::dcee_helper_stage1_fit_nuisance(
        dta = dta,
        id_var = "userid",
        trt_var = "A",
        outcome_var = "Y",
        control_reg_method = "lm",
        control_formula = ~ X + Z, # RHS-only
        cross_fit = TRUE,
        cf_fold = cf_fold
    )

    # Every row should have non-NA mu_hat_0 and mu_hat_1
    expect_true(all(!is.na(out1$dta$mu_hat_0)))
    expect_true(all(!is.na(out1$dta$mu_hat_1)))

    # Last-fold learner objects should exist
    expect_s3_class(out1$regfit_a0, "lm")
    expect_s3_class(out1$regfit_a1, "lm")
})

test_that("dcee() completes with cross-fitting and uneven folds (lm)", {
    skip_on_cran()
    set.seed(456)

    # reuse helper to construct an uneven-fold dataset
    make_toy_distal <- function(n_id = 23, T = 4) {
        ids <- rep(seq_len(n_id), each = T)
        dp <- rep(seq_len(T), times = n_id)
        X <- rnorm(n_id * T)
        Z <- rbinom(n_id * T, 1, 0.4)
        avail <- rbinom(n_id * T, 1, 0.9)
        pA <- plogis(0.2 * X - 0.5 * Z)
        A <- rbinom(n_id * T, 1, ifelse(avail == 1, pA, 0)) * avail
        Y_id <- tapply(X, ids, mean) * 0.8 + tapply(Z, ids, mean) * 0.6 +
            tapply(A, ids, mean) * 1.0 + rnorm(n_id, sd = 0.6)
        Y <- Y_id[as.character(ids)]
        data.frame(
            userid = ids, dp = dp,
            X = X, Z = Z,
            avail = avail, prob_A = pA,
            A = A, Y = as.numeric(Y)
        )
    }

    dta <- make_toy_distal(n_id = 29, T = 4) # 29 not divisible by 6
    cf_fold <- 6
    expect_false(length(unique(dta$userid)) %% cf_fold == 0)

    fit <- dcee(
        data = dta,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = TRUE, cf_fold = cf_fold,
        verbose = FALSE
    )

    # Structure checks
    expect_s3_class(fit, "dcee_fit")
    expect_type(fit$fit$beta_hat, "double")
    expect_true(all(is.finite(fit$fit$beta_hat)))
    expect_true(is.matrix(fit$fit$beta_varcov))
    expect_true(all(is.finite(diag(fit$fit$beta_varcov))))

    # Last-fold models should be present for inspection
    expect_s3_class(fit$fit$regfit_a0, "lm")
    expect_s3_class(fit$fit$regfit_a1, "lm")

    # summary should run and produce a table
    s <- summary(fit)
    expect_s3_class(s, "summary.dcee_fit")
    expect_true(nrow(s$excursion_effect) >= 1)
})

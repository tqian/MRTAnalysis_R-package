test_that("dcee: numeric outputs are stable on a tiny deterministic DGM (lm)", {
    skip_on_cran()

    set.seed(42)

    # ----- Tiny generative model (fast + simple) --------------------------------
    # Structure matches distal format: multiple rows per id, Y is distal (constant within id)
    n <- 60 # subjects
    Tt <- 4 # decision points per subject

    df <- data.frame(
        userid = rep(seq_len(n), each = Tt),
        dp     = rep(seq_len(Tt), times = n)
    )

    # Availability and randomization probability
    df$avail <- 1L
    df$prob_A <- 0.5

    # Time-varying covariates and treatment
    df$X <- rnorm(n * Tt, 0, 1)
    df$Z <- rbinom(n * Tt, 1, 0.4)
    df$A <- rbinom(n * Tt, 1, df$prob_A) * df$avail

    # Distal outcome Y (constant per subject).
    # Simple, linear signal tied to X and A*X so the moderator example is meaningful.
    # Y_i = c0 + c1*mean(X_i.) + c2*mean(A_i.) + c3*mean(A_i.*X_i.) + eps
    c0 <- 0.3
    c1 <- 0.8
    c2 <- 0.5
    c3 <- 1.2
    dY <- aggregate(cbind(X, A) ~ userid, data = df, FUN = mean)
    dAX <- aggregate(I(A * X) ~ userid, data = df, FUN = mean)
    dY <- merge(dY, dAX, by = "userid")
    names(dY) <- c("userid", "mX", "mA", "mAX")
    dY$eps <- rnorm(n, 0, 0.15)
    dY$Y <- with(dY, c0 + c1 * mX + c2 * mA + c3 * mAX + eps)

    df <- merge(df, dY[c("userid", "Y")], by = "userid", all.x = TRUE)

    # ----- Fit 1: Marginal effect (no moderators) -------------------------------
    fit_mar <- dcee(
        data = df,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~1,
        control_formula = ~X, # nuisance model uses X (linear)
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = FALSE,
        verbose = FALSE
    )

    # Current "golden" values: run once, print, and paste below.
    # cat("beta_hat (marginal) =", dput(unclass(fit_mar$fit$beta_hat)))
    # cat("beta_se  (marginal) =", dput(unclass(fit_mar$fit$beta_se)))
    # cat("beta_V   (marginal) =", dput(unclass(fit_mar$fit$beta_varcov)))

    # ---- TODO: paste the numbers you get on first run here ---------------------
    expected_beta_mar <- structure(
        c(
            0.0580563073126 # Intercept
        ),
        .Names = "Intercept"
    )
    expected_se_mar <- structure(
        c(0.0657496189214),
        .Names = "Intercept"
    )
    expected_V_mar <- matrix(
        c(0.00432301238831),
        nrow = 1, dimnames = list("Intercept", "Intercept")
    )
    # ---------------------------------------------------------------------------

    # Check numeric equality with tight tolerances (adjust if needed)
    expect_equal(unclass(fit_mar$fit$beta_hat), expected_beta_mar, tolerance = 1e-8)
    expect_equal(unclass(fit_mar$fit$beta_se), expected_se_mar, tolerance = 1e-8)
    expect_equal(fit_mar$fit$beta_varcov, expected_V_mar, tolerance = 1e-8)

    # Also check CI agrees with beta ± t_crit * se using df from object
    df_mar <- fit_mar$df
    tcrit <- stats::qt(0.975, df = df_mar)
    ci_low <- expected_beta_mar - tcrit * expected_se_mar
    ci_high <- expected_beta_mar + tcrit * expected_se_mar
    got_ci <- fit_mar$fit$conf_int_tquantile
    expect_equal(as.numeric(got_ci[, 1]), as.numeric(ci_low), tolerance = 1e-8)
    expect_equal(as.numeric(got_ci[, 2]), as.numeric(ci_high), tolerance = 1e-8)

    # ----- Fit 2: Moderated effect (moderator = X) ------------------------------
    fit_mod <- dcee(
        data = df,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~X, # estimate intercept + effect moderated by X
        control_formula = ~X, # nuisance model uses X (linear)
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = FALSE,
        verbose = FALSE
    )

    # cat("beta_hat (moderated) =", dput(unclass(fit_mod$fit$beta_hat)))
    # cat("beta_se  (moderated) =", dput(unclass(fit_mod$fit$beta_se)))
    # cat("beta_V   (moderated) =", dput(unclass(fit_mod$fit$beta_varcov)))

    # ---- TODO: paste the numbers you get on first run here ---------------------
    expected_beta_mod <- structure(
        c(
            0.0632080918257, # Intercept
            0.1334112843198
        ), # X
        .Names = c("Intercept", "X")
    )
    expected_se_mod <- structure(
        c(
            0.0646006154111,
            0.0758280636588
        ),
        .Names = c("Intercept", "X")
    )
    expected_V_mod <- matrix(
        c(
            0.004173239511498, 0.000404412691711,
            0.000404412691711, 0.005749895238241
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("Intercept", "X"), c("Intercept", "X"))
    )
    # ---------------------------------------------------------------------------

    expect_equal(unclass(fit_mod$fit$beta_hat), expected_beta_mod, tolerance = 1e-8)
    expect_equal(unclass(fit_mod$fit$beta_se), expected_se_mod, tolerance = 1e-8)
    expect_equal(fit_mod$fit$beta_varcov, expected_V_mod, tolerance = 1e-8)

    # summary() lincomb: test L = [0,1] to extract the X effect
    s <- summary(fit_mod, lincomb = c(0, 1), conf_level = 0.95)
    estL <- as.numeric(expected_beta_mod["X"])
    seL <- as.numeric(expected_se_mod["X"])
    dfL <- fit_mod$df
    tcritL <- stats::qt(0.975, df = dfL)
    expect_equal(as.numeric(s$lincomb$Estimate), estL, tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb$`Std. Error`), seL, tolerance = 1e-8)
    expect_equal(s$lincomb$df, dfL)
    # CIs
    LCL_col <- grep("LCL$", names(s$lincomb), value = TRUE)
    UCL_col <- grep("UCL$", names(s$lincomb), value = TRUE)
    expect_equal(as.numeric(s$lincomb[[LCL_col]]), estL - tcritL * seL, tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb[[UCL_col]]), estL + tcritL * seL, tolerance = 1e-8)
})

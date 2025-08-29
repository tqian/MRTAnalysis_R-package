# DCEE goldens with cross-fitting

.make_tiny_dgm <- function(n = 40, Tt = 4, seed = 123) {
    set.seed(seed)
    df <- data.frame(
        userid = rep(seq_len(n), each = Tt),
        dp     = rep(seq_len(Tt), times = n)
    )
    df$avail <- 1L
    df$prob_A <- 0.5
    df$X <- rnorm(n * Tt, 0, 1)
    df$Z <- rbinom(n * Tt, 1, 0.4)
    df$A <- rbinom(n * Tt, 1, df$prob_A) * df$avail

    # Distal Y
    c0 <- 0.2
    c1 <- 0.7
    c2 <- 0.4
    c3 <- 0.9
    dY <- aggregate(cbind(X, A) ~ userid, data = df, FUN = mean)
    dAX <- aggregate(I(A * X) ~ userid, data = df, FUN = mean)
    dY <- merge(dY, dAX, by = "userid")
    names(dY) <- c("userid", "mX", "mA", "mAX")
    set.seed(seed + 1)
    dY$eps <- rnorm(n, 0, 0.1)
    dY$Y <- with(dY, c0 + c1 * mX + c2 * mA + c3 * mAX + eps)
    merge(df, dY[c("userid", "Y")], by = "userid", all.x = TRUE)
}

test_that("lm with cross-fitting: goldens", {
    skip_on_cran()
    skip_if_not_installed("stats")

    dat <- .make_tiny_dgm()

    set.seed(60601)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "lm",
        cross_fit = TRUE,
        cf_fold = 5,
        verbose = FALSE
    )

    # First run: uncomment to capture goldens
    # cat("lm-xfit beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("lm-xfit beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("lm-xfit beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    expected_beta <- structure(c(0.158722864905109, -0.12200774399407), .Names = c("Intercept", "Z"))
    expected_se <- structure(c(0.0929529031676633, 0.120638079419361), .Names = c("Intercept", "Z"))
    expected_V <- matrix(c(0.00864024220729699, -0.00626751154458443, -0.00626751154458443, 0.014553546205992), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    )

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    s <- summary(fit, lincomb = c(0, 1))
    expect_equal(nrow(s$lincomb), 1L)
})

test_that("gam with cross-fitting: goldens", {
    skip_on_cran()
    skip_if_not_installed("mgcv")

    dat <- .make_tiny_dgm()

    set.seed(70701)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ s(X) + Z,
        availability = "avail",
        control_reg_method = "gam",
        cross_fit = TRUE,
        cf_fold = 5,
        verbose = FALSE
    )

    # First run: uncomment to capture goldens
    # cat("gam-xfit beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("gam-xfit beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("gam-xfit beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    expected_beta <- structure(c(0.149775116348147, -0.0810214070666733), .Names = c("Intercept", "Z"))
    expected_se <- structure(c(0.0951175792493546, 0.123685706719804), .Names = c("Intercept", "Z"))
    expected_V <- matrix(c(0.00904735388225726, -0.00641923562224774, -0.00641923562224774, 0.0152981540467774), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    )

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    s <- summary(fit, lincomb = c(0, 1))
    expect_equal(nrow(s$lincomb), 1L)
})

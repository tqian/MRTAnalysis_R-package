# Shared tiny DGM (fast + deterministic)
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

    # Distal Y (constant per id): linear signal in mean X, mean A, and mean A*X
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

# "DCEE goldens for different learners"

test_that("GAM (mgcv): goldens for beta / se / V and lincomb", {
    skip_on_cran()
    skip_if_not_installed("mgcv")

    dat <- .make_tiny_dgm()

    # Fixed seed for any internal RNG (usually not needed for GAM, but cheap to set)
    set.seed(20201)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z, # => beta: Intercept, Z
        control_formula = ~ s(X) + Z, # spline in nuisance
        availability = "avail",
        control_reg_method = "gam",
        cross_fit = FALSE,
        verbose = FALSE
    )

    # --- First time: uncomment to capture goldens ---
    # cat("GAM beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("GAM beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("GAM beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    # >>>> PASTE YOUR REAL NUMBERS BELOW (placeholders shown) <<<<
    expected_beta <- structure(c(0.151147323591678, -0.102361877063079), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_se <- structure(c(0.0887917218842927, 0.112450883649567), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_V <- matrix(c(0.00788396987517758, -0.00545529153606339, -0.00545529153606339, 0.0126452012335686), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    ) # PLACEHOLDER

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    # lincomb: extract Z effect via [0, 1]
    s <- summary(fit, lincomb = c(0, 1), conf_level = 0.95)
    df <- fit$df
    tcrit <- stats::qt(0.975, df)
    est <- expected_beta["Z"]
    se <- expected_se["Z"]
    LCL <- est - tcrit * se
    UCL <- est + tcrit * se
    LCL_col <- grep("LCL$", names(s$lincomb), value = TRUE)
    UCL_col <- grep("UCL$", names(s$lincomb), value = TRUE)
    expect_equal(as.numeric(s$lincomb$Estimate), as.numeric(est), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb$`Std. Error`), as.numeric(se), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb[[LCL_col]]), as.numeric(LCL), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb[[UCL_col]]), as.numeric(UCL), tolerance = 1e-8)
})

test_that("randomForest: goldens for beta / se / V and lincomb", {
    skip_on_cran()
    skip_if_not_installed("randomForest")

    dat <- .make_tiny_dgm()

    set.seed(30301)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "rf",
        cross_fit = FALSE,
        verbose = FALSE,
        ntree = 300, mtry = 2, nodesize = 5
    )

    # cat("RF beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("RF beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("RF beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    expected_beta <- structure(c(0.150921347945536, -0.144731366766928), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_se <- structure(c(0.0692415601249499, 0.0946123423050273), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_V <- matrix(c(0.00479439364853705, -0.00320667976988224, -0.00320667976988223, 0.00895149531644365), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    ) # PLACEHOLDER

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    s <- summary(fit, lincomb = c(0, 1), conf_level = 0.95)
    df <- fit$df
    tcrit <- stats::qt(0.975, df)
    est <- expected_beta["Z"]
    se <- expected_se["Z"]
    expect_equal(as.numeric(s$lincomb$Estimate), as.numeric(est), tolerance = 1e-8)
    expect_equal(as.numeric(s$lincomb$`Std. Error`), as.numeric(se), tolerance = 1e-8)
})

test_that("ranger: goldens for beta / se / V and lincomb", {
    skip_on_cran()
    skip_if_not_installed("ranger")

    dat <- .make_tiny_dgm()

    set.seed(40401)
    fit <- dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "ranger",
        cross_fit = FALSE,
        verbose = FALSE,
        num.trees = 400, mtry = 2, min.node.size = 5, seed = 40401
    )

    # cat("ranger beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("ranger beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("ranger beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    expected_beta <- structure(c(0.147337094915664, -0.131825565893961), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_se <- structure(c(0.0694862479952816, 0.0954002422648467), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_V <- matrix(c(0.00482833866046178, -0.00330272964784973, -0.00330272964784973, 0.00910120622419144), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    ) # PLACEHOLDER

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    s <- summary(fit, lincomb = c(0, 1), conf_level = 0.95)
    expect_equal(nrow(s$lincomb), 1L)
})

test_that("SuperLearner: goldens for beta / se / V and lincomb", {
    skip_on_cran()
    skip_if_not_installed("SuperLearner")

    suppressMessages(library(SuperLearner))

    dat <- .make_tiny_dgm()

    # Use a deterministic, minimal library; disable shuffling in CV control

    set.seed(50501)
    fit <- suppressPackageStartupMessages(dcee(
        data = dat,
        id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
        moderator_formula = ~Z,
        control_formula = ~ X + Z,
        availability = "avail",
        control_reg_method = "sl",
        cross_fit = FALSE,
        verbose = FALSE,
        # pass through to SuperLearner via ...
        cvControl = SuperLearner::SuperLearner.CV.control(V = 3, shuffle = FALSE)
    ))

    # cat("SL beta_hat = "); dput(unclass(fit$fit$beta_hat)); cat("\n")
    # cat("SL beta_se  = "); dput(unclass(fit$fit$beta_se));  cat("\n")
    # cat("SL beta_V   = "); dput(fit$fit$beta_varcov);       cat("\n")

    expected_beta <- structure(c(0.157693969082668, -0.0905976419823446), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_se <- structure(c(0.089170606422377, 0.109801833613013), .Names = c("Intercept", "Z")) # PLACEHOLDER
    expected_V <- matrix(c(0.00795139704973446, -0.00546213131030219, -0.00546213131030219, 0.0120564426647797), 2,
        byrow = TRUE,
        dimnames = list(c("Intercept", "Z"), c("Intercept", "Z"))
    ) # PLACEHOLDER

    expect_equal(unclass(fit$fit$beta_hat), expected_beta, tolerance = 1e-8)
    expect_equal(unclass(fit$fit$beta_se), expected_se, tolerance = 1e-8)
    expect_equal(fit$fit$beta_varcov, expected_V, tolerance = 1e-8)

    s <- summary(fit, lincomb = c(0, 1))
    expect_equal(nrow(s$lincomb), 1L)
})

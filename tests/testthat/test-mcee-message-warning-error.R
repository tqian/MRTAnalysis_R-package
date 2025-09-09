# mcee family: warnings, errors, and edge cases"

## ---------- Minimal DGM (fast) ----------
make_toy <- function() {
    n <- 6
    Ti <- c(4, 5, 4, 6, 5, 4) # unequal T_i
    id <- rep(seq_len(n), Ti)
    dp <- unlist(lapply(Ti, seq_len))

    set.seed(1)
    I <- rbinom(length(dp), 1, 0.85) # has zeros
    A <- rbinom(length(dp), 1, 0.5)
    M <- rbinom(length(dp), 1, plogis(-0.2 + 0.4 * A + 0.1 * scale(dp)))
    # Y constant in id
    y_lin <- 0.5 * A + 0.6 * M + 0.05 * scale(dp)
    Y <- ave(y_lin, id, FUN = function(v) rep(mean(v), length(v)))
    data.frame(id, dp, I, A, M, Y)
}
dat <- make_toy()

## nontrivial per-row weights (normalized within id)
w_raw <- 0.3 + 0.7 * (dat$dp / ave(dat$dp, dat$id))
w <- ave(w_raw, dat$id, FUN = function(v) v / sum(v))

## ---------- 1) mcee: input validation ----------

test_that("mcee: availability not provided -> message; binary checks; dp & outcome checks", {
    d <- dat[, c("id", "dp", "A", "M", "Y")]
    messages <- capture.output(
        fit <- mcee(
            data = d, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = rep(1, nrow(d)),
            verbose = TRUE
        ),
        type = "message"
    )
    expect_true(any(grepl("assuming all rows available", messages)))

    # treatment not binary
    d2 <- dat
    d2$A <- 2L
    expect_error(
        mcee(
            data = d2, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        "must be coded 0/1"
    )

    # availability not binary
    d3 <- dat
    d3$I <- 7L
    expect_error(
        mcee(
            data = d3, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        "availability.*coded 0/1"
    )

    # outcome not constant within id
    d4 <- dat
    d4$Y[d4$id == 1][1] <- d4$Y[d4$id == 1][1] + 1
    expect_error(
        mcee(
            data = d4, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        "must be constant within each subject"
    )

    # dp not strictly increasing within id
    d5 <- dat
    swap <- which(d5$id == 2)[2:3]
    d5$dp[swap] <- rev(d5$dp[swap])
    expect_error(
        mcee(
            data = d5, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
            availability = "I",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            rand_prob = 0.5,
            weight_per_row = w,
            verbose = FALSE
        ),
        "strictly increasing"
    )
})

test_that("mcee: moderator & control formula checks; rand_prob checks; weights", {
    # moderator must be RHS-only
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = Y ~ dp,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "RHS-only"
    )

    # moderator warns if using vars != dp
    expect_warning(
        suppressMessages(
            mcee(
                dat, "id", "dp", "Y", "A", "M", "I",
                rand_prob = 0.5,
                time_varying_effect_form = ~ dp + M,
                control_formula_with_mediator = ~ dp + M,
                control_reg_method = "glm",
                weight_per_row = w,
                verbose = TRUE
            )
        ),
        regexp = "includes variables beyond 'dp': M",
        fixed = TRUE
    )

    # control formula must be RHS-only and exclude treatment/outcome
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = Y ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "RHS-only"
    )
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M + A,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "must not include the treatment variable"
    )
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M + Y,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "must not include the outcome variable"
    )

    # rand_prob scalar must be in (0,1)
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 1.2,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "in \\(0,1\\)"
    )

    # rand_prob column must exist and be (0,1) where I==1
    d <- dat
    d$p <- 0.5
    d$p[d$I == 1][1] <- 0 # invalid at an available row
    expect_error(
        mcee(d, "id", "dp", "Y", "A", "M", "I",
            rand_prob = "p",
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = w, verbose = FALSE
        ),
        "must be in \\(0,1\\) where availability==1"
    )

    # weights: wrong length, negative, all zeros
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = rep(1, 10), verbose = FALSE
        ),
        "length nrow\\(data\\)"
    )
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = replace(w, 1, -1), verbose = FALSE
        ),
        "nonnegative"
    )
    expect_error(
        mcee(dat, "id", "dp", "Y", "A", "M", "I",
            rand_prob = 0.5,
            time_varying_effect_form = ~1,
            control_formula_with_mediator = ~ dp + M,
            control_reg_method = "glm",
            weight_per_row = rep(0, nrow(dat)), verbose = FALSE
        ),
        "cannot be all zeros"
    )
})

## ---------- 2) mcee_general: config errors & shared checks ----------

test_that("mcee_general: config validation & shared messages", {
    # missing required columns
    expect_error(
        mcee_general(dat[, c("id", "dp", "A", "M")], "id", "dp", "Y", "A", "M",
            availability = NULL,
            time_varying_effect_form = ~1,
            config_p = list(known = 0.5),
            config_q = list(method = "glm", formula = ~ dp + M),
            config_eta = list(method = "glm", formula = ~dp),
            config_mu = list(method = "glm", formula = ~ dp + M),
            config_nu = list(method = "glm", formula = ~dp),
            weight_per_row = rep(1, nrow(dat)), verbose = FALSE
        ),
        "Missing columns in `data`: Y"
    )

    # treatment not binary
    d2 <- dat
    d2$A <- 3L
    expect_error(
        mcee_general(d2, "id", "dp", "Y", "A", "M",
            availability = "I",
            time_varying_effect_form = ~1,
            config_p = list(known = 0.5),
            config_q = list(method = "glm", formula = ~ dp + M),
            config_eta = list(method = "glm", formula = ~dp),
            config_mu = list(method = "glm", formula = ~ dp + M),
            config_nu = list(method = "glm", formula = ~dp),
            weight_per_row = w, verbose = FALSE
        ),
        "treatment.*coded 0/1"
    )

    # availability not binary
    d3 <- dat
    d3$I <- 9L
    expect_error(
        mcee_general(d3, "id", "dp", "Y", "A", "M",
            availability = "I",
            time_varying_effect_form = ~1,
            config_p = list(known = 0.5),
            config_q = list(method = "glm", formula = ~ dp + M),
            config_eta = list(method = "glm", formula = ~dp),
            config_mu = list(method = "glm", formula = ~ dp + M),
            config_nu = list(method = "glm", formula = ~dp),
            weight_per_row = w, verbose = FALSE
        ),
        "availability.*coded 0/1"
    )

    # moderator RHS-only
    expect_error(
        mcee_general(dat, "id", "dp", "Y", "A", "M",
            availability = "I",
            time_varying_effect_form = Y ~ dp,
            config_p = list(known = 0.5),
            config_q = list(method = "glm", formula = ~ dp + M),
            config_eta = list(method = "glm", formula = ~dp),
            config_mu = list(method = "glm", formula = ~ dp + M),
            config_nu = list(method = "glm", formula = ~dp),
            weight_per_row = w, verbose = FALSE
        ),
        "RHS-only"
    )

    # missing method/formula in configs
    expect_error(
        mcee_general(dat, "id", "dp", "Y", "A", "M",
            availability = "I",
            time_varying_effect_form = ~1,
            config_p = list(), # invalid
            config_q = list(method = "glm", formula = ~ dp + M),
            config_eta = list(method = "glm", formula = ~dp),
            config_mu = list(method = "glm", formula = ~ dp + M),
            config_nu = list(method = "glm", formula = ~dp),
            weight_per_row = w, verbose = FALSE
        ),
        "No formula provided|No method provided|known"
    )
})

## ---------- 3) mcee_userfit_nuisance: vector/avail rules & warnings ----------

test_that("mcee_userfit_nuisance: supplied vectors length/values; availability overrides", {
    # Build sane predictions from GLMs for shape/length
    pfit <- glm(A ~ dp, data = dat, family = binomial())
    qfit <- glm(A ~ dp + M, data = dat, family = binomial())
    etaf <- glm(Y ~ A * dp, data = dat)
    muf <- glm(Y ~ A * (dp + M), data = dat)

    p1 <- predict(pfit, type = "response")
    p1[dat$I == 0] <- 1 # enforce p1=1 where I=0
    q1 <- predict(qfit, type = "response")
    q1[dat$I == 0] <- 1 # enforce q1=1 where I=0
    eta1 <- predict(etaf, newdata = transform(dat, A = 1), type = "response")
    eta0 <- predict(etaf, newdata = transform(dat, A = 0), type = "response")
    mu1 <- predict(muf, newdata = transform(dat, A = 1), type = "response")
    mu0 <- predict(muf, newdata = transform(dat, A = 0), type = "response")
    # quick nu regressions
    nu1 <- predict(glm(mu1 ~ dp, data = cbind(dat, mu1 = mu1), subset = (dat$A == 0)), newdata = dat, type = "response")
    nu0 <- predict(glm(mu0 ~ dp, data = cbind(dat, mu0 = mu0), subset = (dat$A == 1)), newdata = dat, type = "response")

    # wrong length
    expect_error(
        mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
            time_varying_effect_form = ~1,
            p1 = p1[-1], q1 = q1,
            eta1 = eta1, eta0 = eta0,
            mu1 = mu1, mu0 = mu0,
            nu1 = nu1, nu0 = nu0,
            weight_per_row = w, verbose = FALSE
        ),
        "length nrow\\(data\\)"
    )

    # invalid p1/q1 where I==1
    badp <- p1
    badp[dat$I == 1][1] <- 0
    expect_error(
        mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
            time_varying_effect_form = ~1,
            p1 = badp, q1 = q1,
            eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
            weight_per_row = w, verbose = FALSE
        ),
        "p1.*\\(0,1\\).*availability==1"
    )
    badq <- q1
    badq[dat$I == 1][2] <- 1
    expect_error(
        mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
            time_varying_effect_form = ~1,
            p1 = p1, q1 = badq,
            eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
            weight_per_row = w, verbose = FALSE
        ),
        "q1.*\\(0,1\\).*availability==1"
    )

    # availability==0 forces p1/q1==1 with warning
    # (manually violate this to trigger the warning)
    p1_bad <- p1
    p1_bad[dat$I == 0] <- 0.7
    q1_bad <- q1
    q1_bad[dat$I == 0] <- 0.2
    expect_warning(
        mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
            time_varying_effect_form = ~1,
            p1 = p1_bad, q1 = q1_bad,
            eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
            weight_per_row = w, verbose = FALSE
        ),
        "*availability==0.*overridden to 1"
    )

    # moderator formula warns if includes vars != dp
    expect_warning(
        suppressMessages(
            mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
                time_varying_effect_form = ~ dp + M,
                p1 = p1, q1 = q1, eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
                weight_per_row = w, verbose = TRUE
            )
        ),
        "includes variables beyond 'dp'"
    )

    # weight_per_row invalid
    expect_error(
        mcee_userfit_nuisance(dat, "id", "dp", "Y", "A", "M", "I",
            time_varying_effect_form = ~1,
            p1 = p1, q1 = q1, eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
            weight_per_row = rep(1, 10), verbose = FALSE
        ),
        "length nrow\\(data\\)"
    )
})

## ---------- 4) summary: nuisance printing for userfit ----------
test_that("summary(..., show_nuisance=TRUE) prints the userfit notice", {
    # A tiny successful userfit run
    set.seed(2)
    pfit <- glm(A ~ dp, data = dat, family = binomial())
    qfit <- glm(A ~ dp + M, data = dat, family = binomial())
    etaf <- glm(Y ~ A * dp, data = dat)
    muf <- glm(Y ~ A * (dp + M), data = dat)

    p1 <- predict(pfit, type = "response")
    p1[dat$I == 0] <- 1 # enforce p1=1 where I=0
    q1 <- predict(qfit, type = "response")
    q1[dat$I == 0] <- 1 # enforce p1=1 where I=0
    eta1 <- predict(etaf, newdata = transform(dat, A = 1), type = "response")
    eta0 <- predict(etaf, newdata = transform(dat, A = 0), type = "response")
    mu1 <- predict(muf, newdata = transform(dat, A = 1), type = "response")
    mu0 <- predict(muf, newdata = transform(dat, A = 0), type = "response")
    nu1 <- predict(glm(mu1 ~ dp, data = cbind(dat, mu1 = mu1), subset = (dat$A == 0)), newdata = dat, type = "response")
    nu0 <- predict(glm(mu0 ~ dp, data = cbind(dat, mu0 = mu0), subset = (dat$A == 1)), newdata = dat, type = "response")

    fit <- mcee_userfit_nuisance(
        dat, "id", "dp", "Y", "A", "M", "I",
        time_varying_effect_form = ~1,
        p1 = p1, q1 = q1, eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0,
        weight_per_row = w, verbose = FALSE
    )
    so <- capture.output(print(summary(fit, show_nuisance = TRUE)))
    expect_true(any(grepl("Fitted values for all nuisance functions were supplied by the user", so)))
})

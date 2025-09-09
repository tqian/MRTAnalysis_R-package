# correct return structure of mcee, mcee_general, mcee_user_fit_nuisance
test_that("mcee, mcee_general, mcee_userfit_nuisance return structure", {
    # ---- small helper to avoid repetition --------------------------------------
    expect_mcee_fit_structure <- function(fit, p_dim, data, dp_name = "dp",
                                          expect_nuisance_models = TRUE,
                                          fn_name = c("mcee", "mcee_general", "mcee_userfit_nuisance")) {
        fn_name <- match.arg(fn_name)

        expect_s3_class(fit, "mcee_fit")
        expect_true(is.list(fit))
        expect_true(all(c("call", "mcee_fit", "nuisance_models", "meta") %in% names(fit)))

        # meta
        expect_true(is.list(fit$meta))
        expect_true(all(c("p_dim", "n_ids", "T_per_id", "weight_per_row_used") %in% names(fit$meta)))
        expect_true(is.numeric(fit$meta$p_dim) && length(fit$meta$p_dim) == 1)
        expect_equal(fit$meta$p_dim, p_dim)
        expect_true(is.numeric(fit$meta$n_ids) && length(fit$meta$n_ids) == 1)
        expect_true(is.integer(fit$meta$T_per_id) && length(fit$meta$T_per_id) == fit$meta$n_ids)
        expect_true(is.numeric(fit$meta$weight_per_row_used) && length(fit$meta$weight_per_row_used) == nrow(data))
        expect_equal(sum(fit$meta$T_per_id), nrow(data))
        expect_equal(fit$meta$T_per_id, as.vector(table(data$id)))

        # stage-2 fit
        expect_true(is.list(fit$mcee_fit))
        expect_true(all(c(
            "alpha_hat", "alpha_se", "beta_hat", "beta_se",
            "varcov", "alpha_varcov", "beta_varcov"
        ) %in% names(fit$mcee_fit)))

        expect_equal(length(fit$mcee_fit$alpha_hat), p_dim)
        expect_equal(length(fit$mcee_fit$alpha_se), p_dim)
        expect_equal(length(fit$mcee_fit$beta_hat), p_dim)
        expect_equal(length(fit$mcee_fit$beta_se), p_dim)

        # basis names from model.matrix(~ dp, data)
        basis_names <- colnames(model.matrix(reformulate(dp_name, response = NULL), data = data))
        expect_identical(names(fit$mcee_fit$alpha_hat), basis_names)
        expect_identical(names(fit$mcee_fit$alpha_se), basis_names)
        expect_identical(names(fit$mcee_fit$beta_hat), basis_names)
        expect_identical(names(fit$mcee_fit$beta_se), basis_names)

        # varcov shapes/names
        expect_true(is.matrix(fit$mcee_fit$varcov))
        expect_equal(nrow(fit$mcee_fit$varcov), 2 * p_dim)
        expect_equal(ncol(fit$mcee_fit$varcov), 2 * p_dim)
        expect_identical(
            rownames(fit$mcee_fit$varcov),
            c(paste0("alpha_", basis_names), paste0("beta_", basis_names))
        )
        expect_identical(
            colnames(fit$mcee_fit$varcov),
            c(paste0("alpha_", basis_names), paste0("beta_", basis_names))
        )

        expect_true(is.matrix(fit$mcee_fit$alpha_varcov))
        expect_equal(dim(fit$mcee_fit$alpha_varcov), c(p_dim, p_dim))
        expect_identical(rownames(fit$mcee_fit$alpha_varcov), basis_names)
        expect_identical(colnames(fit$mcee_fit$alpha_varcov), basis_names)

        expect_true(is.matrix(fit$mcee_fit$beta_varcov))
        expect_equal(dim(fit$mcee_fit$beta_varcov), c(p_dim, p_dim))
        expect_identical(rownames(fit$mcee_fit$beta_varcov), basis_names)
        expect_identical(colnames(fit$mcee_fit$beta_varcov), basis_names)

        # call sanity
        expect_identical(fit$call[[1]], as.name(fn_name))
        expect_identical(fit$call$data, as.name("data"))

        # nuisance_models presence depends on wrapper
        if (isTRUE(expect_nuisance_models)) {
            expect_true(is.list(fit$nuisance_models))
            expect_true(all(c("p", "q", "eta1", "eta0", "mu1", "mu0", "nu1", "nu0")
            %in% names(fit$nuisance_models)))
        } else {
            expect_true(is.null(fit$nuisance_models) || is.character(fit$nuisance_models))
        }
    }


    # ---- generate a simple data ----------------------------------------------

    set.seed(12345)
    n <- 30
    T_val <- 10
    id <- rep(1:n, each = T_val)
    dp <- rep(1:T_val, times = n)
    A <- rbinom(n * T_val, 1, 0.5)
    M <- rbinom(n * T_val, 1, plogis(-0.5 + 0.3 * A + 0.1 * dp))

    # generate Y of length n
    Y_tmp <- rnorm(n * T_val, mean = 1 + 0.5 * A + 0.7 * M + 0.2 * dp, sd = 10)
    # take an average of Y_tmp within each id to create Y
    Y <- ave(Y_tmp, id)
    data <- data.frame(id = id, dp = dp, A = A, M = M, Y = Y)

    # ---- mcee ----------------------------------------------
    fit_mcee <- mcee(
        data = data,
        id = "id",
        dp = "dp",
        outcome = "Y",
        treatment = "A",
        mediator = "M",
        time_varying_effect_form = ~dp,
        control_formula_with_mediator = ~ dp + M,
        control_reg_method = "glm",
        rand_prob = 0.5,
        verbose = FALSE
    )
    expect_mcee_fit_structure(
        fit_mcee,
        p_dim = 2, data = data, dp_name = "dp",
        expect_nuisance_models = TRUE, fn_name = "mcee"
    )

    # a couple of mcee-specific call args
    expect_equal(fit_mcee$call$rand_prob, 0.5)
    expect_equal(fit_mcee$call$control_reg_method, "glm")
    expect_equal(fit_mcee$call$time_varying_effect_form, quote(~dp))
    expect_equal(fit_mcee$call$control_formula_with_mediator, quote(~ dp + M))

    # ---- mcee_general ----------------------------------------------
    config_p <- list(formula = ~dp, method = "glm", family = binomial())
    config_q <- list(formula = ~ dp + M, method = "glm", family = binomial())
    config_eta <- list(formula = ~dp, method = "glm", family = gaussian())
    config_mu <- list(formula = ~ dp + M, method = "glm", family = gaussian())
    config_nu <- list(formula = ~dp, method = "glm", family = gaussian())

    fit_mcee_general <- mcee_general(
        data = data, id = "id", dp = "dp", outcome = "Y",
        treatment = "A", mediator = "M",
        time_varying_effect_form = ~dp,
        config_p = config_p, config_q = config_q,
        config_eta = config_eta, config_mu = config_mu, config_nu = config_nu,
        verbose = FALSE
    )

    expect_mcee_fit_structure(
        fit_mcee_general,
        p_dim = 2, data = data, dp_name = "dp",
        expect_nuisance_models = TRUE, fn_name = "mcee_general"
    )

    # ---- mcee_userfit_nuisance -------------------------------------

    p1_fit <- glm(A ~ dp, data = data, family = binomial())
    p1_vec <- predict(p1_fit, type = "response")
    q1_fit <- glm(A ~ dp + M, data = data, family = binomial())
    q1_vec <- predict(q1_fit, newdata = data, type = "response")
    eta_fit <- glm(Y ~ A * dp, data = data)
    eta1_vec <- predict(eta_fit, newdata = transform(data, A = 1), type = "response")
    eta0_vec <- predict(eta_fit, newdata = transform(data, A = 0), type = "response")
    mu_fit <- glm(Y ~ A * (M + dp), data = data)
    mu1_vec <- predict(mu_fit, newdata = transform(data, A = 1), type = "response")
    mu0_vec <- predict(mu_fit, newdata = transform(data, A = 0), type = "response")
    nu1_fit <- glm(mu1 ~ dp, data = cbind(data, mu1 = mu1_vec), subset = (A == 0))
    nu1_vec <- predict(nu1_fit, newdata = data, type = "response")
    nu0_fit <- glm(mu0 ~ dp, data = cbind(data, mu0 = mu0_vec), subset = (A == 1))
    nu0_vec <- predict(nu0_fit, newdata = data, type = "response")

    fit_mcee_userfit <- mcee_userfit_nuisance(
        data = data, id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
        time_varying_effect_form = ~dp,
        p1 = p1_vec, q1 = q1_vec,
        eta1 = eta1_vec, eta0 = eta0_vec,
        mu1 = mu1_vec, mu0 = mu0_vec,
        nu1 = nu1_vec, nu0 = nu0_vec,
        verbose = FALSE
    )

    expect_mcee_fit_structure(
        fit_mcee_userfit,
        p_dim = 2, data = data, dp_name = "dp",
        expect_nuisance_models = FALSE, fn_name = "mcee_userfit_nuisance"
    )
})

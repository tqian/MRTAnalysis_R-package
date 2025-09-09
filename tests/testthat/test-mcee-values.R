test_that("mcee family: value tests across learners & wrappers", {
    # skip_on_cran()
    set.seed(1234)

    ## ---------- 1) DGM with unequal T_i and availability == 0 ---------
    make_dgm <- function(n = 40, Tmin = 5, Tmax = 9, pA = 0.6) {
        # Per-id length
        Ti <- sample(Tmin:Tmax, size = n, replace = TRUE)
        id <- rep(seq_len(n), times = Ti)
        dp <- unlist(lapply(Ti, function(Ti_i) seq_len(Ti_i)))

        # Availability with some zeros
        I <- rbinom(length(dp), 1, prob = 0.9)
        # Treatment assignment (we'll also use constant rand_prob = pA for mcee)
        A <- rbinom(length(dp), 1, prob = pA)

        # Mediator depends on A and time (mild nonlinearity)
        M_mean <- -0.2 + 0.4 * A + 0.1 * scale(dp)
        M <- rbinom(length(dp), 1, plogis(M_mean))

        # Distal outcome Y: constant within id (use a latent linear index per id)
        base_id <- rnorm(n, 0, 0.3)
        # Average signal across each id's rows (depends on A/M/dp) + id random intercept
        signal <- 0.5 * A + 0.6 * M + 0.08 * scale(dp)
        Y_tmp <- as.numeric(signal)
        Y_bar <- ave(Y_tmp, id, FUN = mean) + base_id[id]
        # Y must be constant within id:
        Y <- ave(Y_bar, id, FUN = function(v) rep(v[1], length(v)))

        data.frame(id, dp, I, A, M, Y, stringsAsFactors = FALSE)
    }

    dat <- make_dgm()

    # Non-trivial per-row weights: emphasize later decision points, normalize within id
    w_raw <- 0.3 + 0.7 * (dat$dp / ave(dat$dp, dat$id, FUN = max))
    w <- ave(w_raw, dat$id, FUN = function(v) v / sum(v)) # normalized within id
    stopifnot(all(abs(tapply(w, dat$id, sum) - 1) < 1e-12))

    # Moderator scenarios
    MODS <- list(
        mod1 = ~1,
        mod2 = ~dp
    )

    # Learners to test (no SL)
    # LEARNERS <- c("glm", "gam", "ranger")
    LEARNERS <- c("glm") # only checks GLM (as others may fail on CRAN due to package versions)

    # Helper: build configs for mcee_general (families auto)
    make_configs <- function(method) {
        list(
            p   = list(formula = ~dp, method = method), # binomial auto
            q   = list(formula = ~ dp + M, method = method), # binomial auto
            eta = list(formula = ~dp, method = method), # gaussian auto
            mu  = list(formula = ~ dp + M, method = method), # gaussian auto
            nu  = list(formula = ~dp, method = method) # gaussian auto
        )
    }

    # Helper to compare two mcee_fit lists numerically
    expect_fit_equal <- function(fa, fb, tol) {
        fields <- c("alpha_hat", "alpha_se", "beta_hat", "beta_se")
        for (nm in fields) {
            expect_equal(unname(fa[[nm]]), unname(fb[[nm]]), tolerance = tol)
        }
        expect_equal(unname(fa$varcov), unname(fb$varcov), tolerance = tol)
    }

    # Tolerances per learner
    TOL <- c(glm = 1e-8, gam = 1e-6, ranger = 1e-5)

    # Storage for first-run values you’ll paste into EXPECTED
    # EXPECTED[[method]][[modkey]] <- list(alpha_hat=..., alpha_se=..., beta_hat=..., beta_se=..., varcov=...)
    EXPECTED <- new.env(parent = emptyenv())

    for (method in LEARNERS) {
        for (modkey in names(MODS)) {
            mod_fml <- MODS[[modkey]]

            ## ---- mcee (MRT; known p) ----
            set.seed(123)
            fit_mrt <- mcee(
                data = dat,
                id = "id",
                dp = "dp",
                outcome = "Y",
                treatment = "A",
                mediator = "M",
                availability = "I",
                rand_prob = 0.6, # constant randomization prob
                time_varying_effect_form = mod_fml,
                control_formula_with_mediator = ~ dp + M,
                control_reg_method = method,
                weight_per_row = w,
                verbose = FALSE
            )

            ## ---- mcee_general (configs) ----
            cfg <- make_configs(method)
            set.seed(123)
            fit_gen <- mcee_general(
                data = dat,
                id = "id", dp = "dp", outcome = "Y",
                treatment = "A", mediator = "M",
                availability = "I",
                time_varying_effect_form = mod_fml,
                config_p = cfg$p, config_q = cfg$q,
                config_eta = cfg$eta, config_mu = cfg$mu, config_nu = cfg$nu,
                weight_per_row = w,
                verbose = FALSE
            )

            cfg$p <- list(known = 0.6)
            set.seed(123)
            fit_gen_constant_rand_prob <- mcee_general(
                data = dat,
                id = "id", dp = "dp", outcome = "Y",
                treatment = "A", mediator = "M",
                availability = "I",
                time_varying_effect_form = mod_fml,
                config_p = cfg$p, config_q = cfg$q,
                config_eta = cfg$eta, config_mu = cfg$mu, config_nu = cfg$nu,
                weight_per_row = w,
                verbose = FALSE
            )

            ## ---- mcee_userfit_nuisance (inject predictions) ----
            set.seed(123)
            fit_usr <- mcee_userfit_nuisance(
                data = dat,
                id = "id", dp = "dp", outcome = "Y",
                treatment = "A", mediator = "M",
                availability = "I",
                time_varying_effect_form = mod_fml,
                p1 = fit_gen$nuisance_fitted$p1,
                q1 = fit_gen$nuisance_fitted$q1,
                eta1 = fit_gen$nuisance_fitted$eta1,
                eta0 = fit_gen$nuisance_fitted$eta0,
                mu1 = fit_gen$nuisance_fitted$mu1,
                mu0 = fit_gen$nuisance_fitted$mu0,
                nu1 = fit_gen$nuisance_fitted$nu1,
                nu0 = fit_gen$nuisance_fitted$nu0,
                weight_per_row = w,
                verbose = FALSE
            )

            # (2) Wrapper equivalence under matched nuisances
            expect_fit_equal(fit_mrt$mcee_fit, fit_gen_constant_rand_prob$mcee_fit, tol = TOL[[method]])
            expect_fit_equal(fit_gen$mcee_fit, fit_usr$mcee_fit, tol = TOL[[method]])

            # (1) Stored numeric expectations for *each* wrapper (paste after your first run)
            key <- paste(method, modkey, sep = "_")

            # --------- HOW TO FILL EXPECTED VALUES ---------
            # 1) Run this test locally with the 3 lines below temporarily uncommented to print:
            # print(list(method=method, modkey=modkey,
            # mrt = fit_mrt$mcee_fit, gen = fit_gen$mcee_fit))
            # 2) Copy the numeric vectors/matrices into the blocks below:
            if (!exists("EXPECTED")) EXPECTED <- list()

            ## ========== glm_mod1 ==========
            EXPECTED[["glm_mod1"]] <- list()

            EXPECTED[["glm_mod1"]]$mrt <- list(
                alpha_hat = c(0.06625261791),
                alpha_se = c(0.06027190583),
                beta_hat = c(0.006380160272),
                beta_se = c(0.01900838516),
                varcov = matrix(c(
                    0.0036327026322, -0.0003823088512,
                    -0.0003823088512, 0.0003613187065
                ), nrow = 2, byrow = TRUE)
            )

            EXPECTED[["glm_mod1"]]$gen <- list(
                alpha_hat = c(0.06507440706),
                alpha_se = c(0.04996718889),
                beta_hat = c(0.005941104296),
                beta_se = c(0.006747936183),
                varcov = matrix(c(
                    0.002496719966, -0.00006425975332,
                    -0.00006425975332, 0.00004553464273
                ), nrow = 2, byrow = TRUE)
            )

            ## ========== glm_mod2 ==========
            EXPECTED[["glm_mod2"]] <- list()

            EXPECTED[["glm_mod2"]]$mrt <- list(
                alpha_hat = c(0.030017337767, 0.007850867826),
                alpha_se = c(0.12801738044, 0.02227183982),
                beta_hat = c(0.041327376655, -0.007571791237),
                beta_se = c(0.038186545710, 0.005928110977),
                varcov = matrix(c(
                    0.0163884496956, -0.0025290626819, -0.0011304776089, 0.0001452819716,
                    -0.0025290626819, 0.0004960348490, 0.0001103524672, -0.00001951057729,
                    -0.0011304776089, 0.0001103524672, 0.0014582122732, -0.0002009640781,
                    0.0001452819716, -0.00001951057729, -0.0002009640781, 0.00003514249975
                ), nrow = 4, byrow = TRUE)
            )

            EXPECTED[["glm_mod2"]]$gen <- list(
                alpha_hat = c(0.033752155192, 0.006786393217),
                alpha_se = c(0.10416435485, 0.01860692425),
                beta_hat = c(0.031295748024, -0.005493429494),
                beta_se = c(0.019296968966, 0.003880748757),
                varcov = matrix(c(
                    0.010850212821, -0.0017057894733, -0.00001932318063, 0.00001730666023,
                    -0.0017057894733, 0.0003462176301, -0.00004637807156, 0.000004525942008,
                    -0.00001932318063, -0.00004637807156, 0.0003723730113, -0.00007019736136,
                    0.00001730666023, 0.000004525942008, -0.00007019736136, 0.00001506021092
                ), nrow = 4, byrow = TRUE)
            )

            ## ========== gam_mod1 ==========
            EXPECTED[["gam_mod1"]] <- list()

            EXPECTED[["gam_mod1"]]$mrt <- list(
                alpha_hat = c(0.06625261787),
                alpha_se = c(0.06027190580),
                beta_hat = c(0.006380160315),
                beta_se = c(0.01900838512),
                varcov = matrix(c(
                    0.0036327026294, -0.0003823088489,
                    -0.0003823088489, 0.0003613187048
                ), nrow = 2, byrow = TRUE)
            )

            EXPECTED[["gam_mod1"]]$gen <- list(
                alpha_hat = c(0.06507440702),
                alpha_se = c(0.04996718887),
                beta_hat = c(0.005941104332),
                beta_se = c(0.006747936188),
                varcov = matrix(c(
                    0.002496719964, -0.00006425975232,
                    -0.00006425975232, 0.00004553464280
                ), nrow = 2, byrow = TRUE)
            )

            ## ========== gam_mod2 ==========
            EXPECTED[["gam_mod2"]] <- list()

            EXPECTED[["gam_mod2"]]$mrt <- list(
                alpha_hat = c(0.030017337812, 0.007850867807),
                alpha_se = c(0.12801738045, 0.02227183982),
                beta_hat = c(0.041327376610, -0.007571791218),
                beta_se = c(0.038186545819, 0.005928110986),
                varcov = matrix(c(
                    0.0163884496965, -0.0025290626818, -0.0011304776136, 0.0001452819728,
                    -0.0025290626818, 0.0004960348488, 0.0001103524672, -0.00001951057723,
                    -0.0011304776136, 0.0001103524672, 0.0014582122816, -0.0002009640794,
                    0.0001452819728, -0.00001951057723, -0.0002009640794, 0.00003514249986
                ), nrow = 4, byrow = TRUE)
            )

            EXPECTED[["gam_mod2"]]$gen <- list(
                alpha_hat = c(0.03375215523, 0.00678639320),
                alpha_se = c(0.10416435485, 0.01860692425),
                beta_hat = c(0.031295747983, -0.005493429477),
                beta_se = c(0.019296969015, 0.003880748759),
                varcov = matrix(c(
                    0.010850212822, -0.0017057894731, -0.00001932318189, 0.00001730666089,
                    -0.0017057894731, 0.0003462176299, -0.00004637807215, 0.000004525942100,
                    -0.00001932318189, -0.00004637807215, 0.0003723730132, -0.00007019736159,
                    0.00001730666089, 0.000004525942100, -0.00007019736159, 0.00001506021093
                ), nrow = 4, byrow = TRUE)
            )

            ## ========== ranger_mod1 ==========
            EXPECTED[["ranger_mod1"]] <- list()

            EXPECTED[["ranger_mod1"]]$mrt <- list(
                alpha_hat = c(0.06130214535),
                alpha_se = c(0.05857425950),
                beta_hat = c(0.01222087510),
                beta_se = c(0.01518365075),
                varcov = matrix(c(
                    0.0034309438754, -0.0002274185378,
                    -0.0002274185378, 0.0002305432500
                ), nrow = 2, byrow = TRUE)
            )

            EXPECTED[["ranger_mod1"]]$gen <- list(
                alpha_hat = c(0.05691244928),
                alpha_se = c(0.05063385737),
                beta_hat = c(0.01411860946),
                beta_se = c(0.008340159273),
                varcov = matrix(c(
                    0.002563787512, -0.00008885874108,
                    -0.00008885874108, 0.00006955825669
                ), nrow = 2, byrow = TRUE)
            )

            ## ========== ranger_mod2 ==========
            EXPECTED[["ranger_mod2"]] <- list()

            EXPECTED[["ranger_mod2"]]$mrt <- list(
                alpha_hat = c(0.029294096898, 0.006934980403),
                alpha_se = c(0.12601744745, 0.02229737329),
                beta_hat = c(0.043447618320, -0.006765699965),
                beta_se = c(0.034762947526, 0.005476128531),
                varcov = matrix(c(
                    0.0158803970621, -0.0024979125103, -0.0008423222012, 0.0001256392113,
                    -0.0024979125103, 0.0004971728557, 0.0001034727549, -0.00002022514408,
                    -0.0008423222012, 0.0001034727549, 0.0012084625207, -0.0001758076070,
                    0.0001256392113, -0.00002022514408, -0.0001758076070, 0.00002998798369
                ), nrow = 4, byrow = TRUE)
            )

            EXPECTED[["ranger_mod2"]]$gen <- list(
                alpha_hat = c(0.029964233772, 0.005838698563),
                alpha_se = c(0.10627577575, 0.01970689175),
                beta_hat = c(0.040487353902, -0.005713148251),
                beta_se = c(0.022361293445, 0.004124347117),
                varcov = matrix(c(
                    0.011294540512, -0.0018436034528, 0.0001528413129, -0.00001538939788,
                    -0.0018436034528, 0.0003883615824, -0.00004995490442, 0.000003119285379,
                    0.0001528413129, -0.00004995490442, 0.0005000274445, -0.00008579069255,
                    -0.00001538939788, 0.000003119285379, -0.00008579069255, 0.00001701023914
                ), nrow = 4, byrow = TRUE)
            )
            # 3) Re-run tests; from then on, CI will check exact values (within tolerance).

            if (exists(key, envir = EXPECTED)) {
                exp <- get(key, envir = EXPECTED)

                # Compare to stored expectations
                expect_fit_equal(fit_mrt$mcee_fit, exp$mrt, tol = TOL[[method]])
                expect_fit_equal(fit_gen$mcee_fit, exp$gen, tol = TOL[[method]])
            } else {
                # No stored expectations yet; provide a gentle hint
                message(sprintf(
                    "[mcee value tests] No EXPECTED values for key='%s'. Run locally, print, and paste results into EXPECTED.",
                    key
                ))
            }
        }
    }
})

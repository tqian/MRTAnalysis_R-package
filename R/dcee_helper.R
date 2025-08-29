#' Stage-1 (nuisance) + Stage-2 (core) driver for DCEE
#'
#' Orchestrates the two-stage distal causal excursion effect (DCEE) estimation:
#' (i) fits nuisance outcome regressions \eqn{\mu_a(H_t)} with a chosen learner
#' (parametric or ML), optionally with cross-fitting; (ii) calls the core
#' estimating-equations routine to obtain \eqn{\beta}, its sandwich variance,
#' and confidence intervals.
#'
#' @param dta A data.frame in long format; one row per decision point.
#' @param id_var Character scalar; column name for subject identifier.
#' @param moderator_formula RHS-only formula for excursion-effect moderators
#'   (e.g., `~ 1`, `~ Z`, `~ Z1 + Z2`). If `~ 1`, the marginal effect is estimated.
#' @param trt_var Character scalar; column name for treatment \code{A_t} (must be 0/1).
#' @param outcome_var Character scalar; column name for (distal) continuous outcome \code{Y}.
#' @param prob_A_var Character scalar; column name for randomization probability
#'   \eqn{P(A_t=1 \mid H_t)}; values must lie strictly in \code{(0, 1)}.
#' @param avail_var Character scalar; column name for availability indicator (0/1).
#' @param control_reg_method Character scalar; nuisance learner to use for Stage 1.
#'   One of \code{"gam"}, \code{"lm"}, \code{"rf"}, \code{"ranger"},
#'   \code{"sl"}, \code{"sl.user-specified-library"}, \code{"set_to_zero"}.
#'   When \code{"gam"}, \code{control_formula} may include smooth terms like \code{s(x)}.
#'   When \code{"sl.user-specified-library"}, supply \code{sl.library} via \code{...}.
#' @param control_formula RHS-only formula specifying covariates used to learn
#'   \eqn{\mu_a(H_t)}. For SuperLearner methods, variables are extracted from this
#'   formula to form the design matrix \code{X}. Use \code{~ 1} only with
#'   \code{control_reg_method = "set_to_zero"}.
#' @param cross_fit Logical; if \code{TRUE}, perform K-fold cross-fitting by subject id.
#' @param cf_fold Integer; number of folds for cross-fitting (default often 10).
#' @param weighting_function_var Optional character scalar; column name giving
#'   decision-point weights \eqn{\omega_t}. If \code{NULL}, defaults to per-timepoint
#'   weight \code{1/T} (with constant \code{T} per subject).
#' @param ... Additional arguments forwarded to the selected learner
#'   (e.g., \code{num.trees}, \code{mtry} for forests; \code{sl.library} for
#'   SuperLearner; \code{family}, \code{knots} for \pkg{mgcv}).
#'
#' @details
#' Stage 1 fits two outcome regressions on the subsets \code{A=0} and \code{A=1}
#' and attaches predicted values \code{mu_hat_0}, \code{mu_hat_1} to \code{dta}.
#' Stage 2 uses these nuisance fits to solve estimating equations for
#' \eqn{\beta} (distal causal excursion effects) and computes a sandwich variance
#' with small-sample \eqn{t}-based confidence intervals.
#'
#' @return A list with at least the following components (used downstream by
#'   \code{summary.dcee}):
#' \itemize{
#'   \item \code{beta_hat} — named numeric vector of \eqn{\beta} estimates
#'   \item \code{beta_varcov} — variance–covariance matrix of \eqn{\beta}
#'   \item \code{beta_se} — standard errors
#'   \item \code{conf_int}, \code{conf_int_tquantile} — Wald and \eqn{t}-based CIs
#'   \item \code{regfit_a0}, \code{regfit_a1} — fitted Stage 1 models (may be \code{NULL} when cross-fitting)
#' }
#'
#' @seealso \code{\link{dcee}} for the user-facing wrapper,
#'   \code{\link{dcee_helper_stage1_fit_nuisance}} and \code{\link{dcee_helper_stage2_estimate_dcee}} for internals.
#'
#' @noRd

dcee_helper_2stage_estimation <- function(
    dta,
    id_var,
    moderator_formula,
    trt_var,
    outcome_var,
    prob_A_var,
    avail_var,
    control_reg_method,
    control_formula,
    cross_fit,
    cf_fold,
    weighting_function_var,
    ...) {
    dta <- as.data.frame(dta)

    # Stage 1: fit nuisance outcome regressions and attach mu_hat_0 / mu_hat_1

    outcome_reg_return <- dcee_helper_stage1_fit_nuisance(
        dta = dta,
        id_var = id_var,
        trt_var = trt_var,
        outcome_var = outcome_var,
        control_reg_method = control_reg_method,
        control_formula = control_formula,
        cross_fit = cross_fit,
        cf_fold = cf_fold,
        ...
    )

    dta <- outcome_reg_return$dta

    # Extract moderator variable names from RHS-only formula (~ 1 -> NULL)
    moderator_var <- all.vars(moderator_formula)
    if (length(moderator_var) == 0L) moderator_var <- NULL

    # Stage 2: estimating equations for beta
    fit <- dcee_helper_stage2_estimate_dcee(
        dta = dta,
        id_var = id_var,
        moderator_var = moderator_var,
        trt_var = trt_var,
        outcome_var = outcome_var,
        mu_hat_1_var = "mu_hat_1",
        mu_hat_0_var = "mu_hat_0",
        prob_A_var = prob_A_var,
        avail_var = avail_var,
        weighting_function_var = weighting_function_var
    )

    # Return with references to first-stage fits (if any)
    c(fit, list(
        regfit_a0 = outcome_reg_return$regfit_a0,
        regfit_a1 = outcome_reg_return$regfit_a1
    ))
}




#' Stage 1: fit nuisance outcome regressions for DCEE
#'
#' Fits the nuisance outcome regressions \eqn{\mu_a(H_t) = E[Y_t \mid H_t, A_t=a]}
#' for \eqn{a \in \{0,1\}} using a chosen learner (parametric or ML), optionally
#' with cross-fitting by subject id. The fitted functions are then used to
#' generate per–decision point predictions \code{mu_hat_0} and \code{mu_hat_1}
#' that are appended to \code{dta} for use in Stage 2.
#'
#' @param dta A data.frame in long format; one row per decision point.
#' @param id_var Character scalar; column name for subject identifier.
#' @param trt_var Character scalar; column name for treatment \code{A_t} (must be 0/1).
#' @param outcome_var Character scalar; column name for the (distal) continuous outcome \code{Y}.
#' @param control_reg_method Character scalar; learner used to model \eqn{\mu_a(H_t)}.
#'   One of \code{"gam"}, \code{"lm"}, \code{"rf"}, \code{"ranger"},
#'   \code{"sl"}, \code{"sl.user-specified-library"}, \code{"set_to_zero"}.
#'   When \code{"gam"}, \code{control_formula} may include smooth terms such as \code{s(x)}.
#'   When \code{"sl.user-specified-library"}, provide \code{sl.library} via \code{...}.
#'   When \code{"set_to_zero"}, both \code{mu_hat_0} and \code{mu_hat_1} are set to 0 and no models are fitted.
#' @param control_formula RHS-only formula specifying covariates used to learn
#'   \eqn{\mu_a(H_t)}. For SuperLearner methods, variables in this formula are
#'   used to build the design matrix \code{X}. Use \code{~ 1} only with
#'   \code{control_reg_method = "set_to_zero"}.
#' @param cross_fit Logical; if \code{TRUE}, perform K-fold cross-fitting by subject id:
#'   at each fold, fit the nuisance models on the training subjects and predict
#'   \code{mu_hat_0}, \code{mu_hat_1} on the held-out subjects.
#' @param cf_fold Integer; number of folds for cross-fitting (default 10).
#' @param ... Additional arguments forwarded to the selected learner
#'   (e.g., \code{num.trees}, \code{mtry} for forests; \code{sl.library} for
#'   SuperLearner; \code{family}, \code{knots} for \pkg{mgcv}).
#'
#' @details
#' The function fits separate models on the subsets \code{A=0} and \code{A=1}.
#' With cross-fitting enabled, subject ids are partitioned into \code{cf_fold}
#' folds; predictions are obtained on held-out subjects to reduce overfitting bias.
#' For SuperLearner-based methods, the package \pkg{SuperLearner} must be available.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{dta}: the input data with two added numeric columns
#'     \code{mu_hat_0} and \code{mu_hat_1} (predicted \eqn{\mu_0(H_t)} and \eqn{\mu_1(H_t)}).
#'   \item \code{regfit_a0}, \code{regfit_a1}: fitted learner objects for \code{A=0} and \code{A=1}
#'     (may be \code{NULL} when \code{cross_fit = TRUE} if per-fold fits are not retained,
#'     or when \code{control_reg_method = "set_to_zero"}).
#' }
#'
#' @seealso \code{\link{dcee_helper_stage2_estimate_dcee}} for Stage 2, and
#'   \code{\link{dcee_helper_2stage_estimation}} for the two-stage driver.
#'
#' @noRd

dcee_helper_stage1_fit_nuisance <- function(
    dta,
    id_var,
    trt_var,
    outcome_var,
    control_reg_method,
    control_formula,
    cross_fit,
    cf_fold,
    ... # forwarded to learners (e.g., num.trees, mtry, sl.library)
    ) {
    # Build outcome formula once; allow s(x) for GAM
    rhs <- as.character(control_formula)[2L]
    control_formula_with_outcome <- stats::as.formula(paste(outcome_var, "~", rhs))

    # Variables-only set for learners that use design matrices (e.g., SuperLearner)
    control_var <- all.vars(control_formula)

    # Ensure s() not used with non-GAM formula learners
    if (control_reg_method %in% c("lm", "rf", "ranger") && grepl("s\\s*\\(", rhs)) {
        stop("control_formula contains s(...), which is only supported when control_reg_method = 'gam'.")
    }

    # Handle SuperLearner library selection from ... when requested
    dots <- list(...)
    if (identical(control_reg_method, "sl")) {
        sl.library <- c("SL.mean", "SL.glm", "SL.earth")
    } else if (identical(control_reg_method, "sl.user-specified-library")) {
        sl.library <- dots$sl.library
        if (is.null(sl.library)) {
            stop("Provide `sl.library` via ... when control_reg_method = 'sl.user-specified-library'.")
        }
    } else {
        sl.library <- NULL
    }

    # set-to-zero shortcut
    if (identical(control_reg_method, "set_to_zero")) {
        dta$mu_hat_0 <- 0
        dta$mu_hat_1 <- 0
        regfit_a0 <- regfit_a1 <- NULL
        return(list(dta = dta, regfit_a0 = regfit_a0, regfit_a1 = regfit_a1))
    }

    # Split by treatment for training
    dta_a0 <- dta[dta[[trt_var]] == 0, , drop = FALSE]
    dta_a1 <- dta[dta[[trt_var]] == 1, , drop = FALSE]

    # Cross-fitting by unique id
    if (isTRUE(cross_fit)) {
        ids <- unique(dta[[id_var]])
        K <- min(cf_fold, length(ids))
        fold_id <- sample(rep_len(seq_len(K), length.out = length(ids)))

        dta_out <- dta[0, , drop = FALSE] # empty like dta

        for (k in seq_len(K)) {
            id_holdout <- ids[fold_id == k]
            id_train <- ids[fold_id != k]

            dta_holdout <- dta[dta[[id_var]] %in% id_holdout, , drop = FALSE]
            dta_train <- dta[dta[[id_var]] %in% id_train, , drop = FALSE]

            dta_train_a0 <- dta_train[dta_train[[trt_var]] == 0, , drop = FALSE]
            dta_train_a1 <- dta_train[dta_train[[trt_var]] == 1, , drop = FALSE]

            # Fit models on train splits
            if (control_reg_method == "gam") {
                regfit_a0 <- mgcv::gam(control_formula_with_outcome, data = dta_train_a0, ...)
                regfit_a1 <- mgcv::gam(control_formula_with_outcome, data = dta_train_a1, ...)
            } else if (control_reg_method == "lm") {
                regfit_a0 <- stats::lm(control_formula_with_outcome, data = dta_train_a0, ...)
                regfit_a1 <- stats::lm(control_formula_with_outcome, data = dta_train_a1, ...)
            } else if (control_reg_method == "rf") {
                regfit_a0 <- randomForest::randomForest(control_formula_with_outcome, data = dta_train_a0, ...)
                regfit_a1 <- randomForest::randomForest(control_formula_with_outcome, data = dta_train_a1, ...)
            } else if (control_reg_method == "ranger") {
                regfit_a0 <- ranger::ranger(formula = control_formula_with_outcome, data = dta_train_a0, ...)
                regfit_a1 <- ranger::ranger(formula = control_formula_with_outcome, data = dta_train_a1, ...)
            } else if (control_reg_method %in% c("sl", "sl.user-specified-library")) {
                if (!requireNamespace("SuperLearner", quietly = TRUE)) {
                    stop("SuperLearner package is required for control_reg_method = 'sl'/'sl.user-specified-library'.")
                }
                Y0 <- dta_train_a0[[outcome_var]]
                X0 <- dta_train_a0[, control_var, drop = FALSE]
                Y1 <- dta_train_a1[[outcome_var]]
                X1 <- dta_train_a1[, control_var, drop = FALSE]
                regfit_a0 <- SuperLearner::SuperLearner(
                    Y = Y0, X = X0, family = stats::gaussian(),
                    verbose = FALSE, SL.library = sl.library, ...
                )
                regfit_a1 <- SuperLearner::SuperLearner(
                    Y = Y1, X = X1, family = stats::gaussian(),
                    verbose = FALSE, SL.library = sl.library, ...
                )
            }

            # Predict for holdout ids
            if (control_reg_method %in% c("sl", "sl.user-specified-library")) {
                newdf <- dta_holdout[, control_var, drop = FALSE]
                pred0 <- stats::predict(regfit_a0, newdata = newdf)$pred
                pred1 <- stats::predict(regfit_a1, newdata = newdf)$pred
            } else if (control_reg_method == "ranger") {
                pred0 <- predict(regfit_a0, data = dta_holdout)$predictions
                pred1 <- predict(regfit_a1, data = dta_holdout)$predictions
            } else if (control_reg_method %in% c("gam", "lm", "rf")) {
                pred0 <- stats::predict(regfit_a0, newdata = dta_holdout, type = "response")
                pred1 <- stats::predict(regfit_a1, newdata = dta_holdout, type = "response")
            }

            dta_holdout$mu_hat_0 <- as.numeric(pred0)
            dta_holdout$mu_hat_1 <- as.numeric(pred1)
            dta_out <- rbind(dta_out, dta_holdout)
        }

        # reorder to original row order just in case
        dta_out <- dta_out[match(seq_len(nrow(dta)), rownames(dta_out)), , drop = FALSE]
        # regfit_a0 <- regfit_a1 <- NULL # comment out if want to return the last-fold fit
        dta <- dta_out
    } else {
        # Fit on full data (no cross-fitting)
        if (control_reg_method == "gam") {
            regfit_a0 <- mgcv::gam(control_formula_with_outcome, data = dta_a0, ...)
            regfit_a1 <- mgcv::gam(control_formula_with_outcome, data = dta_a1, ...)
        } else if (control_reg_method == "lm") {
            regfit_a0 <- stats::lm(control_formula_with_outcome, data = dta_a0, ...)
            regfit_a1 <- stats::lm(control_formula_with_outcome, data = dta_a1, ...)
        } else if (control_reg_method == "rf") {
            regfit_a0 <- randomForest::randomForest(control_formula_with_outcome, data = dta_a0, ...)
            regfit_a1 <- randomForest::randomForest(control_formula_with_outcome, data = dta_a1, ...)
        } else if (control_reg_method == "ranger") {
            regfit_a0 <- ranger::ranger(formula = control_formula_with_outcome, data = dta_a0, ...)
            regfit_a1 <- ranger::ranger(formula = control_formula_with_outcome, data = dta_a1, ...)
        } else if (control_reg_method %in% c("sl", "sl.user-specified-library")) {
            if (!requireNamespace("SuperLearner", quietly = TRUE)) {
                stop("SuperLearner package is required for control_reg_method = 'sl'/'sl.user-specified-library'.")
            }
            Y0 <- dta_a0[[outcome_var]]
            X0 <- dta_a0[, control_var, drop = FALSE]
            Y1 <- dta_a1[[outcome_var]]
            X1 <- dta_a1[, control_var, drop = FALSE]
            regfit_a0 <- SuperLearner::SuperLearner(
                Y = Y0, X = X0, family = stats::gaussian(),
                verbose = FALSE, SL.library = sl.library, ...
            )
            regfit_a1 <- SuperLearner::SuperLearner(
                Y = Y1, X = X1, family = stats::gaussian(),
                verbose = FALSE, SL.library = sl.library, ...
            )
        }

        # Predict on all rows
        if (control_reg_method %in% c("sl", "sl.user-specified-library")) {
            newdf <- dta[, control_var, drop = FALSE]
            pred0 <- stats::predict(regfit_a0, newdata = newdf)$pred
            pred1 <- stats::predict(regfit_a1, newdata = newdf)$pred
        } else if (control_reg_method == "ranger") {
            pred0 <- predict(regfit_a0, data = dta)$predictions
            pred1 <- predict(regfit_a1, data = dta)$predictions
        } else if (control_reg_method %in% c("gam", "lm", "rf")) {
            pred0 <- stats::predict(regfit_a0, newdata = dta, type = "response")
            pred1 <- stats::predict(regfit_a1, newdata = dta, type = "response")
        }

        dta$mu_hat_0 <- as.numeric(pred0)
        dta$mu_hat_1 <- as.numeric(pred1)
    }

    list(dta = dta, regfit_a0 = regfit_a0, regfit_a1 = regfit_a1)
}


#' Stage 2: estimate distal causal excursion effects (DCEE)
#'
#' Solves the estimating equations for the distal causal excursion effect
#' parameters \eqn{\beta} using the Stage-1 nuisance predictions
#' \eqn{\mu_0(H_t)} and \eqn{\mu_1(H_t)}. Inference is based on a
#' cluster-robust (by subject id) sandwich variance.
#'
#' @param dta A data.frame in long format; one row per decision point.
#'   Must contain at least the columns named by the arguments below.
#' @param id_var Character scalar; column name for subject identifier.
#' @param moderator_var Character vector of column names to be used as
#'   moderators in the excursion effect model. The Stage-2 design matrix is
#'   \code{cbind(Intercept = 1, dta[, moderator_var])}. Use \code{NULL} or
#'   \code{character(0)} for a marginal (intercept-only) effect.
#' @param trt_var Character scalar; column name for treatment \code{A_t}
#'   (must be coded \code{0/1}).
#' @param outcome_var Character scalar; column name for the (distal) outcome \code{Y}
#'   (constant within subject in typical DCEE setups).
#' @param mu_hat_1_var,mu_hat_0_var Character scalars; column names for the
#'   Stage-1 predictions \code{mu_hat_1 = \eqn{\mu_1(H_t)}}, \code{mu_hat_0 = \eqn{\mu_0(H_t)}}.
#' @param prob_A_var Character scalar; column name for the randomization probability
#'   \eqn{p_t = P(A_t=1 \mid H_t)}. Values must lie strictly in \code{(0,1)}.
#' @param avail_var Character scalar; column name for availability indicator (0/1).
#' @param weighting_function_var Optional character scalar; column name containing
#'   decision-point weights \eqn{\omega_t}. If \code{NULL}, a default per-timepoint
#'   weight of \code{1/T} is used, where \code{T} is the (assumed constant) number
#'   of decision points per subject.
#'
#' @details
#' Let \eqn{S_t = (1, M_t)^\top} denote the moderator vector (including the
#' intercept), and define \eqn{p_t(1)=p_t}, \eqn{p_t(0)=1-p_t}. The estimating
#' equation is built from the per-timepoint contribution
#' \deqn{
#'   U_t(\beta) \;=\; \omega_t \, S_t \Big[ \,(-1)^{1-A_t}\, \frac{\mathrm{avail}_t}{P(A_t \mid H_t)}\,
#'      \big\{ Y - p_t(0)\mu_1(H_t) - p_t(1)\mu_0(H_t) \big\}
#'      \;-\; S_t^\top \beta \Big],
#' }
#' and \eqn{\hat\beta} solves \eqn{\sum_t U_t(\beta) = 0}. The sandwich variance
#' is computed by clustering at the subject level (summing contributions within
#' id, then taking a cross-product), yielding a standard “bread–meat–bread”
#' estimator. Small-sample \eqn{t}-based intervals can be formed with
#' \code{df = #subjects - #betas} (usually computed by the caller).
#'
#' @return A list with components (consumed by \code{summary.dcee} and the wrapper):
#' \itemize{
#'   \item \code{beta_hat} — named numeric vector of \eqn{\beta} estimates (length \code{1 + length(moderator_var)})
#'   \item \code{beta_varcov} — variance–covariance matrix of \eqn{\beta}
#'   \item \code{beta_se} — standard errors (square roots of diagonal of \code{beta_varcov})
#'   \item \code{conf_int} — large-sample Wald CIs (normal critical value)
#'   \item \code{conf_int_tquantile} — \eqn{t}-based CIs using degrees of freedom supplied by the caller
#' }
#'
#' @seealso \code{\link{dcee_helper_stage1_fit_nuisance}} for Stage 1 nuisance fitting,
#'   and \code{\link{dcee_helper_2stage_estimation}} for the two-stage driver.
#'
#' @noRd


dcee_helper_stage2_estimate_dcee <- function(dta,
                                             id_var,
                                             moderator_var,
                                             trt_var,
                                             outcome_var,
                                             mu_hat_1_var,
                                             mu_hat_0_var,
                                             prob_A_var,
                                             avail_var,
                                             weighting_function_var) {
    total_person_decisionpoint <- nrow(dta)
    sample_size <- length(unique(dta[[id_var]]))
    total_T <- total_person_decisionpoint / sample_size # assuming everyone has the same number of dp

    if (is.null(weighting_function_var)) {
        omega <- rep(1 / total_T, total_person_decisionpoint)
    } else {
        omega <- dta[[weighting_function_var]]
    }

    pt1 <- dta[, prob_A_var] # this is p_t(1|H_t)
    A <- dta[, trt_var]
    Sdm <- as.matrix(cbind(rep(1, nrow(dta)), dta[, moderator_var]))
    Y <- dta[, outcome_var]
    mu_hat_1 <- dta[, mu_hat_1_var]
    mu_hat_0 <- dta[, mu_hat_0_var]
    avail <- dta[, avail_var]

    pt0 <- 1 - pt1 # this is p_t(0|H_t)
    pA <- ifelse(A, pt1, pt0) # this is p_t(A_t|H_t)

    dim_beta <- ncol(Sdm)

    person_first_index <- c(find_change_location(dta[, id_var]), total_person_decisionpoint + 1)
    sample_size <- length(person_first_index) - 1

    #### 1. Calculate beta_hat using its analytic form ####
    # from Goodnotes "2025.01.19 - Asymptotic Variance (updated)"

    ## compute the front inverse term

    sum_SStrans <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (it in 1:total_person_decisionpoint) {
        S_it <- matrix(Sdm[it, ], ncol = 1)
        sum_SStrans <- sum_SStrans + omega[it] * S_it %*% t(S_it)
    }
    avg_SStrans <- sum_SStrans / sample_size

    ## compute the second term (involving a residual)
    residual <- omega * avail * (-1)^(1 - A) / pA * (Y - pt0 * mu_hat_1 - pt1 * mu_hat_0)
    sum_residualS <- matrix(0, nrow = dim_beta, ncol = 1)
    for (it in 1:total_person_decisionpoint) {
        S_it <- matrix(Sdm[it, ], ncol = 1)
        r_it <- as.numeric(residual[it])
        sum_residualS <- sum_residualS + r_it * S_it
    }
    avg_residualS <- sum_residualS / sample_size

    ## compute beta_hat

    beta_hat <- solve(avg_SStrans) %*% avg_residualS # a single column matrix

    #### 2. Calculate standard error of beta_hat ####

    ## The derivative (bread) term is avg_SStrans

    bread <- avg_SStrans

    ## Calculate the cross-product (meat) term

    meat <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (i in 1:sample_size) {
        rowid <- person_first_index[i]:(person_first_index[i + 1] - 1)
        sum_t_ee_it <- numeric(dim_beta)
        for (it in rowid) {
            S_it <- matrix(Sdm[it, ], ncol = 1)
            r_it <- as.numeric(residual[it])
            omega_it <- omega[it]

            ee_it <- r_it * S_it - omega_it * S_it %*% t(S_it) %*% beta_hat
            sum_t_ee_it <- sum_t_ee_it + ee_it
        }
        meat <- meat + sum_t_ee_it %*% t(sum_t_ee_it)
    }
    meat <- meat / sample_size

    ## Calculate the variance of beta_hat using the sandwich estimator
    bread_inv <- solve(bread)
    beta_varcov <- bread_inv %*% meat %*% t(bread_inv) / sample_size
    beta_se <- sqrt(diag(beta_varcov))

    Snames <- c("Intercept", moderator_var)
    beta_hat <- as.vector(beta_hat)
    names(beta_hat) <- names(beta_se) <- Snames
    colnames(beta_varcov) <- rownames(beta_varcov) <- Snames

    #### 3. calculate confidence interval

    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05 / 2, df = sample_size - dim_beta)
    conf_int_tquantile <- cbind(
        beta_hat - c * beta_se,
        beta_hat + c * beta_se
    )
    colnames(conf_int) <- colnames(conf_int_tquantile) <- c("2.5 %", "97.5 %")

    return(list(
        beta_hat = beta_hat,
        beta_se = beta_se,
        conf_int = conf_int,
        conf_int_tquantile = conf_int_tquantile,
        beta_varcov = beta_varcov
    ))
}

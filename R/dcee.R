#' Distal Causal Excursion Effect (DCEE) Estimation
#'
#' @description
#' Fits distal causal excursion effects in micro-randomized trials using a
#' **two-stage** estimator: (i) learn nuisance outcome regressions
#' \eqn{\mu_a(H_t)} with a specified learner (parametric/ML), optionally with
#' cross-fitting; (ii) solve estimating equations for the distal excursion
#' effect parameters (\eqn{\beta}).
#'
#' This wrapper standardizes inputs and delegates computation to
#' [dcee_helper_2stage_estimation()].
#'
#' @param data A data.frame in long format.
#' @param id Character scalar: column name for subject identifier.
#' @param outcome Character scalar: column name for proximal/distal outcome.
#' @param treatment Character scalar: column name for binary treatment \{0,1\}.
#' @param rand_prob Character scalar: column name for randomization probability
#'   giving \eqn{P(A_t=1\mid H_t)} (must lie in (0,1)).
#' @param moderator_formula RHS-only formula of moderators of the excursion effect
#'   (e.g., `~ 1`, `~ Z`, or `~ Z1 + Z2`).
#' @param control_formula RHS-only formula of covariates for learning nuisance
#'   outcome regressions. When `control_reg_method = "gam"`, `s(x)` terms are
#'   allowed (e.g., `~ x1 + s(x2)`). For SuperLearner methods, variables are
#'   extracted from this formula to build the design matrix `X`.
#' @param availability Optional character scalar: column name for availability
#'   indicator (0/1). If `NULL`, availability is taken as 1 for all rows.
#' @param control_reg_method One of `"gam"`, `"lm"`, `"rf"`, `"ranger"`,
#'   `"sl"`, `"sl.user-specified-library"`, `"set_to_zero"`. See Details.
#' @param cross_fit Logical; if `TRUE`, perform K-fold cross-fitting by subject id.
#' @param cf_fold Integer; number of folds if `cross_fit = TRUE` (default 10).
#' @param weighting_function Either a single numeric constant applied to all
#'   rows, or a character column name in `data` giving decision-point weights
#'   \eqn{\omega_t}.
#' @param verbose Logical; print minimal preprocessing messages (default `TRUE`).
#' @param ... Additional arguments passed through to the chosen learner
#'   (e.g., `num.trees`, `mtry` for random forests; `sl.library` when
#'   `control_reg_method = "sl.user-specified-library"`).
#'
#' @return An object of class `"dcee_fit"` with components:
#' \describe{
#'   \item{call}{The matched call to \code{dcee()}.}
#'
#'   \item{fit}{A list returned by the two–stage helper with elements:
#'     \describe{
#'       \item{\code{beta_hat}}{Named numeric vector of distal causal excursion
#'         effect estimates \eqn{\beta}. Names are \code{"Intercept"} and the
#'         moderator names (if any) from \code{moderator_formula}.}
#'       \item{\code{beta_se}}{Named numeric vector of standard errors for
#'         \code{beta_hat} (same order/names).}
#'       \item{\code{beta_varcov}}{Variance–covariance matrix of \code{beta_hat}
#'         (square matrix; row/column names match \code{names(beta_hat)}).}
#'       \item{\code{conf_int}}{Matrix of large-sample (normal) Wald
#'         95\% confidence intervals for \code{beta_hat};
#'         columns are \code{"2.5 \%"} and \code{"97.5 \%"}.}
#'       \item{\code{conf_int_tquantile}}{Matrix of small-sample
#'         (t-quantile) 95\% confidence intervals for \code{beta_hat};
#'         columns are \code{"2.5 \%"} and \code{"97.5 \%"}; degrees of freedom
#'         are provided in \code{$df} of the \code{"dcee_fit"} object.}
#'       \item{\code{regfit_a0}}{Stage-1 nuisance regression fit for
#'         \eqn{\mu_0(H_t)} (outcome model among \code{A=0}), or \code{NULL}
#'         when \code{control_reg_method = "set_to_zero"}. \strong{Note:} when
#'         \code{cross_fit = TRUE}, this is the learner object from the
#'         \emph{last} fold and is provided for inspection only (do not use for
#'         out-of-fold prediction).}
#'       \item{\code{regfit_a1}}{Stage-1 nuisance regression fit for
#'         \eqn{\mu_1(H_t)} (outcome model among \code{A=1}); same caveats as
#'         \code{regfit_a0} regarding \code{cross_fit}.}
#'     }
#'   }
#'
#'   \item{df}{Small-sample degrees of freedom used for t-based intervals:
#'     number of unique subjects minus \code{length(fit$beta_hat)}.}
#' }
#'
#' @details
#' **Learners.**
#' - `gam` uses \pkg{mgcv} and supports `s(.)` terms in `control_formula`.
#' - `lm` uses base \code{stats::lm}.
#' - `rf` uses \pkg{randomForest}; `ranger` uses \pkg{ranger}.
#' - `sl` / `sl.user-specified-library` use \pkg{SuperLearner}. For the former,
#'   `sl.library = c("SL.mean", "SL.glm", "SL.earth")` are used. For the latter,
#'   please provide `sl.library = c("SL.mean", ...)` via `...`.
#'
#' **Notes.**
#' - Treatment must be coded 0/1; `rand_prob` must lie strictly in (0,1).
#' - `control_formula = ~ 1` is only valid with `control_reg_method = "set_to_zero"`.
#'
#' @examples
#' data(data_distal_continuous, package = "MRTAnalysis")
#'
#' ## Fast example: marginal effect with linear nuisance (CRAN-friendly)
#' fit_lm <- dcee(
#'     data = data_distal_continuous,
#'     id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
#'     moderator_formula = ~1, # marginal (no moderators)
#'     control_formula = ~X, # simple linear nuisance
#'     availability = "avail",
#'     control_reg_method = "lm",
#'     cross_fit = FALSE
#' )
#' summary(fit_lm)
#' summary(fit_lm, show_control_fit = TRUE) # show Stage-1 fit info
#'
#' ## Moderated effect with GAM nuisance (allows smooth terms); may be slower
#' \donttest{
#' fit_gam <- dcee(
#'     data = data_distal_continuous,
#'     id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
#'     moderator_formula = ~Z, # test moderation by Z
#'     control_formula = ~ s(X) + Z, # smooth in nuisance via mgcv::gam
#'     availability = "avail",
#'     control_reg_method = "gam",
#'     cross_fit = TRUE, cf_fold = 5
#' )
#' summary(fit_gam, lincomb = c(0, 1)) # linear combo: the Z coefficient
#' summary(fit_gam, show_control_fit = TRUE) # show Stage-1 fit info
#' }
#'
#' ## Optional: SuperLearner (runs only if installed)
#' \donttest{
#' library(SuperLearner)
#' fit_sl <- dcee(
#'     data = data_distal_continuous,
#'     id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
#'     moderator_formula = ~1,
#'     control_formula = ~ X + Z,
#'     availability = "avail",
#'     control_reg_method = "sl",
#'     cross_fit = FALSE
#' )
#' summary(fit_sl)
#' }
#' @export

dcee <- function(
    data,
    id,
    outcome,
    treatment,
    rand_prob,
    moderator_formula,
    control_formula,
    availability = NULL,
    control_reg_method = c("gam", "lm", "rf", "ranger", "sl", "sl.user-specified-library", "set_to_zero"),
    cross_fit = FALSE,
    cf_fold = 10,
    weighting_function = NULL,
    verbose = TRUE,
    ...) {
    # ---- small helpers ---------------------------------------------------------
    as_name <- function(x, arg) {
        if (!is.character(x) || length(x) != 1) stop(sprintf("`%s` must be a single character column name.", arg))
        x
    }
    rhs_vars <- function(fml, arg) {
        if (!inherits(fml, "formula")) stop(sprintf("`%s` must be a RHS-only formula like `~ x1 + x2`.", arg))
        tl <- attr(stats::terms(fml), "term.labels")
        if (length(tl) == 0L) character(0) else tl
    }

    # ---- normalize & validate inputs ------------------------------------------
    stopifnot(is.data.frame(data))
    id <- as_name(id, "id")
    outcome <- as_name(outcome, "outcome")
    treatment <- as_name(treatment, "treatment")
    rand_prob <- as_name(rand_prob, "rand_prob")
    if (!is.null(availability)) availability <- as_name(availability, "availability")

    control_reg_method <- match.arg(control_reg_method)

    # required columns present
    need <- c(id, outcome, treatment, rand_prob)
    if (!is.null(availability)) need <- c(need, availability)
    miss <- setdiff(need, names(data))
    if (length(miss)) stop("Missing columns in `data`: ", paste(miss, collapse = ", "))

    # ---- basic checks (treatment), then set up local data & availability -------
    A <- data[[treatment]]
    if (!all(A %in% c(0, 1))) stop("`treatment` must be coded 0/1.")

    # work on a local copy to avoid side-effects
    data_local <- data

    # availability: use user-supplied column if given, otherwise add an all-ones temp
    avail_var <- availability
    if (is.null(avail_var)) {
        tmp <- "..__avail__.."
        while (tmp %in% names(data_local)) tmp <- paste0(tmp, "_x")
        data_local[[tmp]] <- 1L
        avail_var <- tmp
    } else {
        # if user provided availability, ensure it is 0/1
        avail_vals <- data_local[[avail_var]]
        if (!all(avail_vals %in% c(0, 1))) {
            stop("`availability` must be coded 0/1 when provided.")
        }
    }

    # ---- rand_prob validation conditional on availability == 1 -----------------
    pA1 <- data_local[[rand_prob]]
    avail_vec <- data_local[[avail_var]]
    bad <- (avail_vec == 1) & (!is.finite(pA1) | pA1 <= 0 | pA1 >= 1)
    if (any(bad, na.rm = TRUE)) {
        stop("`rand_prob` must lie strictly in (0,1) whenever availability = 1.")
    }

    # parse formulas (moderators may be empty for marginal effect)
    moderator_vars <- rhs_vars(moderator_formula, "moderator_formula")
    control_vars <- rhs_vars(control_formula, "control_formula")

    # learner-specific guardrails
    if (identical(control_reg_method, "set_to_zero") && length(control_vars) > 0L) {
        if (verbose) message("[dcee] control_reg_method='set_to_zero' - ignoring control_formula terms.")
    }
    cf_str <- paste(deparse(control_formula), collapse = " ")
    if (!identical(control_reg_method, "gam") && grepl("s\\s*\\(", cf_str)) {
        stop("`s( )` terms are only supported for control_reg_method = 'gam'.")
    }
    if (length(control_vars) == 0L &&
        !(control_reg_method %in% c("set_to_zero", "lm", "gam"))) {
        stop("`control_formula = ~ 1` is only allowed with control_reg_method in {'set_to_zero','lm','gam'}.")
    }

    # handle weighting_function: numeric constant or column name (local-only)
    weighting_function_var <- NULL
    if (!is.null(weighting_function)) {
        if (is.character(weighting_function) && length(weighting_function) == 1) {
            if (!weighting_function %in% names(data_local)) stop("`weighting_function` column not found in `data`.")
            weighting_function_var <- weighting_function
        } else if (is.numeric(weighting_function) && length(weighting_function) == 1) {
            tmpw <- "..__omega__.."
            while (tmpw %in% names(data_local)) tmpw <- paste0(tmpw, "_x")
            data_local[[tmpw]] <- as.numeric(weighting_function)
            weighting_function_var <- tmpw
        } else {
            stop("`weighting_function` must be a single numeric or a single column name.")
        }
    }

    if (isTRUE(verbose)) {
        msg_m <- if (length(moderator_vars)) paste(moderator_vars, collapse = ", ") else "<marginal>"
        msg_c <- if (length(control_vars)) paste(control_vars, collapse = ", ") else "<none>"
        message(sprintf(
            "[dcee] Learner: %s | Moderators: %s | Controls: %s",
            control_reg_method, msg_m, msg_c
        ))
        if (isTRUE(cross_fit)) message(sprintf("[dcee] Cross-fitting enabled with %d folds", as.integer(cf_fold)))
    }

    # ---- delegate to helper (no modifications to its API) ----------------------
    fit <- dcee_helper_2stage_estimation(
        dta = data_local,
        id_var = id,
        moderator_formula = moderator_formula,
        trt_var = treatment,
        outcome_var = outcome,
        prob_A_var = rand_prob,
        avail_var = avail_var,
        control_reg_method = control_reg_method,
        control_formula = control_formula,
        cross_fit = cross_fit,
        cf_fold = as.integer(cf_fold),
        weighting_function_var = weighting_function_var,
        ...
    )

    # ---- assemble S3 object matching emee-style shape --------------------------
    out <- list(fit = fit)
    out$call <- match.call()
    n_id <- length(unique(data[[id]])) # use original data for df tally
    p <- length(fit$beta_hat) # number of beta parameters
    out$df <- n_id - p
    out <- out[c("call", "fit", "df")]
    class(out) <- "dcee_fit"
    out
}

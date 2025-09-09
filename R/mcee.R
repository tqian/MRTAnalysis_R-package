#' Mediated Causal Excursion Effects for MRTs (streamlined)
#'
#' Estimates the Natural Direct Excursion Effect (NDEE; \eqn{\alpha}) and
#' Natural Indirect Excursion Effect (NIEE; \eqn{\beta}) for distal outcomes
#' in micro-randomized trials (MRTs). Assumes the randomization probability
#' is known (via \code{rand_prob}) and fits all nuisance functions using the
#' same learner specified by \code{control_reg_method}.
#'
#' @param data A data.frame in long format (one row per id-by-decision point).
#' @param id Character. Column name for subject identifier.
#' @param dp Character. Column name for decision point index (must increase strictly within subject).
#' @param outcome Character. Column name for distal outcome (constant within subject).
#' @param treatment Character. Column name for treatment (coded 0/1).
#' @param mediator Character. Column name for mediator.
#' @param availability Optional character. Column name for availability (0/1). If \code{NULL}, all rows are treated as available.
#' @param rand_prob Either a column name in \code{data} or a scalar giving the known randomization probability \eqn{P(A_t=1 \mid H_t)}.  (Technically,
#'   this is \eqn{P(A_t=I_t\mid H_t)}, but the user is allowed to input \eqn{P(A_t=1\mid H_t)}
#'   and the function will automatically correct it by setting \code{p1 = 1} when \eqn{I_t = 0}.)
#' @param time_varying_effect_form RHS-only formula for the basis \eqn{f(t)} (e.g., \code{~ 1}, \code{~ dp}, \code{~ poly(dp,2)}).
#' @param control_formula_with_mediator RHS-only formula for control variables used in nuisance models that may include the mediator (the wrapper will drop the mediator internally for nuisances that must exclude it).
#' @param control_reg_method Learner for nuisance fits: one of \code{"glm"}, \code{"gam"}, \code{"rf"}, \code{"ranger"}, \code{"sl"}.
#' @param weight_per_row Optional numeric vector of row weights (nonnegative, length \code{nrow(data)}). If \code{NULL}, uniform within-id weights are used.
#' @param specific_dp_only Optional numeric vector of decision points to target; internally converted to \code{weight_per_row}.
#' @param verbose Logical; print progress messages.
#' @param SL.library Optional character vector of SuperLearner libraries (used when \code{control_reg_method = "sl"}).
#'
#' @details
#' Requirements: rows grouped by subject, strictly increasing \code{dp} within subject,
#' no missing (\code{NA}/\code{NaN}/\code{Inf}) in relevant variables. If \code{availability}
#' is supplied, the wrapper enforces at \eqn{I=0}: \eqn{p_1=q_1=1} in the nuisances.
#'
#' @return An object of class \code{"mcee_fit"} with elements:
#' \itemize{
#'   \item \code{mcee_fit}: list with \code{alpha_hat}, \code{beta_hat}, \code{alpha_se}, \code{beta_se},
#'         \code{varcov}, \code{alpha_varcov}, \code{beta_varcov}.
#'   \item \code{nuisance_models}: fitted Stage-1 models for \code{p,q,eta,mu,nu}.
#'   \item \code{nuisance_fitted}: per-row fitted values for the nuisance functions.
#'   \item \code{meta}: list with basis dimension, number of ids, per-id lengths, weights used.
#'   \item \code{call}: the matched call.
#' }
#'
#' @seealso \code{\link{summary.mcee_fit}}, \code{\link{mcee_general}}, \code{\link{mcee_userfit_nuisance}}
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' T <- 4
#' id <- rep(1:n, each = T)
#' dp <- rep(1:T, times = n)
#' A <- rbinom(n * T, 1, 0.5)
#' M <- rbinom(n * T, 1, plogis(-0.2 + 0.3 * A + 0.1 * dp))
#' Y <- ave(0.5 * A + 0.6 * M + 0.1 * dp + rnorm(n * T), id)
#' dat <- data.frame(id, dp, A, M, Y)
#'
#' fit <- mcee(dat, "id", "dp", "Y", "A", "M",
#'     time_varying_effect_form = ~1,
#'     control_formula_with_mediator = ~ dp + M,
#'     control_reg_method = "glm",
#'     rand_prob = 0.5, verbose = TRUE
#' )
#' summary(fit)
#'
#' @export
mcee <- function(
    data,
    id,
    dp,
    outcome,
    treatment,
    mediator,
    availability = NULL,
    rand_prob,
    time_varying_effect_form,
    control_formula_with_mediator,
    control_reg_method = c("glm", "gam", "rf", "ranger", "sl"),
    weight_per_row = NULL,
    specific_dp_only = NULL,
    verbose = TRUE,
    SL.library = NULL) {
    cl <- match.call()

    # ---- Basic checks ----------------------------------------------------------
    .mcee_assert_df(data)
    .mcee_require_cols(data, c(id, dp, outcome, treatment, mediator))
    if (!is.null(availability)) .mcee_require_cols(data, availability)

    # Check for missingness
    vars_core <- c(id, dp, outcome, treatment, mediator, if (!is.null(availability)) availability)
    vars_ctrl <- .mcee_vars_in_rhs(control_formula_with_mediator)
    vars_mod <- .mcee_vars_in_rhs(time_varying_effect_form)

    .mcee_check_no_missing_vars(
        data,
        vars = unique(c(vars_core, vars_ctrl, vars_mod)),
        where = "mcee inputs & control_formula_with_mediator & time_varying_effect_form"
    )

    # Check legal values of treatment, availability, outcome const within id
    .mcee_check_binary_col(data, treatment, allow_all1 = TRUE, label = "treatment")
    .mcee_check_binary_col(data, availability, allow_all1 = TRUE, label = "availability")
    .mcee_message_if_no_availability_provided(availability, verbose)
    .mcee_check_outcome_constant_within_id(data, id, outcome)

    # check for each id, the rows are grouped together, and that dp strictly increasing within id
    .mcee_check_id_rows_grouped(data, id)
    .mcee_check_dp_strictly_increasing(data, id, dp)

    # Known randomization prob p_t(1|H_t): column or scalar; validate (0,1) where avail==1
    p_vec <- .mcee_resolve_rand_prob(data, rand_prob, availability)

    # Basis (f-matrix) & moderator formula sanity
    .mcee_check_time_varying_effect_form(time_varying_effect_form, dp)
    f_nrows <- .mcee_build_f_matrix(time_varying_effect_form, data)

    # Row-weights omega(i,t)
    omega_nrows <- .mcee_build_weights(
        data, id, dp,
        weight_per_row = weight_per_row,
        specific_dp_only = specific_dp_only,
        verbose = verbose
    )

    # Control formula checks
    .mcee_check_control_formula(control_formula_with_mediator, treatment, outcome, dp,
        label = "control_formula_with_mediator"
    )
    control_reg_method <- match.arg(control_reg_method)

    # Drop mediator for eta/nu
    control_formula_mediator_removed <-
        .mcee_drop_var_from_rhs(control_formula_with_mediator, mediator)

    # ---- Build nuisance configs via mcee_config_maker --------------------------
    config_p <- mcee_config_maker(target = "p", known = p_vec)
    config_q <- mcee_config_maker(
        target = "q", method = control_reg_method,
        formula = control_formula_with_mediator,
        SL.library = SL.library
    )
    config_eta <- mcee_config_maker(
        target = "eta", method = control_reg_method,
        formula = control_formula_mediator_removed,
        SL.library = SL.library
    )
    config_mu <- mcee_config_maker(
        target = "mu", method = control_reg_method,
        formula = control_formula_with_mediator,
        SL.library = SL.library
    )
    config_nu <- mcee_config_maker(
        target = "nu", method = control_reg_method,
        formula = control_formula_mediator_removed,
        SL.library = SL.library
    )

    if (isTRUE(verbose)) {
        message(
            "[mcee] Learners: p=known; q/eta/mu/nu='", control_reg_method,
            "' | basis dim p=", ncol(f_nrows)
        )
    }

    # ---- Delegate to mcee_general ---------------------------------------------
    fit <- mcee_helper_2stage_estimation(
        data          = data,
        id_var        = id,
        dp_var        = dp,
        outcome_var   = outcome,
        treatment_var = treatment,
        mediator_var  = mediator,
        avail_var     = availability,
        config_p      = config_p,
        config_q      = config_q,
        config_eta    = config_eta,
        config_mu     = config_mu,
        config_nu     = config_nu,
        omega_nrows   = omega_nrows,
        f_nrows       = f_nrows
    )

    # ---- Assemble return --------------------------------------------------------
    T_per_id <- as.integer(table(data[[id]]))
    out <- list(
        call = cl,
        mcee_fit = fit$mcee_fit,
        nuisance_models = fit$nuisance_models,
        nuisance_fitted = fit$nuisance_fitted,
        meta = list(
            p_dim               = ncol(f_nrows),
            n_ids               = length(unique(data[[id]])),
            T_per_id            = T_per_id,
            weight_per_row_used = omega_nrows
        )
    )
    class(out) <- "mcee_fit"
    out
}

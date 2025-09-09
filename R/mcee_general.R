#' Mediated Causal Excursion Effects (configurable nuisance models)
#'
#' Like \code{\link{mcee}}, but each nuisance function is configured explicitly
#' via \code{config_*} objects (formula/method/family or known).
#'
#' @inheritParams mcee
#' @param config_p,config_q,config_eta,config_mu,config_nu Lists created by
#'   \code{mcee_config_maker()} or helpers (e.g., \code{mcee_config_glm()},
#'   \code{mcee_config_known()}, \code{mcee_config_ranger()}, etc.).
#'
#' @details
#' Use this wrapper for observational studies (estimate \code{p}) or when you want
#' different learners per nuisance. The same data requirements as \code{\link{mcee}} apply.
#'
#' @return An \code{"mcee_fit"} object; see \code{\link{mcee}}.
#' @seealso \code{\link{mcee}}, \code{\link{mcee_userfit_nuisance}}, \code{\link{mcee_config_maker}}
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
#' cfg <- list(
#'     p   = mcee_config_known("p", 0.5),
#'     q   = mcee_config_glm("q", ~ dp + M),
#'     eta = mcee_config_glm("eta", ~dp),
#'     mu  = mcee_config_glm("mu", ~ dp + M),
#'     nu  = mcee_config_glm("nu", ~dp)
#' )
#' fit_gen <- mcee_general(dat, "id","dp","Y","A","M",
#'     time_varying_effect_form = ~ dp,
#'     config_p=cfg$p, config_q=cfg$q, config_eta=cfg$eta, config_mu=cfg$mu, config_nu=cfg$nu)
#' @export
mcee_general <- function(
    data,
    id,
    dp,
    outcome,
    treatment,
    mediator,
    availability = NULL,
    time_varying_effect_form,
    config_p,
    config_q,
    config_eta,
    config_mu,
    config_nu,
    weight_per_row = NULL,
    verbose = TRUE) {
    cl <- match.call()

    # ---- Basic checks ----------------------------------------------------------
    .mcee_assert_df(data)
    .mcee_require_cols(data, c(id, dp, outcome, treatment, mediator))
    if (!is.null(availability)) .mcee_require_cols(data, availability)

    # Check for missingness
    vars_core <- c(id, dp, outcome, treatment, mediator, if (!is.null(availability)) availability)
    vars_cfg <- unique(c(
        .mcee_vars_in_config(config_p),
        .mcee_vars_in_config(config_q),
        .mcee_vars_in_config(config_eta),
        .mcee_vars_in_config(config_mu),
        .mcee_vars_in_config(config_nu)
    ))
    vars_mod <- .mcee_vars_in_rhs(time_varying_effect_form)

    .mcee_check_no_missing_vars(
        data,
        vars = unique(c(vars_core, vars_cfg, vars_mod)),
        where = "mcee_general inputs & config formulas & time_varying_effect_form"
    )

    # Check legal values of treatment, availability, outcome const within id
    .mcee_check_binary_col(data, treatment, allow_all1 = TRUE, label = "treatment")
    .mcee_check_binary_col(data, availability, allow_all1 = TRUE, label = "availability")
    .mcee_message_if_no_availability_provided(availability, verbose)
    .mcee_check_outcome_constant_within_id(data, id, outcome)

    .mcee_check_all_formulas(
        config_p = config_p,
        config_q = config_q,
        config_eta = config_eta,
        config_mu = config_mu,
        config_nu = config_nu,
        mediator = mediator
    )

    # check for each id, the rows are grouped together, and that dp strictly increasing within id
    .mcee_check_id_rows_grouped(data, id)
    .mcee_check_dp_strictly_increasing(data, id, dp)

    # Basis (f-matrix) & moderator formula sanity
    .mcee_check_time_varying_effect_form(time_varying_effect_form, dp)
    f_nrows <- .mcee_build_f_matrix(time_varying_effect_form, data)

    # Row-weights omega(i,t) (no specific_dp_only here -- pass exactly what user gave)
    omega_nrows <- .mcee_build_weights(
        data, id, dp,
        weight_per_row = weight_per_row,
        specific_dp_only = NULL,
        verbose = verbose
    )

    # ---- Two-stage estimation ---------------------------------------------------
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

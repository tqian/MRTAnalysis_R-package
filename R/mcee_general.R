# ------------------------------------------------------------------------------
# mcee_general(): Flexible wrapper (user supplies nuisance configs)
# ------------------------------------------------------------------------------

#' Mediated Causal Excursion Effects (MCEE) — General (Config-Driven)
#'
#' @description
#' General wrapper where users **provide five nuisance configs** (`config_p`,
#' `config_q`, `config_eta`, `config_mu`, `config_nu`). The function handles
#' data checks, sorting, basis construction, weights, and calls the two-stage
#' helper to estimate \eqn{\alpha}, \eqn{\beta} and their joint variance.
#'
#' @param data,id,dp,outcome,treatment,mediator,availability See \code{mcee()}.
#' @param time_varying_effect_form RHS-only formula to build the Stage-2 basis \(f(t)\)
#'   (functions of the decision point).
#' @param config_p,config_q,config_eta,config_mu,config_nu Nuisance model
#'   configuration lists (see `mcee_config_maker()` and the package docs).
#' @param weight_per_row Optional numeric vector of per-row weights \(\omega(i,t)\).
#'   If `NULL`, all ones are used.
#' @param verbose Logical; print brief diagnostics.
#'
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
    verbose = TRUE
) {
  cl <- match.call()

  # ---- Basic checks ----------------------------------------------------------
  .mcee_assert_df(data)
  .mcee_require_cols(data, c(id, dp, outcome, treatment, mediator))
  if (!is.null(availability)) .mcee_require_cols(data, availability)

  # Check for missingness
  vars_core <- c(id, dp, outcome, treatment, mediator, if (!is.null(availability)) availability)
  vars_cfg  <- unique(c(
    .mcee_vars_in_config(config_p),
    .mcee_vars_in_config(config_q),
    .mcee_vars_in_config(config_eta),
    .mcee_vars_in_config(config_mu),
    .mcee_vars_in_config(config_nu)
  ))
  vars_mod  <- .mcee_vars_in_rhs(time_varying_effect_form)

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

  .mcee_check_all_formulas(config_p = config_p,
                           config_q = config_q,
                           config_eta = config_eta,
                           config_mu = config_mu,
                           config_nu = config_nu,
                           mediator = mediator)

  # check for each id, the rows are grouped together, and that dp strictly increasing within id
  .mcee_check_id_rows_grouped
  .mcee_check_dp_strictly_increasing(data, id, dp)

  # Basis (f-matrix) & moderator formula sanity
  .mcee_check_time_varying_effect_form(time_varying_effect_form, dp)
  f_nrows <- .mcee_build_f_matrix(time_varying_effect_form, data)

  # Row-weights ω(i,t) (no specific_dp_only here—pass exactly what user gave)
  omega_nrows <- .mcee_build_weights(
    data, id, dp,
    weight_per_row = weight_per_row,
    specific_dp_only = NULL,
    verbose = verbose
  )

  # ---- Two-stage estimation ---------------------------------------------------
  fit2 <- mcee_helper_2stage_estimation(
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
    call            = cl,
    mcee_fit        = fit2$mcee_fit,
    nuisance_models = fit2$nuisance_models,
    nuisance_fitted = fit2$nuisance_fitted,
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

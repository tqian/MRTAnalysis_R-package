# ------------------------------------------------------------------------------
# mcee(): Streamlined MRT wrapper (known randomization probability)
# ------------------------------------------------------------------------------

#' Mediated Causal Excursion Effects (MCEE) for Distal Outcomes in MRTs (Streamlined)
#'
#' @description
#' Estimates **natural direct** (NDEE; \eqn{\alpha}) and **natural indirect**
#' (NIEE; \eqn{\beta}) excursion effects via a two-stage procedure.
#' Stage 1 fits nuisance models; Stage 2 solves estimating equations with
#' robust sandwich variance and small-sample \eqn{t}-inference.
#'
#' This streamlined wrapper assumes **known MRT randomization probabilities**
#' (`rand_prob`) and builds nuisance configs via `mcee_config_maker()`.
#'
#' @inheritParams mcee_general
#' @param rand_prob Character column name or numeric scalar for known
#'   randomization probability \eqn{p_t(1\mid H_t)}.
#' @param control_formula_with_mediator RHS-only formula for nuisance learning
#'   that **may include the mediator** (used for `q` and `mu`; automatically
#'   dropped for `eta` and `nu` internally).
#' @param control_reg_method One of `"glm"`, `"gam"`, `"rf"`, `"ranger"`, `"sl"`.
#' @param specific_dp_only Optional numeric vector of decision-point values to
#'   target (rows with `dp` not in this set receive weight 0).
#' @param SL.library Optional character vector of SuperLearner base learners
#'   (used when `control_reg_method="sl"`).
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
    SL.library = NULL
) {
  cl <- match.call()

  # ---- Basic checks ----------------------------------------------------------
  .mcee_assert_df(data)
  .mcee_require_cols(data, c(id, dp, outcome, treatment, mediator))
  if (!is.null(availability)) .mcee_require_cols(data, availability)

  # Check for missingness
  vars_core <- c(id, dp, outcome, treatment, mediator, if (!is.null(availability)) availability)
  vars_ctrl <- .mcee_vars_in_rhs(control_formula_with_mediator)
  vars_mod  <- .mcee_vars_in_rhs(time_varying_effect_form)

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
  .mcee_check_id_rows_grouped
  .mcee_check_dp_strictly_increasing(data, id, dp)

  # Known randomization prob p_t(1|H_t): column or scalar; validate (0,1) where avail==1
  p_vec <- .mcee_resolve_rand_prob(data, rand_prob, availability)

  # Basis (f-matrix) & moderator formula sanity
  .mcee_check_time_varying_effect_form(time_varying_effect_form, dp)
  f_nrows <- .mcee_build_f_matrix(time_varying_effect_form, data)

  # Row-weights ω(i,t)
  omega_nrows <- .mcee_build_weights(
    data, id, dp,
    weight_per_row = weight_per_row,
    specific_dp_only = specific_dp_only,
    verbose = verbose
  )

  # Control formula checks
  .mcee_check_control_formula(control_formula_with_mediator, treatment, outcome, dp,
                              label = "control_formula_with_mediator")
  control_reg_method <- match.arg(control_reg_method)

  # Drop mediator for eta/nu
  control_formula_mediator_removed <-
    .mcee_drop_var_from_rhs(control_formula_with_mediator, mediator)

  # ---- Build nuisance configs via mcee_config_maker --------------------------
  config_p   <- mcee_config_maker(target = "p",   known   = p_vec)
  config_q   <- mcee_config_maker(target = "q",   method  = control_reg_method,
                                  formula = control_formula_with_mediator,
                                  SL.library = SL.library)
  config_eta <- mcee_config_maker(target = "eta", method  = control_reg_method,
                                  formula = control_formula_mediator_removed,
                                  SL.library = SL.library)
  config_mu  <- mcee_config_maker(target = "mu",  method  = control_reg_method,
                                  formula = control_formula_with_mediator,
                                  SL.library = SL.library)
  config_nu  <- mcee_config_maker(target = "nu",  method  = control_reg_method,
                                  formula = control_formula_mediator_removed,
                                  SL.library = SL.library)

  if (isTRUE(verbose)) {
    message("[mcee] Learners: p=known; q/eta/mu/nu='", control_reg_method,
            "' | basis dim p=", ncol(f_nrows))
  }

  # ---- Delegate to mcee_general ---------------------------------------------
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

# ------------------------------------------------------------------------------
# MCEE two-stage helper: Stage 1 (nuisance fits) + Stage 2 (alpha/beta & variance)
# ------------------------------------------------------------------------------

#' Two-stage helper for mediated causal excursion effects (MCEE)
#'
#' Fits all nuisance components (Stage 1) and then computes the MCEE parameters
#' (Stage 2) and their sandwich variance. This is a low-level driver used by the
#' high-level wrapper; it assumes `omega_nrows` and `f_nrows` are already aligned
#' to the rows of `data`.
#'
#' @param data          A long-format `data.frame` (one row per subject-by-decision point).
#' @param id_var        Character scalar. Name of the subject ID column.
#' @param dp_var        Character scalar. Name of the decision point column
#'   (values need not be consecutive; they may vary in count across subjects).
#' @param outcome_var   Character scalar. Name of the distal outcome column.
#' @param treatment_var Character scalar. Name of the binary treatment column (coded 0/1).
#' @param mediator_var  Character scalar. Name of the mediator column.
#' @param avail_var     Character scalar or `NULL`. Name of the availability column
#'   (1 = available, 0 = unavailable). If `NULL`, availability is treated as all 1.
#'
#' @param config_p      Configuration for \(p_t(a\mid H_t)\) (propensity). A **list**
#'   using one of the following schemas:
#'   \itemize{
#'     \item \emph{Known constant(s)} (skips fitting):
#'       \code{list(known = c(...))} or arm-specific \code{known_a1}/\code{known_a0}.
#'     \item \emph{Model fit}: \code{list(formula = ~ rhs, method = m, ...)} where
#'       \code{method} is one of \code{"glm"}, \code{"gam"}, \code{"rf"},
#'       \code{"ranger"}, \code{"sl"}, \code{"sl.user-specified-library"}.
#'       Optional fields:
#'       \itemize{
#'         \item \code{family}: a GLM/GAM family. If omitted, **auto-detected** as
#'               \code{binomial()} for \eqn{p} and \eqn{q}, otherwise \code{gaussian()}.
#'         \item \code{clipping}: numeric length-2 \code{c(lo, hi)} to clip probabilities
#'               into \code{[lo, hi]} (useful for stability).
#'         \item For \code{method = "sl"}: \code{SL.library} (character vector of learners);
#'         if omitted, a simple default library is used: c("SL.mean", "SL.glm", "SL.gam").
#'       }
#'   }
#' @param config_q      Configuration for \(q_t(a\mid H_t, M_t)\). Same schema as \code{config_p}.
#' @param config_eta    Configuration for \(\eta_t(a, H_t)\) (outcome given \(A,H\)).
#'   Same schema as \code{config_p}; default family auto-detected to \code{gaussian()} if omitted.
#' @param config_mu     Configuration for \(\mu_t(a, H_t, M_t)\) (outcome given \(A,H,M\)).
#'   Same schema as \code{config_p}; default family auto-detected to \code{gaussian()} if omitted.
#' @param config_nu     Configuration for \(\nu_t(a, H_t)\) (cross-world ICE based on \(\mu\)).
#'   Same schema as \code{config_p}; default family auto-detected to \code{gaussian()} if omitted.
#'
#' @param omega_nrows   Numeric vector of length \code{nrow(data)}. Per-row weights
#'   \(\omega(i,t) \ge 0\). Rows are aligned with \code{data} (no reordering).
#' @param f_nrows       Numeric matrix with \code{nrow(data)} rows and \code{p} columns.
#'   Row \code{r} contains \(f(t_r)^\top\) (the basis evaluated at the decision point
#'   of row \code{r}). The same basis is used for both \(\alpha\) (NDEE) and \(\beta\) (NIEE).
#'
#' @details
#' \strong{Availability handling:}
#' When \code{avail_var} exists and equals 0, Stage 1 sets the working probabilities
#' to 1 for that row (e.g., \eqn{\hat p_t(1\mid H_t)=1}, \eqn{\hat p_t(0\mid H_t)=1}, similarly
#' for \eqn{\hat q_t}). This prevents division-by-zero in the estimating equations.
#'
#' \strong{Auto-family rules:}
#' If \code{family} is omitted in a GLM/GAM config, it defaults to \code{binomial()}
#' for \code{config_p} and \code{config_q}, and to \code{gaussian()} for
#' \code{config_eta}, \code{config_mu}, and \code{config_nu}.
#'
#' \strong{Learners:}
#' \itemize{
#'   \item \code{"glm"}: uses \code{stats::glm()}.
#'   \item \code{"gam"}: uses \code{mgcv::gam()} (supports \code{s()} smooths).
#'   \item \code{"rf"}: uses \code{randomForest::randomForest()}.
#'   \item \code{"ranger"}: uses \code{ranger::ranger()}.
#'   \item \code{"sl"}: uses \code{SuperLearner::SuperLearner()}.
#'         If \code{SL.library} is not given, a simple default library is used:
#'         c("SL.mean", "SL.glm", "SL.gam").
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{fit}}{A list with entries
#'     \code{alpha_hat}, \code{alpha_se}, \code{beta_hat}, \code{beta_se},
#'     and \code{varcov} (the \(2p\times 2p\) sandwich variance for \((\alpha^\top,\beta^\top)^\top\)).}
#'   \item{\code{nuisance_models}}{A list of fitted Stage-1 objects:
#'     \code{p}, \code{q}, \code{eta1}, \code{eta0}, \code{mu1}, \code{mu0}, \code{nu1}, \code{nu0}.
#'     (For \code{known}/\code{known_a0}/\code{known_a1}, a small descriptor list is returned.)}
#' }
#'
#' @seealso \code{\link{mcee_general}} for a high-level wrapper that constructs
#' \code{omega_nrows} and \code{f_nrows} from user-friendly arguments.
#'
#' @examples
#' \dontrun{
#' # Minimal sketch (assuming `df` has columns id, t, A, M, Y, I):
#' fit <- mcee_helper_2stage_estimation(
#'   data = df,
#'   id_var = "id", dp_var = "t", outcome_var = "Y",
#'   treatment_var = "A", mediator_var = "M", avail_var = "I",
#'   config_p   = list(formula = ~ t + M, method = "glm"),     # binomial auto
#'   config_q   = list(formula = ~ t + M + A, method = "glm"), # binomial auto
#'   config_eta = list(formula = ~ t, method = "gam"),         # gaussian auto
#'   config_mu  = list(formula = ~ t + s(M), method = "gam"),  # gaussian auto
#'   config_nu  = list(formula = ~ t, method = "glm"),         # gaussian auto
#'   omega_nrows = rep(1, nrow(df)),
#'   f_nrows     = cbind(1)  # marginal (p = 1)
#' )
#' fit$fit$alpha_hat
#' fit$fit$beta_hat
#' }
mcee_helper_2stage_estimation <- function(
    data,
    id_var,
    dp_var,
    outcome_var,
    treatment_var,
    mediator_var,
    avail_var = NULL,
    config_p,
    config_q,
    config_eta,
    config_mu,
    config_nu,
    omega_nrows,
    f_nrows
) {

  # ---- Stage 1: Fit nuisance parameters (and keep the models) ----------------
  stage1 <- mcee_helper_stage1_fit_nuisance(
    data            = data,
    id_var          = id_var,
    dp_var          = dp_var,
    outcome_var     = outcome_var,
    treatment_var   = treatment_var,
    mediator_var    = mediator_var,
    avail_var       = avail_var,
    config_p        = config_p,
    config_q        = config_q,
    config_eta      = config_eta,
    config_mu       = config_mu,
    config_nu       = config_nu
  )

  nuisance_models <- stage1$nuisance_models
  nuisance_fitted <- stage1$nuisance_fitted

  # ---- Stage 2: Estimate (alpha, beta) and variance --------------------------
  fit <- mcee_helper_stage2_estimate_mcee(
    data         = data,
    id_var       = id_var,
    dp_var       = dp_var,
    outcome_var  = outcome_var,
    treatment_var= treatment_var,
    avail_var    = avail_var,
    p1           = nuisance_fitted$p1,
    p0           = nuisance_fitted$p0,
    q1           = nuisance_fitted$q1,
    q0           = nuisance_fitted$q0,
    eta1         = nuisance_fitted$eta1,
    eta0         = nuisance_fitted$eta0,
    mu1          = nuisance_fitted$mu1,
    mu0          = nuisance_fitted$mu0,
    nu1          = nuisance_fitted$nu1,
    nu0          = nuisance_fitted$nu0,
    omega_nrows  = omega_nrows,   # length nrow(data)
    f_nrows      = f_nrows       # nrow(data) × p
  )

  # Return Stage-2 results + Stage-1 fitted model objects
  list(
    mcee_fit = fit,
    nuisance_models = nuisance_models,
    nuisance_fitted = nuisance_fitted
  )
}


# ------------------------------------------------------------------------------
# Stage 1: fit nuisance models and return both predictions & fitted objects
# ------------------------------------------------------------------------------

#' Fit nuisance models for MCEE and return predictions + model objects
#'
#' @return list with:
#'   - nuisance_fitted: per-row vectors (p1, p0, q1, q0, eta1, eta0, mu1, mu0, nu1, nu0)
#'   - nuisance_models: fitted objects (or "known" descriptors) for each component
mcee_helper_stage1_fit_nuisance <- function(
    data,
    id_var,
    dp_var,
    outcome_var,
    treatment_var,
    mediator_var,
    avail_var,
    config_p,
    config_q,
    config_eta,
    config_mu,
    config_nu
) {
  nrows <- nrow(data)

  Y <- data[[outcome_var]]
  A <- data[[treatment_var]]
  # M <- data[[mediator_var]]  # extracted automatically via config formulas
  I <- if (is.null(avail_var)) rep(1, nrows) else data[[avail_var]]

  # ---- p_t(1 | H_t) ----------------------------------------------------------
  res_p <- .mcee_fit_nuisance(
    config               = config_p,
    data_for_fitting     = data[I == 1, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = treatment_var,
    param_name           = "p_t(1|H_t)",
    data_for_fitting_name = "all available decision points"
  )
  p1 <- res_p$pred
  p1[I == 0] <- 1
  p0 <- 1 - p1
  p0[I == 0] <- 1
  model_p <- res_p$model

  # ---- q_t(1 | H_t, M_t) -----------------------------------------------------
  res_q <- .mcee_fit_nuisance(
    config               = config_q,
    data_for_fitting     = data[I == 1, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = treatment_var,
    param_name           = "q_t(1|H_t, M_t)",
    data_for_fitting_name = "all available decision points"
  )
  q1 <- res_q$pred
  q1[I == 0] <- 1
  q0 <- 1 - q1
  q0[I == 0] <- 1
  model_q <- res_q$model

  # ---- eta_t(a, H_t): outcome regression without M_t -------------------------
  # A=1 fit (plus I==0 rows to stabilize when always-available rows exist)
  res_eta1 <- .mcee_fit_nuisance(
    config               = config_eta,
    data_for_fitting     = data[A == 1 | I == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = outcome_var,
    param_name           = "eta_t(1, H_t)",
    data_for_fitting_name = "all treated or unavailable decision points"
  )
  # A=0 fit
  res_eta0 <- .mcee_fit_nuisance(
    config               = config_eta,
    data_for_fitting     = data[A == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = outcome_var,
    param_name           = "eta_t(0, H_t)",
    data_for_fitting_name = "all untreated decision points"
  )
  eta1 <- res_eta1$pred
  eta0 <- res_eta0$pred
  model_eta1 <- res_eta1$model
  model_eta0 <- res_eta0$model

  # ---- mu_t(a, H_t, M_t): outcome regression with M_t ------------------------
  res_mu1 <- .mcee_fit_nuisance(
    config               = config_mu,
    data_for_fitting     = data[A == 1 | I == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = outcome_var,
    param_name           = "mu_t(1, H_t, M_t)",
    data_for_fitting_name = "all treated or unavailable decision points"
  )
  res_mu0 <- .mcee_fit_nuisance(
    config               = config_mu,
    data_for_fitting     = data[A == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = outcome_var,
    param_name           = "mu_t(0, H_t, M_t)",
    data_for_fitting_name = "all untreated decision points"
  )
  mu1 <- res_mu1$pred
  mu0 <- res_mu0$pred
  model_mu1 <- res_mu1$model
  model_mu0 <- res_mu0$model

  # ---- nu_t(a, H_t): cross-world ICE based on mu_t ---------------------------
  # Fit nu_t(1, H_t) with lhs = mu1_predicted among A==0 rows
  res_nu1 <- .mcee_fit_nuisance(
    config               = config_nu,
    data_for_fitting     = cbind(data, mu1_predicted = mu1)[A == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = "mu1_predicted",
    param_name           = "nu_t(1, H_t)",
    data_for_fitting_name = "all untreated decision points"
  )
  # Fit nu_t(0, H_t) with lhs = mu0_predicted among A==1 or I==0 rows
  res_nu0 <- .mcee_fit_nuisance(
    config               = config_nu,
    data_for_fitting     = cbind(data, mu0_predicted = mu0)[A == 1 | I == 0, , drop = FALSE],
    data_for_predicting  = data,
    lhs_var              = "mu0_predicted",
    param_name           = "nu_t(0, H_t)",
    data_for_fitting_name = "all treated or unavailable decision points"
  )
  nu1 <- res_nu1$pred
  nu0 <- res_nu0$pred
  model_nu1 <- res_nu1$model
  model_nu0 <- res_nu0$model

  list(
    nuisance_fitted = list(
      p1 = p1, p0 = p0, q1 = q1, q0 = q0,
      eta1 = eta1, eta0 = eta0, mu1 = mu1, mu0 = mu0, nu1 = nu1, nu0 = nu0
    ),
    nuisance_models = list(
      p = model_p, q = model_q,
      eta1 = model_eta1, eta0 = model_eta0,
      mu1  = model_mu1,  mu0 = model_mu0,
      nu1  = model_nu1,  nu0 = model_nu0
    )
  )
}


# ------------------------------------------------------------------------------
# Stage 2: build φ’s per row, then call numeric core
# ------------------------------------------------------------------------------

#' Stage-2 estimator for MCEE given Stage-1 predictions
#'
#' @param omega_nrows Numeric vector length nrow(data): ω(i,t).
#' @param f_nrows     Numeric matrix nrow(data) × p: per-row basis rows f(t)^T.
#'
#' @return list(alpha_hat, alpha_se, beta_hat, beta_se, varcov, alpha_varcov, beta_varcov)
mcee_helper_stage2_estimate_mcee <- function(
    data,
    id_var,
    dp_var,
    outcome_var,
    treatment_var,
    avail_var = NULL,
    p1,
    p0,
    q1,
    q0,
    eta1,
    eta0,
    mu1,
    mu0,
    nu1,
    nu0,
    omega_nrows,
    f_nrows
) {
  stopifnot(is.data.frame(data))
  nrows <- nrow(data)

  # Validate ω and f
  if (!is.numeric(omega_nrows) || length(omega_nrows) != nrows)
    stop("`omega_nrows` must be a numeric vector of length nrow(data).")
  if (any(!is.finite(omega_nrows)) || any(omega_nrows < 0))
    stop("`omega_nrows` must be nonnegative and finite.")
  if (!is.matrix(f_nrows) || nrow(f_nrows) != nrows)
    stop("`f_nrows` must be a numeric matrix with nrow(f_nrows) == nrow(data).")

  # Subject index 1..n (stable by first appearance)
  id_chr    <- as.character(data[[id_var]])
  id_levels <- unique(id_chr)
  i_index   <- match(id_chr, id_levels)
  n         <- length(id_levels)

  # Build indicators (A==1, A==0); allow I==0 to force I1=1 per your construction
  Y <- data[[outcome_var]]
  A <- data[[treatment_var]]
  I <- if (is.null(avail_var)) rep(1, nrows) else data[[avail_var]]

  indic_A_equal_dt1 <- as.numeric(A == 1 | I == 0)
  indic_A_equal_dt0 <- as.numeric(A == 0)

  # φ’s per row (vectorized); Stage-1 should ensure safe denominators at I==0
  phi11_vec <- indic_A_equal_dt1 / p1 * Y - (indic_A_equal_dt1 - p1) / p1 * eta1
  phi00_vec <- indic_A_equal_dt0 / p0 * Y - (indic_A_equal_dt0 - p0) / p0 * eta0

  phi10_vec <- indic_A_equal_dt1 * q0 * (Y - mu1) / (p0 * q1) +
    indic_A_equal_dt0 * (mu1 - nu1) / p0 +
    nu1

  # Call the numerical core that implements the formulas exactly
  core <- .mcee_core_rows(
    n            = n,
    f_nrows      = f_nrows,
    omega_nrows  = omega_nrows,
    i_index      = as.integer(i_index),
    phi11_vec    = phi11_vec,
    phi10_vec    = phi10_vec,
    phi00_vec    = phi00_vec
  )

  # Give names to coefficients, se, and varcov
  coef_names <- colnames(f_nrows)
  names(core$alpha_hat) <- names(core$alpha_se) <- coef_names
  names(core$beta_hat) <- names(core$beta_se) <- coef_names
  rownames(core$alpha_varcov) <- colnames(core$alpha_varcov) <- coef_names
  rownames(core$beta_varcov) <- colnames(core$beta_varcov) <- coef_names
  rownames(core$varcov) <- colnames(core$varcov) <-
    c(paste0("alpha_", coef_names), paste0("beta_", coef_names))

  core
}


# ------------------------------------------------------------------------------
# Fit-one-nuisance helper: supports 'glm', 'lm', 'gam', 'rf', 'ranger', 'sl',
# and 'sl.user-specified-library', or known constants.
# ------------------------------------------------------------------------------

#' Fit a single nuisance component and predict on all rows
#'
#' @param config A list describing how to fit the nuisance:
#'   - One of:
#'       * known = <numeric scalar>              # same value for all rows
#'       * known_a1 / known_a0 = <numeric>       # arm-specific constant
#'       * formula = ~ x1 + s(x2) + ...          # RHS-only formula
#'         and method ∈ {"glm","lm","gam","rf","ranger","sl","sl.user-specified-library"}
#'   - Optional common fields:
#'       * family = stats::binomial() or stats::gaussian()   # for glm/gam/SL; inferred if absent
#'   - Optional method-specific fields:
#'       * glm.args / gam.args / rf.args / ranger.args : named lists passed to the fitter
#'       * SL.library (optional library for "sl")
#'       * SL.control (optional control list for "sl")
#'
#' @param data_for_fitting data.frame used to fit the model
#' @param data_for_predicting data.frame used to generate predictions
#' @param lhs_var character, name of response to model
#' @param param_name character, label for messages/errors
#' @param data_for_fitting_name character, for printing the data set used in fitting
#' the nuisance parameter.
#'
#' @return list(pred = numeric, model = fitted object or "known" descriptor)
.mcee_fit_nuisance <- function(config, data_for_fitting, data_for_predicting,
                               lhs_var, param_name, data_for_fitting_name) {
  out <- list(pred = NULL, model = NULL)

  # ----- helpers ---------------------------------------------------------------
  normalize_formula <- function(rhs_like, lhs_name) {
    rhs_txt <- if (inherits(rhs_like, "formula")) {
      deparse(rhs_like[[2L]])
    } else {
      deparse(stats::as.formula(rhs_like)[[2L]])
    }
    stats::as.formula(paste0(lhs_name, " ~ ", rhs_txt))
  }
  detect_binomial <- function(cfg, y) {
    if (!is.null(cfg$family)) {
      fam <- cfg$family
      if (inherits(fam, "family")) return(identical(fam$family, "binomial"))
      # if user passed e.g. "binomial"
      if (is.character(fam)) return(tolower(fam) == "binomial")
    }
    # fallback: infer from y being 0/1 only
    uy <- unique(na.omit(as.numeric(y)))
    all(uy %in% c(0, 1))
  }
  # Build model matrix X for SuperLearner from a RHS-only formula (drop intercept)
  build_X_from_rhs <- function(rhs_like, data) {
    rhs_txt <- if (inherits(rhs_like, "formula")) {
      deparse(rhs_like[[2L]])
    } else {
      deparse(stats::as.formula(rhs_like)[[2L]])
    }
    mm <- stats::model.matrix(stats::as.formula(paste0("~", rhs_txt)), data)
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
    mm
  }

  # ----- known constants -------------------------------------------------------
  if (!is.null(config$known)) {
    out$pred  <- rep_len(config$known, nrow(data_for_predicting))
    out$model <- list(type = "known", value = config$known, for_param = param_name)
    return(out)
  }
  if (!is.null(config$known_a1) && grepl("1", param_name, fixed = TRUE)) {
    out$pred  <- rep_len(config$known_a1, nrow(data_for_predicting))
    out$model <- list(type = "known", value = config$known_a1, for_param = param_name)
    return(out)
  }
  if (!is.null(config$known_a0) && grepl("0", param_name, fixed = TRUE)) {
    out$pred  <- rep_len(config$known_a0, nrow(data_for_predicting))
    out$model <- list(type = "known", value = config$known_a0, for_param = param_name)
    return(out)
  }

  # ----- model-based paths -----------------------------------------------------
  if (is.null(config$formula)) {
    stop("No formula provided for nuisance parameter: ", param_name)
  }
  if (is.null(config$method)) {
    stop("No method provided for nuisance parameter: ", param_name,
         ". Supported: 'glm','lm','gam','rf','ranger','sl','sl.user-specified-library'.")
  }

  method <- match.arg(config$method,
                      c("glm","lm","gam","rf","ranger","sl","sl.user-specified-library"))
  form <- normalize_formula(config$formula, lhs_var)

  # family & task type
  y_fit <- data_for_fitting[[lhs_var]]
  is_binom <- detect_binomial(config, y_fit)

  # helper: choose a symbolic family when using defaults
  .family_sym <- function(config_family, is_binom) {
    if (!is.null(config_family)) return(config_family)           # pass through user-supplied
    if (is_binom) quote(stats::binomial) else quote(stats::gaussian)
  }

  # ----- GLM -------------------------------------------------------------------
  if (identical(method, "glm")) {
    family_sym <- .family_sym(config$family, is_binom)
    args <- c(
      list(
        formula = form,
        data    = quote(data_for_fitting),   # keep name in recorded call
        family  = family_sym                 # keep name (binomial/gaussian) in call
      ),
      config$glm.args %||% list()
    )
    fit <- do.call("glm", args)
    fam <- fit$family$family      # "gaussian", "binomial", ...
    lnk <- fit$family$link        # "identity", "logit", ...
    defaults <- list(
      binomial = "logit", gaussian = "identity", Gamma = "inverse",
      poisson  = "log", quasibinomial = "logit", quasipoisson = "log",
      "inverse.gaussian" = "1/mu^2", quasi = "identity"
    )
    fit$call$family <- if (!is.null(defaults[[fam]]) && identical(lnk, defaults[[fam]])) {
      as.name(fam)                        # e.g., gaussian
    } else {
      as.call(list(as.name(fam), link = lnk))   # e.g., binomial(link = "probit")
    }
    fit$call$data   <- as.name(data_for_fitting_name)
    pred <- stats::predict(fit, newdata = data_for_predicting, type = "response")
    out$pred  <- as.numeric(pred)
    out$model <- fit
    return(out)
  }

  # ----- LM --------------------------------------------------------------------
  if (identical(method, "lm")) {
    if (is_binom) {
      warning("`lm` used with an apparent binary outcome in ", param_name,
              "; consider `glm(family=binomial)`.")
    }
    args <- c(
      list(
        formula = form,
        data    = quote(data_for_fitting)
      ),
      config$lm.args %||% list()
    )
    fit  <- do.call("lm", args)
    fit$call$data   <- as.name(data_for_fitting_name)
    pred <- stats::predict(fit, newdata = data_for_predicting)
    out$pred  <- as.numeric(pred)
    out$model <- fit
    return(out)
  }

  # ----- GAM (mgcv) ------------------------------------------------------------
  if (identical(method, "gam")) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("mgcv is required for method='gam'.")
    }
    family_sym <- .family_sym(config$family, is_binom)
    args <- c(
      list(
        formula = form,
        data    = quote(data_for_fitting),
        family  = family_sym
      ),
      config$gam.args %||% list()
    )
    # call mgcv::gam (namespace-loaded but not attached)
    fit <- do.call(getExportedValue("mgcv", "gam"), args)
    # make the printed call pretty
    fit$call[[1]] <- as.name("gam")
    fit$call$data   <- as.name(data_for_fitting_name)
    pred <- stats::predict(fit, newdata = data_for_predicting, type = "response")
    out$pred  <- as.numeric(pred)
    out$model <- fit
    return(out)
  }

  # ----- RandomForest ----------------------------------------------------------
  if (identical(method, "rf")) {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("randomForest is required for method='rf'.")
    }
    if (is_binom) {
      data_for_fitting[[lhs_var]] <- factor(data_for_fitting[[lhs_var]], levels = c(0, 1))
    }
    args <- c(
      list(
        formula = form,
        data    = quote(data_for_fitting)
      ),
      config$rf.args %||% list()
    )
    fit <- do.call(getExportedValue("randomForest", "randomForest"), args)
    fit$call[[1]] <- as.name("randomForest")  # pretty call
    fit$call$data   <- as.name(data_for_fitting_name)

    if (is_binom) {
      pr <- stats::predict(fit, newdata = data_for_predicting, type = "prob")
      col1 <- if ("1" %in% colnames(pr)) "1" else colnames(pr)[2L]
      pred <- pr[, col1]
    } else {
      pred <- stats::predict(fit, newdata = data_for_predicting)
    }
    out$pred  <- as.numeric(pred)
    out$model <- fit
    return(out)
  }

  # ----- ranger ----------------------------------------------------------------
  if (identical(method, "ranger")) {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("ranger is required for method='ranger'.")
    }
    dep_name <- lhs_var
    if (is_binom) {
      data_for_fitting[[dep_name]] <- factor(data_for_fitting[[dep_name]], levels = c(0, 1))
    }
    args <- c(
      list(
        formula = form,
        data    = quote(data_for_fitting)
      ),
      config$ranger.args %||% list()
    )
    if (is_binom && is.null(args$probability)) args$probability <- TRUE

    fit <- do.call(getExportedValue("ranger", "ranger"), args)
    fit$call[[1]] <- as.name("ranger")  # pretty call
    fit$call$data   <- as.name(data_for_fitting_name)

    pr <- stats::predict(fit, data = data_for_predicting)
    if (is_binom) {
      mat <- pr$predictions
      pred <- if (is.vector(mat)) as.numeric(mat) else {
        col1 <- if ("1" %in% colnames(mat)) "1" else colnames(mat)[2L]
        mat[, col1]
      }
    } else {
      pred <- as.numeric(pr$predictions)
    }
    out$pred  <- as.numeric(pred)
    out$model <- fit
    return(out)
  }

  # ----- SuperLearner -----------------------------------------------------------
  if (identical(method, "sl")) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("SuperLearner is not installed but method='", method, "' was requested.")
    }

    SL.lib <- config$SL.library
    if (identical(method, "sl") && is.null(SL.lib)) SL.lib <- c("SL.mean", "SL.glm")

    X_fit <- build_X_from_rhs(config$formula, data_for_fitting)
    X_new <- build_X_from_rhs(config$formula, data_for_predicting)
    Y_fit <- data_for_fitting[[lhs_var]]
    family_val <- if (!is.null(config$family)) config$family else
      if (is_binom) stats::binomial() else stats::gaussian()

    SL.ctrl <- config$SL.control %||% list()

    fit <- SuperLearner::SuperLearner(
      Y = Y_fit,
      X = as.data.frame(X_fit),
      family = family_val,
      SL.library = SL.lib,
      verbose = FALSE,
      control = SL.ctrl
    )
    fit$call$family <- as.name(family_val$family)

    pr <- stats::predict(fit, newdata = as.data.frame(X_new))$pred
    out$pred  <- as.numeric(pr)
    out$model <- fit
    return(out)
  }

  stop("Unsupported method for nuisance parameter: ", param_name)
}

# little infix helper for optional lists
`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# Numeric core (row-wise): implements the math literally and supports T_i unequal
# ------------------------------------------------------------------------------

#' Numeric core for MCEE (row-wise; supports unequal T_i)
#'
#' Implements:
#'   α̂  = S^{-1} (1/n) Σ_{i,t} ω(i,t){φ^{10}-φ^{00}} f(t)
#'   β̂  = S^{-1} (1/n) Σ_{i,t} ω(i,t){φ^{11}-φ^{10}} f(t)
#' where S = (1/n) Σ_{i,t} ω(i,t) f(t)f(t)^T.
#'
#' Var((α̂,β̂)) = Bread^{-1} Meat Bread^{-T} / n,
#' with Bread = blockdiag(S, S), and
#' Meat = (1/n) Σ_i u_i u_i^T, u_i = Σ_t ω(i,t) [ {φ^{10}-φ^{00}-f^Tα̂}f ; {φ^{11}-φ^{10}-f^Tβ̂}f ].
#'
#' @return list(alpha_hat, alpha_se, beta_hat, beta_se, varcov, alpha_varcov, beta_varcov)
.mcee_core_rows <- function(
    n,
    f_nrows,
    omega_nrows,
    i_index,
    phi11_vec, phi10_vec, phi00_vec
) {
  # ---- Checks ----------------------------------------------------------------
  if (!is.matrix(f_nrows)) stop("`f_nrows` must be a numeric matrix (nrows × p).")
  nrows <- nrow(f_nrows)
  p     <- ncol(f_nrows)

  if (!is.numeric(omega_nrows) || length(omega_nrows) != nrows)
    stop("`omega_nrows` must be a numeric vector of length nrow(f_nrows).")
  if (any(!is.finite(omega_nrows)) || any(omega_nrows < 0))
    stop("`omega_nrows` must be nonnegative and finite.")

  if (length(i_index) != nrows ||
      length(phi11_vec) != nrows ||
      length(phi10_vec) != nrows ||
      length(phi00_vec) != nrows)
    stop("Lengths of i_index and φ-vectors must equal nrow(f_nrows).")

  if (!is.integer(i_index)) i_index <- as.integer(i_index)
  if (any(i_index < 1L | i_index > n))
    stop("`i_index` contains values outside 1..n.")

  # ---- Precompute row-wise pieces -------------------------------------------
  d_alpha_row <- (phi10_vec - phi00_vec) * omega_nrows   # ω * {φ10 - φ00}
  d_beta_row  <- (phi11_vec - phi10_vec) * omega_nrows   # ω * {φ11 - φ10}

  # S = (1/n) Σ ω f f^T
  S <- (t(f_nrows) %*% (f_nrows * omega_nrows)) / n
  qrS <- qr(S)
  if (qrS$rank < p) stop("Bread block S is singular; try adjusting f/omega.")
  S_inv <- solve(S)

  # α̂, β̂
  alpha_rhs <- as.vector(colSums(f_nrows * d_alpha_row) / n)
  beta_rhs  <- as.vector(colSums(f_nrows * d_beta_row)  / n)
  alpha_hat <- as.vector(S_inv %*% alpha_rhs)
  beta_hat  <- as.vector(S_inv %*% beta_rhs)

  # Meat = (1/n) Σ_i u_i u_i^T, with
  # u_i blocks built from residuals times f(t), summed per subject.
  f_alpha <- as.numeric(f_nrows %*% alpha_hat)
  f_beta  <- as.numeric(f_nrows %*% beta_hat)
  res_alpha_row <- (phi10_vec - phi00_vec) - f_alpha
  res_beta_row  <- (phi11_vec - phi10_vec) - f_beta

  wres_alpha <- omega_nrows * res_alpha_row
  wres_beta  <- omega_nrows * res_beta_row

  tmp1 <- matrix(0, nrow = n, ncol = p)  # Σ_t ω * residual_alpha * f(t)
  tmp2 <- matrix(0, nrow = n, ncol = p)  # Σ_t ω * residual_beta  * f(t)

  for (r in seq_len(nrows)) {
    i <- i_index[r]
    if (wres_alpha[r] != 0) tmp1[i, ] <- tmp1[i, ] + wres_alpha[r] * f_nrows[r, ]
    if (wres_beta[r]  != 0) tmp2[i, ] <- tmp2[i, ] + wres_beta[r]  * f_nrows[r, ]
  }

  U    <- cbind(tmp1, tmp2)          # n × (2p)
  Meat <- crossprod(U) / n           # (2p × 2p)

  # Bread^{-1} = blockdiag(S^{-1}, S^{-1})
  Binv <- matrix(0, nrow = 2 * p, ncol = 2 * p)
  Binv[1:p, 1:p] <- S_inv
  Binv[(p + 1):(2 * p), (p + 1):(2 * p)] <- S_inv

  varcov <- Binv %*% Meat %*% t(Binv) / n  # (2p × 2p)

  alpha_varcov <- varcov[1:p, 1:p, drop = FALSE]
  beta_varcov  <- varcov[(p + 1):(2 * p), (p + 1):(2 * p), drop = FALSE]
  alpha_se     <- sqrt(if (p == 1) alpha_varcov else diag(alpha_varcov))
  beta_se      <- sqrt(if (p == 1) beta_varcov  else diag(beta_varcov))

  list(
    alpha_hat    = as.vector(alpha_hat),
    alpha_se     = as.vector(alpha_se),
    beta_hat     = as.vector(beta_hat),
    beta_se      = as.vector(beta_se),
    varcov       = varcov,
    alpha_varcov = alpha_varcov,
    beta_varcov  = beta_varcov
  )
}

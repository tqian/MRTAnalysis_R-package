# ------------------------------------------------------------
# Helpers to build config lists for MCEE nuisance components
# ------------------------------------------------------------

# Internal: pick default family if omitted
.mcee_default_family <- function(target, method) {
  # target in {"p","q","eta","mu","nu"}
  if (method %in% c("glm", "gam")) {
    if (target %in% c("p", "q")) stats::binomial() else stats::gaussian()
  } else {
    NULL
  }
}

# Internal: validate method & clipping
.mcee_validate_method <- function(method) {
  ok <- c("glm", "gam", "lm", "rf", "ranger", "sl", "sl.user-specified-library")
  if (is.null(method)) return(invisible(TRUE))
  if (!method %in% ok) {
    stop("Unsupported method '", method,
         "'. Choose one of: ", paste(ok, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.mcee_validate_clipping <- function(clipping) {
  if (is.null(clipping)) return(invisible(TRUE))
  if (!is.numeric(clipping) || length(clipping) != 2L) {
    stop("`clipping` must be numeric length-2 (e.g., c(0.01, 0.99)).", call. = FALSE)
  }
  lo <- clipping[1]; hi <- clipping[2]
  if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi >= 1 || lo >= hi) {
    stop("`clipping` must satisfy 0 < lo < hi < 1. Got c(", lo, ", ", hi, ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' General factory for MCEE nuisance configs
#'
#' @param target One of `"p"`, `"q"`, `"eta"`, `"mu"`, `"nu"`.
#' @param method One of `"glm"`, `"gam"`, `"lm"`, `"rf"`, `"ranger"`, `"sl"`, `"sl.user-specified-library"`.
#'   Omit `method` when using `known` or `known_a*`.
#' @param formula RHS-only formula (e.g., `~ x + s(t)`). Ignored if `known*` is used.
#' @param family GLM/GAM family. If omitted, auto-detected: `binomial()` for `p`/`q`, else `gaussian()`.
#' @param known,known_a1,known_a0 Known value(s) to bypass fitting
#'   (vector or scalar). If set, `method/formula` are ignored.
#' @param clipping Optional numeric length-2 to bound probabilities for stability
#'   (mainly for `p`/`q`): e.g., `c(0.01, 0.99)`.
#' @param SL.library Character vector of SuperLearner wrappers for `method="sl"`.
#'If omitted for `"sl"`, default c("SL.mean", "SL.glm", "SL.gam") is used.
#' @param ... Reserved for future extensions (silently stored but not used).
#'
#' @return A list suitable for `config_*` in `mcee_general()` / `mcee()`.
mcee_config_maker <- function(target,
                             method = NULL,
                             formula = NULL,
                             family = NULL,
                             known = NULL,
                             known_a1 = NULL,
                             known_a0 = NULL,
                             clipping = NULL,
                             SL.library = NULL,
                             ...) {
  target <- match.arg(tolower(target), c("p","q","eta","mu","nu"))
  # Known overrides: ignore method/formula/family
  if (!is.null(known) || !is.null(known_a1) || !is.null(known_a0)) {
    return(list(
      nuisance_parameter = target,
      known = known,
      known_a1 = known_a1,
      known_a0 = known_a0,
      clipping = clipping  # harmless if present
    ))
  }

  # Validate method & clipping
  .mcee_validate_method(method)
  .mcee_validate_clipping(clipping)

  # Family default if needed (GLM/GAM only)
  if (is.null(family)) {
    family <- .mcee_default_family(target, method)
  }

  # Default SL library
  if (identical(method, "sl") && is.null(SL.library)) {
    SL.library <- c("SL.mean", "SL.glm", "SL.gam")
  }

  out <- list(
    nuisance_parameter = target,
    method = method,
    formula = formula,
    family  = family,
    clipping = clipping
  )
  if (!is.null(SL.library)) out$SL.library <- SL.library

  out
}

# ------------------------------
# Convenience thin wrappers
# ------------------------------

mcee_config_known <- function(target, value = NULL, a1 = NULL, a0 = NULL) {
  # if both scalar & arm-specific provided, arm-specific wins
  list(known = value, known_a1 = a1, known_a0 = a0)
}

mcee_config_glm <- function(target, formula, family = NULL, clipping = NULL) {
  mcee_config_maker(target = target, method = "glm",
                   formula = formula, family = family, clipping = clipping)
}

mcee_config_gam <- function(target, formula, family = NULL, clipping = NULL) {
  mcee_config_maker(target = target, method = "gam",
                   formula = formula, family = family, clipping = clipping)
}

mcee_config_lm <- function(target, formula) {
  mcee_config_maker(target = target, method = "lm", formula = formula)
}

mcee_config_rf <- function(target, formula) {
  mcee_config_maker(target = target, method = "rf", formula = formula)
}

mcee_config_ranger <- function(target, formula) {
  mcee_config_maker(target = target, method = "ranger", formula = formula)
}

mcee_config_sl <- function(target, formula, SL.library = NULL, clipping = NULL) {
  mcee_config_maker(target = target, method = "sl",
                   formula = formula, SL.library = SL.library, clipping = clipping)
}

mcee_config_sl_user <- function(target, formula, SL.library, clipping = NULL) {
  mcee_config_maker(target = target, method = "sl.user-specified-library",
                   formula = formula, SL.library = SL.library, clipping = clipping)
}

# ------------------------------
# Usage examples (sketch)
# ------------------------------
# config_p   <- mcee_config_glm("p",   ~ s(t) + X + M_prev, clipping = c(0.01, 0.99))
# config_q   <- mcee_config_gam("q",   ~ s(t) + X + s(M),   clipping = c(0.01, 0.99))
# config_eta <- mcee_config_lm ("eta", ~ t + X)
# config_mu  <- mcee_config_gam("mu",  ~ s(t) + X + s(M))
# config_nu  <- mcee_config_glm("nu",  ~ t + X)
# # Known propensity (e.g., constant 0.6):
# config_p_known <- mcee_config_known(value = 0.6)
# # SuperLearner with explicit library:
# config_mu_sl <- mcee_config_sl_user("mu", ~ t + X + M,
#                                     SL.library = c("SL.glm", "SL.gam", "SL.mean"))



# ------------------------------------------------------------------------------
# Helpers to lint nuisance model formulas
# ------------------------------------------------------------------------------

#' Check config formula for inclusion/exclusion of mediator
#'
#' @param config A nuisance config list (may or may not contain `formula`).
#' @param target Character scalar, one of "p", "q", "eta", "mu", "nu".
#' @param mediator Character scalar: mediator variable name.
#'
#' @return Invisibly TRUE. Warnings are produced if formulas look suspicious.
.mcee_check_formula_mediator <- function(config, target, mediator) {
  # If no formula, nothing to check
  if (is.null(config$formula)) return(invisible(TRUE))

  # Extract variables from formula RHS
  vars <- all.vars(config$formula)

  if (target %in% c("p", "eta", "nu")) {
    # Mediator should NOT be included
    if (mediator %in% vars) {
      warning(sprintf(
        "[mcee] For nuisance parameter '%s': mediator '%s' appears in formula, but it should be excluded.",
        target, mediator
      ), call. = FALSE)
    }
  } else if (target %in% c("q", "mu")) {
    # Mediator SHOULD be included
    if (!(mediator %in% vars)) {
      warning(sprintf(
        "[mcee] For nuisance parameter '%s': mediator '%s' not found in formula, but it is usually recommended to include it.",
        target, mediator
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Vectorized helper to check all configs at once
# ------------------------------------------------------------------------------

.mcee_check_all_formulas <- function(config_p, config_q, config_eta, config_mu, config_nu,
                                     mediator) {
  .mcee_check_formula_mediator(config_p,   "p",   mediator)
  .mcee_check_formula_mediator(config_q,   "q",   mediator)
  .mcee_check_formula_mediator(config_eta, "eta", mediator)
  .mcee_check_formula_mediator(config_mu,  "mu",  mediator)
  .mcee_check_formula_mediator(config_nu,  "nu",  mediator)
  invisible(TRUE)
}



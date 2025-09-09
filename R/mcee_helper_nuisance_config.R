#' Build a nuisance-configuration object for \code{mcee_general()}
#'
#' @description
#' Creates a configuration list describing **how to obtain a nuisance function**
#' used by \code{\link{mcee_general}}. You may either:
#' \enumerate{
#'   \item supply **known values** (bypasses learning), or
#'   \item specify a **learning method** (e.g., GLM/GAM/RF/Ranger/SL) with a formula.
#' }
#'
#' @param target Character; which nuisance to configure. One of
#'   \code{"p"}, \code{"q"}, \code{"eta"}, \code{"mu"}, \code{"nu"}.
#' @param method Optional character learner name when *not* using known values.
#'   Supported:\cr
#'   \itemize{
#'     \item \code{"glm"}, \code{"gam"}, \code{"lm"} (formula-based);
#'     \item \code{"rf"} (randomForest), \code{"ranger"};
#'     \item \code{"sl"} (SuperLearner).
#'   }
#'   Ignored if any of \code{known}, \code{known_a1}, \code{known_a0} is provided.
#' @param formula RHS-only formula describing predictors for the learner
#'   (used when \code{method} is formula-based; ignored for \code{known}).
#'   For \code{method = "gam"}, \code{s()} terms are allowed (via \pkg{mgcv}).
#' @param family Optional GLM/GAM family. If \code{NULL}, a default is chosen
#'   based on \code{target} and \code{method} (typically \code{binomial()}
#'   for \code{p}/\code{q}, \code{gaussian()} for \code{eta}/\code{mu}/\code{nu}).
#' @param known Optional numeric scalar/vector of **known values** for the nuisance.
#'   Commonly used for \code{target = "p"} in MRTs when the randomization
#'   probability is known.
#' @param known_a1,known_a0 Optional numeric scalar/vector providing known values
#'   for the \emph{treatment-specific} versions of a nuisance (e.g., \code{eta1}/\code{eta0},
#'   \code{mu1}/\code{mu0}, \code{nu1}/\code{nu0}). If supplied, the function
#'   returns a "known" config and \code{method}/\code{formula}/\code{family} are ignored.
#' @param clipping Optional numeric vector of length 2, \code{c(lower, upper)},
#'   to truncate predictions (e.g., probabilities into \eqn{[\epsilon, 1-\epsilon]}).
#'   Validated but not required.
#' @param SL.library Character vector of SuperLearner libraries (only used
#'   when \code{method = "sl"}). If \code{NULL}, defaults to
#'   \code{c("SL.mean", "SL.glm", "SL.gam")}.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @details
#' If any of \code{known}, \code{known_a1}, or \code{known_a0} is provided,
#' the returned configuration is of type “known” and **no learner will be fit**.
#' Otherwise, the configuration records the requested learner, formula, family,
#' optional clipping, and (for SL) the library.
#'
#' Internally, helper validators ensure \code{method} is supported and
#' \code{clipping} (if provided) is sane. Family defaults are chosen when
#' \code{family = NULL} for GLM/GAM methods.
#'
#' @return A named \code{list} describing the configuration. For known configs:
#' \preformatted{
#' list(
#'   nuisance_parameter = <target>,
#'   known   = <numeric or NULL>,
#'   known_a1 = <numeric or NULL>,
#'   known_a0 = <numeric or NULL>,
#'   clipping = <numeric length-2 or NULL>
#' )
#' }
#' For learner configs:
#' \preformatted{
#' list(
#'   nuisance_parameter = <target>,
#'   method  = <character>,
#'   formula = <formula or NULL>,
#'   family  = <family or NULL>,
#'   clipping = <numeric length-2 or NULL>,
#'   SL.library = <character vector; only when method == "sl">
#' )
#' }
#'
#' @seealso
#' \code{\link{mcee_general}},
#' helper constructors like
#' \code{mcee_config_known()},
#' \code{mcee_config_glm()},
#' \code{mcee_config_gam()},
#' \code{mcee_config_lm()},
#' \code{mcee_config_rf()},
#' \code{mcee_config_ranger()},
#' \code{mcee_config_sl()},
#' \code{mcee_config_sl_user()}.
#'
#' @examples
#' # Known p (MRT randomization), GLM for other nuisances
#' cfg_p <- mcee_config_maker("p", known = 0.5)
#' cfg_q <- mcee_config_maker("q", method = "glm", formula = ~ dp + M)
#' cfg_eta <- mcee_config_maker("eta", method = "glm", formula = ~dp)
#' cfg_mu <- mcee_config_maker("mu", method = "glm", formula = ~ dp + M)
#' cfg_nu <- mcee_config_maker("nu", method = "glm", formula = ~dp)
#'
#' # SuperLearner with default library (set explicitly if you prefer)
#' # cfg_q_sl <- mcee_config_maker("q", method = "sl", formula = ~ dp + M,
#' #                               SL.library = c("SL.mean","SL.glm","SL.ranger"))
#'
#' # Known treatment-specific outcome regressions (e.g., from external source)
#' # cfg_eta_known <- mcee_config_maker("eta", known_a1 = rep(1, 100),
#' #                                         known_a0 = rep(0, 100))
#'
#' @export

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
    target <- match.arg(tolower(target), c("p", "q", "eta", "mu", "nu"))
    # Known overrides: ignore method/formula/family
    if (!is.null(known) || !is.null(known_a1) || !is.null(known_a0)) {
        return(list(
            nuisance_parameter = target,
            known = known,
            known_a1 = known_a1,
            known_a0 = known_a0,
            clipping = clipping # harmless if present
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
        family = family,
        clipping = clipping
    )
    if (!is.null(SL.library)) out$SL.library <- SL.library

    out
}

# ------------------------------
# Convenience thin wrappers
# ------------------------------

#' Configure known constant values for MCEE nuisance parameters
#'
#' Creates a configuration for nuisance parameters with known constant values,
#' bypassing model fitting. Useful for known randomization probabilities in MRTs.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param value Numeric scalar. Single constant value for all observations.
#' @param a1,a0 Numeric scalars. Arm-specific constants for A=1 and A=0 conditions.
#'   If provided, these override \code{value}.
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # Known randomization probability
#' cfg_p <- mcee_config_known("p", 0.6)
#'
#' # Arm-specific known values
#' cfg_eta <- mcee_config_known("eta", a1 = 0.8, a0 = 0.2)
#' @export
mcee_config_known <- function(target, value = NULL, a1 = NULL, a0 = NULL) {
    # if both scalar & arm-specific provided, arm-specific wins
    list(known = value, known_a1 = a1, known_a0 = a0)
}

#' Configure GLM for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using generalized linear models
#' via \code{stats::glm()}.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + poly(time, 2)}).
#' @param family Optional GLM family. Defaults to \code{binomial()} for "p"/"q",
#'   \code{gaussian()} for "eta"/"mu"/"nu".
#' @param clipping Optional numeric vector \code{c(lo, hi)} to clip predictions
#'   into [lo, hi] for numerical stability.
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # Binary outcome model for propensity
#' cfg_q <- mcee_config_glm("q", ~ dp + M, family = binomial())
#'
#' # Gaussian outcome model
#' cfg_eta <- mcee_config_glm("eta", ~ dp + X1)
#' @export
mcee_config_glm <- function(target, formula, family = NULL, clipping = NULL) {
    mcee_config_maker(
        target = target, method = "glm",
        formula = formula, family = family, clipping = clipping
    )
}

#' Configure GAM for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using generalized additive models
#' via \code{mgcv::gam()}. Supports smooth terms like \code{s()}.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + s(time) + s(X2, k=5)}).
#' @param family Optional GLM family. Defaults to \code{binomial()} for "p"/"q",
#'   \code{gaussian()} for "eta"/"mu"/"nu".
#' @param clipping Optional numeric vector \code{c(lo, hi)} to clip predictions
#'   into [lo, hi] for numerical stability.
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # GAM with smooth time effect
#' cfg_eta <- mcee_config_gam("eta", ~ X1 + s(dp, k = 4))
#'
#' # GAM with multiple smooths
#' cfg_mu <- mcee_config_gam("mu", ~ s(dp) + s(M, X1, k = 10))
#' @export
mcee_config_gam <- function(target, formula, family = NULL, clipping = NULL) {
    mcee_config_maker(
        target = target, method = "gam",
        formula = formula, family = family, clipping = clipping
    )
}

#' Configure linear model for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using linear models
#' via \code{stats::lm()}. Only appropriate for continuous outcomes.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + poly(dp, 2)}).
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # Linear model for continuous outcome
#' cfg_eta <- mcee_config_lm("eta", ~ dp + X1 + X2)
#' @export
mcee_config_lm <- function(target, formula) {
    mcee_config_maker(target = target, method = "lm", formula = formula)
}

#' Configure Random Forest for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using random forests
#' via \code{randomForest::randomForest()}. Good for nonlinear patterns.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + dp}).
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # Random forest for complex propensity model
#' cfg_q <- mcee_config_rf("q", ~ dp + M + X1 + X2)
#' @export
mcee_config_rf <- function(target, formula) {
    mcee_config_maker(target = target, method = "rf", formula = formula)
}

#' Configure Ranger Random Forest for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using ranger random forests
#' via \code{ranger::ranger()}. Faster alternative to \code{randomForest}.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + dp}).
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # Ranger random forest for outcome model
#' cfg_eta <- mcee_config_ranger("eta", ~ dp + X1 + X2 + X3)
#' @export
mcee_config_ranger <- function(target, formula) {
    mcee_config_maker(target = target, method = "ranger", formula = formula)
}

#' Configure SuperLearner for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using SuperLearner
#' via \code{SuperLearner::SuperLearner()}. Automatically selects among
#' multiple learning algorithms.
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + dp}).
#' @param SL.library Optional character vector of learner names. If \code{NULL},
#'   uses default library: \code{c("SL.mean", "SL.glm", "SL.gam")}.
#' @param clipping Optional numeric vector \code{c(lo, hi)} to clip predictions
#'   into [lo, hi] for numerical stability.
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # SuperLearner with default library
#' cfg_q <- mcee_config_sl("q", ~ dp + M + X1)
#'
#' # SuperLearner with custom library
#' cfg_eta <- mcee_config_sl("eta", ~ dp + X1,
#'     SL.library = c("SL.glm", "SL.rf", "SL.ranger")
#' )
#' @export
mcee_config_sl <- function(target, formula, SL.library = NULL, clipping = NULL) {
    mcee_config_maker(
        target = target, method = "sl",
        formula = formula, SL.library = SL.library, clipping = clipping
    )
}

#' Configure SuperLearner with user-specified library for MCEE nuisance parameters
#'
#' Creates a configuration to fit nuisance parameters using SuperLearner
#' with a user-specified library (required parameter).
#'
#' @param target Character. Nuisance parameter name ("p", "q", "eta", "mu", "nu").
#' @param formula RHS-only formula (e.g., \code{~ X1 + X2 + dp}).
#' @param SL.library Character vector of learner names (required).
#' @param clipping Optional numeric vector \code{c(lo, hi)} to clip predictions
#'   into [lo, hi] for numerical stability.
#'
#' @return A configuration list for use with \code{\link{mcee_general}}.
#'
#' @examples
#' # SuperLearner with specific library
#' cfg_mu <- mcee_config_sl_user("mu", ~ dp + M + X1,
#'     SL.library = c("SL.glm", "SL.earth", "SL.nnet")
#' )
#' @export
mcee_config_sl_user <- function(target, formula, SL.library, clipping = NULL) {
    mcee_config_maker(
        target = target, method = "sl.user-specified-library",
        formula = formula, SL.library = SL.library, clipping = clipping
    )
}


#' Select default GLM family based on nuisance parameter type
#' @param target Nuisance parameter name ("p", "q", "eta", "mu", "nu")
#' @param method Learning method name
#' @return GLM family object or NULL for non-GLM methods
.mcee_default_family <- function(target, method) {
    # target in {"p","q","eta","mu","nu"}
    if (method %in% c("glm", "gam")) {
        if (target %in% c("p", "q")) stats::binomial() else stats::gaussian()
    } else {
        NULL
    }
}

#' Validate that learning method is supported
#' @param method Method name to validate
#' @return Invisibly TRUE if valid, otherwise stops with error
.mcee_validate_method <- function(method) {
    ok <- c("glm", "gam", "lm", "rf", "ranger", "sl", "sl.user-specified-library")
    if (is.null(method)) {
        return(invisible(TRUE))
    }
    if (!method %in% ok) {
        stop("Unsupported method '", method,
            "'. Choose one of: ", paste(ok, collapse = ", "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate clipping bounds for probability predictions
#' @param clipping Numeric vector of length 2 with lower and upper bounds
#' @return Invisibly TRUE if valid, otherwise stops with error
.mcee_validate_clipping <- function(clipping) {
    if (is.null(clipping)) {
        return(invisible(TRUE))
    }
    if (!is.numeric(clipping) || length(clipping) != 2L) {
        stop("`clipping` must be numeric length-2 (e.g., c(0.01, 0.99)).", call. = FALSE)
    }
    lo <- clipping[1]
    hi <- clipping[2]
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi >= 1 || lo >= hi) {
        stop("`clipping` must satisfy 0 < lo < hi < 1. Got c(", lo, ", ", hi, ").", call. = FALSE)
    }
    invisible(TRUE)
}


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
    if (is.null(config$formula)) {
        return(invisible(TRUE))
    }

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
    .mcee_check_formula_mediator(config_p, "p", mediator)
    .mcee_check_formula_mediator(config_q, "q", mediator)
    .mcee_check_formula_mediator(config_eta, "eta", mediator)
    .mcee_check_formula_mediator(config_mu, "mu", mediator)
    .mcee_check_formula_mediator(config_nu, "nu", mediator)
    invisible(TRUE)
}

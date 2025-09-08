#' Summary for MCEE fits
#'
#' @description
#' Produce inference tables for mediated causal excursion effects from an
#' [mcee()] or [mcee_general()] fit. Inference uses small-sample t-tests with
#' `df = number of subjects − 2*p`, where `p` is the dimension of the basis
#' for \eqn{f(t)} (i.e., length of \code{alpha_hat} / \code{beta_hat}).
#'
#' @param object An object of class \code{"mcee_fit"} returned by [mcee()] or
#'   [mcee_general()].
#' @param lincomb_alpha Optional numeric vector or matrix specifying linear
#'   combinations \(L_\alpha \alpha\). If a vector of length `p`, computes a
#'   single combination; if a matrix, must have `p` columns.
#' @param lincomb_beta Optional numeric vector or matrix specifying linear
#'   combinations \(L_\beta \beta\). Same rules as \code{lincomb_alpha}.
#' @param lincomb_joint Optional numeric vector or matrix specifying linear
#'   combinations \(L \theta\) where \(\theta = (\alpha^\top,\beta^\top)^\top\).
#'   If a vector of length `2p`, computes a single combination; if a matrix,
#'   must have `2p` columns.
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param show_nuisance Logical; if `TRUE`, prints compact information about
#'   Stage-1 nuisance fits saved in \code{object$nuisance_models}.
#' @param ... Currently ignored.
#'
#' @return An object of class \code{"summary.mcee_fit"} with components:
#' \itemize{
#'   \item \code{call} — original call
#'   \item \code{df} — degrees of freedom for t-tests
#'   \item \code{conf_level} — requested confidence level
#'   \item \code{alpha} — data.frame with estimates/SE/t/CI for \eqn{\alpha}
#'   \item \code{beta} — data.frame with estimates/SE/t/CI for \eqn{\beta}
#'   \item \code{lincomb_alpha}, \code{lincomb_beta}, \code{lincomb_joint} — optional
#'         data.frames with linear-combination results (present only if requested)
#' }
#' If \code{show_nuisance=TRUE}, the print method appends a compact section
#' describing nuisance fits.
#'
#' @export
summary.mcee_fit <- function(object,
                             lincomb_alpha = NULL,
                             lincomb_beta  = NULL,
                             lincomb_joint = NULL,
                             conf_level = 0.95,
                             show_nuisance = FALSE,
                             ...) {
  stopifnot(inherits(object, "mcee_fit"))

  fit <- object$mcee_fit
  if (is.null(fit) || is.null(fit$alpha_hat) || is.null(fit$beta_hat) || is.null(fit$varcov)) {
    stop("`mcee_fit` object is missing required components (`alpha_hat`, `beta_hat`, `varcov`).")
  }

  # Keep names on the estimates
  alpha <- fit$alpha_hat
  beta  <- fit$beta_hat
  Vfull <- fit$varcov

  # Basic checks
  p <- length(alpha)
  if (!is.matrix(Vfull) || nrow(Vfull) != 2 * p || ncol(Vfull) != 2 * p) {
    stop("`varcov` must be a (2p x 2p) matrix with p = length(alpha_hat).")
  }

  # ------- Coefficient names (do not overwrite informative names) -------------
  basis_names <- NULL

  if (!is.null(names(alpha)) && length(names(alpha)) == p && all(nzchar(names(alpha)))) {
    basis_names <- names(alpha)
  } else if (!is.null(fit$f) && is.matrix(fit$f) && ncol(fit$f) == p && !is.null(colnames(fit$f))) {
    basis_names <- colnames(fit$f)
  } else if (!is.null(object$meta$basis_names) && length(object$meta$basis_names) == p) {
    basis_names <- object$meta$basis_names
  } else {
    basis_names <- paste0("f", seq_len(p))
  }

  # make sure both alpha and beta carry the same (resolved) names
  names(alpha) <- basis_names
  names(beta)  <- basis_names

  # df = n_ids − 2p if available
  n_ids <- if (!is.null(object$meta$n_ids)) object$meta$n_ids else NA_real_
  df <- if (is.finite(n_ids)) as.numeric(n_ids) - 2 * p else NA_real_

  # Convenience function to build a coefficient table (preserves rownames)
  build_coef_table <- function(est_named, Vblock, rn, df, conf_level) {
    est <- as.numeric(est_named)
    names(est) <- rn
    se   <- sqrt(diag(Vblock))
    tval <- est / se
    pval <- 2 * stats::pt(abs(tval), df = df, lower.tail = FALSE)
    alpha_lvl <- 1 - conf_level
    tcrit <- stats::qt(1 - alpha_lvl / 2, df = df)
    lcl <- est - tcrit * se
    ucl <- est + tcrit * se

    tab <- data.frame(
      Estimate     = as.numeric(est),
      `Std. Error` = as.numeric(se),
      `t value`    = as.numeric(tval),
      df           = rep_len(df, length(tval)),
      `Pr(>|t|)`   = as.numeric(pval),
      check.names  = FALSE,
      row.names    = rn
    )
    tab[[sprintf("%g%% LCL", 100 * conf_level)]] <- as.numeric(lcl)
    tab[[sprintf("%g%% UCL", 100 * conf_level)]] <- as.numeric(ucl)
    tab
  }

  # Extract alpha/beta blocks
  Valpha <- Vfull[seq_len(p), seq_len(p), drop = FALSE]
  Vbeta  <- Vfull[(p + 1):(2 * p), (p + 1):(2 * p), drop = FALSE]

  alpha_tab <- build_coef_table(alpha, Valpha, basis_names, df, conf_level)
  beta_tab  <- build_coef_table(beta,  Vbeta,  basis_names, df, conf_level)

  # Linear combinations helper
  lincomb_table <- function(L, est_named, V, rn_prefix, df, conf_level) {
    est <- as.numeric(est_named)
    # Normalize L into a matrix
    if (is.numeric(L) && is.null(dim(L))) {
      if (length(L) != length(est)) stop("`lincomb` has wrong length.")
      L <- matrix(L, nrow = 1)
      rownames(L) <- rn_prefix
    } else if (is.matrix(L)) {
      if (ncol(L) != length(est)) stop("`lincomb` must have ", length(est), " columns.")
      if (is.null(rownames(L))) rownames(L) <- paste0(rn_prefix, seq_len(nrow(L)))
    } else {
      stop("`lincomb` must be a numeric vector or matrix.")
    }

    estL <- as.vector(L %*% est)
    varL <- diag(L %*% V %*% t(L))
    seL  <- sqrt(varL)
    tL   <- estL / seL
    pL   <- 2 * stats::pt(abs(tL), df = df, lower.tail = FALSE)
    alpha_lvl <- 1 - conf_level
    tcrit <- stats::qt(1 - alpha_lvl / 2, df = df)
    lclL <- estL - tcrit * seL
    uclL <- estL + tcrit * seL

    tab <- data.frame(
      Estimate     = estL,
      `Std. Error` = seL,
      `t value`    = tL,
      df           = rep_len(df, length(tL)),
      `Pr(>|t|)`   = pL,
      check.names  = FALSE
    )
    tab[[sprintf("%g%% LCL", 100 * conf_level)]] <- lclL
    tab[[sprintf("%g%% UCL", 100 * conf_level)]] <- uclL
    rownames(tab) <- rownames(L)
    tab
  }

  lin_alpha_tab <- if (!is.null(lincomb_alpha)) {
    lincomb_table(lincomb_alpha, alpha, Valpha, "Lα", df, conf_level)
  } else NULL

  lin_beta_tab <- if (!is.null(lincomb_beta)) {
    lincomb_table(lincomb_beta, beta, Vbeta, "Lβ", df, conf_level)
  } else NULL

  lin_joint_tab <- NULL
  if (!is.null(lincomb_joint)) {
    theta <- c(alpha, beta)
    rn_theta <- c(paste0("α:", basis_names), paste0("β:", basis_names))
    names(theta) <- rn_theta
    lin_joint_tab <- lincomb_table(lincomb_joint, theta, Vfull, "Lθ", df, conf_level)
  }

  out <- list(
    call          = object$call,
    df            = df,
    conf_level    = conf_level,
    alpha         = alpha_tab,
    beta          = beta_tab,
    lincomb_alpha = lin_alpha_tab,
    lincomb_beta  = lin_beta_tab,
    lincomb_joint = lin_joint_tab,
    .nuisance     = if (isTRUE(show_nuisance)) object$nuisance_models else NULL
  )
  class(out) <- "summary.mcee_fit"
  out
}

#' @export
print.summary.mcee_fit <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)

  cat(sprintf("\nInference: small-sample t; df = %s\n",
              ifelse(is.na(x$df), "NA", format(x$df, digits = 6))))
  cat(sprintf("Confidence level: %g%%\n", 100 * x$conf_level))

  cat("\nNatural Direct Excursion Effect (alpha):\n")
  .mcee_print_coef_table(x$alpha)

  cat("\nNatural Indirect Excursion Effect (beta):\n")
  .mcee_print_coef_table(x$beta)

  if (!is.null(x$lincomb_alpha)) {
    cat("\nLinear combinations of alpha (L * alpha):\n")
    .mcee_print_coef_table(x$lincomb_alpha)
  }
  if (!is.null(x$lincomb_beta)) {
    cat("\nLinear combinations of beta (L * beta):\n")
    .mcee_print_coef_table(x$lincomb_beta)
  }
  if (!is.null(x$lincomb_joint)) {
    cat("\nJoint linear combinations (L * (alpha, beta)):\n")
    .mcee_print_coef_table(x$lincomb_joint)
  }

  if (!is.null(x$.nuisance)) {
    cat("\nStage-1 nuisance fits:\n")

    if (is.character(x$.nuisance) && length(x$.nuisance) == 1L) {
      # e.g., from mcee_userfit_nuisance:
      # "Fitted values for all nuisance functions were supplied by the user."
      cat("  ", x$.nuisance, "\n", sep = "")
    } else if (is.list(x$.nuisance)) {
      nm <- names(x$.nuisance)
      if (is.null(nm)) nm <- as.character(seq_along(x$.nuisance))
      for (k in nm) {
        cat(sprintf("  [%s] %s\n", k, .mcee_compact_model_info(x$.nuisance[[k]])))
      }
      cat("  Note: For full details, inspect `$nuisance_models` directly.\n")
    } else {
      cat("  (Unrecognized `nuisance_models` structure.)\n")
    }
  }

  invisible(x)
}

# --- small internal printers --------------------------------------------

.mcee_print_coef_table <- function(tab) {
  # Reorder columns for nicer print: Estimate, LCL, UCL, StdErr, t, df, p
  nm <- colnames(tab)
  lcl <- grep("LCL$", nm)
  ucl <- grep("UCL$", nm)
  ord <- c(
    match("Estimate", nm),
    lcl, ucl,
    match("Std. Error", nm),
    match("t value", nm),
    match("df", nm),
    match("Pr(>|t|)", nm)
  )
  ord <- ord[!is.na(ord)]
  tab <- tab[, ord, drop = FALSE]
  printCoefmat(as.matrix(tab), P.values = TRUE, has.Pvalue = TRUE)
}

# Try to produce a one-line description of a nuisance object
.mcee_compact_model_info <- function(obj) {
  # Known constants (our helper stores a small list)
  if (is.list(obj) && !is.null(obj$type) && identical(obj$type, "known")) {
    val <- tryCatch(as.numeric(obj$value), error = function(e) NA_real_)
    if (length(unique(val)) == 1) {
      # if val only have a single level, print out the single level
      return(paste0("known constant: ", unique(val)))
    } else if (length(unique(val)) >= 2) {
      # val has multiple values, print out a summary table
      return(paste0("known constant: multiple values (e.g., ",
                    paste(head(val, 3), collapse = ", "),
                    if (length(val) > 3) ", ..." else "",
                    ")"))
    } else {
      stop("Length of unique values in known constant is 0.")
    }
  }
  # Standard model objects
  cl <- class(obj)[1]
  # Try to extract a compact formula RHS
  rhs <- NULL
  fm <- tryCatch(obj$call$formula, error = function(e) NULL)
  if (!is.null(fm)) {
    rhs <- tryCatch(as.character(stats::as.formula(fm))[3L], error = function(e) NULL)
  }
  if (!is.null(rhs)) {
    return(sprintf("%s: ~ %s", cl, rhs))
  }
  sprintf("%s object", cl)
}

# Safe "NULL coalesce"
`%||%` <- function(a, b) if (!is.null(a)) a else b

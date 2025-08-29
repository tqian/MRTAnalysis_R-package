#' Summary for DCEE fits
#'
#' Produce inference tables for distal causal excursion effects from a
#' [dcee()] model. By default uses small-sample \eqn{t}-tests with
#' \code{df = object$df} (subjects minus number of betas). If \code{df}
#' is missing or nonpositive, falls back to large-sample normal (z) inference.
#'
#' @param object An object of class \code{"dcee_fit"} returned by [dcee()].
#' @param lincomb Optional numeric vector or matrix specifying linear
#'   combinations \eqn{L \beta}. If a vector of length \eqn{p} (number of betas),
#'   a single linear combination is evaluated. If a matrix, it must have \eqn{p}
#'   columns; each row defines one combination. Row names (if present) are used
#'   as labels.
#' @param conf_level Confidence level for intervals (default \code{0.95}).
#' @param show_control_fit Logical; if \code{TRUE}, include compact information about
#'   the Stage-1 nuisance regressions (if available). When \code{cross_fit = TRUE}
#'   in [dcee()], \code{regfit_a0}/\code{regfit_a1} refer to the \emph{last fold}
#'   fit and are provided for inspection only.
#' @param ... Currently ignored.
#'
#' @return A list of class \code{"summary.dcee_fit"} with components:
#' \itemize{
#'   \item \code{call} — the original call
#'   \item \code{df} — degrees of freedom used for t-tests (may be \code{NA})
#'   \item \code{conf_level} — the confidence level
#'   \item \code{excursion_effect} — data frame with coefficient table for \eqn{\beta}
#'   \item \code{lincomb} — optional data frame with linear-combination results
#'   \item \code{control_fit} — optional list describing Stage-1 fits (only if \code{show_control_fit})
#' }
#'
#' @export
#' @method summary dcee_fit
summary.dcee_fit <- function(object,
                             lincomb = NULL,
                             conf_level = 0.95,
                             show_control_fit = FALSE,
                             ...) {
    stopifnot(inherits(object, "dcee_fit"))
    fit <- object$fit
    if (is.null(fit$beta_hat) || is.null(fit$beta_varcov)) {
        stop("dcee_fit$fit must contain `beta_hat` and `beta_varcov`.")
    }

    beta <- fit$beta_hat
    V <- fit$beta_varcov
    p <- length(beta)

    # df & distribution
    df <- object$df
    use_t <- is.finite(df) && df > 0
    alpha <- 1 - conf_level
    crit <- if (use_t) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)

    se <- sqrt(diag(V))
    stat <- beta / se
    pval <- if (use_t) {
        2 * stats::pt(abs(stat), df = df, lower.tail = FALSE)
    } else {
        2 * stats::pnorm(abs(stat), lower.tail = FALSE)
    }
    lcl <- beta - crit * se
    ucl <- beta + crit * se

    coef_tab <- data.frame(
        Estimate = as.numeric(beta),
        `Std. Error` = as.numeric(se),
        `t value` = as.numeric(stat),
        df = rep_len(if (use_t) df else NA_real_, length(stat)),
        `Pr(>|t value|)` = as.numeric(pval),
        check.names = FALSE
    )
    coef_tab[[sprintf("%g%% LCL", 100 * conf_level)]] <- as.numeric(lcl)
    coef_tab[[sprintf("%g%% UCL", 100 * conf_level)]] <- as.numeric(ucl)
    rownames(coef_tab) <- names(beta)

    # Linear combinations
    lincomb_tab <- NULL
    if (!is.null(lincomb)) {
        L <- lincomb
        if (is.numeric(L) && is.null(dim(L))) {
            if (length(L) != p) stop("`lincomb` vector must have length ", p, ".")
            lbl <- if (!is.null(names(L))) paste0("L: ", paste(names(L)[L != 0], collapse = " + ")) else "L1"
            L <- matrix(L, nrow = 1)
            rownames(L) <- lbl
        } else if (is.matrix(L)) {
            if (ncol(L) != p) {
                if (nrow(L) == p) L <- t(L) else stop("`lincomb` matrix must have ", p, " columns.")
            }
            if (is.null(rownames(L))) rownames(L) <- paste0("L", seq_len(nrow(L)))
        } else {
            stop("`lincomb` must be a numeric vector or matrix.")
        }

        estL <- as.vector(L %*% beta)
        varL <- diag(L %*% V %*% t(L))
        seL <- sqrt(varL)
        statL <- estL / seL
        pL <- if (use_t) {
            2 * stats::pt(abs(statL), df = df, lower.tail = FALSE)
        } else {
            2 * stats::pnorm(abs(statL), lower.tail = FALSE)
        }
        lclL <- estL - crit * seL
        uclL <- estL + crit * seL

        lincomb_tab <- data.frame(
            Estimate = estL,
            `Std. Error` = seL,
            `t value` = statL,
            df = rep_len(if (use_t) df else NA_real_, length(statL)),
            `Pr(>|t value|)` = pL,
            check.names = FALSE
        )
        lincomb_tab[[sprintf("%g%% LCL", 100 * conf_level)]] <- lclL
        lincomb_tab[[sprintf("%g%% UCL", 100 * conf_level)]] <- uclL
        rownames(lincomb_tab) <- rownames(L)
    }

    out <- list(
        call = object$call,
        df = if (use_t) df else NA_real_,
        conf_level = conf_level,
        excursion_effect = coef_tab,
        lincomb = lincomb_tab
    )

    if (isTRUE(show_control_fit)) {
        # Store actual fit objects; print method will render proper summaries.
        ctrl <- list()
        if (!is.null(fit$regfit_a0)) ctrl$A0 <- fit$regfit_a0
        if (!is.null(fit$regfit_a1)) ctrl$A1 <- fit$regfit_a1
        if (length(ctrl)) {
            ctrl$note <- "If cross_fit = TRUE in dcee(), these are last-fold models for inspection only."
            out$control_fit <- ctrl
        }
    }

    class(out) <- "summary.dcee_fit"
    out
}

#' @export
#' @method print summary.dcee_fit
print.summary.dcee_fit <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    if (is.na(x$df)) {
        cat("\nInference: large-sample normal (z); df not available/positive.\n")
    } else {
        cat(sprintf("\nInference: small-sample t; df = %d\n", x$df))
    }
    cat(sprintf("Confidence level: %g%%\n", 100 * x$conf_level))

    cat("\nDistal causal excursion effect (beta):\n")
    .dcee_print_coef_table(x$excursion_effect)

    if (!is.null(x$lincomb)) {
        cat("\nLinear combinations (L * beta):\n")
        .dcee_print_coef_table(x$lincomb)
    }

    if (!is.null(x$control_fit)) {
        cat("\nStage-1 nuisance fits:\n")

        # Recover labels from the original dcee() call (best-effort)
        oc <- x$call
        outcome_lbl <- tryCatch(as.character(oc$outcome), error = function(e) NULL)
        rhs_lbl <- tryCatch(as.character(stats::as.formula(oc$control_formula))[2L],
            error = function(e) NULL
        )
        trt_lbl <- tryCatch(as.character(oc$treatment), error = function(e) NULL)

        # [A = 0]
        cat("\n  [A = 0] Outcome model\n")
        .dcee_render_stage1(
            obj = x$control_fit$A0,
            subset_label = paste0(trt_lbl %||% "A", " == 0"),
            outcome = outcome_lbl,
            rhs = rhs_lbl,
            indent = 2
        )

        # [A = 1]
        cat("\n  [A = 1] Outcome model\n")
        .dcee_render_stage1(
            obj = x$control_fit$A1,
            subset_label = paste0(trt_lbl %||% "A", " == 1"),
            outcome = outcome_lbl,
            rhs = rhs_lbl,
            indent = 2
        )

        # Cross-fitting caveat (if present)
        if (!is.null(x$control_fit$note)) {
            cat(sprintf("\n  Note: %s\n", x$control_fit$note))
        }

        # Always add this final hint
        cat("\n  For full details, inspect $fit$regfit_a0 and $fit$regfit_a1.\n")
    }

    invisible(x)
}

# ---- helpers ---------------------------------------------------------------

.dcee_print_coef_table <- function(tab) {
    nm <- colnames(tab)
    lcl <- grep(" LCL$", nm)
    ucl <- grep(" UCL$", nm)
    ord <- c(
        match("Estimate", nm), lcl, ucl,
        match("Std. Error", nm), match("t value", nm),
        match("df", nm), match("Pr(>|t value|)", nm)
    )
    ord <- ord[!is.na(ord)]
    tab <- tab[, ord, drop = FALSE]
    printCoefmat(as.matrix(tab), P.values = TRUE, has.Pvalue = TRUE)
}

# Safely print a model summary, with indentation and fallbacks.
.dcee_render_model_summary <- function(obj, indent = 0) {
    pad <- if (indent > 0) paste(rep(" ", indent), collapse = "") else ""
    pr <- function(lines) cat(paste0(pad, lines), sep = "\n")

    # Try summary() first
    ok <- FALSE
    so <- tryCatch(summary(obj), error = function(e) e)
    if (!inherits(so, "error")) {
        out <- utils::capture.output(print(so))
        pr(out)
        ok <- TRUE
    }

    # If summary() failed or was too terse, try print(obj)
    if (!ok) {
        out <- tryCatch(utils::capture.output(print(obj)), error = function(e) NULL)
        if (!is.null(out)) {
            pr(out)
            ok <- TRUE
        }
    }

    # Last resort: compact info
    if (!ok) {
        pr(sprintf("<%s> (no printable summary available)", paste(class(obj), collapse = "/")))
    }
}


`%||%` <- function(a, b) if (is.null(a)) b else a

# Primary method label based on class
.dcee_primary_method <- function(obj) {
    if (is.null(obj)) {
        return("set_to_zero")
    }
    cl <- class(obj)
    prefs <- c("SuperLearner", "gam", "ranger", "randomForest", "glm", "lm")
    hit <- prefs[prefs %in% cl][1]
    if (!is.na(hit)) {
        if (hit == "gam") {
            return("mgcv::gam")
        }
        return(hit)
    }
    cl[1]
}

# Make header for lm/gam branch
.dcee_make_stage1_header <- function(obj, outcome, rhs, subset) {
    method <- .dcee_primary_method(obj)
    n <- tryCatch(length(stats::fitted(obj)), error = function(e) NA_integer_)
    fml <- if (!is.null(outcome) && !is.null(rhs)) sprintf("%s ~ %s", outcome, rhs) else NULL
    pieces <- c(
        sprintf("Method: %s", method),
        if (!is.null(fml)) sprintf("Formula: %s", fml) else NULL,
        sprintf("Subset: %s", subset),
        if (is.finite(n)) sprintf("#person-decision-points = %d", n) else NULL
    )
    paste(pieces, collapse = " | ")
}

# Render a stage-1 model according to your rules
.dcee_render_stage1 <- function(obj, subset_label, outcome, rhs, indent = 0) {
    pad <- if (indent > 0) paste(rep(" ", indent), collapse = "") else ""
    pr <- function(lines) cat(paste0(pad, lines), sep = "\n")

    method <- .dcee_primary_method(obj)

    # set_to_zero case (obj is NULL)
    if (identical(method, "set_to_zero") || is.null(obj)) {
        pr("Nuisance outcome model set to zero (control_reg_method = 'set_to_zero').")
        return(invisible())
    }

    # Branch: lm / gam → header + summary without 'Call:'
    if (method %in% c("lm", "mgcv::gam", "glm")) {
        pr(.dcee_make_stage1_header(obj, outcome = outcome, rhs = rhs, subset = subset_label))
        # capture summary() and strip 'Call:' block
        s <- tryCatch(summary(obj), error = function(e) e)
        if (!inherits(s, "error")) {
            lines <- utils::capture.output(print(s))
            lines <- .dcee_drop_call_block(lines)
            pr(lines)
            return(invisible())
        }
        # fallback to print(obj)
        out <- tryCatch(utils::capture.output(print(obj)), error = function(e) NULL)
        if (!is.null(out)) pr(out)
        return(invisible())
    }

    # Branch: rf / ranger / SuperLearner → keep original call & default print()
    if (method %in% c("randomForest", "ranger", "SuperLearner")) {
        # Just default-print, which typically includes the original call
        out <- tryCatch(utils::capture.output(print(obj)), error = function(e) NULL)
        if (!is.null(out)) {
            pr(out)
        } else {
            pr(sprintf("<%s> (no printable output)", paste(class(obj), collapse = "/")))
        }
        return(invisible())
    }

    # Unknown: try summary then print
    out <- tryCatch(utils::capture.output(print(summary(obj))), error = function(e) NULL)
    if (!is.null(out)) {
        out <- .dcee_drop_call_block(out)
        pr(out)
        return(invisible())
    }
    out <- tryCatch(utils::capture.output(print(obj)), error = function(e) NULL)
    if (!is.null(out)) pr(out) else pr(sprintf("<%s> (no printable output)", paste(class(obj), collapse = "/")))
}

# Strip the 'Call:' block from captured output
.dcee_drop_call_block <- function(xs) {
    i <- which(trimws(xs) == "Call:")
    if (!length(i)) {
        return(xs)
    }
    # find first empty line after 'Call:'
    jrel <- which(xs[(i + 1):length(xs)] == "")[1]
    if (is.na(jrel)) {
        return(xs[-c(i:length(xs))])
    }
    xs[-c(i:(i + jrel - 1))]
}

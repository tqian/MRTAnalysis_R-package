#' Summarize Causal Excursion Effect Fits for MRT with Continuous Outcomes
#'
#' \code{summary} method for class "wcls_fit".
#'
#' @param object An object of class "wcls_fit".
#' @param lincomb A vector of length p (p is the number of moderators including
#' intercept) or a matrix with p columns. When not set to `NULL`,
#' the summary will include the specified linear combinations of the causal excursion
#' effect coefficients and the corresponding confidence interval, standard error,
#' and p-value.
#' @param conf_level A numeric value indicating the confidence level for
#' confidence intervals. Default to 0.95.
#' @param show_control_fit A logical value of whether the fitted coefficients
#' for the control variables will be printed in the summary. Default to FALSE.
#' (Interpreting the fitted coefficients for control variables is not recommended.)
#' @param ... Further arguments passed to or from other methods.
#'
#' @return the original function call and the estimated causal excursion effect
#' coefficients, 95% confidence interval, Standard Error, Hotelling's T-statistic
#' value or Wald-statistic value (depending on whether sample size is < 50),
#' degrees of freedom, and p-value.
#' @import stats
#' @export
#'
#' @examples fit <- wcls(
#'     data = data_mimicHeartSteps,
#'     id = "userid",
#'     outcome = "logstep_30min",
#'     treatment = "intervention",
#'     rand_prob = 0.6,
#'     moderator_formula = ~1,
#'     control_formula = ~logstep_pre30min,
#'     availability = "avail",
#'     numerator_prob = 0.6
#' )
#' summary(fit)
summary.wcls_fit <- function(
    object,
    lincomb = NULL,
    conf_level = 0.95,
    show_control_fit = FALSE,
    ...) {
    # default value from Boruvka's original implementation
    combos <- NULL
    omnibus <- FALSE
    null <- 0
    small <- TRUE
    conf.int <- 0.95
    normal <- FALSE

    p <- length(object$moderator_names)
    q <- length(object$control_names)

    if (!is.null(lincomb)) {
        stopifnot(is.numeric(lincomb))
        if (is.vector(lincomb)) {
            stopifnot(length(lincomb) == p)
            lincomb <- matrix(lincomb, nrow = 1)
        } else if (is.matrix(lincomb)) {
            stopifnot(ncol(lincomb) == p)
        } else {
            stop(paste0(
                "lincomb needs to be a vector of length p or a matrix with p columns. ",
                "(p is the number of moderators including intercept.)"
            ))
        }
    }

    if (is.null(combos)) {
        # combos is used to extract all the coefficients (alpha and beta)
        combos <- diag(length(coef(object)))
        rownames(combos) <- names(coef(object))
        omnibus <- FALSE
    }

    if (!is.null(lincomb)) {
        # add lincomb (for beta only) to combos

        # construct variable names for linear combination terms
        rowname_lincomb <- c()
        rowname_out_beta <- object$moderator_names
        if (!is.null(lincomb)) {
            for (irow in 1:nrow(lincomb)) {
                this_name <- ""
                for (icol in 1:ncol(lincomb)) {
                    lincomb_entry <- lincomb[irow, icol]
                    if (lincomb_entry == 1) {
                        if (this_name != "") {
                            this_name <- paste0(this_name, " + ", rowname_out_beta[icol])
                        } else {
                            this_name <- paste0(rowname_out_beta[icol])
                        }
                    } else if (lincomb_entry > 0) {
                        if (this_name != "") {
                            this_name <- paste0(
                                this_name, " + ", lincomb_entry, "*",
                                rowname_out_beta[icol]
                            )
                        } else {
                            this_name <- paste0(lincomb_entry, "*", rowname_out_beta[icol])
                        }
                    } else if (lincomb_entry == -1) {
                        if (this_name != "") {
                            this_name <- paste0(this_name, " - ", rowname_out_beta[icol])
                        } else {
                            this_name <- paste0("-", rowname_out_beta[icol])
                        }
                    } else if (lincomb_entry < 0) {
                        if (this_name != "") {
                            this_name <- paste0(
                                this_name, " - ", abs(lincomb_entry), "*",
                                rowname_out_beta[icol]
                            )
                        } else {
                            this_name <- paste0(
                                "- ", abs(lincomb_entry), "*",
                                rowname_out_beta[icol]
                            )
                        }
                    }
                }
                rowname_lincomb <- c(rowname_lincomb, this_name)
            }
        }
        lincomb <- cbind(matrix(0, nrow = nrow(lincomb), ncol = q), lincomb)
        rownames(lincomb) <- rowname_lincomb
        combos <- rbind(combos, lincomb)
    }

    est <- combos %*% coef(object)
    if (nrow(est) != length(null)) null <- rep(null[1], nrow(est))
    ## apply Mancl and DeRouen's (2001) small sample correction
    if (is.logical(small)) small <- small * 50
    n <- cluster.number(object, overall = FALSE)
    d1 <- if (omnibus) {
        nrow(combos)
    } else {
        apply(combos != 0, 1, sum)
    }
    d2 <- n - length(coef(object))
    ## apply Hotelling's T-squared test, following Liao et al. (2016)
    if (n <= small & !normal) {
        type <- "Hotelling"
        adj <- d1 * (d1 + d2 - 1) / d2
        qfun <- function(p) mapply(qf, p = p, df1 = d1, df2 = d2) / adj
        pfun <- function(q) 1 - mapply(pf, q = q * adj, df1 = d1, df2 = d2)
    } else {
        type <- "Wald"
        qfun <- if (normal) {
            function(p) qnorm((1 + p) / 2)
        } else {
            function(p) mapply(qf, p = p, df1 = d1, df2 = d2)
        }
        pfun <- if (normal) {
            function(q) 1 - mapply(pchisq, q = q, df = d1)
        } else {
            function(q) 1 - mapply(pf, q = q, df1 = d1, df2 = d2)
        }
    }
    var.est <- combos %*% vcov_geeglm(object, small = small) %*% t(combos)
    se.est <- sqrt(diag(var.est))
    crit <- sqrt(qfun(conf.int))
    lcl <- est - se.est * crit
    ucl <- est + se.est * crit
    stat <- if (omnibus) {
        rep(t(est - null) %*% solve(var.est) %*% (est - null), d1)
    } else {
        (est - null)^2 / diag(var.est)
    }
    pvalue <- pfun(stat)
    out_alpha_beta_lincomb <- cbind(est, lcl, ucl, se.est, stat, d1, d2, pvalue)
    rownames(out_alpha_beta_lincomb) <- rownames(combos)
    colnames(out_alpha_beta_lincomb) <- c(
        "Estimate",
        paste0(round(conf.int * 100), "% ", c("LCL", "UCL")),
        "StdErr", type, "df1", "df2", "p-value"
    )

    if (nrow(out_alpha_beta_lincomb) - q == 1) {
        out_beta <- matrix(out_alpha_beta_lincomb[(q + 1):nrow(out_alpha_beta_lincomb), ], nrow = 1)
    } else {
        out_beta <- out_alpha_beta_lincomb[(q + 1):nrow(out_alpha_beta_lincomb), ]
    }
    colnames(out_beta) <- colnames(out_alpha_beta_lincomb)
    if (is.null(lincomb)) {
        rownames(out_beta) <- object$moderator_names
    } else {
        rownames(out_beta) <- c(object$moderator_names, rowname_lincomb)
    }

    res <- list(
        call = object$call,
        causal_excursion_effect = out_beta
    )

    if (show_control_fit) {
        if (q == 1) {
            out_alpha <- matrix(out_alpha_beta_lincomb[1:q, ], nrow = 1)
        } else if (q > 1) {
            out_alpha <- out_alpha_beta_lincomb[1:q, ]
        }
        rownames(out_alpha) <- object$control_names
        colnames(out_alpha) <- colnames(out_alpha_beta_lincomb)

        res <- c(res, list(control_variables = out_alpha))
        message("Interpreting the fitted coefficients for control variables is not recommended.")
    }

    res
}

# print.wcls_fit <- function(object) {
#     summary(object)
# }

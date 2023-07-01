#' Summarize Causal Excursion Effect Fits for MRT with Binary Outcomes
#'
#' \code{summary} method for class "emee_fit".
#'
#' @param object An object of class "emee_fit".
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
#' coefficients, confidence interval with conf_level, standard error, t-statistic value,
#' degrees of freedom, and p-value.
#' @import stats
#' @export
#' @examples fit <- emee(
#'     data = data_binary,
#'     id = "userid",
#'     outcome = "Y",
#'     treatment = "A",
#'     rand_prob = "rand_prob",
#'     moderator_formula = ~time_var1,
#'     control_formula = ~ time_var1 + time_var2,
#'     availability = "avail",
#'     numerator_prob = 0.5,
#'     start = NULL
#' )
#' summary(fit)
summary.emee_fit <- function(
    object,
    lincomb = NULL,
    conf_level = 0.95,
    show_control_fit = FALSE,
    ...) {
    p <- length(object$fit$beta_hat)
    q <- length(object$fit$alpha_hat)

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

    # output table for beta
    est <- object$fit$beta_hat
    se <- object$fit$beta_se_adjusted

    if (!is.null(lincomb)) {
        est_lincomb <- as.vector(lincomb %*% est)
        varcov_beta <- object$fit$varcov_adjusted[(q + 1):(q + p), (q + 1):(q + p)]
        se_lincomb <- sqrt(diag(lincomb %*% varcov_beta %*% t(lincomb)))

        est <- c(est, est_lincomb)
        se <- c(se, se_lincomb)
    }

    crit <- qt(1 - (1 - conf_level) / 2, df = object$df)
    lcl <- est - se * crit
    ucl <- est + se * crit
    t_value <- est / se
    p_value <- 2 * pt(-abs(t_value), df = object$df)

    out_beta <- cbind(est, lcl, ucl, se, t_value, object$df, p_value)
    colnames(out_beta) <- c(
        "Estimate",
        paste0(round(conf_level * 100), "% ", c("LCL", "UCL")),
        "StdErr", "t_value", "df", "p-value"
    )

    rowname_out_beta <- names(object$fit$beta_hat)

    # construct variable names for linear combination terms
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
            rowname_out_beta <- c(rowname_out_beta, this_name)
        }
    }
    rownames(out_beta) <- rowname_out_beta

    res <- list(
        call = object$call,
        causal_excursion_effect = out_beta
    )

    if (show_control_fit) {
        # output table for beta
        est <- object$fit$alpha_hat
        se <- object$fit$alpha_se_adjusted
        crit <- qt(1 - (1 - conf_level) / 2, df = object$df)
        lcl <- est - se * crit
        ucl <- est + se * crit
        t_value <- est / se
        p_value <- 2 * pt(-abs(t_value), df = object$df)

        out_alpha <- cbind(est, lcl, ucl, se, t_value, object$df, p_value)
        colnames(out_alpha) <- c(
            "Estimate",
            paste0(round(conf_level * 100), "% ", c("LCL", "UCL")),
            "StdErr", "t_value", "df", "p-value"
        )

        rownames(out_alpha) <- names(object$fit$alpha_hat)

        res <- c(res, list(control_variables = out_alpha))
        message("Interpreting the fitted coefficients for control variables is not recommended.")
    }

    res
}

# print.emee_fit <- function(object) {
#     summary(object)
# }

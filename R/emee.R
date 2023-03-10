#' Estimates the causal excursion effect for binary outcome MRT
#'
#' Returns the estimated causal excursion effect (on log relative risk scale) and the estimated standard error.
#' Small sample correction using the "Hat" matrix in the variance estimate is implemented.
#' All variables should correspond to columns in data and should not be in quotation marks.
#'
#' @param data A data set in long format.
#' @param id The subject id variable.
#' @param outcome The outcome variable.
#' @param treatment The binary treatment assignment variable.
#' @param rand_prob The randomization probability variable.
#' @param moderator_formula A formula for the moderator variables. This should
#' start with ~ followed by the moderator variables. When set to \code{~ 1}, a fully
#' marginal excursion effect (no moderators) is estimated.
#' @param control_formula A formula for the control variables. This should
#' start with ~ followed by the control variables. When set to \code{~ 1}, only an
#' intercept is included as the control variable.
#' @param availability The availability variable. Use the default value (\code{NULL})
#' if your MRT doesn't have availability considerations.
#' @param numerator_prob Either a number between 0 and 1, or a variable name for
#' a column in data. If you are not sure what this is, use the default value (\code{NULL}).
#' @param start A vector of the initial value of the estimators used in the numerical
#' solver. If using default value (\code{NULL}), a vector of 0 will be used internally.
#' If specifying a non-default value, this needs to be a numeric vector of length
#' (number of moderator variables including the intercept) +
#' (number of control variables including the intercept).
#' @param verbose If default (`TRUE`), additional messages will be printed
#' during data preprocessing.
#'
#' @return An object of type "emee_fit"
# The estimated beta with its intercept and alpha with its intercept,
#         the standard error of the estimated beta with its intercept and alpha with its intercept,
#         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
#         the dimension of the moderated variables, and the dimension of the control variables,
#         the value of the estimating function at the estimated beta and alpha
#' @import rootSolve stats
#' @export
#'
#' @examples
#'
#' ## estimating the fully marginal excursion effect by setting
#' ## moderator_formula = ~ 1
#' emee(
#'     data = data_binary,
#'     id = userid,
#'     outcome = Y,
#'     treatment = A,
#'     rand_prob = rand_prob,
#'     moderator_formula = ~1,
#'     control_formula = ~ time_var1 + time_var2,
#'     availability = avail
#' )
#'
#' ## estimating the causal excursion effect moderated by time_var1
#' ## by setting moderator_formula = ~ time_var1
#' emee(
#'     data = data_binary,
#'     id = userid,
#'     outcome = Y,
#'     treatment = A,
#'     rand_prob = rand_prob,
#'     moderator_formula = ~time_var1,
#'     control_formula = ~ time_var1 + time_var2,
#'     availability = avail
#' )
emee <- function(data,
                 id,
                 outcome,
                 treatment,
                 rand_prob,
                 moderator_formula,
                 control_formula,
                 availability = NULL,
                 numerator_prob = NULL,
                 start = NULL,
                 verbose = TRUE) {
    # Save input to `mf`
    mf <- match.call(expand.dots = FALSE)

    # preprocess input, including checking types and missingness
    preprocessed_input <- preprocess_input(data, mf,
        outcome_type = "binary",
        verbose = verbose
    )
    id <- preprocessed_input$id
    outcome <- preprocessed_input$outcome
    treatment <- preprocessed_input$treatment
    rand_prob <- preprocessed_input$rand_prob
    moderator_matrix <- preprocessed_input$moderator_matrix
    control_matrix <- preprocessed_input$control_matrix
    availability <- preprocessed_input$availability
    numerator_prob <- preprocessed_input$numerator_prob

    # model fitting
    fit <- compute_emee(
        id = as.numeric(as.factor(id)), # convert id into consecutive integers
        A = treatment,
        Y = outcome,
        p_t = rand_prob,
        avail = availability,
        p_t_tilde = numerator_prob,
        moderator_matrix = moderator_matrix,
        control_matrix = control_matrix,
        estimator_initial_value = mf$start
    )

    # output call, estimators, df
    output <- list(fit)
    names(output) <- c("fit")
    output$call <- match.call()
    output$df <- length(unique(id)) - output$fit$dims$p - output$fit$dims$q
    ord <- c("call", "fit", "df")
    output <- output[ord]
    class(output) <- "emee_fit"

    output
}

#' Estimates the marginal excursion effect for binary outcome MRT
#'
#' Returns the estimates for the marginal excursion effect estimator
#' and provides the estimated variance and standard error for the estimators,
#' with small sample correction for an MRT with binary outcome using the "Hat" matrix
#' in the variance estimate and t-distribution or F-distribution critical value with
#' corrected degrees of freedom.
#'
#' @param id a vector which identifies the subject IDs. The length of ‘id’ should be the same as the number of observations
#' @param A the variable that specifies the assigned treatments for subjects
#' @param Y the variable that specifies the outcomes for subjects
#' @param p_t the variable that specifies the treatment randomizing probability in the data set
#' @param avail the variable that specifies the availability of the subjects, default to be always available at any decision points using NULL
#' @param p_t_tilde a number between 0 and 1, default to 0.5 if NULL
#' @param moderator_matrix a matrix of the effect modifiers, this could be NULL
#' @param control_matrix a matrix of the control modifiers, this could be NULL
#' @param estimator_initial_value a numeric vector of the initial value for the estimator,
#'                                its length should be the sum of the length of control and
#'                                moderator variables plus 2
#'                                default to be all 0's using NULL
#'
#' @return the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
#'         the dimension of the moderated variables, and the dimension of the control variables,
#'         the value of the estimating function at the estimated beta and alpha
#' @import rootSolve stats
#' @noRd
#'
#' @examples compute_emee(
#'     id = data_binary$userid,
#'     A = data_binary$A,
#'     Y = data_binary$Y,
#'     p_t = data_binary$rand_prob,
#'     avail = data_binary$avail,
#'     p_t_tilde = rep(0.5, 3000),
#'     moderator_matrix = model.matrix(as.formula("~time_var1"), data = data_binary),
#'     control_matrix = model.matrix(as.formula("~time_var1 + time_var2"), data = data_binary),
#'     estimator_initial_value = NULL
#' )
compute_emee <- function(
    id,
    A,
    Y,
    p_t,
    avail,
    p_t_tilde,
    moderator_matrix,
    control_matrix,
    estimator_initial_value = NULL) {
    Xdm <- moderator_matrix
    Xnames <- colnames(Xdm)
    Zdm <- control_matrix
    Znames <- colnames(Zdm)

    # Estimation --------------------------------------------------------------

    ### 1. preparation ###

    sample_size <- length(unique(id))
    total_person_decisionpoint <- length(id)
    cA <- A - p_t # centered A
    cA_tilde <- A - p_t_tilde

    WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))

    p <- dim(moderator_matrix)[2] # dimension of beta
    q <- dim(control_matrix)[2] # dimension of alpha

    ### 2. estimation ###

    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q + 1):(q + p)])

        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))

        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- exp_AXdm_beta^(-1)

        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum(weight * residual * avail * WCLS_weight * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum(weight * residual * avail * WCLS_weight * cA_tilde * Xdm[, i])
        }

        ef <- ef / sample_size
        return(ef)
    }

    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }

    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside estimator_EMEE():")
            message("\nThe program cannot find a numerical solution to the estimating eqaution.")
            message(cond)
            return(list(
                root = rep(NaN, p + q), msg = cond,
                f.root = rep(NaN, p + q)
            ))
        }
    )

    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)

    ### 3. asymptotic variance ###

    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###

    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p + q, p + q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction

    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p + q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p + q)

    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }

        pre_multiplier <- exp(-A[it] * Xbeta) * WCLS_weight[it]

        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- 0
        partialD_partialtheta[1:q, (q + 1):(q + p)] <- -pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q + 1):(q + p), 1:q] <- 0
        partialD_partialtheta[(q + 1):(q + p), (q + 1):(q + p)] <- -pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])

        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term

        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term

        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- -exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta

        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2, 3), sum) / sample_size
    Mn_inv <- solve(Mn)

    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###

    Sigman_summand <- array(NA, dim = c(sample_size, p + q, p + q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction

    person_first_index <- c(find_change_location(id), total_person_decisionpoint + 1)

    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i]:(person_first_index[i + 1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i]:(person_first_index[i + 1] - 1)], ncol = 1)

        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2, 3), sum) / sample_size

    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q + 1):(q + p)])


    ### 4. small sample correction ###

    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i]:(person_first_index[i + 1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i]:(person_first_index[i + 1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i]:(person_first_index[i + 1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)

        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size

    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q + 1):(q + p)])


    ### 6. return the result with variable names ###

    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames

    ### 7. calculate confidence interval

    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05 / 2, df = sample_size - p - q)
    conf_int_adjusted <- cbind(
        beta_hat - c * beta_se_adjusted,
        beta_hat + c * beta_se_adjusted
    )
    colnames(conf_int) <- colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")

    return(list(
        beta_hat = beta_hat, alpha_hat = alpha_hat,
        beta_se = beta_se, alpha_se = alpha_se,
        beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
        varcov = varcov,
        varcov_adjusted = varcov_adjusted,
        conf_int = conf_int, conf_int_adjusted = conf_int_adjusted,
        dims = list(p = p, q = q),
        f.root = solution$f.root
    ))
}

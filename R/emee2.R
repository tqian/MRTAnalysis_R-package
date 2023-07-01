#' Estimates the causal excursion effect for binary outcome MRT
#'
#' Returns the estimated causal excursion effect (on log relative risk scale) and the estimated standard error.
#' Small sample correction using the "Hat" matrix in the variance estimate is implemented.
#' This is a slightly altered version of \code{emee()}, where the treatment
#' assignment indicator is also centered in the residual term. It would have
#' similar (but not exactly the same) numerical output as \code{emee()}. This
#' is the estimator based on which the sample size calculator for binary outcome
#' MRT is developed. (See R package \code{MRTSampleSizeBinary}.)
#'
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
# Returns the estimated beta with its intercept and alpha with its intercept,
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
#' ## estimating the fully marginal excursion effect by setting
#' ## moderator_formula = ~ 1
#' emee2(
#'     data = data_binary,
#'     id = "userid",
#'     outcome = "Y",
#'     treatment = "A",
#'     rand_prob = "rand_prob",
#'     moderator_formula = ~1,
#'     control_formula = ~ time_var1 + time_var2,
#'     availability = "avail"
#' )
#'
#' ## estimating the causal excursion effect moderated by time_var1
#' ## by setting moderator_formula = ~ time_var1
#' emee2(
#'     data = data_binary,
#'     id = "userid",
#'     outcome = "Y",
#'     treatment = "A",
#'     rand_prob = "rand_prob",
#'     moderator_formula = ~time_var1,
#'     control_formula = ~ time_var1 + time_var2,
#'     availability = "avail"
#' )
emee2 <- function(data,
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
    fit <- compute_emee2(
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
#' where the treatment indicator is always centered
#'
#' This estimator assumes that the randomization probability does not
#' depend on H_t, and that the treatment indicator A is always centered.
#' This estimator returns the estimates for the marginal excursion effect estimator
#' with the above two assumptions and provides the estimated variance and
#' standard error for the estimators, with small sample correction for
#' an MRT with binary outcome using the "Hat" matrix in the variance estimate
#' and t-distribution or F-distribution critical value with corrected degrees of freedom.
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
#' @return Returns the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its
#'         intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha
#'         with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha
#'         with its intercept for small sample,
#'         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence
#'         interval for beta_hat,
#'         the dimension of the moderated variables, and the dimension of the control variables,
#'         the value of the estimating function at the estimated beta and alpha
#' @import rootSolve stats
#' @noRd
#'
#' @examples compute_emee2(
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
compute_emee2 <- function(
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

    # Estimating --------------------------------------------------------------
    ### 1. preparation ###

    sample_size <- length(unique(id))
    total_person_decisionpoint <- length(id)
    cA <- A - p_t # centered A

    p <- dim(moderator_matrix)[2] # dimension of beta
    q <- dim(control_matrix)[2] # dimension of alpha

    ### 2. estimation ###

    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q + 1):(q + p)])

        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_negcAXdm_beta <- exp(-cA * (Xdm %*% beta))

        residual <- exp_negcAXdm_beta * Y - exp_Zdm_alpha

        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum(residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum(residual * avail * cA * Xdm[, i])
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
            message("\nCatched error in multiroot inside estimator_EMEE_alwaysCenterA():")
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

    ### 3. asymptotic variance and small sample correction ###

    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###

    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p + q, p + q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction

    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
        if (p == 1) {
            Xbeta_it <- Xdm[it, ] * beta_hat # a scalar
        } else {
            Xbeta_it <- as.numeric(Xdm[it, ] %*% beta_hat) # a scalar
        }
        if (q == 1) {
            Zalpha_it <- Zdm[it, ] * alpha_hat # a scalar
        } else {
            Zalpha_it <- as.numeric(Zdm[it, ] %*% alpha_hat) # a scalar
        }

        Zdm_it <- as.numeric(Zdm[it, ])
        Xdm_it <- as.numeric(Xdm[it, ])

        stopifnot(class(Zalpha_it) == "numeric")
        stopifnot(class(Zdm_it) == "numeric")
        stopifnot(length(Zdm_it) == q)

        stopifnot(class(Xbeta_it) == "numeric")
        stopifnot(class(Xdm_it) == "numeric")
        stopifnot(length(Xdm_it) == p)

        Mn_summand[it, 1:q, 1:q] <- -avail[it] * exp(Zalpha_it) * (Zdm_it %o% Zdm_it)
        Mn_summand[it, 1:q, (q + 1):(q + p)] <- -avail[it] * exp(-cA[it] * Xbeta_it) * Y[it] * cA[it] * (Zdm_it %o% Xdm_it)
        Mn_summand[it, (q + 1):(q + p), 1:q] <- -avail[it] * exp(Zalpha_it) * cA[it] * (Xdm_it %o% Zdm_it)
        Mn_summand[it, (q + 1):(q + p), (q + 1):(q + p)] <- -avail[it] * exp(-cA[it] * Xbeta_it) * Y[it] * cA[it]^2 * (Xdm_it %o% Xdm_it)
    }
    Mn <- apply(Mn_summand, c(2, 3), sum) / sample_size
    Mn_inv <- solve(Mn)

    ### 3.2 Compute \Sigma_n matrix and \tilde{\Sigma}_n matrix ###

    Sigman_summand <- array(NA, dim = c(sample_size, p + q, p + q))
    Sigman_tilde_summand <- array(NA, dim = c(sample_size, p + q, p + q))

    person_first_index <- c(find_change_location(id), total_person_decisionpoint + 1)

    for (i in 1:sample_size) {
        T_i <- person_first_index[i + 1] - person_first_index[i] # number of time points for individual i
        index_i <- person_first_index[i]:(person_first_index[i + 1] - 1) # row number of individual i's data in dta

        Zdm_i <- Zdm[index_i, ] # each row in Zdm_i corresponds to a time point
        Xdm_i <- Xdm[index_i, ] # each row in Xdm_i corresponds to a time point
        cA_i <- cA[index_i] # each entry in cA_i corresponds to a time point
        Y_i <- Y[index_i]

        if (p == 1) {
            Xdm_i <- matrix(Xdm_i, ncol = 1)
        }
        if (q == 1) {
            Zdm_i <- matrix(Zdm_i, ncol = 1)
        }

        stopifnot(class(Zdm_i)[1] == "matrix")
        stopifnot(class(Xdm_i)[1] == "matrix")
        stopifnot(nrow(Zdm_i) == T_i & ncol(Zdm_i) == q)
        stopifnot(nrow(Xdm_i) == T_i & ncol(Xdm_i) == p)

        D_i <- cbind(Zdm_i, cA_i * Xdm_i)
        r_i <- matrix(exp(-cA_i * Xdm_i %*% beta_hat) * Y_i - exp(Zdm_i %*% alpha_hat), nrow = T_i, ncol = 1)
        I_i <- diag(avail[index_i])

        stopifnot(nrow(D_i) == T_i & ncol(D_i) == (q + p))
        stopifnot(nrow(r_i) == T_i & ncol(r_i) == 1)
        stopifnot(nrow(I_i) == T_i & ncol(I_i) == T_i)

        Sigman_summand[i, , ] <- t(D_i) %*% I_i %*% r_i %*% t(r_i) %*% t(I_i) %*% D_i

        # deriv_r_i is \partial r(\theta) / \partial \theta^T for the i-th individual
        deriv_r_i <- cbind(
            -as.numeric(exp(Zdm_i %*% alpha_hat)) * Zdm_i,
            -as.numeric(exp(-cA_i * Xdm_i %*% beta_hat)) * Y_i * cA_i * Xdm_i
        )
        stopifnot(nrow(deriv_r_i) == T_i & ncol(deriv_r_i) == (q + p))

        H_ii <- deriv_r_i %*% Mn_inv %*% t(D_i) / sample_size
        stopifnot(nrow(H_ii) == T_i & ncol(H_ii) == T_i)

        I_minus_H_i <- diag(1, nrow = T_i, ncol = T_i)
        I_minus_H_i_inv <- solve(I_minus_H_i - H_ii)

        Sigman_tilde_summand[i, , ] <- t(D_i) %*% I_i %*% I_minus_H_i_inv %*% r_i %*%
            t(r_i) %*% t(I_minus_H_i_inv) %*% t(I_i) %*% D_i
    }
    Sigman <- apply(Sigman_summand, c(2, 3), sum) / sample_size

    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q + 1):(q + p)])

    Sigman_tilde <- apply(Sigman_tilde_summand, c(2, 3), sum) / sample_size

    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q + 1):(q + p)])

    ### 4. return the result with variable names ###

    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames

    ### 5. calculate confidence interval

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

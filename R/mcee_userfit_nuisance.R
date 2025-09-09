#' Mediated Causal Excursion Effects with user-supplied nuisance predictions
#'
#' Skips Stage-1 model fitting and uses user-provided nuisance predictions.
#'
#' @inheritParams mcee
#' @param p1,q1,eta1,eta0,mu1,mu0,nu1,nu0 Numeric vectors (or column names) of
#'   per-row predictions aligned with \code{data}. See Details for definitions.
#'
#' @details
#' Nuisance definitions:
#' \itemize{
#'   \item \code{p1}: \eqn{P(A_t=1\mid H_t)} (known in MRTs). (Technically,
#'   this is \eqn{P(A_t=I_t\mid H_t)}, but the user is allowed to input \eqn{P(A_t=1\mid H_t)}
#'   and the function will automatically correct it by setting \code{p1 = 1} when \eqn{I_t = 0}.)
#'   \item \code{q1}: \eqn{P(A_t=1\mid H_t, M_t)}. (Technically,
#'   this is \eqn{P(A_t=I_t\mid H_t, M_t)}, but the user is allowed to input \eqn{P(A_t=1\mid H_t, M_t)}
#'   and the function will automatically correct it by setting \code{q1 = 1} when \eqn{I_t = 0}.)
#'   \item \code{eta1}, \code{eta0}: \eqn{E(Y\mid H_t, A_t=1)} and \eqn{E(Y\mid H_t, A_t=0)}.
#'   \item \code{mu1}, \code{mu0}: \eqn{E(Y\mid H_t, A_t=1, M_t)} and \eqn{E(Y\mid H_t, A_t=0, M_t)}.
#'   \item \code{nu1}, \code{nu0}: cross-world regressions; see vignette and paper for definitions.
#' }
#' If \code{availability} is provided, rows with \eqn{I=0} are coerced to \code{p1=q1=1}
#' (and hence \code{p0=q0=1}); a warning is emitted if overrides occur.
#'
#' @return An \code{"mcee_fit"} object; see \code{\link{mcee}}.
#' @seealso \code{\link{mcee}}, \code{\link{mcee_general}}
#' @examples
#' set.seed(1)
#' n <- 10
#' T <- 4
#' id <- rep(1:n, each = T)
#' dp <- rep(1:T, times = n)
#' A <- rbinom(n * T, 1, 0.5)
#' M <- rbinom(n * T, 1, plogis(-0.2 + 0.3 * A + 0.1 * dp))
#' Y <- ave(0.5 * A + 0.6 * M + 0.1 * dp + rnorm(n * T), id)
#' dat <- data.frame(id, dp, A, M, Y)
#'
#' fit_usr <- mcee_userfit_nuisance(dat, "id","dp","Y","A","M",
#'     time_varying_effect_form = ~ dp,
#'     p1 = rep(0.5, nrow(dat)),
#'     q1 = runif(nrow(dat),.3,.7),
#'     eta1 = rnorm(nrow(dat)), eta0 = rnorm(nrow(dat)),
#'     mu1 = rnorm(nrow(dat)),  mu0 = rnorm(nrow(dat)),
#'     nu1 = rnorm(nrow(dat)),  nu0 = rnorm(nrow(dat)))
#' @export
mcee_userfit_nuisance <- function(
    data,
    id,
    dp,
    outcome,
    treatment,
    mediator,
    availability = NULL,
    time_varying_effect_form,
    p1,
    q1,
    eta1,
    eta0,
    mu1,
    mu0,
    nu1,
    nu0,
    weight_per_row = NULL,
    verbose = TRUE) {
    cl <- match.call()

    # ---- Basic checks ----------------------------------------------------------
    .mcee_assert_df(data)
    .mcee_require_cols(data, c(id, dp, outcome, treatment, mediator))
    if (!is.null(availability)) .mcee_require_cols(data, availability)

    # Check for missingness
    vars_core <- c(id, dp, outcome, treatment, mediator, if (!is.null(availability)) availability)
    vars_mod <- .mcee_vars_in_rhs(time_varying_effect_form)
    .mcee_check_no_missing_vars(
        data,
        vars = unique(c(vars_core, vars_mod)),
        where = "mcee_userfit_nuisance inputs & time_varying_effect_form"
    )

    # Nuisance predictions
    .mcee_check_no_missing_vec(p1, "p1")
    .mcee_check_no_missing_vec(q1, "q1")
    .mcee_check_no_missing_vec(eta1, "eta1")
    .mcee_check_no_missing_vec(eta0, "eta0")
    .mcee_check_no_missing_vec(mu1, "mu1")
    .mcee_check_no_missing_vec(mu0, "mu0")
    .mcee_check_no_missing_vec(nu1, "nu1")
    .mcee_check_no_missing_vec(nu0, "nu0")

    # Check legal values of treatment, availability, outcome const within id
    .mcee_check_binary_col(data, treatment, allow_all1 = TRUE, label = "treatment")
    .mcee_check_binary_col(data, availability, allow_all1 = TRUE, label = "availability")
    .mcee_message_if_no_availability_provided(availability, verbose)
    .mcee_check_outcome_constant_within_id(data, id, outcome)

    # check for each id, the rows are grouped together, and that dp strictly increasing within id
    .mcee_check_id_rows_grouped(data, id)
    .mcee_check_dp_strictly_increasing(data, id, dp)

    # Resolve column-or-vector helper
    resolve_nv <- function(x, nm) {
        if (is.character(x) && length(x) == 1L) {
            .mcee_require_cols(data, x)
            as.numeric(data[[x]])
        } else if (is.numeric(x) && length(x) == nrow(data)) {
            as.numeric(x)
        } else {
            stop(nm, " must be a column name or numeric vector of length nrow(data).", call. = FALSE)
        }
    }

    # Pull nuisance vectors
    p1_vec <- resolve_nv(p1, "p1")
    q1_vec <- resolve_nv(q1, "q1")
    eta1_vec <- resolve_nv(eta1, "eta1")
    eta0_vec <- resolve_nv(eta0, "eta0")
    mu1_vec <- resolve_nv(mu1, "mu1")
    mu0_vec <- resolve_nv(mu0, "mu0")
    nu1_vec <- resolve_nv(nu1, "nu1")
    nu0_vec <- resolve_nv(nu0, "nu0")

    # Availability & bounds for p1, q1 (only when I==1)
    I <- if (is.null(availability)) rep(1, nrow(data)) else as.numeric(data[[availability]])
    if (any(I == 1 & (p1_vec <= 0 | p1_vec >= 1))) {
        stop("`p1` must lie in (0,1) when availability==1.", call. = FALSE)
    }
    if (any(I == 1 & (q1_vec <= 0 | q1_vec >= 1))) {
        stop("`q1` must lie in (0,1) when availability==1.", call. = FALSE)
    }

    # Enforce I==0 => p1=p0=q1=q0=1 (warn if overridden)
    p1_adj <- p1_vec
    q1_adj <- q1_vec
    if (!is.null(availability)) {
        changed <- which(I == 0 & ((p1_vec != 1) | (q1_vec != 1)))
        if (length(changed) > 0) {
            warning(sprintf(
                "For %d row(s) with availability==0, p1/q1 were overridden to 1 (and p0/q0 to 1).",
                length(changed)
            ))
        }
        p1_adj[I == 0] <- 1
        q1_adj[I == 0] <- 1
    }
    p0_adj <- 1 - p1_adj
    q0_adj <- 1 - q1_adj
    if (!is.null(availability)) {
        p0_adj[I == 0] <- 1
        q0_adj[I == 0] <- 1
    }

    # Basis and weights
    .mcee_check_time_varying_effect_form(time_varying_effect_form, dp)
    f_nrows <- .mcee_build_f_matrix(time_varying_effect_form, data)
    omega_nrows <- .mcee_build_weights(
        data, id, dp,
        weight_per_row = weight_per_row,
        specific_dp_only = NULL,
        verbose = verbose
    )

    if (isTRUE(verbose)) {
        message(
            "[mcee_userfit_nuisance] Stage-2 with p-dim = ", ncol(f_nrows),
            ", subjects = ", length(unique(data[[id]])), "."
        )
    }

    # ---- Stage-2 estimation ----------------------------------------------------
    core <- mcee_helper_stage2_estimate_mcee(
        data = data,
        id_var = id,
        dp_var = dp,
        outcome_var = outcome,
        treatment_var = treatment,
        avail_var = availability,
        p1 = p1_adj, p0 = p0_adj,
        q1 = q1_adj, q0 = q0_adj,
        eta1 = eta1_vec, eta0 = eta0_vec,
        mu1 = mu1_vec, mu0 = mu0_vec,
        nu1 = nu1_vec, nu0 = nu0_vec,
        omega_nrows = omega_nrows,
        f_nrows = f_nrows
    )

    # ---- Assemble return -------------------------------------------------------
    T_per_id <- as.integer(table(data[[id]]))
    out <- list(
        call = cl,
        mcee_fit = core,
        nuisance_models = "Fitted values for all nuisance functions were supplied by the user.", # no Stage-1 fits returned in this wrapper
        nuisance_fitted = list(
            p1 = p1_adj,
            p0 = p0_adj,
            q1 = q1_adj,
            q0 = q0_adj,
            eta1 = eta1_vec,
            eta0 = eta0_vec,
            mu1 = mu1_vec,
            mu0 = mu0_vec,
            nu1 = nu1_vec,
            nu0 = nu0_vec
        ),
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

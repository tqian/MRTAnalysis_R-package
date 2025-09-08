# code to prepare data_time_varying_mediator_distal_outcome goes here

dgm_time_varying_mediator_distal_outcome <- function(
        n,
        T_val,
        intervention_time = NULL,
        set_At_for_Y      = NULL,
        set_At_for_M      = NULL
        # intervention_time, set_At_for_Y, set_At_for_M were arguments
        # from the simulation study, where we set A_t to a fixed value at a certain time point
        # separately for Y and for M_t. They are not used in the R package.
) {
    expit <- function(x) 1 / (1 + exp(-x))

    # ---- Parameters & helper functions ----
    rho_XX <- 0.3; rho_XA <- 0.2; rho_XM <- 0.2; sigma_X <- 1
    rho_Iint <- 1.5; rho_IA <- -0.3; rho_IM <- -0.3; rho_IX <- 0.3
    rho_Aint <- 0; rho_AA <- 0.2; rho_AM <- 0.2; rho_AX <- 0.3
    rho_MM <- 0.4; rho_MX <- 0.3; rho_MA <- 0.6; rho_MAlag1 <- 0.4; sigma_M <- 1
    # rho_XX <- 0; rho_XA <- 0; rho_XM <- 0; sigma_X <- 1
    # rho_Iint <- 99; rho_IA <- 0; rho_IM <- 0; rho_IX <- 0
    # rho_Aint <- 0; rho_AA <- 0; rho_AM <- 0; rho_AX <- 0
    # rho_MM <- 0; rho_MX <- 0; rho_MA <- 10; rho_MAlag1 <- 0; sigma_M <- 1

    rho_YX  <- rep(0.3, T_val)
    rho_YM  <- rep(0.4, T_val)
    rho_YA  <- rep(0.2, T_val)
    rho_YAM <- rep(0.1, T_val)
    sigma_Y <- 1

    # New helper functions for A and M
    h1 <- function(t_val, z) {
        t_transformed <- 3 * (2 * t_val - T_val) / T_val # [-3, 3]
        tanh(t_transformed) + sin(z)
        # z_transformed <- plogis(z)
        # 0.5 * dbeta(t_transformed, shape1 = 2, shape2 = 5) + 0.5 * dbeta(z_transformed, shape1 = 2, shape2 = 5)
    }
    h2 <- function(t_val, z) {
        t_transformed <- 3 * (2 * t_val - T_val) / T_val # [-3, 3]
        tanh(t_transformed) + sin(z)
        # z_transformed <- plogis(z)
        # 0.5 * dbeta(t_transformed, shape1 = 5, shape2 = 2) + 0.5 * dbeta(z_transformed, shape1 = 5, shape2 = 2)
    }

    # h1 <- function(t_val, z) z
    # h2 <- function(t_val, z) z
    h3 <- h2
    h4 <- h1

    # Outcome g-functions
    g1 <- h1
    g2 <- h2
    g3 <- h1

    # ---- Storage ----
    X_mat     <- matrix(NA_real_, n, T_val)
    I_mat     <- matrix(NA_integer_, n, T_val)
    A_mat     <- matrix(NA_integer_, n, T_val)
    M_mat     <- matrix(NA_real_, n, T_val)

    Xprev_mat <- matrix(NA_real_, n, T_val)
    Iprev_mat <- matrix(NA_integer_, n, T_val)
    Aprev_mat <- matrix(NA_integer_, n, T_val)
    Mprev_mat <- matrix(NA_real_, n, T_val)

    muX_mat   <- matrix(NA_real_, n, T_val)
    pI_mat    <- matrix(NA_real_, n, T_val)
    pA_mat    <- matrix(NA_real_, n, T_val)
    muM_mat   <- matrix(NA_real_, n, T_val)

    mu_Y      <- numeric(n)

    # ---- Initial previous values ----
    X_prev <- rnorm(n)
    A_prev <- integer(n)
    M_prev <- numeric(n)
    I_prev <- integer(n)

    # ---- Time loop ----
    for (t_val in seq_len(T_val)) {
        # Store previous values
        Xprev_mat[, t_val] <- X_prev
        Iprev_mat[, t_val] <- I_prev
        Aprev_mat[, t_val] <- A_prev
        Mprev_mat[, t_val] <- M_prev

        # X_t
        mX   <- rho_XX * X_prev + rho_XA * A_prev + rho_XM * M_prev
        mX    <- round(mX, 2)
        X_t  <- rnorm(n, mean = mX, sd = sigma_X)
        X_t   <- round(X_t, 2)
        X_mat[, t_val]   <- X_t
        muX_mat[, t_val] <- mX

        # I_t
        pI   <- expit(rho_Iint + rho_IA * A_prev + rho_IM * M_prev + rho_IX * X_t)
        pI    <- round(pI, 2)
        I_t  <- rbinom(n, 1, pI)
        I_mat[, t_val]   <- I_t
        pI_mat[, t_val]  <- pI

        # Determine A_t for M and for Y
        # Natural A
        pA   <- expit(rho_Aint + rho_AA * A_prev + rho_AM * h1(t_val, M_prev) + rho_AX * h2(t_val, X_t))
        pA    <- round(pA, 2)
        A_nat <- ifelse(I_t == 1, rbinom(n, 1, pA), 0L)
        # Intervention override for Y
        A_Y   <- if (!is.null(intervention_time) && t_val == intervention_time && !is.null(set_At_for_Y)) {
            ifelse(I_t == 1, set_At_for_Y, 0L)
        } else {
            A_nat
        }
        # Intervention override for M
        A_M   <- if (!is.null(intervention_time) && t_val == intervention_time && !is.null(set_At_for_M)) {
            ifelse(I_t == 1, set_At_for_M, 0L)
        } else {
            A_Y
        }

        A_mat[, t_val]   <- A_Y
        pA_mat[, t_val]  <- pA

        # M_t via h3, h4
        mM   <- rho_MM * h3(t_val, M_prev) +
            rho_MX * h4(t_val, X_t) +
            rho_MA * A_M + rho_MAlag1 * A_prev
        mM    <- round(mM, 2)
        M_t  <- rnorm(n, mean = mM, sd = sigma_M)
        M_t   <- round(M_t, 2)
        M_mat[, t_val]   <- M_t
        muM_mat[, t_val] <- mM

        # Accumulate pre-noise outcome
        mu_Y <- mu_Y +
            rho_YX[t_val]  * g1(t_val, X_t) +
            rho_YM[t_val]  * g2(t_val, M_t) +
            rho_YA[t_val]  * A_Y +
            rho_YAM[t_val] * A_Y * g3(t_val, M_t)
        mu_Y  <- round(mu_Y, 2)

        # Update previous for next iteration
        X_prev <- X_t;  A_prev <- A_Y;  M_prev <- M_t;  I_prev <- I_t
    }

    # ---- Final outcome + noise ----
    Y <- mu_Y + rnorm(n, mean = 0, sd = sigma_Y)
    Y <- round(Y, 2)

    # ---- Build long-format data.frame ----
    data.frame(
        id      = rep(seq_len(n),       each = T_val),
        dp   = rep(seq_len(T_val),   times = n),

        X_prev  = c(t(Xprev_mat)),
        I_prev  = c(t(Iprev_mat)),
        A_prev  = c(t(Aprev_mat)),
        M_prev  = c(t(Mprev_mat)),

        X       = c(t(X_mat)),
        I       = c(t(I_mat)),
        A       = c(t(A_mat)),
        M       = c(t(M_mat)),

        mu_X    = c(t(muX_mat)),
        p_I     = c(t(pI_mat)),
        p_A     = c(t(pA_mat)),
        mu_M    = c(t(muM_mat)),

        mu_Y    = rep(mu_Y, each = T_val),
        Y       = rep(Y,    each = T_val)
    )

}

set.seed(123)
data_time_varying_mediator_distal_outcome <- dgm_time_varying_mediator_distal_outcome(
  n = 50,
  T_val = 10
)
usethis::use_data(data_time_varying_mediator_distal_outcome, overwrite = TRUE)

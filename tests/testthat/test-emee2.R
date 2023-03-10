# # Generate Data -----------------------------------------------------------
# set.seed(123)
#
# dgm_demo <- function(sample_size, total_T) {
#     alpha_0 <- -1.5
#     alpha_1 <- 0.5
#     alpha_2 <- 0.3
#
#     beta_0 <- 0.1
#     beta_1 <- 0.3
#
#     # With the above parameter specification, the range of success probability of Y would be
#     # [exp(-1.5), exp(-1.5 + 0.5 + 0.3 + 0.1 + 0.3)] = [0.223, 0.741]
#
#     df_names <- c("userid", "time", "time_var1", "time_var2", "Y", "A", "avail", "prob_Y", "prob_Y_A0", "rand_prob")
#     # time_var1 is time / total_T
#     # time_var2 is 1(time > total_T/2)
#
#     dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
#     names(dta) <- df_names
#
#     dta$userid <- rep(1:sample_size, each = total_T)
#     dta$time <- rep(1:total_T, times = sample_size)
#     dta$time_var1 <- dta$time / total_T
#     dta$time_var2 <- as.numeric(dta$time > (total_T / 2))
#
#     for (t in 1:total_T) {
#         # row index for the rows corresponding to time t for every subject
#         row_index <- seq(from = t, by = total_T, length = sample_size)
#
#         dta$avail[row_index] <- rbinom(sample_size, 1, 0.8) # 0.8 probability to be available
#
#         dta$rand_prob[row_index] <- ifelse(t %% 3 == 1, 0.3, ifelse(t %% 3 == 2, 0.5, 0.7))
#         dta$A[row_index] <- rbinom(sample_size, 1, dta$rand_prob[row_index]) * dta$avail[row_index] # A can only be 1 if avail = 1
#         # We keep rand_prob as-is for those observations with avail = 0. It's OK because those rand_prob won't be used in the estimation.
#
#         dta$prob_Y_A0[row_index] <- exp(alpha_0 + alpha_1 * dta$time_var1[row_index] + alpha_2 * dta$time_var2[row_index])
#         dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * (beta_0 + beta_1 * dta$time_var1[row_index]))
#         dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
#     }
#
#     return(dta)
# }
#
# data_binary <- dgm_demo(sample_size = 100, total_T = 30)


# Tests -------------------------------------------------------------------

test_that(
    "check beta_hat",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$beta_hat),
            as.vector(c(0.08502509, 0.41260611)),
            tolerance = 1e-7
        )
    }
)

test_that(
    "check alpha_hat",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$alpha_hat),
            as.vector(c(-1.4358563, 0.6935913, 0.2647090)),
            tolerance = 1e-6
        )
    }
)

test_that(
    "check beta_se",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$beta_se),
            as.vector(c(0.1250711, 0.1810747)),
            tolerance = 1e-6
        )
    }
)

test_that(
    "check alpha_se",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$alpha_se),
            as.vector(c(0.06213613, 0.17279523, 0.10357601)),
            tolerance = 1e-6
        )
    }
)

test_that(
    "check beta_se_adjusted",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$beta_se_adjusted),
            as.vector(c(0.1266788, 0.1834912)),
            tolerance = 1e-6
        )
    }
)

test_that(
    "check alpha_se_adjusted",
    {
        expect_equal(
            as.numeric(emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$alpha_se_adjusted),
            as.vector(c(0.06290973, 0.17500825, 0.10478227)),
            tolerance = 1e-6
        )
    }
)

test_that(
    "check varcov",
    {
        expect_equal(
            emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$varcov,
            matrix(
                c(
                    0.0038608991, -0.007512966, 0.0020372384, -0.0005068623, 0.001975262,
                    -0.0075129657, 0.029858190, -0.0156664102, 0.0026684541, -0.008951825,
                    0.0020372384, -0.015666410, 0.0107279889, -0.0008563799, 0.003400086,
                    -0.0005068623, 0.002668454, -0.0008563799, 0.0156427857, -0.020762192,
                    0.0019752617, -0.008951825, 0.0034000860, -0.0207621925, 0.032788043
                ),
                nrow = 5, ncol = 5, byrow = TRUE
            ),
            tolerance = 1e-7
        )
    }
)

test_that(
    "check varcov_adjusted",
    {
        expect_equal(
            emee2(
                outcome = Y,
                treatment = A,
                moderator_formula = ~time_var1,
                control_formula = ~ time_var1 + time_var2,
                availability = avail,
                rand_prob = rand_prob,
                numerator_prob = 0.5,
                data = data_binary,
                id = userid,
                start = NULL,
                verbose = FALSE
            )$fit$varcov_adjusted,
            matrix(
                c(
                    0.0039576339, -0.0077190622, 0.0021002006, -0.0005141107, 0.0020213942,
                    -0.0077190622, 0.0306278873, -0.0160546827, 0.0027628732, -0.0092120177,
                    0.0021002006, -0.0160546827, 0.0109793233, -0.0009021267, 0.0035113704,
                    -0.0005141107, 0.0027628732, -0.0009021267, 0.0160475226, -0.0213114237,
                    0.0020213942, -0.0092120177, 0.0035113704, -0.0213114237, 0.0336690206
                ),
                nrow = 5, ncol = 5, byrow = TRUE
            ),
            tolerance = 1e-7
        )
    }
)

# extract the confidence interval from the output and drop its column and row names

conf_int <- emee2(
    outcome = Y,
    treatment = A,
    moderator_formula = ~time_var1,
    control_formula = ~ time_var1 + time_var2,
    availability = avail,
    rand_prob = rand_prob,
    numerator_prob = 0.5,
    data = data_binary,
    id = userid,
    start = NULL,
    verbose = FALSE
)$fit$conf_int
dimnames(conf_int) <- c()

test_that(
    "check conf_int",
    {
        expect_equal(conf_int,
            matrix(c(-0.16011431, 0.3301645, 0.05769972, 0.7675125),
                nrow = 2, ncol = 2, byrow = TRUE
            ),
            tolerance = 1e-6
        )
    }
)

# extract the adjusted confidence interval from the output and drop its column and row names

conf_int_adjusted <- emee2(
    outcome = Y,
    treatment = A,
    moderator_formula = ~time_var1,
    control_formula = ~ time_var1 + time_var2,
    availability = avail,
    rand_prob = rand_prob,
    numerator_prob = 0.5,
    data = data_binary,
    id = userid,
    start = NULL,
    verbose = FALSE
)$fit$conf_int_adjusted
dimnames(conf_int_adjusted) <- c()

test_that(
    "check conf_int_adjusted",
    {
        expect_equal(conf_int_adjusted,
            matrix(c(-0.16646416, 0.3365143, 0.04833002, 0.7768822),
                nrow = 2, ncol = 2, byrow = TRUE
            ),
            tolerance = 1e-6
        )
    }
)

# test_that(
#     "check error when moderator_formula has Y in formula",
#     {
#         expect_error(
#             emee2(
#                 outcome = Y,
#                 treatment = A,
#                 moderator_formula = Y ~ time_var1,
#                 control_formula = ~ time_var1 + time_var2,
#                 availability = avail,
#                 rand_prob = rand_prob,
#                 numerator_prob = 0.5,
#                 data = data_binary,
#                 id = userid,
#                 start = NULL,
verbose <- FALSE
#             ),
#             "It seems like you included variables to left of ~. moderator_formula should look like ~1 or ~ mod_var1 + mod_var2."
#         )
#     }
# )

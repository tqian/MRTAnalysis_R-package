# create equations of alphas and betas
f1 <- function(x) {
    c(
        F1 = x[1] + x[2] - 1,
        F2 = x[1] - x[2] - 7
    )
}
f2 <- function(y) {
    c(
        F1 = y[1] + y[2] - 1,
        F2 = y[1] - y[2] - 7,
        F3 = y[3] + y[4] - 1,
        F4 = y[3] - y[4] - 5,
        F5 = y[5] + y[6] - 1,
        F6 = y[5] - y[6] - 3
    )
}
# generate solutions using the multiroot function
solution1 <- rootSolve::multiroot(f = f1, start = c(1, 1))
solution2 <- rootSolve::multiroot(f = f2, start = rep(1, 6))

# tests
test_that("check alpha when p = 1 and q = 1", {
    expect_equal(
        get_alpha_beta_from_multiroot_result(
            root = solution1,
            p = 1,
            q = 1
        )$alpha,
        4
    )
})

test_that("check beta when p = 1 and q = 1", {
    expect_equal(
        get_alpha_beta_from_multiroot_result(
            root = solution1,
            p = 1,
            q = 1
        )$beta,
        -3
    )
})
test_that("check beta when p > 1 and q > 1", {
    expect_equal(
        get_alpha_beta_from_multiroot_result(
            root = solution2,
            p = 3,
            q = 3
        )$alpha,
        as.matrix(c(4, -3, 3))
    )
})

test_that("check beta when p > 1 and q > 1", {
    expect_equal(
        get_alpha_beta_from_multiroot_result(
            root = solution2,
            p = 3,
            q = 3
        )$beta,
        as.matrix(c(-2, 2, -1))
    )
})

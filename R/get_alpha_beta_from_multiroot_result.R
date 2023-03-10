#' Extract roots of alpha and beta from output of
#' the multiroot function
#'
#' @param root a list of function output from the multiroot function
#' @param p a numeric value that specifies number of variables in beta
#' @param q a numeric value that specifies number of variables in alpha
#'
#' @return return a vector of roots for alpha,
#'         and a vector of roots for beta
#' @noRd
#' @examples f <- function(x) {
#'     c(
#'         F1 = x[1] + x[2] - 1,
#'         F2 = x[1] - x[2] - 7
#'     )
#' }
#' solution <- rootSolve::multiroot(f = f, start = c(1, 1))
#' get_alpha_beta_from_multiroot_result(
#'     root = solution,
#'     p = 1,
#'     q = 1
#' )
get_alpha_beta_from_multiroot_result <- function(root, p, q) {
    if (p == 1) {
        beta_root <- root$root[q + 1]
    } else {
        beta_root <- as.matrix(root$root[(q + 1):(q + p)])
    }
    if (q == 1) {
        alpha_root <- root$root[1]
    } else {
        alpha_root <- as.matrix(root$root[1:q])
    }
    return(list(alpha = alpha_root, beta = beta_root))
}

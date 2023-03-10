#' Find the locations at which the value changes from the previous one
#'
#' @param v a vector that specifies different values
#'
#' @return Return a vector that specifies the indexes of locations
#'         at which the value changes from the previous one
#'
#' @noRd
#' @examples v <- c("a", "a", "b", "c", "c")
#' find_change_location(v)
find_change_location <- function(v) {
    n <- length(v)
    if (n <= 1) {
        stop("The vector need to have length > 1.")
    }
    return(c(1, 1 + which(v[1:(n - 1)] != v[2:n])))
}

#' Preprocess the input to `wcls()`, `emee()`, or `emee2()`
#'
#' @param data A data set in long format.
#' @param mf An object constructed using `match.call()`.
#' @param outcome_type If `"binary"`, an additional check will be performed to
#' ensure the outcome only takes value 0 or 1.
#' @param verbose If default (`TRUE`), additional messages will be printed
#' during data preprocessing.
#'
#' @return Returns a list of preprocessed variables: id, outcome, treatment,
#' moderator_matrix, control_matrix,availability, numerator_prob.
#'
#' @noRd

preprocess_input <- function(data, mf, outcome_type, verbose = TRUE) {
    if ("tibble" %in% class(data)) {
        warning("The function may not work for data in tibble format. Use as.data.frame() to convert into data frame.")
    }

    # creating id, outcome, treatment
    id <- data[, as.character(mf$id)]
    outcome <- as.numeric(data[, as.character(mf$outcome)])
    treatment <- as.numeric(data[, as.character(mf$treatment)])

    # creating availability
    if (is.null(mf$availability)) {
        availability <- rep(1, nrow(data))
        if (verbose) {
            message("availability = NULL: defaulting availability to always available.")
        }
    } else {
        availability <- as.numeric(data[, as.character(mf$availability)])
    }

    # creating rand_prob (randomization probability)
    if (is.numeric(mf$rand_prob)) {
        # user-input constant randomization probability
        rand_prob <- rep(mf$rand_prob, nrow(data))
        if (verbose) {
            message(paste0("Constant randomization probability ", mf$rand_prob, " is used."))
        }
    } else {
        rand_prob <- as.numeric(data[, as.character(mf$rand_prob)])
    }

    # creating and checking moderator variables
    moderator_formula <- as.formula(mf$moderator_formula)
    parse.fm_mdr <- as.character(mf$moderator_formula)
    if (length(parse.fm_mdr) == 2) {
        # the formula looks like "~x"
    } else if (length(parse.fm_mdr) == 3) {
        # the formula looks like "y ~ x"
        stop("It seems like you included variables to left of ~. moderator_formula should look like ~1 or ~ mod_var1 + mod_var2.")
    } else {
        stop("Unknown moderator_formula pattern! moderator_formula should look like ~1 or ~ mod_var1 + mod_var2.")
    }
    moderator_matrix <- model.matrix(as.formula(mf$moderator_formula), data = data)

    # creating and checking control variables
    control_formula <- as.formula(mf$control_formula)
    parse.fm_ctl <- as.character(mf$control_formula)
    if (length(parse.fm_ctl) == 2) {
        # the formula looks like "~x"
    } else if (length(parse.fm_ctl) == 3) {
        # the formula looks like "y ~ x"
        stop("It seems like you included variables to left of ~. control_formula should look like ~1 or ~ ctrl_var1 + ctrl_var2.")
    } else {
        stop("Unknown control_formula pattern! control_formula should look like ~1 or ~ ctrl_var1 + ctrl_var2.")
    }
    control_matrix <- model.matrix(as.formula(mf$control_formula), data = data)

    # creating and checking numerator_prob (numerator probability)
    if (is.null(mf$numerator_prob)) {
        if (is.numeric(mf$rand_prob)) {
            numerator_prob <- rep(mf$rand_prob, nrow(data))
        } else if (length(unique(rand_prob)) == 1) {
            numerator_prob <- rep(rand_prob[1], nrow(data))
        } else {
            numerator_prob <- rep(0.5, nrow(data))
        }
    } else if (is.numeric(mf$numerator_prob)) {
        numerator_prob <- rep(mf$numerator_prob, nrow(data))
        if (verbose) {
            message(paste0("Constant numerator probability ", mf$numerator_prob, " is used."))
        }
    } else {
        numerator_prob <- as.numeric(data[, as.character(mf$numerator_prob)])
        if (length(unique(numerator_prob)) > 1) {
            if (verbose) {
                message(paste0(
                    "Non-constant numerator probability is used. ",
                    "Make sure you know how to appropriately choose a non-constant numerator probability ",
                    "(that it can only depend on the moderators). ",
                    "Otherwise the estimated causal effect can be biased."
                ))
            }
        }
    }

    # checking variable types
    stopifnot(is.numeric(outcome))
    stopifnot(is.numeric(treatment))
    stopifnot(is.numeric(rand_prob))
    stopifnot(is.numeric(moderator_matrix))
    stopifnot(is.numeric(control_matrix))
    stopifnot(is.numeric(availability))
    stopifnot(is.numeric(numerator_prob))

    # more checking for availability
    if (any(is.na(availability))) {
        stop("NA in availability variable. This package is unable to handle this.")
    }
    if (!all(availability %in% 0:1)) {
        stop("availability variable should only contain 0 or 1.")
    }

    # checking NA in other variables among availability time points
    index_avail <- which(availability == 1)
    if (any(is.na(outcome[index_avail]))) {
        stop("NA in outcome. This package is unable to handle this.")
    }
    if (any(is.na(treatment[index_avail]))) {
        stop("NA in treatment. This package is unable to handle this.")
    }
    if (any(is.na(rand_prob[index_avail]))) {
        stop("NA in rand_prob. This package is unable to handle this.")
    }
    if (any(is.na(moderator_matrix[index_avail, ]))) {
        stop("NA in moderator variables. This package is unable to handle this.")
    }
    if (any(is.na(control_matrix[index_avail, ]))) {
        stop("NA in control variables. This package is unable to handle this.")
    }
    if (any(is.na(numerator_prob[index_avail]))) {
        stop("NA in numerator_prob. This package is unable to handle this.")
    }

    # checking variables that should be binary among availability time points
    if (!all(treatment[index_avail] %in% 0:1)) {
        stop("treatment variable should only contain 0 or 1.")
    }
    if (any(rand_prob[index_avail] <= 0) | any(rand_prob[index_avail] >= 1)) {
        stop("rand_prob should be greater than 0 and less than 1.")
    }
    if (any(numerator_prob[index_avail] <= 0) | any(numerator_prob[index_avail] >= 1)) {
        stop("rand_prob should be greater than 0 and less than 1.")
    }

    if (outcome_type == "binary") {
        if (!all(outcome[index_avail] %in% 0:1)) {
            stop("outcome variable should only contain 0 or 1.")
        }
    }

    list(
        id = id,
        outcome = outcome,
        treatment = treatment,
        rand_prob = rand_prob,
        moderator_matrix = moderator_matrix,
        control_matrix = control_matrix,
        availability = availability,
        numerator_prob = numerator_prob
    )
}

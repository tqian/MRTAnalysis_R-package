# -----------------------------------
# Internal utilities for MCEE wrappers
# -----------------------------------

.mcee_assert_df <- function(data) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  invisible(TRUE)
}

.mcee_require_cols <- function(data, cols, where = "data") {
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Missing columns in `", where, "`: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Generic binary-column checker (optionally forbid all-1)
.mcee_check_binary_col <- function(data, col, allow_all1 = TRUE, label = NULL) {
  if (is.null(col)) return(invisible(TRUE))
  vals <- unique(data[[col]])
  if (!all(vals %in% c(0, 1))) {
    nm <- if (is.null(label)) col else paste0(label, " ('", col, "')")
    stop(nm, " must be coded 0/1.", call. = FALSE)
  }
  if (!allow_all1 && setequal(vals, 1)) {
    nm <- if (is.null(label)) col else paste0(label, " ('", col, "')")
    stop(nm, " cannot be all ones.", call. = FALSE)
  }
  invisible(TRUE)
}

.mcee_check_outcome_constant_within_id <- function(data, id, outcome) {
  k <- tapply(data[[outcome]], data[[id]], function(x) length(unique(x)))
  if (any(k != 1)) {
    stop("`outcome` column '", outcome, "' must be constant within each subject (id='", id, "').",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Check that rows for each id appear in a single contiguous block
.mcee_check_id_rows_grouped <- function(data, id, max_show = 5) {
  id_vec <- data[[id]]
  # contiguity check: each id value should appear in exactly one run
  runs <- rle(as.character(id_vec))$values
  dup_runs <- names(which(table(runs) > 1L))

  if (length(dup_runs) > 0L) {
    offenders <- paste(utils::head(dup_runs, max_show), collapse = ", ")
    more <- if (length(dup_runs) > max_show) ", ..." else ""
    stop("Rows for id column '", id,
         "' must be contiguous (all rows for a given id must appear together). ",
         "Offending id(s): ", offenders, more, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.mcee_check_dp_strictly_increasing <- function(data, id, dp) {
  bad <- tapply(data[[dp]], data[[id]], function(v) any(diff(as.numeric(v)) <= 0, na.rm = TRUE))
  if (any(bad)) {
    offenders <- names(bad)[bad]
    stop("`dp` column '", dp, "' must be strictly increasing within each subject. Offending id(s): ",
         paste(head(offenders, 5), collapse = ", "),
         if (length(offenders) > 5) ", ...",
         call. = FALSE)
  }
  invisible(TRUE)
}

.mcee_message_if_no_availability_provided <- function(availability, verbose) {
  if (is.null(availability) && isTRUE(verbose)) {
    message("`availability` not provided; assuming all rows available.")
  }
  invisible(TRUE)
}

.mcee_check_time_varying_effect_form <- function(time_varying_effect_form, dp) {
  if (!inherits(time_varying_effect_form, "formula") || length(time_varying_effect_form) != 2L) {
    stop("`time_varying_effect_form` must be RHS-only (e.g., ~ 1, ~ poly(", dp, ", 2)).", call. = FALSE)
  }
  vars <- all.vars(time_varying_effect_form)
  extra <- setdiff(vars, dp)
  if (length(extra)) {
    warning("`time_varying_effect_form` includes variables beyond '", dp, "': ",
            paste(extra, collapse = ", "),
            ". Only functions of the decision point are intended; ",
            "precomputed basis columns are allowed.")
  }
  invisible(TRUE)
}

.mcee_build_f_matrix <- function(time_varying_effect_form, data) {
  stats::model.matrix(time_varying_effect_form, data = data)
}

.mcee_check_control_formula <- function(control_formula, treatment, outcome, dp, label) {
  if (!inherits(control_formula, "formula") || length(control_formula) != 2L) {
    stop("`", label, "` must be RHS-only (e.g., ~ X1 + X2 + splines::ns(", dp, ", 3)).",
         call. = FALSE)
  }
  vars <- all.vars(control_formula)
  if (treatment %in% vars) {
    stop("`", label, "` must not include the treatment variable '", treatment, "'.", call. = FALSE)
  }
  if (outcome %in% vars) {
    stop("`", label, "` must not include the outcome variable '", outcome, "'.", call. = FALSE)
  }
  invisible(TRUE)
}

# Drop one variable from RHS-only formula
.mcee_drop_var_from_rhs <- function(rhs_only_formula, var) {
  tt <- stats::terms(rhs_only_formula)
  labs <- attr(tt, "term.labels")
  labs <- setdiff(labs, var)
  if (length(labs) == 0) stats::as.formula("~ 1") else
    stats::as.formula(paste("~", paste(labs, collapse = " + ")))
}

# Build row weights ω(i,t)
.mcee_build_weights <- function(data, id, dp,
                                weight_per_row = NULL,
                                specific_dp_only = NULL,
                                verbose = TRUE) {
  if (!is.null(specific_dp_only)) {
    if (!is.numeric(specific_dp_only)) {
      stop("`specific_dp_only` must be numeric dp value(s).", call. = FALSE)
    }
    dps <- unique(data[[dp]])
    if (!all(specific_dp_only %in% dps)) {
      miss <- setdiff(specific_dp_only, dps)
      stop("All `specific_dp_only` values must appear in '", dp, "'. Missing: ",
           paste(miss, collapse = ", "),
           call. = FALSE)
    }
    w <- as.numeric(data[[dp]] %in% specific_dp_only)
  } else if (is.null(weight_per_row)) {
    if (isTRUE(verbose)) {
      message("`weight_per_row` not provided; using uniform weights (all ones).")
    }
    w <- rep(1, nrow(data))
  } else {
    w <- as.numeric(weight_per_row)
  }

  if (length(w) != nrow(data) || any(!is.finite(w)) || any(w < 0)) {
    stop("`weight_per_row` must be a nonnegative numeric vector of length nrow(data).", call. = FALSE)
  }
  if (all(w == 0)) stop("`weight_per_row` cannot be all zeros.", call. = FALSE)

  w
}

# Resolve known randomization prob as column or scalar; validate (0,1) when available
.mcee_resolve_rand_prob <- function(data, rand_prob, availability = NULL) {
  if (is.character(rand_prob)) {
    .mcee_require_cols(data, rand_prob)
    p_vec <- as.numeric(data[[rand_prob]])
  } else if (is.numeric(rand_prob) && length(rand_prob) == 1L) {
    p_vec <- rep_len(as.numeric(rand_prob), nrow(data))
  } else {
    stop("`rand_prob` must be a column name or a numeric scalar in (0,1).", call. = FALSE)
  }

  if (!is.null(availability)) {
    I <- as.numeric(data[[availability]])
    ok <- (I == 0) | (is.finite(p_vec) & p_vec > 0 & p_vec < 1)
    if (!all(ok)) stop("`rand_prob` must be in (0,1) where availability==1.", call. = FALSE)
  } else {
    if (any(!is.finite(p_vec) | p_vec <= 0 | p_vec >= 1))
      stop("`rand_prob` must be in (0,1).", call. = FALSE)
  }

  p_vec
}


# -------------------------------------------------------------------
# Detect missing / NaN / Inf in data columns
# -------------------------------------------------------------------
.mcee_check_no_missing_vars <- function(data, vars, where = NULL, max_show = 5) {
  stopifnot(is.data.frame(data))
  vars <- unique(vars)
  vars <- vars[vars %in% names(data)]  # silently ignore unknown names

  offenders <- list()

  for (v in vars) {
    x <- data[[v]]
    bad <- if (is.numeric(x)) {
      which(is.na(x) | !is.finite(x))  # catches NA, NaN, Inf, -Inf
    } else {
      which(is.na(x))                  # for non-numeric: NA only
    }
    if (length(bad)) offenders[[v]] <- utils::head(bad, max_show)
  }

  if (length(offenders)) {
    pieces <- vapply(names(offenders), function(v) {
      paste0(v, " at rows ", paste(offenders[[v]], collapse = ", "))
    }, character(1))
    ctx <- if (is.null(where)) "" else paste0(" (", where, ")")
    stop(
      "Missing/NaN/Inf detected in the following variable(s)", ctx, ":\n  - ",
      paste(pieces, collapse = "\n  - "),
      "\nThe software currently does not support handling missing data. ",
      "Please remove or impute missing values before calling this function.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------
# Detect missing / NaN / Inf in a supplied numeric vector (nuisance preds)
# -------------------------------------------------------------------
.mcee_check_no_missing_vec <- function(vec, name, max_show = 5) {
  bad <- which(is.na(vec) | !is.finite(vec))  # catches NA, NaN, Inf, -Inf
  if (length(bad)) {
    stop(
      "Missing/NaN/Inf detected in '", name, "' at rows ",
      paste(utils::head(bad, max_show), collapse = ", "), ". ",
      "The software currently does not support handling missing data. ",
      "Consider imputing missing values before calling this function.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# -------------------------------------------------------------------
# Extract RHS variable names from a RHS-only formula
# -------------------------------------------------------------------
.mcee_vars_in_rhs <- function(rhs_only_formula) {
  if (inherits(rhs_only_formula, "formula") && length(rhs_only_formula) == 2L) {
    all.vars(rhs_only_formula)
  } else {
    character(0)
  }
}

# -------------------------------------------------------------------
# Collect RHS variables from a single config (if it has a formula)
#   config: list(method=..., formula=~ ..., ...)
# -------------------------------------------------------------------
.mcee_vars_in_config <- function(cfg) {
  if (is.list(cfg) && !is.null(cfg$formula)) .mcee_vars_in_rhs(cfg$formula) else character(0)
}

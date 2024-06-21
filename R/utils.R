`%||%` <- function(x, y) if (is.null(x)) y else x

is_scalar <- function(x) {
    length(x) == 1L
}

is_scalar_numeric <- function(x) {
    is_scalar(x) && is.numeric(x)
}

is_number <- function(x) is_scalar_numeric(x) && !is.na(x)

pkg_nm <- function() utils::packageName(topenv(environment()))

# mimic rlang::exprs, but will omit `{`
exprs <- function(...) {
    dots <- rlang::exprs(...)
    if (...length() == 1L && identical(dots[[1L]][[1L]], as.symbol("{"))) {
        dots <- as.list(dots[[1L]][-1L])
    }
    dots
}

fclass <- function(x) class(x)[1L]

on_exit <- function(expr, add = FALSE, after = TRUE, envir = parent.frame()) {
    expr <- rlang::enquo(expr)
    expr <- rlang::expr(on.exit(expr = !!expr, add = !!add, after = !!after))
    rlang::eval_tidy(expr, env = envir)
}

set_seed <- function(seed = NULL, add = FALSE, envir = parent.frame()) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        oseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        oseed <- NULL
    }
    on_exit(restore_rng(oseed), add = add, envir = envir)
    seed <- seed %||% random_seed(1L)
    set.seed(seed)
}

restore_rng <- function(oseed) {
    if (is.null(oseed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
    } else {
        assign(".Random.seed", oseed, envir = .GlobalEnv, inherits = FALSE)
    }
}

random_seed <- function(n) sample.int(1e6L, n)

use_names_to_integer_indices <- function(use, names,
                                         arg = rlang::caller_arg(use),
                                         call = rlang::caller_env()) {
    force(arg)
    if (anyNA(use)) {
        rlang::abort(
            sprintf("%s cannot contain `NA`", style_arg(arg)),
            call = call
        )
    }
    if (isTRUE(use)) {
        use <- seq_along(names)
    } else if (isFALSE(use)) {
        use <- integer(0L)
    } else if (is.character(use)) {
        index <- match(use, names)
        if (anyNA(index)) {
            rlang::abort(sprintf(
                "%s contains invalid values (%s)",
                style_arg(arg), style_val(use[is.na(index)])
            ), call = call)
        }
        use <- index
    } else if (is.numeric(use)) {
        use <- as.integer(use)
        if (any(use < 1L) || any(use > length(names))) {
            rlang::abort(sprintf(
                "%s contains out-of-bounds indices", style_arg(arg)
            ), call = call)
        }
    } else {
        rlang::abort(
            sprintf(
                "%s must be a bool or an atomic numeric/character",
                style_arg(arg)
            ),
            call = call
        )
    }
    use
}

get_attrs <- function(data, attrs,
                      arg = rlang::caller_arg(attrs),
                      call = rlang::caller_call()) {
    force(arg)
    if (anyNA(attrs)) {
        rlang::abort(
            sprintf("%s cannot contain `NA`", style_arg(arg)),
            call = call
        )
    }
    if (is.character(attrs)) {
        missing <- setdiff(attrs, names(data))
        if (length(missing)) {
            cli::cli_abort(
                "Cannot find {missing} attributes in {.arg {arg}}",
                call = call
            )
        }
        attrs <- .subset(data, attrs)
    } else if (isTRUE(attrs)) {
        attrs <- data
    } else if (isFALSE(attrs)) {
        attrs <- NULL
    } else if (!is.null(attrs) && !is.list(attrs)) {
        cli::cli_abort(
            "{.arg {arg}} must be a bool or a character or a list"
        )
    }
    attrs
}

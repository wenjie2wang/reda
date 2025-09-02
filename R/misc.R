##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2025
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


### some trivial internal functions ============================================
## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list()) {
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## warning if x contains NA (or NaN)
na_warning <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        warning(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## stop if x contains NA (or NaN)
na_stop <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        stop(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## is x a numeric matrix (optionally of nRow rows and nCol columns)
isNumMatrix <- function(x, nRow = NULL, nCol = NULL,
                        warn_na = FALSE, error_na = FALSE,
                        sub_env = "parent", ...)
{
    out <- is.numeric(x) && is.matrix(x)
    if (out) {
        nDim <- dim(x)
        if (! is.null(nRow)) out <- out && nDim[1L] == nRow
        if (! is.null(nCol)) out <- out && nDim[2L] == nCol
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}

## is x a numeric vector
isNumVector <- function(x, warn_na = FALSE, error_na = FALSE,
                        sub_env = "parent", ...)
{
    out <- is.numeric(x) && is.vector(x)
    if (out) {
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}

## is x a numeric value
isNumOne <- function(x, sub_env = "grandparent", ...)
{
    isNumVector(x, sub_env = sub_env, ...) && length(x) == 1L
}

## is x a character vector
isCharVector <- function(x, warn_na = FALSE, error_na = FALSE,
                         sub_env = "parent", ...)
{
    out <- is.character(x) && is.vector(x)
    if (out) {
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}

## is x a character value
isCharOne <- function(x, sub_env = "grandparent", ...)
{
    isCharVector(x, sub_env = sub_env, ...) && length(x) == 1L
}

## is x a logical vector
isLogicalVector <- function(x, warn_na = FALSE, error_na = FALSE,
                            sub_env = "parent", ...)
{
    out <- is.logical(x) && is.vector(x)
    if (out) {
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}

## is x a logical value
isLogicalOne <- function(x, sub_env = "grandparent", ...) {
    isLogicalVector(x, sub_env = sub_env, ...) && length(x) == 1L
}

## is `x` object of class `foo`?
## is x a Survr object
is.Survr <- function(x) {
    is(x, "Survr")
}
## is x a mcf.sample object
is.mcf.formula <- function(x) {
    is(x, "mcf.formula")
}
## is x a rateReg object
is.rateReg <- function(x) {
    is(x, "rateReg")
}

## throw warnings if `...` is specified by mistake
warn_dots <- function(...) {
    dotsList <- list(...)
    if (length(dotsList) > 0) {
        .fun_name <- as.character(sys.call(- 1L)[[1L]])
        list_names <- names(dotsList)
        if (is.null(list_names)) {
            warning(wrapMessages(
                sprintf(paste("Some invalid argument(s) went into `...`",
                              "of %s()"),
                        .fun_name)
            ), call. = FALSE)
        } else {
            list_names <- list_names[list_names != ""]
            if (length(list_names) > 2) {
                all_names <- paste(sprintf("'%s'", list_names), collapse = ", ")
                all_names <- gsub("(.+), (.+)$", "\\1, and \\2", all_names)
            } else {
                all_names <- paste(sprintf("'%s'", list_names),
                                   collapse = " and ")
            }
            warning(wrapMessages(
                sprintf("Invalid argument %s went into `...` of %s()",
                        all_names, .fun_name)
            ), call. = FALSE)
        }
    }
    invisible(NULL)
}

## simplified version of utils::modifyList with keep.null = TRUE
modify_list <- function (x, val) {
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    for (v in vnames) {
        x[v] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
                    list(modify_list(x[[v]], val[[v]]))
                else val[v]
    }
    x
}

## check if the suggested package is available
suggest_pkg <- function(pkg_name)
{
    if (! requireNamespace(pkg_name, quietly = TRUE)) {
        stop(sprintf("The package '%s' is required for this function.",
                     pkg_name), call. = FALSE)
    }
    invisible()
}

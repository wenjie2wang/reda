### some internal functions taken from formatR v1.7: usage.R
formatR_env <- new.env()
with(formatR_env, {
    deparse_collapse = function(x) {
        d = deparse(x)
        if (length(d) > 1L) {
            paste(trimws(d, which = 'both'), collapse = ' ')
        } else {
            d
        }
    }

    count_tokens = function(.call) {
        if (length(.call) == 1L) {
            ## +2 for '()'
            return(nchar(.call) + 2L)
        }
        ## +1 for value-delimiting '(', ',', or ')'
        cnt_val = nchar(vapply(.call, deparse_collapse, character(1L))) + 1L
        nms = names(.call[-1L])
        if (is.null(nms)) nms = character(length(.call[-1L]))
        ## nchar() of argument names
        cnt_nm = nchar(nms)
        ## +3 for ' = ', for argument-value pairs
        cnt_nm[cnt_nm != 0L] = cnt_nm[cnt_nm != 0L] + 3L
        ## +1 for space before name, beyond the first argument
        cnt_nm[-1L] = cnt_nm[-1L] + 1L
        ## function itself is not a named component
        cnt_nm = c(0L, cnt_nm)
        cumsum(cnt_val + cnt_nm)
    }
    ## counts is a strictly increasing, positive integer vector
    find_breaks = function(counts, width, indent, track, counted = 0L) {
        if (!length(counts)) {
            return(list(breaks = NULL, overflow = NULL))
        }
        overflow = NULL
        shift = if (counted == 0L) 0L else indent
        fits = counts - counted + shift <= width
        i = which.min(fits) - 1L
        if (i == 0L) {
            if (fits[[1L]]) {
                ## all components of fits_on_line are TRUE
                i = length(counts)
            } else {
                ## all components of fits_on_line are FALSE
                overflow = track(counted, counts[1L], shift)
                i = 1L
            }
        }
        post_space = if (i == 1L && counted == 0L) 0L else 1L
        rest = Recall(counts[-(1L:i)], width, indent, track,
                      counts[i] + post_space)
        list(
            breaks   = c(counts[i], rest$breaks),
            overflow = c(overflow, rest$overflow)
        )
    }

    overflow_message = function(overflow, width, indent, text) {
        header = sprintf(
            'Could not fit all lines to width %s (with indent %s):',
            width, indent
        )
        idxs = seq_along(overflow)
        args = vapply(idxs[idxs %% 3L == 1L], function(i) {
            l = paste(c(rep(' ', overflow[i + 2L]),
                        trimws(substr(text, overflow[i] + 1L, overflow[i + 1L]),
                               which = 'left')),
                      collapse = '')
            sprintf('(%s) \"%s\"', nchar(l), l)
        }, character(1L))
        paste(c(header, args), collapse = '\n')
    }

    tidy_usage = function(nm, usg, width, indent, fail) {
        text = paste(trimws(usg, which = 'both'), collapse = ' ')
        text = sub(sprintf('^%s\\s*', nm), nm, text)
        expr = parse(text = text)[[1L]]
        track_overflow = if (fail == 'none') function(...) NULL else base::c
        breaks = find_breaks(count_tokens(expr), width, indent, track_overflow)
        if (length(breaks$overflow)) {
            signal = switch(fail, stop = 'stop', warn = 'warning')
            msg = overflow_message(breaks$overflow, width, indent, text)
            getFromNamespace(signal, 'base')(msg, call. = FALSE)
        }
        breaks = c(0L, breaks$breaks)
        newline = paste(c('\n', character(indent)), collapse = ' ')
        paste(
            vapply(1L:(length(breaks) - 1L), function(i) {
                trimws(substr(text, breaks[i] + 1L, breaks[i + 1L]),
                       which = 'left')
            }, character(1L)),
            collapse = newline
        )
    }
})


### re-format usage section following the style of base R
format_rd_usage <- function(rd_file) {
    fun_doc <- readLines(rd_file)
    ## locate usage section
    has_usage <- grepl("\\usage{", fun_doc, fixed = TRUE)
    if (! any(has_usage)) {
        return(invisible(NULL))
    }
    usage_begin <- which(has_usage)
    usage_end <- which(grepl("\\arguments{", fun_doc, fixed = TRUE)) - 1L
    usage_idx <- seq.int(usage_begin, usage_end)
    usage_doc <- fun_doc[usage_idx]
    ## main function
    format_usage <- function(usage_doc) {
        ## locate all left and right parentheses
        left_m <- gregexpr("(", usage_doc, fixed = TRUE)
        right_m <- gregexpr(")", usage_doc, fixed = TRUE)
        foo <- function(x) {
            do.call(c, lapply(seq_along(x), function(i) {
                if (any(x[[i]] > 0)) {
                    rep(i, length(x[[i]]))
                } else {
                    NULL
                }
            }))
        }
        left_par <- foo(left_m)
        right_par <- foo(right_m)
        if (! length(left_par))
            return(usage_doc[- c(1, length(usage_doc))])
        ## get the pair of function left and right parentheses
        drop_idx <- left_par[- 1L] <= right_par[- length(right_par)]
        if (any(drop_idx)) {
            drop_which <- which(drop_idx)
            left_par <- left_par[- (drop_which + 1)]
            right_par <- right_par[- drop_which]
        }
        ## get function name
        get_fun_name <- function(x) {
            if (grepl("{", x, fixed = TRUE)) {
                m <- regexpr("\\{.*?\\}", x)
                nm <- regmatches(x, m)
                gsub("^\\{|\\}$", "", nm)
            } else {
                gsub("\\(.*$", "", x)
            }
        }
        out <- lapply(seq_along(left_par), function(i) {
            fun_name <- get_fun_name(usage_doc[left_par[i]])
            num_exdent <- nchar(fun_name) + 1L
            fun_doc <- gsub(
                "[ ]+", " ",
                paste(usage_doc[seq.int(left_par[i], right_par[i])],
                      sep = "", collapse = "")
            )
            fun_name_doc <- gsub("(^.*?)\\(.*\\)$", "\\1", fun_doc)
            fun_args <- trimws(gsub("^.*?\\((.*)\\)$", "\\1", fun_doc))
            for (k in seq_len(4L)) {
                width_set <- 72L - 5L + (k - 1) * 5L
                res <- tryCatch(
                    formatR_env$tidy_usage(
                                    fun_name,
                                    c(fun_name, sprintf("(%s)", fun_args)),
                                    width = width_set,
                                    indent = num_exdent,
                                    fail = "stop"
                                ),
                    error = function(e) e
                )
                if (! "error" %in% class(res)) {
                    break
                }
                if (width_set > 80) {
                    res <- sprintf("%s(%s)", fun_name, fun_args)
                }
            }
            res <- sprintf(gsub(
                "^.*?\\((.*)\\)$", "%s(\\1)", res
            ), fun_name_doc)
            c(do.call(c, strsplit(res, split = "[ ]?\\n")), "")
        })
        out <- do.call(c, out)
        out[- length(out)]
    }
    new_rd <- c(fun_doc[seq.int(usage_begin)],
                format_usage(usage_doc),
                fun_doc[seq.int(usage_end, length(fun_doc))])
    writeLines(new_rd, rd_file)
}

## overwrite rd files in man
rd_dir <- "man"
rd_files <- list.files(path = rd_dir, pattern = "\\.Rd$")
for (i in rd_files) {
    format_rd_usage(file.path(rd_dir, i))
}

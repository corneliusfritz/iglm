#' @importFrom Rcpp compileAttributes
#' @importFrom RcppArmadillo fastLm
#' @importFrom Matrix bdiag sparseMatrix
# Utility: deparse to a single string
.deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

.split_plus <- function(expr) {
  out <- list()
  rec <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      rec(e[[2L]])
      rec(e[[3L]])
    } else {
      out[[length(out) + 1L]] <<- e
    }
  }
  rec(expr)
  out
}

find_ranks <- function(x) {
  match(x, sort(unique(x)))
}

get_rank <- function(new_vals, ref_data) {
  ref_sorted <- sort(ref_data)
  n_less <- findInterval(new_vals, ref_sorted, left.open = TRUE)
  n_less_equal <- findInterval(new_vals, ref_sorted)
  ((n_less + 1) + (n_less_equal + 1)) / 2
}

update_formula_remove_terms <- function(formula, terms_to_remove) {
  rhs_terms <- attr(terms(formula), "term.labels")
  rhs_terms_updated <- rhs_terms[!rhs_terms %in% terms_to_remove]
  new_formula <- as.formula(paste(
    deparse(formula[[2]]), "~",
    paste(rhs_terms_updated, collapse = " + ")
  ))
  environment(new_formula) <- environment(formula)
  new_formula
}

.make_unique_name <- function(nm, existing) {
  if (!(nm %in% existing)) {
    return(nm)
  }
  k <- 2L
  while (paste0(nm, "_", k) %in% existing) k <- k + 1L
  paste0(nm, "_", k)
}

eval_change <- function(formula, additional_args = NULL, object) {
  term_labels <- attr(terms(formula), "term.labels")
  # Initialize a list to hold all results
  all_results <- list()
  # Iterate over each term label
  for (i in seq_along(term_labels)) {
    term_string <- term_labels[i]
    # Convert the string term back into an expression object
    expression_obj <- parse(text = term_string)[[1]]
    # Convert the expression into a list of its components
    # The first element is the function name (a symbol), the rest are arguments
    call_list <- as.list(expression_obj)
    # Function name is the first element, converted to a character string
    func_name <- as.character(call_list[[1]])
    # Arguments are the remaining elements
    args <- call_list[-1]
    if (!is.null(additional_args)) {
      args[[length(args) + 1]] <- additional_args[[i]]
      names(args)[length(args)] <- names(additional_args)[i]
    }
    args$plot <- FALSE
    # The function/method to call is the named item inside the R6 object
    method_to_call <- object[[func_name]]
    if (is.function(method_to_call)) {
      # Execute the function call using the
      # iglm.data object's environment (self)
      result <- do.call(method_to_call, args)

      # Store the result
      all_results[[term_string]] <- result
    } else {
      warning(paste("Method not found for term:", func_name))
    }
  }
  return(all_results)
}

rhs_terms_as_list <- function(formula, env = NULL, evaluate_calls = FALSE) {
  formula <- as.formula(formula)
  if (is.null(env)) {
    env <- environment(formula)
    if (is.null(env)) env <- parent.frame()
  }
  rhs_expr <- formula[[3L]]
  terms_exprs <- .split_plus(rhs_expr)

  out <- list()
  taken_names <- character(0L)

  for (term_expr in terms_exprs) {
    if (is.symbol(term_expr)) {
      base_name <- as.character(term_expr)
      if (base_name %in% taken_names) {
        next
      }
      taken_names <- c(taken_names, base_name)

      out[[base_name]] <- list(
        label = .deparse1(term_expr),
        base_name = base_name
      )
    } else if (is.call(term_expr)) {
      fun_sym <- term_expr[[1L]]
      base_name <- if (is.symbol(fun_sym)) as.character(fun_sym) else .deparse1(fun_sym)

      # raw argument expressions
      arg_exprs <- as.list(term_expr)[-1L]
      arg_names <- names(arg_exprs)
      if (is.null(arg_names)) arg_names <- rep("", length(arg_exprs))

      # prepare container with named, evaluated arguments
      # unnamed arguments get positional names ..1, ..2, ...
      pos_names <- ifelse(arg_names == "", paste0("..", seq_along(arg_exprs)),
        arg_names
      )
      arg_vals <- list()
      for (i in seq_along(arg_exprs)) {
        val_expr <- arg_exprs[[i]]
        # If it's a character literal, keep it as character
        if (is.character(val_expr)) {
          val <- val_expr
        } else {
          val <- try(eval(val_expr, envir = env), silent = TRUE)
          if (inherits(val, "try-error")) {
            cond <- attr(val, "condition")
            err_msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(val)
            stop(sprintf(
              "Could not evaluate argument '%s' in term '%s': %s",
              pos_names[i], base_name, trimws(err_msg)
            ), call. = FALSE)
          }
        }

        nm <- arg_names[i]
        if (nm == "") {
          arg_vals[[paste0("..", i)]] <- val
        } else {
          arg_vals[[nm]] <- val
        }
      }

      # optionally evaluate the whole call
      evaluated <- NULL
      if (evaluate_calls) {
        tmp <- try(eval(term_expr, envir = env), silent = TRUE)
        if (!inherits(tmp, "try-error")) evaluated <- tmp
      }

      entry <- c(
        list(label = gsub(pattern = '\\\"', replacement = "'", x = .deparse1(term_expr))),
        arg_vals,
        list(base_name = base_name)
      )
      if (!is.null(evaluated)) entry$.evaluated <- evaluated

      # For the elt_name (key in the output list), we still want the full name
      # but we MUST NOT change base_name in the entry list
      name_addon <- ""
      if (!is.null(arg_exprs$type)) name_addon <- paste0(name_addon, "_", .deparse1(arg_exprs$type))
      if (!is.null(arg_exprs$data)) name_addon <- paste0(name_addon, "_", .deparse1(arg_exprs$data))
      if (!is.null(arg_exprs$mode)) name_addon <- paste0(name_addon, "_", .deparse1(arg_exprs$mode))
      if (!is.null(arg_exprs$variant)) name_addon <- paste0(name_addon, "_", .deparse1(arg_exprs$variant))

      elt_name <- paste0(base_name, name_addon)
      # Ensure uniqueness
      if (elt_name %in% names(out)) {
        suffix <- 2
        while (paste0(elt_name, ".", suffix) %in% names(out)) suffix <- suffix + 1
        elt_name <- paste0(elt_name, ".", suffix)
      }
      out[[elt_name]] <- entry
    }
  }
  class(out) <- "iglm.formulainfo"
  out
}

#' @export
#' @method print iglm.formulainfo
print.iglm.formulainfo <- function(x, ..., max_items = 5) {
  n_terms <- length(x)
  cat(
    "<iglm.formulainfo> object with", n_terms,
    if (n_terms == 1L) "term" else "terms", "\n\n"
  )

  for (nm in names(x)) {
    term <- x[[nm]]

    cat("$", nm, " (", term$type, ")\n", sep = "")
    cat("  label: ", term$label, "\n", sep = "")

    if (term$type == "symbol") {
      # show value class or small preview
      val <- term$value
      if (is.null(val)) {
        cat("  value: NULL\n\n")
      } else {
        cat("  value: <", class(val)[1L], ">",
          if (is.atomic(val) && length(val) <= max_items) {
            paste0(" ", toString(val))
          },
          if (length(val) > max_items) " ...", "\n\n",
          sep = ""
        )
      }
    } else if (term$type == "call") {
      # show each argument entry
      arg_names <- setdiff(names(term), c("name", "type", "label", ".evaluated"))
      if (length(arg_names) == 0L) {
        cat("  (no arguments)\n\n")
      } else {
        for (an in arg_names) {
          val <- term[[an]]
          cat("  ", an, " = ", sep = "")
          if (is.null(val)) {
            cat("NULL\n")
          } else if (is.atomic(val) && length(val) <= max_items) {
            cat(toString(val), "\n")
          } else {
            cat("<", class(val)[1L], ">", sep = "")
            if (is.data.frame(val)) cat(" [", nrow(val), "x", ncol(val), "]", sep = "")
            cat("\n")
          }
        }
        cat("\n")
      }
    }
  }
  invisible(x)
}


map_to_mat <- function(map, n_actor) {
  # Generate empty network
  mat <- matrix(0, nrow = n_actor, ncol = n_actor)
  for (i in 1:n_actor) {
    if (length(map[[i + 1]]) > 0) {
      mat[i, map[[i + 1]]] <- 1
    }
  }
  return(mat)
}

set_to_vec <- function(set, n_actor) {
  # Generate empty vector
  vec <- numeric(length = n_actor)
  vec[set] <- 1
  return(vec)
}

XZ_to_R <- function(x_attribute, z_network, n_actor) {
  x_attribute <- set_to_vec(set = x_attribute, n_actor = n_actor)
  z_network <- map_to_mat(map = z_network, n_actor = n_actor)
  return(list(x_attribute = x_attribute, z_network = z_network))
}

XYZ_to_R <- function(x_attribute, y_attribute, z_network, n_actor, return_adj_mat) {
  # x_attribute = set_to_vec(set = x_attribute,n_actor = n_actor)
  # y_attribute = set_to_vec(set = y_attribute,n_actor = n_actor)
  if (return_adj_mat) {
    z_network_tmp <- map_to_mat(map = z_network, n_actor = n_actor)
  } else {
    z_network_tmp <- do.call(rbind, lapply(1:n_actor, FUN = function(x) {
      tmp <- z_network[[x + 1]]
      if (length(tmp) == 0) {
        return(NA)
      } else {
        return(cbind(x, tmp))
      }
    }))
    z_network_tmp <- z_network_tmp[!is.na(z_network_tmp[, 1]), ]
    if (length(z_network_tmp) == 0) {
      z_network_tmp <- matrix(numeric(0), ncol = 2)
    }
    colnames(z_network_tmp) <- c("from", "to")
  }
  return(list(
    x_attribute = x_attribute, y_attribute = y_attribute,
    z_network = z_network_tmp
  ))
}

check_overlap <- function(mat_1, mat_2) {
  colnames(mat_1) <- colnames(mat_2)
  combined <- rbind(mat_1, mat_2)
  return(duplicated(combined, fromLast = TRUE)[seq_len(nrow(mat_1))])
}

iglm.data.neighborhood <- function(neighborhood, directed = NA, n_actor = NA) {
  if (!is.matrix(neighborhood) && !is.data.frame(neighborhood)) {
    if (length(neighborhood) == 0) {
      neighborhood <- matrix(numeric(0), nrow = 0, ncol = 2)
    } else {
      stop("`neighborhood` must be a matrix or data frame.", call. = FALSE)
    }
  }

  if (nrow(neighborhood) == 0) {
    if (ncol(neighborhood) != 0 && ncol(neighborhood) != 2) {
      stop("Empty neighborhood matrix must have 0 or 2 columns (edgelist).", call. = FALSE)
    }
    res <- list(
      neighborhood = matrix(numeric(0), nrow = 0, ncol = 2),
      overlap = matrix(numeric(0), nrow = 0, ncol = 2)
    )
    class(res) <- "iglm.data.neighborhood"
    return(res)
  }
  if (is.na(n_actor)) {
    if (ncol(neighborhood) > 2) {
      n_actor <- nrow(neighborhood)
    } else {
      n_actor <- max(neighborhood)
    }
  }
  if (is.na(n_actor)) {
    stop("n_actor could not be inferred. Please provide n_actor.")
  }
  if (ncol(neighborhood) == 2) {
    sp_nb <- spMatrix(
      nrow = n_actor, ncol = n_actor,
      i = neighborhood[, 1], j = neighborhood[, 2], x = rep(1, length(neighborhood[, 2]))
    )
    sp_nb_trans <- sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)

    overlap <- sp_nb %*% sp_nb_trans
    overlap <- as(overlap, "TsparseMatrix")
    overlap <- cbind(
      overlap@i + 1,
      overlap@j + 1
    )
    overlap <- overlap[overlap[, 1] != overlap[, 2], ]
  } else {
    positions <- which(neighborhood == 1, arr.ind = T)
    sp_nb <- spMatrix(
      nrow = ncol(neighborhood), ncol = ncol(neighborhood),
      i = positions[, 1], j = positions[, 2], x = rep(1, length(positions[, 2]))
    )
    sp_nb_trans <- sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)

    overlap <- as.matrix(sp_nb %*% sp_nb_trans > 0)
    diag(overlap) <- 0
    overlap <- which(overlap == 1, arr.ind = T)
    neighborhood <- which(neighborhood == 1, arr.ind = T)
  }
  res <- list(
    neighborhood = neighborhood,
    overlap = overlap
  )
  class(res) <- "iglm.data.neighborhood"
  return(res)
}

get_i <- function(x, i) {
  stopifnot(is.list(x), length(i) == 1L, is.numeric(i), i >= 1L)
  k <- 1L
  for (el in x) {
    if (k == i) {
      return(el)
    }
    k <- k + 1L
  }
  stop("subscript out of bounds")
}

#' @export
#' @method [[ iglm.data.list
`[[.iglm.data.list` <- function(x, i, ...) {
  item <- get_i(x, i)
  # browser()
  item$set_neighborhood_overlap(
    attr(x, "neighborhood")$neighborhood,
    attr(x, "neighborhood")$overlap
  )
  item
}

append_iglm.data <- function(x, y) {
  tmp <- c(x, y)
  attributes(tmp) <- attributes(x)
  class(tmp) <- "iglm.data.list"
  return(tmp)
}
#' @export
#' @method [ iglm.data.list
`[.iglm.data.list` <- function(x, i, ...) {
  # browser()
  res <- list()
  k <- 1
  for (j in i) {
    item <- get_i(x, j)
    item$set_neighborhood_overlap(
      attr(x, "neighborhood")$neighborhood,
      attr(x, "neighborhood")$overlap
    )
    res[[k]] <- item
    names(res)[k] <- j
    k <- k + 1
  }
  res
}

#' @export
#' @method print iglm.data.list
print.iglm.data.list <- function(x, ...) {
  # Header
  n_items <- length(x)
  cat(
    "List of iglm.data object with", n_items,
    if (n_items == 1L) "entry\n" else "entries\n"
  )

  # Summarize elements
  if (n_items == 0L) {
    cat("(empty list)\n")
    return(invisible(x))
  }

  nm <- names(x)
  if (is.null(nm)) nm <- paste0("[[", seq_len(n_items), "]]\n")

  for (i in seq_len(n_items)) {
    el <- x[[i]]
    name <- nm[i]
    cat(name, sep = "")
    print(el)
    cat("\n")
  }

  invisible(x)
}

formula_preprocess <- function(formula) {
  data_object <- eval(formula[[2]], envir = environment(formula))
  if (!inherits(data_object, "iglm.data")) {
    stop("The response in the formula must be an iglm.data object.")
  }

  includes_degrees <- "degrees" %in% all.vars(formula)
  formula <- stats::update(formula, . ~ . - degrees)
  formula_info <- rhs_terms_as_list(formula)
  if (length(formula_info) == 0 && !includes_degrees) {
    stop("Formula must contain at least one term on the right-hand side.", call. = FALSE)
  }

  term_names <- character(length(formula_info))
  data_list <- list()
  type_list <- numeric(length(formula_info))
  coef_names <- character(length(formula_info))

  for (i in seq_along(formula_info)) {
    arglist <- formula_info[[i]]
    # Call the modular initialization system
    init <- InitIglmTerm(data_object = data_object, arglist = arglist)

    term_names[i] <- init$term_name
    data_list[[i]] <- if (is.null(init$data)) matrix(1) else init$data
    type_list[i] <- if (is.null(init$type)) 1L else init$type
    coef_names[i] <- init$coef_name
  }

  return(list(
    data_object = data_object,
    data_list = data_list,
    type_list = type_list,
    coef_names = coef_names,
    term_names = term_names,
    includes_degrees = includes_degrees
  ))
}

is_string_a_function_execution <- function(s) {
  obj <- try(str2lang(s), silent = TRUE)

  if (inherits(obj, "try-error")) {
    return(FALSE)
  }

  if (!is.call(obj)) {
    return(FALSE)
  }
  head_obj <- obj[[1]]

  if (is.symbol(head_obj)) {
    func_name <- as.character(head_obj)

    # List of special operators that are 'calls' but not 'executions'.
    special_operators <- c(
      "(", "[", "[[", "{", "$",
      "+", "-", "*", "/", "^", "%%", "%/%", "%*%",
      "<", "<=", "==", "!=", ">=", ">",
      "&", "&&", "|", "||", "!",
      "~", ":", "=", "<-", "<<-", "->", "->>"
    )
    if (func_name %in% special_operators) {
      return(FALSE)
    }

    return(TRUE)
  }

  return(TRUE)
}

#' @noRd
draw_simulation_envelope <- function(x, sim_matrix, color, alpha = 0.4, lwd_mean = 2) {
  x_poly <- c(x, rev(x))
  y_poly <- c(colMins(sim_matrix), rev(colMaxs(sim_matrix)))
  graphics::polygon(x_poly, y_poly, col = add_alpha(color, alpha_level = alpha), border = NA)
  graphics::lines(x, colMeans(sim_matrix), type = "l", col = color, lwd = lwd_mean)
  graphics::lines(x, colMins(sim_matrix), type = "l", col = color, lwd = 1)
  graphics::lines(x, colMaxs(sim_matrix), type = "l", col = color, lwd = 1)
}

#' @noRd
extract_assessment_matrix <- function(sim_list, metric, subkey = NULL) {
  sims <- lapply(sim_list, function(x) {
    val <- x[[metric]]
    if (!is.null(subkey)) val <- val[[subkey]]
    val
  })
  do.call("rbind", sims)
}

#' @noRd
plot_assessment_multi <- function(observed, sim_main, sim_dots = list(),
                                  model_names, colors, xlab, ylab = "Percentage",
                                  x_at = NULL, x_labels = NULL, x_positions = NULL,
                                  x_margin = 0.3, lwd_mean = 2) {
  all_sims <- c(list(sim_main), sim_dots)
  
  # Check if observed or simulation matrices have column names to align on
  all_names <- names(observed)
  for (s in all_sims) {
    if (!is.null(colnames(s))) all_names <- union(all_names, colnames(s))
  }
  
  if (!is.null(all_names) && is.null(x_positions)) {
    # If all names are numeric, sort them numerically
    num_names <- suppressWarnings(as.numeric(all_names))
    if (!any(is.na(num_names))) {
      all_names <- as.character(sort(num_names))
      x <- as.numeric(all_names)
    } else {
      x <- seq_along(all_names)
    }
    
    obs_aligned <- stats::setNames(rep(0, length(all_names)), all_names)
    common_obs <- intersect(names(observed), all_names)
    obs_aligned[common_obs] <- as.numeric(observed[common_obs])
    observed <- obs_aligned
    
    all_sims <- lapply(all_sims, function(sim_mat) {
      if (!is.null(colnames(sim_mat))) {
        aligned <- matrix(0, nrow = nrow(sim_mat), ncol = length(all_names), dimnames = list(NULL, all_names))
        common <- intersect(colnames(sim_mat), all_names)
        aligned[, common] <- sim_mat[, common]
        aligned
      } else {
        sim_mat
      }
    })
  } else {
    x <- if (!is.null(x_positions)) {
      as.numeric(x_positions)
    } else if (!is.null(names(observed))) {
      as.numeric(names(observed))
    } else {
      seq_along(observed)
    }
  }
  
  xlim <- c(min(x) - x_margin, max(x) + x_margin)
  all_vals <- c(as.numeric(observed), unlist(all_sims))
  ylim <- range(all_vals, na.rm = TRUE)
  
  plot(x, as.vector(observed),
       type = "n", xlab = xlab, ylab = ylab,
       xlim = xlim, ylim = ylim, las = 1, axes = FALSE
  )
  
  if (is.null(x_at)) {
    x_at <- pretty(range(x), n = 10)
  }
  
  if (is.null(x_labels)) {
    axis(side = 1, at = x_at, lwd = 0, lwd.ticks = 1)
  } else {
    axis(side = 1, at = x_at, labels = x_labels, lwd = 0, lwd.ticks = 1)
  }
  axis(side = 2, las = 1, lwd = 0, lwd.ticks = 1)
  
  box(bty = "l", lwd = 1)
  
  draw_simulation_envelope(x, all_sims[[1]], color = colors[1], alpha = 0.1, lwd_mean = lwd_mean)
  
  if (length(all_sims) > 1) {
    for (m in 2:length(all_sims)) {
      draw_simulation_envelope(x, all_sims[[m]], color = colors[m], alpha = 0.1, lwd_mean = lwd_mean)
    }
  }
  lines(x, as.vector(observed), type = "l", col = "black", lwd = 2)
  
  legend("topright",
         legend = c("Observed", model_names),
         col = c("black", colors),
         lty = c(1, 1, rep(1, max(0, length(model_names) - 1))),
         lwd = c(2, rep(lwd_mean, length(model_names))),
         bty = "n"
  )
}

#' @noRd
plot_assessment_single <- function(observed, sim_matrix, xlab, ylab = "Percentage",
                                   box_col = "#87CEEB80", line_col = "#D55E00",
                                   x_at = NULL, x_labels = NULL, x_positions = NULL,
                                   x_margin = 0.3, use_envelope = FALSE) {
  ylim <- range(c(sim_matrix, as.numeric(observed)), na.rm = TRUE)
  
  if (is.null(x_positions)) {
    x <- if (!is.null(names(observed))) as.numeric(names(observed)) else seq_along(observed)
  } else {
    x <- as.numeric(x_positions)
  }
  
  xlim <- c(min(x) - x_margin, max(x) + x_margin)
  
  # 1. Initialize plot (axes = FALSE suppresses default axes and box)
  plot(x, as.vector(observed),
       type = "n", xlab = xlab, ylab = ylab,
       xlim = xlim, ylim = ylim, las = 1, axes = FALSE
  )
  
  if (is.null(x_at)) {
    x_at <- pretty(range(x), n = 10)
  }
  
  # 2. Draw ticks and labels directly on the plot boundary (suppress axis lines)
  if (is.null(x_labels)) {
    axis(side = 1, at = x_at, lwd = 0, lwd.ticks = 1)
  } else {
    axis(side = 1, at = x_at, labels = x_labels, lwd = 0, lwd.ticks = 1)
  }
  axis(side = 2, las = 1, lwd = 0, lwd.ticks = 1)
  
  # 3. Draw the exact continuous L-shaped corner along the plot boundary
  box(bty = "l", lwd = 1)
  
  # 4. Plot data layers
  if (use_envelope) {
    draw_simulation_envelope(x, sim_matrix, color = box_col, alpha = 0.4, lwd_mean = 2)
  } else {
    boxplot(sim_matrix,
            at = x, boxwex = 0.5,
            add = TRUE, col = box_col, axes = FALSE, las = 1
    )
  }
  
  lines(x, as.vector(observed), type = "l", col = line_col, lwd = 2)
}

#' @noRd
plot_multitrace <- function(mat, xlab = "Iteration", ylab = "Coefficients", las = 1, bty = "l", ...) {
  mat <- as.matrix(mat)
  plot(NA,
    xlim = c(1, nrow(mat)), ylim = range(mat, na.rm = TRUE),
    xlab = xlab, ylab = ylab, las = las, bty = bty, ...
  )
  for (tmp in seq_len(ncol(mat))) {
    lines(y = mat[, tmp], x = seq_len(nrow(mat)), col = tmp)
  }
}

#' @noRd
get_assessment_constraint_suffix <- function(name, base_name) {
  raw_suffix <- sub(paste0("^", base_name, "_?"), "", name)
  if (!nzchar(raw_suffix)) return("")
  matches <- gregexpr("(x_i|x_j|y_i|y_j|mode)_([^_]+)", raw_suffix)
  reg_matches <- regmatches(raw_suffix, matches)[[1]]
  if (length(reg_matches) > 0) {
    formatted <- sub("^([a-z_]+)_(.+)$", "\\1 = \\2", reg_matches)
    return(paste0(" (", paste(formatted, collapse = ", "), ")"))
  } else {
    return(paste0(" (", gsub("_", " ", raw_suffix), ")"))
  }
}


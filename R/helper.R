#' @importFrom Rcpp compileAttributes 
#' @importFrom RcppArmadillo fastLm
#' @importFrom Matrix bdiag sparseMatrix
# Utility: deparse to a single string
.deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

.split_plus <- function(expr) {
  out <- list()
  rec <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      rec(e[[2L]]); rec(e[[3L]])
    } else {
      out[[length(out) + 1L]] <<- e
    }
  }
  rec(expr)
  out
}

find_ranks <- function(x){
  match(x, sort(unique(x)))  
}

get_rank <- function(new_vals, ref_data) {
  ref_sorted <- sort(ref_data)
  n_less <- findInterval(new_vals, ref_sorted, left.open = TRUE)
  n_less_equal <- findInterval(new_vals, ref_sorted)
  avg_rank <- ( (n_less + 1) + (n_less_equal + 1) ) / 2
  
  return(avg_rank)
}

update_formula_remove_terms <- function(formula, terms_to_remove){
  rhs_terms <- attr(terms(formula), "term.labels")
  rhs_terms_updated <- rhs_terms[!rhs_terms %in% terms_to_remove]
  new_formula <- as.formula(paste(deparse(formula[[2]]), "~", 
                                  paste(rhs_terms_updated, collapse = " + ")))
  environment(new_formula) <- environment(formula)
  return(new_formula)
}

.make_unique_name <- function(nm, existing) {
  if (!(nm %in% existing)) return(nm)
  k <- 2L
  while (paste0(nm, "_", k) %in% existing) k <- k + 1L
  paste0(nm, "_", k)
}

eval_change <- function(formula, additional_args = NULL,object) {
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
    if(!is.null(additional_args)){
      args[[length(args)+1]] <- additional_args[[i]]
      names(args)[length(args)] <- names(additional_args)[i]
    }
    # The function/method to call is the named item inside the R6 object
    method_to_call <- object[[func_name]]
    if (is.function(method_to_call)) {
      
      # Execute the function call using the iglm.data object's environment (self)
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
  # browser()
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
      # cat("Symbol")
      base_name <- as.character(term_expr)
      if(base_name %in% taken_names){
        next
      }
      taken_names <- c(taken_names, base_name)
      
      out[[base_name]] <- list( 
        label = .deparse1(term_expr),
        base_name = base_name
      )
      
    } else if (is.call(term_expr)) {
      # cat("call term")
      fun_sym  <- term_expr[[1L]]
      base_name <- if (is.symbol(fun_sym)) as.character(fun_sym) else .deparse1(fun_sym)
      
      # raw argument expressions
      arg_exprs  <- as.list(term_expr)[-1L]
      arg_names  <- names(arg_exprs)
      if (is.null(arg_names)) arg_names <- rep("", length(arg_exprs))
      
      # prepare container with named, evaluated arguments
      # unnamed arguments get positional names ..1, ..2, ...
      pos_names <- ifelse(arg_names == "", paste0("..", seq_along(arg_exprs)), arg_names)
      arg_vals  <- vector("list", length(arg_exprs))
      names(arg_vals) <- pos_names
      
      for (i in seq_along(arg_exprs)) {
        v <- try(eval(arg_exprs[[i]], envir = env), silent = TRUE)
        arg_vals[[i]] <- if (inherits(v, "try-error")) NULL else v
      }
      
      # optionally evaluate the whole call
      evaluated <- NULL
      if (evaluate_calls) {
        tmp <- try(eval(term_expr, envir = env), silent = TRUE)
        if (!inherits(tmp, "try-error")) evaluated <- tmp
      }
      
      
      entry <- c(
        label = gsub(pattern = '\\\"', replacement = "'", x = .deparse1(term_expr)),
        arg_vals, 
        base_name = base_name
      )
      if (!is.null(evaluated)) entry$.evaluated <- evaluated
      name_addon <- ""
      if(!is.null(arg_exprs$type)){
        name_addon <- paste0(name_addon,.deparse1(arg_exprs$type))
      }
      if(!is.null(arg_exprs$data)){
        name_addon <- paste0(name_addon,.deparse1(arg_exprs$data))
      }
      if(!is.null(arg_exprs$mode)){
        name_addon <- paste0(name_addon,arg_exprs$mode)
      }
      if(!is.null(arg_exprs$variant)){
        name_addon <- paste0(name_addon,arg_exprs$variant)
      }
      elt_name  <- paste0(base_name,"_",name_addon)
      taken_names <- c(taken_names, elt_name)
      out[[elt_name]] <- entry
      
    }
  }
  class(out) <- "iglm.formulainfo"
  out
}

#' @method print iglm_formulainfo 
print.iglm_formulainfo <- function(x, ..., max_items = 5) {
  n_terms <- length(x)
  cat("<iglm.formulainfo> object with", n_terms, 
      if (n_terms == 1L) "term" else "terms", "\n\n")
  
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
            if (is.atomic(val) && length(val) <= max_items)
              paste0(" ", toString(val)),
            if (length(val) > max_items) " ...", "\n\n", sep = "")
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



map_to_mat = function(map,n_actor){
  # Generate empty network
  mat = matrix(0,nrow = n_actor,ncol = n_actor)
  for(i in 1:n_actor){
    mat[as.numeric(names(map)[i]),map[[i]]] = 1
  }
  return(mat)
}

set_to_vec = function(set,n_actor){
  # Generate empty vector
  vec = numeric(length = n_actor)
  vec[set] = 1
  return(vec)
}

XZ_to_R = function(x_attribute,z_network, n_actor) {
  x_attribute = set_to_vec(set = x_attribute,n_actor = n_actor)
  z_network = map_to_mat(map = z_network,n_actor = n_actor)
  return(list(x_attribute = x_attribute, z_network = z_network))
}

XYZ_to_R = function(x_attribute,y_attribute ,z_network,n_actor, return_adj_mat) {
  # x_attribute = set_to_vec(set = x_attribute,n_actor = n_actor)
  # y_attribute = set_to_vec(set = y_attribute,n_actor = n_actor)
  if(return_adj_mat){
    z_network_tmp = map_to_mat(map = z_network,n_actor = n_actor)  
  } 
  else {
    z_network = z_network[order(as.numeric(names(z_network)))]
    z_network_tmp = do.call(rbind, lapply(1:n_actor, FUN = function(x) {
      tmp = z_network[[x]]
      if(length(tmp) == 0){
        return(NA)
      } else {
        return(cbind(x, tmp))
      }
    }))
    z_network_tmp = z_network_tmp[!is.na(z_network_tmp[,1]),]
    if(length(z_network_tmp) == 0){
      z_network_tmp = matrix(numeric(0), ncol = 2)
    }
    colnames(z_network_tmp) = c("from", "to")  
  }
  return(list(x_attribute = x_attribute, y_attribute = y_attribute, 
              z_network = z_network_tmp))
}

check_overlap <- function(mat_1, mat_2){
  colnames(mat_1) <- colnames(mat_2)
  combined <- rbind(mat_1, mat_2)
  return(duplicated(combined, fromLast = TRUE)[1:nrow(mat_1)])
}

iglm.data.neighborhood = function(neighborhood, directed = NA, n_actor = NA){
  
  if(is.na(n_actor)){
    if(ncol(neighborhood)>2){
      n_actor = nrow(neighborhood)  
    } else {
      n_actor = max(neighborhood)
    }
  }
  if(is.na(n_actor)){
    stop("n_actor could not be inferred. Please provide n_actor.")
  }
  if(ncol(neighborhood) == 2){
    sp_nb = spMatrix(nrow = n_actor,ncol = n_actor, 
                     i = neighborhood[,1], j = neighborhood[,2],x = rep(1,length(neighborhood[,2])))
    sp_nb_trans = sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)
    
    overlap =  sp_nb %*%sp_nb_trans
    overlap = as(overlap, "TsparseMatrix")
    overlap = cbind(overlap@i + 1,
                    overlap@j+1)
    overlap = overlap[overlap[,1] != overlap[,2],]
  } else {
    positions = which(neighborhood == 1, arr.ind = T)
    sp_nb = spMatrix(nrow = ncol(neighborhood),ncol = ncol(neighborhood), 
                     i = positions[,1], j = positions[,2],x = rep(1,length(positions[,2])))
    sp_nb_trans = sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)
    
    overlap =  as.matrix(sp_nb %*%sp_nb_trans>0)
    diag(overlap) = 0
    overlap = which(overlap == 1, arr.ind = T)
    neighborhood = which(neighborhood == 1, arr.ind = T)
  }
  res = list(neighborhood = neighborhood,  
             overlap = overlap)
  class(res) = "iglm.data.neighborhood"
  return(res)
}

get_i <- function(x, i) {
  stopifnot(is.list(x), length(i) == 1L, is.numeric(i), i >= 1L)
  k <- 1L
  for (el in x) {
    if (k == i) return(el)
    k <- k + 1L
  }
  stop("subscript out of bounds")
}

#' @export
#' @method [[ iglm.data.list
`[[.iglm.data.list` <- function(x, i, ...) {
  item <- get_i(x,i)
  # browser()
  item$set_neighborhood_overlap(attr(x,"neighborhood")$neighborhood, 
                                attr(x,"neighborhood")$overlap)
  item
}

append_iglm.data <- function(x,y){
  tmp <- c(x,y)
  attributes(tmp)<- attributes(x)
  class(tmp) <- "iglm.data.list"
  return(tmp)
}
#' @export
#' @method [ iglm.data.list
`[.iglm.data.list` <- function(x, i, ...) {
  # browser()
  res <- list()
  k <- 1
  for(j in i){
    item <- get_i(x,j)
    item$set_neighborhood_overlap(attr(x,"neighborhood")$neighborhood, 
                                  attr(x,"neighborhood")$overlap)
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
  cat("List of iglm.data object with", n_items,
      if (n_items == 1L) "entry\n" else "entries\n")
  
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
    cat(name,sep = "")
    print(el)
    cat("\n")
  }
  
  invisible(x)
}

formula_preprocess = function(formula){
  data_object = eval(formula[[2]],envir = environment(formula))
  includes_degrees <- "degrees" %in% all.vars(formula)
  formula <- stats::update(formula, .~ .-degrees)
  # debugonce(rhs_terms_as_list)
  formula_info <- rhs_terms_as_list(formula)
  
  term_per_term <-  unlist(lapply(formula_info, function(x) x$base_name))
  special_terms <- c("spillover_yy_scaled","spillover_xx_scaled","spillover_yx_scaled",
                     "spillover_xy_scaled","gwdegree","gwdsp","gwesp", "edges", 
                     "mutual", "cov_z", "cov_z_out", "cov_z_in",
                     "inedges_y", "outedges_y","attribute_xy", "inedges_x", "outedges_x")
  is_special <- term_per_term %in%special_terms
  
  mode_per_term <- lapply(formula_info, function(x) x$mode)
  mode_per_term <- unlist(lapply(mode_per_term, function(x){
    if(is.null(x)){"global"
    } else{
        if(!x %in% c("global","local", "alocal")){
          stop(paste0("Mode '", x, "' not recognized. Please use 'global', 'local' or 'alocal'."))
        } else x
      }
  }))
  term_per_term[is_special] <- paste0(term_per_term[is_special], "_", mode_per_term[is_special], sep = "")
  
  is_very_special <- term_per_term %in% c("gwdsp_global","gwdsp_local","gwesp_global", "gwesp_local")
  for(m in which(is_very_special)){
    formula_info[[m]]$data <- matrix(formula_info[[m]]$decay)
  }
  variant_per_term <- lapply(formula_info, function(x) x$variant)
  variant_per_term <- unlist(lapply(variant_per_term, function(x){
    if(is.null(x)){"OSP"
    } else{
      if(!x %in% c("ITP","ISP", "OTP", "OSP")){
        stop(paste0("Mode '", x, "' not recognized. Please use 'ITP','ISP', 'OTP' or 'OSP'."))
      } else x
    }
  }))
  term_per_term[is_very_special] <- paste0(term_per_term[is_very_special], "_", variant_per_term[is_very_special], sep = "")
  
  type_per_term <- lapply(formula_info, function(x) x$type)
  type_per_term <- unlist(lapply(type_per_term, function(x) if(is.null(x)){1} else{x}))
  
  data_per_term <-  lapply(formula_info, function(x) x$data)
  data_per_term <- lapply(data_per_term, function(x) if(is.null(x)){matrix(1)} else{x})
  name_per_term <- unlist(lapply(formula_info, function(x) x$label))
  return(list(data_object = data_object, 
              data_list = data_per_term,
              type_list = type_per_term,
              coef_names = name_per_term,
              term_names = term_per_term, 
              includes_degrees = includes_degrees))
  
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

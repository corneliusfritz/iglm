#' Compute Statistics
#'
#' @description
#' Computes statistics.
#'
#' @param formula A model `formula` object. The left-hand side should be the
#'   name of a \code{\link{iglm.data}} object available in the calling environment. 
#'   Alternatively, the left-hand side can be a \code{iglm.data.list} object to compute statistics
#'   for multiple \code{\link{iglm.data}} objects at once (is, e.g., the normal outcome of all simulations).
#'   See \code{\link{model.terms}} for details on specifying the right-hand side terms.

#'
#' @return A named numeric vector. Each element corresponds to a term in the
#'   `formula`, and its value is the calculated observed feature
#'   for that term based on the data in the \code{\link{iglm.data}} object. The names of the
#'   vector match the coefficient names derived from the formula terms.
#'
#' @examples
#' # Create a iglm.data object
#' n_actors = 10
#' neighborhood = matrix(1, nrow = n_actors, ncol = n_actors)
#' type_x <- "binomial"
#' type_y <- "binomial"
#' x_attr_data <- rbinom(n_actors, 1, 0.5)
#' y_attr_data <- rbinom(n_actors, 1, 0.5)
#' z_net_data  <- matrix(0, nrow = n_actors, ncol = n_actors) 
#' object = iglm.data(z_network = z_net_data, x_attribute = x_attr_data,
#'  y_attribute = y_attr_data, neighborhood = neighborhood,
#'  directed = FALSE,type_x = type_x,type_y = type_y)
#' count_statistics(object ~ edges(mode = "local") + attribute_y + attribute_x)
#' @export
count_statistics = function(formula) {
  tmp_obj <- eval(formula[[2]],envir = environment(formula))
  if(inherits(tmp_obj, "iglm.data.list")) {
    k <- length(tmp_obj)
    lhs_call <- formula[[2]]
    rhs_part <- formula[[3]]
    f_env <- environment(formula)
    formula_list <- lapply(1:k, function(i) {
      new_lhs <- call("[[", lhs_call, i)
      new_formula_call <- call("~", new_lhs, rhs_part)
      as.formula(new_formula_call, env = f_env)
    })
    counts <- lapply(formula_list, function(x){
      count_statistics(x)
    })
    counts <- do.call(counts,what = "rbind")
  } else{
    preprocessed = formula_preprocess(formula)
    counts = as.vector(xyz_count_statistics(preprocessed))
    names(counts) = preprocessed$coef_names
  }
  return(counts)
}

# xz_count_statistics = function(preprocessed, ...) {
#   n_actor = length(preprocessed$data_object$x_attribute)
#   xz_count_global(network = preprocessed$data_object$z_network,
#                   attribute =preprocessed$data_object$x_attribute,
#                   neighborhood = preprocessed$data_object$neighborhood,
#                   overlap=preprocessed$data_object$overlap,
#                   directed = preprocessed$data_object$directed,
#                   terms = preprocessed$term_names,
#                   type = preprocessed$data_object$type,
#                   data_list = preprocessed$data_list,
#                   type_list = preprocessed$type_list,
#                   n_actors = n_actors)
#   
#   
# }

xyz_count_statistics = function(preprocessed, ...) {
  n_actor_tmp = length(preprocessed$data_object$x_attribute)
  # browser()
  xyz_count_global(z_network =  preprocessed$data_object$z_network,
                   x_attribute = preprocessed$data_object$x_attribute,
                   y_attribute = preprocessed$data_object$y_attribute,
                   type_x = preprocessed$data_object$type_x,
                   type_y = preprocessed$data_object$type_y,
                   attr_x_scale = preprocessed$data_object$scale_x,
                   attr_y_scale = preprocessed$data_object$scale_y,
                   neighborhood = preprocessed$data_object$neighborhood,
                   overlap=preprocessed$data_object$overlap,
                   directed = preprocessed$data_object$directed,
                   terms = preprocessed$term_names,
                   data_list = preprocessed$data_list,
                   type_list = preprocessed$type_list,
                   n_actor = n_actor_tmp)
}

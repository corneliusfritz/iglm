#' Simulate responses and connections
#'
#' @description
#' Simulate responses and connections.
#'
#' @param formula A model `formula` object. The left-hand side should be the
#'   name of a `iglm.data` object available in the calling environment. 
#'   See \code{\link{model_terms}} for details on specifying the right-hand side terms.
#' @param coef Numeric vector containing the coefficient values for the structural
#'   (non-popularity) terms defined in the `formula`.
#' @param coef_popularity Numeric vector specifying the popularity coefficient
#'   values (expansiveness/attractiveness). This is required \strong{only if} the
#'   `formula` includes popularity terms. Its length must be `n_actor` (for
#'   undirected networks) or `2 * n_actor` (for directed networks), where
#'   `n_actor` is determined from the `iglm.data` object in the formula. 
#' @param sampler An object of class `sampler_iglm` (created by
#'   `sampler.iglm()`) specifying the MCMC sampling parameters. This includes
#'   the number of simulations (`n_simulation`), burn-in iterations (`n_burn_in`),
#'   initialization settings (`init_empty`), and component sampler settings
#'   (`sampler.x`, `sampler.y`, etc.).
#'   If `NULL` (default), default settings from `sampler.iglm()` are used.
#' @param only_stats (logical). If \code{TRUE} (default, consistent with the usage signature), the
#'   function returns only the matrix of features calculated
#'   for each simulation. The full simulated \code{iglm.data} objects are discarded to minimize
#'   memory usage. If \code{FALSE}, the complete simulated \code{iglm.data} objects are created
#'   and returned within the \code{samples} component of the output list.
#' @param display_progress Logical. If `TRUE`, progress messages or a progress
#'   bar (depending on the backend implementation) are displayed during simulation.
#'   Default is `FALSE`.
#' @param offset_nonoverlap Numeric scalar value passed to the C++ simulator.
#'   This value is typically added to the linear predictor for dyads that are
#'   \strong{not} part of the 'overlap' set defined in the `iglm.data` object, potentially
#'   modifying tie probabilities outside the primary neighborhood. Default is `0`.
#' @param cluster Optional parallel cluster object created, for example, by
#'   ``parallel::makeCluster``. If provided and valid, the function performs a
#'   single burn-in simulation on the main R process, then distributes the
#'   remaining `n_simulation` tasks across the cluster workers using
#'   ``parallel::parLapply``. Seeds for component samplers are offset
#'   for each worker to ensure different random streams. If `NULL` (default),
#'   all simulations are run sequentially in the main R process.
#' @param fix_x Logical. If `TRUE`, the simulation holds the `x_attribute` fixed
#'   at its initial state (from the \code{\link{iglm.data}} object) and only simulates the
#'   `y_attribute` and `z_network`. If `FALSE` (default), all components (x, y, z)
#'   are simulated according to the model and sampler settings. 
#'
#' @details
#'
#' \strong{Parallel Execution:} When a `cluster` object is provided, the simulation
#' process is adapted:
#' \enumerate{
#'   \item A single simulation run (including burn-in specified by `sampler$n_burn_in`)
#'     is performed on the master node to obtain a starting state for the parallel chains.
#'   \item The total number of requested simulations (`sampler$n_simulation`) is divided
#'     among the cluster workers.
#'   \item ``parallel::parLapply`` is used to run simulations on each worker.
#'     Each worker starts from the state obtained after the initial burn-in, performs
#'     \strong{zero} additional burn-in (`n_burn_in = 0` passed to workers), and generates
#'     its assigned share of the simulations. Component sampler seeds are offset
#'     based on the worker ID to ensure pseudo-independent random number streams.
#'   \item Results (simulated objects or statistics) from all workers are collected
#'     and combined.
#' }
#' This approach ensures that the initial burn-in phase happens only once, saving time.
#'
#' @return A list containing two components:
#' \describe{
#'   \item{`samples`}{If `only_stats = FALSE`, this is a list of length
#'     `sampler$n_simulation` where each element is a `iglm.data` object
#'     representing one simulated draw from the model. The list has the S3 class
#'     `"iglm.data.list"`. If `only_stats = TRUE`, this is typically an empty list.}
#'   \item{`stats`}{A numeric matrix with `sampler$n_simulation` rows and
#'     `length(coef)` columns. Each row contains the features
#'     (corresponding to the model terms in `formula`) calculated for one
#'     simulation draw. Column names are set to match the term names.}
#' }
#'
#' @section Errors:
#' The function stops with an error if:
#' \itemize{
#'   \item The length of `coef` does not match the number of terms derived from `formula`.
#'   \item `formula_preprocess` fails.
#'   \item The `sampler` object is not of class `sampler.iglm`.
#'   \item The C++ backend `xyz_simulate_cpp` encounters an error.
#'   \item Helper functions like `XYZ_to_R` or `is_cluster_active` are not found.
#' }
#' Warnings may be issued if default sampler settings are used.
#'
#' @seealso \code{iglm} for creating the model object,
#'   \code{sampler.iglm} for creating the sampler object,
#'   \code{iglm.data} for the data object structure. 
#'
#' @export
#' @importFrom parallel parLapply
simulate_iglm = function(formula,coef,coef_popularity = NULL, 
                        sampler = NULL, 
                        only_stats = T, 
                        display_progress = F, 
                        offset_nonoverlap = 0,
                        cluster = NULL, 
                        fix_x = F) {
  if(is.null(sampler)){
    sampler= sampler.iglm()
    # if no specifications of the sampler are provided use the default one
  }
  if("sampler.iglm" %in% class(sampler)){
    sampler= sampler.iglm()
  }
  # Search for all things in the environment of the formula
  # attr(formula, ".Environment")
  preprocessed = formula_preprocess(formula) 
  popularity <- preprocessed$includes_popularity
  n_actor = length(preprocessed$data_object$x_attribute)
  if(length(coef) != length(preprocessed$term_names)){
    return("Wrong number of coefficients for the wanted terms.")
  }
  
  if(!is_cluster_active(cluster)){
    cluster = NULL
  }
  if(is.null(cluster)){
    res = xyz_simulate_cpp(coef = coef,coef_popularity = coef_popularity,
                           terms= preprocessed$term_names,
                           n_actor = n_actor,
                           x_attribute=preprocessed$data_object$x_attribute,
                           y_attribute=preprocessed$data_object$y_attribute,
                           z_network=preprocessed$data_object$z_network,
                           type_x =  preprocessed$data_object$type_x,
                           type_y =  preprocessed$data_object$type_y,
                             attr_x_scale =  preprocessed$data_object$scale_x,
                           attr_y_scale =  preprocessed$data_object$scale_y,
                           init_empty = sampler$init_empty,
                           neighborhood=preprocessed$data_object$neighborhood,
                           overlap=preprocessed$data_object$overlap,
                           directed =preprocessed$data_object$directed,
                           data_list = preprocessed$data_list,
                           type_list = preprocessed$type_list,
                           n_burn_in = sampler$n_burn_in,
                           seed_x = sampler$sampler.x$seed,
                           n_proposals_x = sampler$sampler.x$n_proposals,
                           seed_y = sampler$sampler.y$seed,
                           n_proposals_y = sampler$sampler.y$n_proposals,
                           seed_z = sampler$sampler.z$seed,
                           n_proposals_z = sampler$sampler.z$n_proposals,
                           n_simulation = sampler$n_simulation,
                           only_stats =only_stats,
                           display_progress = display_progress, 
                           popularity = popularity, 
                           offset_nonoverlap = offset_nonoverlap,
                           fix_x = fix_x)  
  } else {
    if(display_progress){
      cat("Starting with burn-in\n")
    }
    
    res_burn_in = xyz_simulate_cpp(coef = coef,
                                   coef_popularity = coef_popularity,
                                   terms= preprocessed$term_names,
                                   n_actor = n_actor,
                                   x_attribute=preprocessed$data_object$x_attribute,
                                   y_attribute=preprocessed$data_object$y_attribute,
                                   z_network=preprocessed$data_object$z_network,
                                   init_empty = sampler$init_empty,
                                   neighborhood= preprocessed$data_object$neighborhood,
                                   type_x = preprocessed$data_object$type_x, 
                                   type_y = preprocessed$data_object$type_y, 
                                   attr_x_scale =  preprocessed$data_object$scale_x,
                                   attr_y_scale =  preprocessed$data_object$scale_y,
                                   overlap=preprocessed$data_object$overlap,
                                   directed =preprocessed$data_object$directed,
                                   data_list = preprocessed$data_list,
                                   type_list = preprocessed$type_list,
                                   n_burn_in = sampler$n_burn_in,
                                   seed_x = sampler$sampler.x$seed,
                                   n_proposals_x = sampler$sampler.x$n_proposals,
                                   seed_y = sampler$sampler.y$seed,
                                   n_proposals_y = sampler$sampler.y$n_proposals,
                                   seed_z = sampler$sampler.z$seed,
                                   n_proposals_z = sampler$sampler.z$n_proposals,
                                   n_simulation = 1,
                                   only_stats =F,
                                   display_progress = display_progress, 
                                   popularity = popularity, 
                                   offset_nonoverlap = offset_nonoverlap, fix_x = fix_x)
    res_burnin = XYZ_to_R(x_attribute = res_burn_in$simulation_attributes_x[[1]], 
                          y_attribute = res_burn_in$simulation_attributes_y[[1]], 
                          z_network = res_burn_in$simulation_networks_z[[1]],
                          n_actor = n_actor, return_adj_mat = FALSE)
    # tmp_split = split(1:(round(sampler$n_simulation/length(cluster))*length(cluster)), 1:length(cluster))
    tmp_split = suppressWarnings(split(1:sampler$n_simulation, 1:length(cluster)))
    
    if(display_progress){
      cat("Starting with parallel simulations")
    }
    
    res_parallel = parLapply(cl = cluster, X = tmp_split, fun = function(x, preprocessed, n_actor, coef, 
                                                                         coef_popularity, popularity, sampler, 
                                                                         res_burnin, offset_nonoverlap){
      xyz_simulate_cpp(coef = coef,coef_popularity = coef_popularity,
                       terms= preprocessed$term_names,
                       n_actor = n_actor,
                       type_x = preprocessed$data_object$type_x, 
                       type_y = preprocessed$data_object$type_y, 
                       attr_x_scale =  preprocessed$data_object$scale_x,
                       attr_y_scale =  preprocessed$data_object$scale_y,
                       x_attribute=res_burnin$x_attribute,
                       y_attribute=res_burnin$y_attribute,
                       z_network=res_burnin$z_network,
                       init_empty = sampler$init_empty,
                       neighborhood=preprocessed$data_object$neighborhood,
                       overlap=preprocessed$data_object$overlap,
                       directed =preprocessed$data_object$directed,
                       data_list = preprocessed$data_list,
                       type_list = preprocessed$type_list,
                       n_burn_in = sampler$n_burn_in,
                       seed_x = sampler$sampler.x$seed + min(x),
                       n_proposals_x = sampler$sampler.x$n_proposals,
                       seed_y = sampler$sampler.y$seed + 2*min(x),
                       n_proposals_y = sampler$sampler.y$n_proposals,
                       seed_z = sampler$sampler.z$seed +  3*min(x),
                       n_proposals_z = sampler$sampler.z$n_proposals,
                       n_simulation = length(x),
                       only_stats =only_stats,
                       display_progress = F, 
                       popularity = popularity, fix_x = fix_x,
                       offset_nonoverlap = offset_nonoverlap)
    },preprocessed = preprocessed, n_actor = n_actor, coef = coef, 
    coef_popularity = coef_popularity, popularity = popularity, 
    sampler = sampler, res_burnin = res_burnin, offset_nonoverlap = offset_nonoverlap)
    
    res = list()
    res$simulation_attributes_x = unlist(lapply(res_parallel, function(x){x$simulation_attributes_x}),recursive = F)
    res$simulation_attributes_y = unlist(lapply(res_parallel, function(x){x$simulation_attributes_y}),recursive = F)
    res$simulation_networks = unlist(lapply(res_parallel, function(x){x$simulation_networks}),recursive = F)
    res$stats = do.call(rbind,lapply(res_parallel, function(x){x$stats}))
  }
  tmp = lapply(1:length(res$simulation_networks),
               function(x){XYZ_to_R(x_attribute = res$simulation_attributes_x[[x]], 
                                    y_attribute = res$simulation_attributes_y[[x]],
                                    z_network = res$simulation_networks[[x]],
                                    n_actor = n_actor, return_adj_mat = FALSE)})
  tmp = lapply(1:length(res$simulation_networks),
               function(x){iglm.data(x_attribute = tmp[[x]]$x_attribute,
                                   y_attribute = tmp[[x]]$y_attribute,
                                   z_network = tmp[[x]]$z_network, 
                                   n_actor = length(tmp[[x]]$x_attribute),
                                   type_x = preprocessed$data_object$type_x, 
                                   type_y = preprocessed$data_object$type_y, 
                                   scale_x = preprocessed$data_object$scale_x, 
                                   scale_y = preprocessed$data_object$scale_y, 
                                   return_neighborhood = FALSE)})
  # browser()
  class(tmp) <- "iglm.data.list"
  attr(tmp, "neighborhood") <- iglm.data.neighborhood(preprocessed$data_object$neighborhood)
  colnames(res$stats) <- preprocessed$coef_names
  return(list(samples = tmp, stats = res$stats))
}

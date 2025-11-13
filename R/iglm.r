#' @docType class
#' @title An R6 class for Network GLM (Generalized Linear Model) Objects
#' @description
#' The `iglm_object` class encapsulates all components required to define,
#' estimate, and simulate from a network generalized linear model. This includes
#' the model formula, coefficients, the underlying network and attribute data
#' (via a `iglm.data` object), sampler controls, estimation controls, and storage
#' for results.
#' @references 
#' Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A Regression Framework for Studying Relationships among Attributes under Network Interference. Journal of the American Statistical Association, to appear.
#' 
#' Stewart, J. R. and M. Schweinberger (2025). Pseudo-Likelihood-Based M-Estimation of Random Graphs with Dependent Edges and Parameter Vectors of Increasing Dimension. Annals of Statistics, to appear. 
#' 
#' Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.
#' 
#' @importFrom R6 R6Class
#' @importFrom methods as
#' @importFrom stats as.formula pnorm terms quantile update printCoefmat na.omit sd
#' @importFrom utils modifyList
#' @export
iglm_object_generator <- R6::R6Class("iglm_object",
                                     private = list(
                                       .formula = NULL,
                                       .iglm.data = NULL,
                                       .coef = NULL,
                                       .coef_popularity = NULL,
                                       .coef_popularity_internal = NULL,
                                       .sampler = NULL,
                                       .control = NULL,
                                       .preprocess = NULL,
                                       .time_estimation = NULL,
                                       .count_statistics = NULL,
                                       .results = list(),
                                       #' @description
                                       #' Internal method to calculate the observed count statistics based on the
                                       #' model formula and the data in the `iglm.data` object. Populates the
                                       #' `private$.count_statistics` field.
                                       .calc_count_statistics = function() {
                                         counts <- as.vector(xyz_count_global(z_network =  private$.iglm.data$z_network,
                                                                              x_attribute =private$.iglm.data$x_attribute,
                                                                              y_attribute = private$.iglm.data$y_attribute,
                                                                              type_x = private$.iglm.data$type_x,
                                                                              type_y = private$.iglm.data$type_y,
                                                                              attr_x_scale = private$.iglm.data$scale_x,
                                                                              attr_y_scale = private$.iglm.data$scale_y,
                                                                              neighborhood = private$.iglm.data$neighborhood,
                                                                              overlap=private$.iglm.data$overlap,
                                                                              directed = private$.iglm.data$directed,
                                                                              terms = private$.preprocess$term_names,
                                                                              data_list = private$.preprocess$data_list,
                                                                              type_list = private$.preprocess$type_list,
                                                                              n_actor = private$.iglm.data$n_actor))
                                         names(counts) <- private$.preprocess$coef_names
                                         private$.count_statistics <- counts
                                       },
                                       #' @description
                                       #' Internal validation method. Checks the consistency and validity of
                                       #' all components of the `iglm_object`. Stops with an error if any
                                       #' check fails.
                                       .validate = function() {
                                         # browser()
                                         if(is.null(private$.formula) || !inherits(private$.formula, "formula")) {
                                           stop("Invalid formula in iglm_object.", call. = FALSE)
                                         }
                                         if(is.null(private$.coef) || !is.numeric(private$.coef)) {
                                           stop("Invalid coef in iglm_object.", call. = FALSE)
                                         }
                                         if(length(private$.coef) != length(private$.preprocess$coef_names)) {
                                           stop("Length of coef does not match number of terms in formula.", call. = FALSE)
                                         }
                                         if(!inherits(private$.sampler, "sampler_iglm")) {
                                           stop("Invalid sampler in iglm_object.", call. = FALSE)
                                         }
                                         
                                         if(!is.null(private$.coef_popularity)){
                                           if(!is.numeric(private$.coef_popularity)) {
                                             stop("Invalid coef_popularity in iglm_object.", call. = FALSE)
                                           }
                                           expected_length <- private$.iglm.data$n_actor +
                                             private$.iglm.data$directed*private$.iglm.data$n_actor
                                           if(length(private$.coef_popularity) != expected_length) {
                                             stop("Length of coef_popularity does not match number of actors in data object.", call. = FALSE)
                                           }
                                           if(length(private$.coef_popularity_internal) != expected_length) {
                                             stop("Length of coef_popularity does not match number of actors in data object.", call. = FALSE)
                                           }
                                         }
                                         if(!is.null(private$.results)){
                                           if(!inherits(private$.results, "results")) {
                                             stop("Invalid results object in iglm_object.", call. = FALSE)
                                           }
                                         }
                                         if(!inherits(private$.iglm.data, "iglm.data")) {
                                           stop("Invalid iglm.data object in iglm_object.", call. = FALSE)
                                         }
                                         if(!inherits(private$.control, "control.iglm")) {
                                           stop("Invalid control object in iglm_object.", call. = FALSE)
                                         }
                                         
                                         if(is.null(private$.count_statistics)){
                                           stop("Count statistics have not been computed yet.", call. = FALSE)
                                         }
                                       }
                                     ),
                                     public = list(
                                       #' @description
                                       #' Creates a new `iglm_object`. This involves parsing the formula,
                                       #' linking the data object, initializing coefficients, setting up sampler
                                       #' and control objects, calculating initial statistics, and validating.
                                       #' @param formula A model `formula` object. The left-hand side should be the
                                       #'   name of a \code{\link{iglm.data}} object available in the calling environment. 
                                       #'   See \code{\link{model_terms}} for details on specifying the right-hand side terms.
                                       #' @param coef A numeric vector of initial coefficients for the terms in
                                       #'   the formula (excluding popularity). If `NULL`, coefficients are
                                       #'   initialized to zero.
                                       #' @param coef_popularity An optional numeric vector of initial popularity
                                       #'   coefficients. Should be `NULL` if the formula does
                                       #'   not include popularity terms.
                                       #' @param sampler A \code{\link{sampler_iglm}} object specifying the MCMC sampler
                                       #'   settings. If `NULL`, default settings are used. 
                                       #' @param control A \code{\link{control.iglm}} object specifying estimation control
                                       #'   parameters. If `NULL`, default settings are used.
                                       #' @param file (character or `NULL`) If provided, loads the sampler state from
                                       #'  the specified .rds file instead of initializing from parameters.
                                       #' @return A new `iglm_object`.
                                       initialize = function(formula = NULL, coef = NULL, 
                                                             coef_popularity = NULL, 
                                                             sampler = NULL, 
                                                             control = NULL, file = NULL) {
                                         if(!is.null(file)){
                                           if (missing(file) || !is.character(file) || length(file) != 1) {
                                             stop("A valid 'file' (character string) must be provided.", call. = FALSE)
                                           }
                                           data_loaded <- readRDS(file)
                                           required_fields <- c("formula", "preprocess", "iglm.data", "coef", 
                                                                "coef_popularity", "time_estimation", 
                                                                "count_statistics", "results", 
                                                                "control", "sampler")
                                           if (!is.list(data_loaded) || !all(required_fields %in% names(data_loaded))) {
                                             stop("File does not contain a valid iglm state.", call. = FALSE)
                                           }
                                           private$.formula <- data_loaded$formula
                                           private$.preprocess <- data_loaded$preprocess
                                           private$.iglm.data <- iglm.data(x_attribute = data_loaded$iglm.data$x_attribute,
                                                                           y_attribute = data_loaded$iglm.data$y_attribute,
                                                                           z_network = data_loaded$iglm.data$z_network,
                                                                           n_actor = data_loaded$iglm.data$n_actor,
                                                                           type_x = data_loaded$iglm.data$type_x,
                                                                           type_y = data_loaded$iglm.data$type_y,
                                                                           scale_x = data_loaded$iglm.data$scale_x,
                                                                           scale_y = data_loaded$iglm.data$scale_y,
                                                                           neighborhood = data_loaded$iglm.data$neighborhood,
                                                                           directed = data_loaded$iglm.data$directed)
                                           private$.coef <- data_loaded$coef
                                           private$.coef_popularity <- data_loaded$coef_popularity
                                           private$.coef_popularity_internal <- private$.coef_popularity
                                           private$.time_estimation <- data_loaded$time_estimation
                                           private$.count_statistics <- data_loaded$count_statistics
                                           
                                           private$.results <- results(size_coef = length(private$.coef), 
                                                                       size_coef_popularity =length(private$.coef_popularity)*
                                                                         private$.preprocess$includes_popularity)
                                           private$.results$update(coefficients_path = data_loaded$results$coefficients_path, 
                                                                   var = data_loaded$results$var,
                                                                   fisher_popularity = data_loaded$results$fisher_popularity,
                                                                   fisher_nonpopularity = data_loaded$results$fisher_nonpopularity,
                                                                   score_popularity = data_loaded$results$score_popularity,
                                                                   score_nonpopularity = data_loaded$results$score_nonpopularity,
                                                                   llh = data_loaded$results$llh,
                                                                   estimated = data_loaded$results$estimated)
                                           if(length(data_loaded$results$samples) > 0){
                                             private$.results$update(stats = data_loaded$results$stats,
                                                                     samples = data_loaded$results$samples)
                                           }
                                           if(length(data_loaded$results$model_assessment) > 0){
                                             private$.results$set_model_assessment(res = 
                                                                                     data_loaded$results$model_assessment)  
                                           }
                                           private$.control <- data_loaded$control
                                           private$.sampler <- data_loaded$sampler
                                         } else{
                                           if (!inherits(formula, "formula")) {
                                             stop("`formula` must be a formula object.")
                                           }
                                           private$.formula <- formula
                                           private$.preprocess <- formula_preprocess(formula)
                                           private$.iglm.data <- private$.preprocess$data_object
                                           private$.preprocess$data_object <- NULL
                                           
                                           if(is.null(coef)){
                                             private$.coef <- rep(0, length(private$.preprocess$coef_names))
                                           } else {
                                             private$.coef <- coef
                                           }
                                           if(is.null(coef_popularity)){
                                             private$.coef_popularity <- rep(0, 
                                                                             private$.iglm.data$n_actor +
                                                                               private$.iglm.data$directed*private$.iglm.data$n_actor)
                                           } else {
                                             if(private$.preprocess$includes_popularity == FALSE){
                                               stop("The formula does not include popularity terms, so `coef_popularity` should be NULL.")
                                             }
                                             private$.coef_popularity <- coef_popularity
                                           }
                                           private$.coef_popularity_internal <- private$.coef_popularity
                                           private$.results <- results(size_coef = length(private$.coef), 
                                                                       size_coef_popularity =length(private$.coef_popularity)*
                                                                         private$.preprocess$includes_popularity)
                                           
                                           if (is.null(control)) {
                                             private$.control <- control.iglm() 
                                           } else {
                                             private$.control <- control
                                           }
                                           # --- Handle sampler ---
                                           if (is.null(sampler)) {
                                             private$.sampler <- sampler.iglm() 
                                           } else {
                                             if (!inherits(sampler, "sampler_iglm")) {
                                               stop("`sampler` must be a 'sampler_iglm' object created with sampler.iglm().")
                                             }
                                             private$.sampler <- sampler
                                           }
                                           private$.calc_count_statistics()
                                         }
                                         private$.validate()
                                         invisible(self)
                                       },
                                       #' @description
                                       #' Performs model assessment by calculating specified network statistics
                                       #' on the observed network and comparing their distribution to the
                                       #' distribution obtained from simulated networks based on the current
                                       #' model parameters. Requires simulations to have been run first (via
                                       #'  \code{iglm_object$simulate} or \code{iglm_object_generator$estimate}).
                                       #'
                                       #' @param formula A formula specifying the network statistics to assess
                                       #'   (e.g., `~ degree_distribution() + geodesic_distances_distribution()`).
                                       #'   The terms should correspond to methods available in the \code{\link{iglm.data}} 
                                       #'   object that end with `distributions`. 
                                       #'   If the term mcmc_diagnostics is included, MCMC diagnostics will also be computed.
                                       #' @return An object of class `iglm_model_assessment` containing the
                                       #'   observed statistics and the distribution of simulated statistics.
                                       #'   The result is also stored internally.
                                       model_assessment = function(formula){
                                         if(length(private$.results$samples) ==0){
                                           self$simulate()
                                         }
                                         names_tmp <- attr(terms(formula),"term.labels")
                                         if("mcmc_diagnostics" %in% names_tmp){
                                           names_tmp <- names_tmp[names_tmp != "mcmc_diagnostics"]
                                           formula <- update(formula, . ~ . - mcmc_diagnostics)
                                           count_statistics <- private$.count_statistics
                                           include_mcmc <- TRUE
                                         } else {
                                           include_mcmc <- FALSE
                                           count_statistics <- NULL
                                         }
                                         names_tmp <- gsub("\"", "", names_tmp)
                                         names_tmp <- gsub("\\(", "_", names_tmp)
                                         names_tmp <- gsub("\\)", "", names_tmp)
                                         names_tmp <- gsub("=", "_", names_tmp)
                                         names_tmp <- gsub(" ", "", names_tmp)
                                         observed <- eval_change(formula = formula,object = private$.iglm.data)
                                         names(observed) = names_tmp
                                         base_name <- unlist(lapply(rhs_terms_as_list(update(formula, a ~ .)), function(x){
                                           x$base_name
                                         }))
                                         ranges_tmp <- lapply(observed, function(x) {
                                           if(is.numeric(x)){
                                             return(range(as.numeric(names(x))[is.finite(as.numeric(names(x)))]))
                                           } else if(is.list(x)){
                                             return(range(as.numeric(unlist(lapply(x, function(y){names(y)})))))
                                           } else {
                                             return(NULL)
                                           }
                                         })
                                         names(ranges_tmp) <- rep("value_range", length(names_tmp))
                                         simulated <- lapply(private$.results$samples, 
                                                             function(object, info_tmp, names_tmp) {
                                                               res <- eval_change(formula = formula,object = object, additional_args = ranges_tmp)
                                                               names(res) = names_tmp
                                                               return(res)
                                                             }, names_tmp = names_tmp, info_tmp = info_tmp
                                         )
                                         
                                         res <- list(observed = observed, 
                                                     simulated = simulated, 
                                                     base_name = base_name,
                                                     names = names_tmp, 
                                                     include_mcmc = include_mcmc, 
                                                     count_statistics = count_statistics)
                                         class(res) <- "iglm_model_assessment"
                                         private$.results$set_model_assessment(res)
                                         invisible(res)
                                       }, 
                                       #' @description
                                       #' Print a summary of the `iglm_object`. If estimation results are
                                       #' available, they are printed in a standard coefficient table format.
                                       #' @param digits (integer) Number of digits for rounding numeric output.
                                       #' @param ... Additional arguments (not used).
                                       print = function(digits = 4, ...) {
                                         cat("iglm object\n")
                                         cat(strrep("-", 50), "\n", sep = "")
                                         
                                         # --- Display formula
                                         cat("Formula:\n  ", deparse(private$.formula), "\n", sep = "")
                                         # --- Results 
                                         if(nrow(private$.results$coefficients_path) >0){
                                           # cat("\n")
                                           cat("Results: \n\n")
                                           names = rownames(private$.coef)
                                           est = as.vector(private$.coef)
                                           stderr <- sqrt(diag(private$.results$var))
                                           tvalue <- est / stderr
                                           pvalue <- 2 * pnorm(-abs(tvalue))
                                           
                                           coef_table <- cbind(est, stderr, tvalue, pvalue, private$.count_statistics)
                                           # Check if we can add space between the columns
                                           colnames(coef_table) <- c("Estimate", "SE", "t-value", "Pr(>|t|)", "Count")
                                           rownames(coef_table) <- names
                                           coef_table <- round(coef_table, digits)
                                           eps_threshold <- 10^(-digits)
                                           which_wrong <- coef_table == 0
                                           which_wrong[,4] <- FALSE
                                           coef_table[coef_table == 0] <- 
                                             paste0("< ", format(eps_threshold,scientific = FALSE))
                                           print(coef_table, quote = FALSE, right = TRUE)
                                           
                                           cat(paste("\nTime for estimation: ",
                                                     round(as.numeric(private$.time_estimation),3),
                                                     attr(private$.time_estimation, "units")))
                                         } else {
                                           cat("\n")
                                           cat("Observed Sufficient Statistics:\n")
                                           print(round(private$.count_statistics, digits))
                                         }
                                       }, 
                                       #' @description
                                       #' Plot the estimation results, including coefficient convergence
                                       #' paths and model assessment diagnostics if available.
                                       #' @param stats (logical) If `TRUE`, plot the observed vs. simulated
                                       #'  statistics from model assessment. Default is `FALSE`.
                                       #' @param trace (logical) If `TRUE`, plot the coefficient convergence
                                       #'  paths. Default is `FALSE`.
                                       #' @param model_assessment (logical) If `TRUE`, plot diagnostics from the
                                       #'  model assessment (if already carried out). Default is `FALSE`.
                                       plot = function(stats = FALSE, trace = FALSE, model_assessment = FALSE) {
                                         private$.results$plot(stats = stats, trace = trace,model_assessment =  model_assessment)
                                       },
                                       #' @description
                                       #' Gathers all components of the \code{\link{iglm_object}} into a single list for
                                       #' easy saving or inspection.
                                       #' @return A list containing all key components of the \code{\link{iglm_object}}.
                                       #'   This includes the formula, coefficients, sampler, control settings,
                                       #'   preprocessing info, time taken for estimation, count statistics,
                                       #'   results, and the underlying \code{\link{iglm.data}} data object.
                                       gather = function(){
                                         list(formula = private$.formula,
                                              coef = private$.coef,
                                              coef_popularity = private$.coef_popularity,
                                              sampler = private$.sampler,
                                              control = private$.control,
                                              preprocess = private$.preprocess,
                                              time_estimation = private$.time_estimation,
                                              count_statistics = private$.count_statistics,
                                              results = private$.results$gather(),
                                              iglm.data = private$.iglm.data$gather())
                                       }, 
                                       #' @description
                                       #' Save the \code{\link{iglm_object}} to a file in RDS format.
                                       #' @param file (character) File path to save the object to.
                                       #' @return Invisibly returns `NULL`.
                                       save = function(file = NULL) {
                                         if (missing(file) || !is.character(file) || length(file) != 1) {
                                           stop("A valid 'file' (character string) must be provided.", call. = FALSE)
                                         }
                                         data_to_save <- self$gather()
                                         saveRDS(data_to_save, file = file)
                                         message(paste("Object state saved to:", file))
                                       },
                                       #' @description
                                       #' Estimate the model parameters using the specified control settings.
                                       #' Stores the results internally and updates the coefficient fields.
                                       #' @return If the no preprocessing should be returned 
                                       #'   (as per control settings), this function returns a list containing detailed estimation results, invisibly.
                                       #'   Includes final coefficients, variance-covariance matrix, convergence
                                       #'   path, Fisher information, score vector, log-likelihood, and any
                                       #'   simulations performed during estimation.
                                       #'   Else, the function returns a list of the desired preprocessed data (as a data.frame) and needed time. 
                                       estimate = function(){
                                         # private$.coef_popularity[is.na(private$.coef_popularity)] <- -50
                                         # browser()
                                         now <- Sys.time()
                                         info <-  estimate_xyz(formula = private$.formula,
                                                               preprocessed = private$.preprocess, 
                                                               control = private$.control,
                                                               sampler = private$.sampler,
                                                               beg_coef =  private$.coef,
                                                               beg_coef_popularity = private$.coef_popularity_internal, 
                                                               data_object = private$.iglm.data, 
                                                               start = nrow(private$.results$coefficients_path))
                                         
                                         private$.time_estimation <- Sys.time() - now
                                         if(private$.control$estimate_model){
                                           if(private$.preprocess$includes_popularity){
                                             coefficients_popularity_internal <- info$coefficients_popularity
                                             tmp <- private$.iglm.data$degree()
                                             if(private$.iglm.data$directed){
                                               coefficients_popularity_internal[which(tmp$in_degree_seq ==0) + private$.iglm.data$n_actor] = NA
                                               coefficients_popularity_internal[which(tmp$out_degree_seq ==0)] = NA
                                             } else {
                                               coefficients_popularity_internal[which(tmp$degree_seq ==0)] = NA
                                             }
                                           }
                                           
                                           if(is.null(info$var)){
                                             if(private$.preprocess$includes_popularity){
                                               info$var <- solve(info$fisher_nonpopularity)  
                                             } else {
                                               info$var <- solve(info$fisher)
                                             }
                                           }
                                           if(length(info$where_wrong)>0){
                                             if(private$.preprocess$includes_popularity){
                                               private$.results$resize(size_coef = length(info$coefficients_nonpopularity), 
                                                                       size_coef_popularity = 
                                                                         (private$.iglm.data$n_actor + private$.iglm.data$n_actor* private$.iglm.data$directed)*
                                                                         private$.preprocess$includes_popularity)
                                               
                                             } else {
                                               private$.results$resize(size_coef = length(info$coefficients), 
                                                                       size_coef_popularity = 
                                                                         (private$.iglm.data$n_actor + private$.iglm.data$n_actor* private$.iglm.data$directed)*
                                                                         private$.preprocess$includes_popularity)
                                             }
                                             private$.formula <- update_formula_remove_terms(private$.formula, 
                                                                                             private$.preprocess$coef_names[info$where_wrong+1])
                                             private$.preprocess$coef_names <- private$.preprocess$coef_names[-(info$where_wrong+1)]
                                             
                                             
                                           }
                                           
                                           if(length(info$simulations) > 0){
                                             tmp = lapply(1:length(info$simulations),
                                                          function(x){iglm.data(x_attribute = info$simulations[[x]]$x_attribute,
                                                                                y_attribute = info$simulations[[x]]$y_attribute,
                                                                                z_network = info$simulations[[x]]$z_network, 
                                                                                n_actor = private$.iglm.data$n_actor,
                                                                                return_neighborhood = FALSE,
                                                                                directed = private$.iglm.data$directed,
                                                                                type_x = private$.iglm.data$type_x, 
                                                                                type_y = private$.iglm.data$type_y, 
                                                                                scale_x = private$.iglm.data$scale_x, 
                                                                                scale_y = private$.iglm.data$scale_y)})
                                             class(tmp) <- "iglm.data.list"
                                             attr(tmp, "neighborhood") <- iglm.data.neighborhood(private$.iglm.data$neighborhood)
                                             info$simulations <- tmp
                                             colnames(info$stats) <- private$.preprocess$coef_names
                                             
                                             
                                           }
                                           
                                           # Remove samples because they would not anymore be consistent with the currect estimates
                                           private$.results$remove_samples()
                                           # Update the internal results object
                                           if(private$.preprocess$includes_popularity){
                                             private$.results$update(samples = info$simulations, 
                                                                     var = info$var, 
                                                                     coefficients_path = info$coefficients_path, 
                                                                     fisher_popularity = info$fisher_popularity,
                                                                     fisher_nonpopularity = info$fisher_nonpopularity,
                                                                     score_popularity = info$score_popularity,
                                                                     score_nonpopularity = info$score_nonpopularity,
                                                                     llh = info$llh, 
                                                                     stats = info$stats, 
                                                                     estimated = TRUE)
                                             
                                           } else {
                                             private$.results$update(samples = info$simulations, 
                                                                     var = info$var, 
                                                                     coefficients_path = info$coefficients_path, 
                                                                     fisher_nonpopularity = info$fisher,
                                                                     score_nonpopularity = info$score,
                                                                     llh = info$llh, 
                                                                     stats = info$stats, 
                                                                     estimated = TRUE)
                                             
                                           }
                                           
                                           if(private$.preprocess$includes_popularity){
                                             private$.coef <- info$coefficients_nonpopularity
                                             private$.coef_popularity <- coefficients_popularity_internal
                                             private$.coef_popularity_internal <- info$coefficients_popularity  
                                           } else {
                                             private$.coef <- info$coefficients  
                                           }
                                           private$.validate()
                                           
                                           if(private$.control$display_progress){
                                             cat("\nResults: \n\n")
                                             if(private$.preprocess$includes_popularity){
                                               names = rownames(info$coefficients_nonpopularity)
                                               est = as.vector(info$coefficients_nonpopularity)
                                               stderr <- sqrt(diag(info$var))
                                               coef_table <- cbind(est, stderr)
                                               colnames(coef_table) <- c("Estimate", "Std. Error")
                                               rownames(coef_table) <- names
                                               print(round(coef_table, 3), row.names = TRUE)
                                             } else {
                                               names = rownames(info$coefficients)
                                               est = as.vector(info$coefficients)
                                               stderr <- sqrt(diag(info$var))
                                               coef_table <- cbind(est, stderr)
                                               colnames(coef_table) <- c("Estimate", "Std. Error")
                                               rownames(coef_table) <- names
                                               print(round(coef_table, 3), row.names = TRUE)
                                             }
                                           }
                                         } else {
                                           if(private$.control$display_progress){
                                             cat("\nEstimation skipped as per control settings.\n")
                                           }
                                         }
                                         if(!is.null(info$preprocess)){
                                           if(private$.control$display_progress){
                                             cat("\nPreprocessing is returned from estimation.\n")
                                           }
                                           return(info$preprocess)
                                         }
                                         
                                         invisible(info)
                                       },
                                       #' @description
                                       #' Provides a summary of the estimation results.
                                       #' Requires the model to have been estimated first.
                                       #' @param digits (integer) Number of digits for rounding numeric output.
                                       #' @return Prints the summary to the console and returns `NULL` invisibly.
                                       summary = function(digits = 3){
                                         if(nrow(private$.results$coefficients_path) ==0){
                                           stop("No estimation results available. Please run `estimate()` first.", call. = FALSE)
                                         }
                                         cat("Results: \n\n")
                                         names = rownames(private$.coef)
                                         est = as.vector(private$.coef)
                                         
                                         if(is.null(private$.results$var)){
                                           stderr <- sqrt(diag(solve(private$.results$fisher_nonpopularity)))
                                         } else {
                                           stderr <- sqrt(diag(private$.results$var))
                                         }
                                         tvalue <- est / stderr
                                         pvalue <- 2 * pnorm(-abs(tvalue))
                                         
                                         coef_table <- cbind(est, stderr, tvalue, pvalue, private$.count_statistics)
                                         # Check if we can add space between the columns
                                         colnames(coef_table) <- c("Estimate", "SE", "t-value", "Pr(>|t|)", "Count")
                                         rownames(coef_table) <- names
                                         coef_table <- round(coef_table, digits)
                                         eps_threshold <- 10^(-digits)
                                         which_wrong <- coef_table == 0
                                         which_wrong[,4] <- FALSE
                                         coef_table[coef_table == 0] <- 
                                           paste0("< ", format(eps_threshold,scientific = FALSE))
                                         print(coef_table, quote = FALSE, right = TRUE)
                                         
                                         cat(paste("\nTime for estimation: ",
                                                   round(as.numeric(private$.time_estimation),3),
                                                   attr(private$.time_estimation, "units")))
                                       },
                                       #' @description
                                       #' Simulate networks from the fitted model or a specified model. Stores
                                       #' the simulations and/or summary statistics internally. The simulation 
                                       #' is carried out using the internal MCMC sampler described in \code{\link{simulate_iglm}}.
                                       #' @param nsim (integer) Number of networks to simulate. Default is 1.
                                       #' @param only_stats (logical) If `TRUE`, only calculate and store summary
                                       #'   statistics for each simulation, discarding the network object itself.
                                       #'   Default is `FALSE`.
                                       #' @param display_progress (logical) If `TRUE` (default), display a
                                       #'   progress bar during simulation.
                                       #' @param offset_nonoverlap (numeric) Offset to apply for non-overlapping
                                       #'   dyads during simulation (if applicable to the sampler). This option 
                                       #'   is useful if the sparsity of edges of units with non-overlapping 
                                       #'   neighborhoods is known. Default is 0.
                                       #' @return A list containing the simulated networks (`samples`, as a
                                       #'   `iglm.data.list` if `only_stats = FALSE`) and/or their summary
                                       #'   statistics (`stats`), invisibly.
                                       simulate = function (nsim = 1, only_stats = FALSE, display_progress=TRUE,
                                                            offset_nonoverlap= 0) {
                                         info <- simulate_iglm(formula = private$.formula, coef = private$.coef, 
                                                               coef_popularity = private$.coef_popularity,
                                                               sampler = private$.sampler, 
                                                               only_stats = only_stats, 
                                                               fix_x = private$.control$fix_x,
                                                               display_progress = display_progress, 
                                                               offset_nonoverlap = offset_nonoverlap, 
                                                               cluster = private$.sampler$cluster)
                                         
                                         private$.results$update(samples = info$samples,
                                                                 stats = info$stats,
                                                                 estimated = length(private$.coef)>0 )
                                         private$.validate()
                                         invisible(private$.results$simulation)
                                       }, 
                                       #' @description
                                       #' Retrieve the simulated networks stored in the object.
                                       #' Requires \code{simulate} or \code{estimate} to have been run first.
                                       #' @return A list of \code{\link{iglm.data}} objects representing
                                       #'   the simulated networks, invisibly. Returns an error if no samples
                                       #'   are available.
                                       get_samples = function() {
                                         if(is.null(private$.results$samples)) {
                                           stop("No samples available. Please run `simulate()` first.", call. = FALSE)
                                         } else {
                                           invisible(private$.results$samples)  
                                         }
                                         
                                       },
                                       #' @description
                                       #' Replace the internal MCMC sampler with a new one.
                                       #' This is useful for changing the sampling scheme without
                                       #' redefining the entire model.
                                       #' @param sampler A \code{\link{sampler_iglm}} object.
                                       #'  @return The \code{\link{iglm_object}} itself, invisibly.
                                       set_sampler = function(sampler) {
                                         if (!inherits(sampler, "sampler_iglm")) {
                                           stop("`sampler` must be a 'sampler_iglm' object created with sampler.iglm().")
                                         }
                                         private$.sampler <- sampler
                                         if(private$.control$display_progress){
                                           cat("Sampler has been set successfully.\n")
                                         }
                                         invisible(self)
                                       },
                                       #' @description
                                       #' Replace the internal `iglm.data` data object with a new one. This is
                                       #' useful for applying a fitted model to new observed data. Recalculates
                                       #' count statistics and re-validates the object.
                                       #' @param x A \code{\link{iglm.data}} `` object containing the new observed data.
                                       #' @return The \code{\link{iglm_object}} itself, invisibly.
                                       set_target = function(x) {
                                         if(!"iglm.data" %in% class(x)) {
                                           stop("The target object must be of class 'iglm.data'.", call. = FALSE)
                                         } else {
                                           private$.iglm.data <- x
                                           # Reset the results object
                                           private$.results <- results(size_coef = length(private$.coef), 
                                                                       size_coef_popularity =length(private$.coef_popularity)*
                                                                         private$.preprocess$includes_popularity)
                                           
                                         }
                                         if(private$.control$display_progress){
                                           cat("Target iglm.data object has been set successfully.\n")
                                         }
                                         private$.calc_count_statistics()
                                         private$.validate()
                                       }
                                     ),
                                     
                                     active = list(
                                       #' @field formula (`formula`) Read-only. The model formula specifying terms and data object.
                                       formula = function(value) { if(missing(value)) private$.formula else stop("`formula` is read-only.", call. = FALSE) },
                                       
                                       #' @field coef (`numeric`) Read-only. The current vector of non-popularity coefficient estimates or initial values.
                                       coef = function(value) { if(missing(value)) private$.coef else stop("`coef` is read-only.", call. = FALSE) },
                                       
                                       #' @field coef_popularity (`numeric` or `NULL`) Read-only. The current vector of popularity coefficient estimates or initial values, or `NULL` if not applicable.
                                       coef_popularity = function(value) { if(missing(value)) private$.coef_popularity else stop("`coef_popularity` is read-only.", call. = FALSE) },
                                       
                                       #' @field results (`results`) Read-only. The \code{\link{results}} R6 object containing all estimation and simulation outputs.
                                       results = function(value) { if(missing(value)) private$.results else stop("`results` is read-only.", call. = FALSE) },
                                       
                                       #' @field iglm.data (`iglm.data`) Read-only. The associated \code{\link{iglm.data}} R6 object containing the network and attribute data.
                                       iglm.data = function(value) { if(missing(value)) private$.iglm.data else stop("`iglm.data` is read-only.", call. = FALSE) },
                                       
                                       #' @field control (`control.iglm`) Read-only. The \code{\link{control.iglm}} object specifying estimation parameters.
                                       control = function(value) { if(missing(value)) private$.control else stop("`control` is read-only.", call. = FALSE) },
                                       
                                       #' @field sampler (`sampler_iglm`) Read-only. The  \code{\link{sampler_iglm}} object specifying MCMC sampling parameters.
                                       sampler = function(value) { if(missing(value)) private$.sampler else stop("`sampler` is read-only.", call. = FALSE) },
                                       
                                       #' @field count_statistics (`numeric`) Read-only. A named vector of the observed network statistics corresponding to the model terms, calculated on the current `iglm.data` data.
                                       count_statistics = function(value) { if(missing(value)) private$.count_statistics else stop("`count_statistics` is read-only.", call. = FALSE) }
                                     )
)

#' @title Construct a iglm Model Specification Object
#' @description 
#' The \code{iglm} package implements a comprehensive regression framework  introduced in Fritz et al. (2025) for 
#' studying relationships among attributes \eqn{(X, Y)} under network interference \eqn{(Z)}. 
#' It is based on a joint probability model for dependent 
#' outcomes (\eqn{Y}) and network connections \eqn{(Z)}, conditional on a fixed set of 
#' predictors (X). This approach generalizes standard 
#' Generalized Linear Models (GLMs) to settings where the responses and 
#' connections of units are interdependent. The framework is
#' designed to be interpretable by representing conditional distributions 
#' as GLMs, scalable to large networks via pseudo-likelihood 
#' and convex optimization, and provides insight into 
#' outcome-connection dependencies (i.e., spillover effects) that are 
#' missed by conditional models.
#'
#' The joint probability density is specified as an exponential-family model
#' of the form:
#' \deqn{f_{\theta}(y,z,x) \propto \Big[\prod_{i=1}^{N} a_y(y_i) \exp(\theta_g^T g_i(x_i, y_i^*)) \Big] \times  
#' \Big[\prod_{i \ne j} a_z(z_{i,j}) \exp(\theta_h^T h_{i,j}(x, y_i^*, y_j^*, z)) \Big],}
#' which is defined by two distinct sets of user-specified features:
#' \itemize{
#'   \item \strong{\eqn{g_i(x,y,z)}}: A vector of actor-level functions (or "g-terms")
#'     that describe the relationship between an individual actor \eqn{i}'s
#'     predictors (\eqn{x_i}) and their own response (\eqn{y_i}).
#'   \item \strong{\eqn{h_{i,j}(x,y,z)}}: A vector of pair-level functions (or "h-terms")
#'     that specify how the connections (\eqn{z}) and responses (\eqn{y_i, y_j})
#'     of a pair of units \eqn{\{i,j\}} depend on each other and the wider 
#'     network structure.
#' }
#' This separation allows the model to simultaneously capture individual-level
#' behavior (via \eqn{g_i}) and dyadic, network-based dependencies (via \eqn{h_{i,j}}), 
#' including local dependence limited to overlapping neighborhoods (see, Fritz et al., 2025). 
#' This help page documents the various statistics available in 'iglm',
#' corresponding to the \eqn{g_i} (attribute-level) and \eqn{h_{i,j}} (pair-level)
#' components of the joint model.
#' This is a user-facing constructor for creating a \code{\link{iglm_object}}. This \code{R6} object 
#' encompasses the complete model specification, linking the formula, data (\code{\link{iglm.data}} object), 
#' initial coefficients, MCMC sampler settings, and estimation controls.
#' It serves as the primary input for subsequent methods like \code{$estimate()} and \code{$simulate()}.
#' @return An object of class \code{\link{iglm_object}}.
#' 
#'
#' @param formula A model `formula` object. The left-hand side should be the
#'   name of a `iglm.data` object available in the calling environment. 
#'   See \code{\link{model_terms}} for details on specifying the right-hand side terms.
#' @param coef Optional numeric vector of initial coefficients for the structural
#'   (non-popularity) terms in `formula`. If `NULL`, coefficients are
#'   initialized to zero. Length must match the number of terms.
#' @param coef_popularity Optional numeric vector specifying the initial popularity
#'   coefficients. Required if `formula` includes popularity terms, otherwise
#'   should be `NULL`. Length must match `n_actor` (for undirected) or
#'   `2 * n_actor` (for directed).
#' @param sampler An object of class \code{\link{sampler_iglm}}, controlling the MCMC sampling scheme. If `NULL`,
#'   default sampler settings will be used.
#' @param control An object of class \code{\link{control.iglm}}, specifying parameters for the estimation algorithm.
#'   If `NULL`, default control settings will be used.
#' @param file Optional character string specifying a file path to load a
#'  previously saved  \code{\link{iglm_object}} from disk (in RDS format). If provided,
#'  other arguments are ignored and the object is loaded from the file.
#' @aliases iglm_object
#' @examples
#' # Example usage:
#' library(iglm)
#' # Create a iglm.data data object (example)
#' n_actors <- 50
#' neighborhood <- matrix(1, nrow = n_actors, ncol = n_actors)
#' xyz_obj <- iglm.data(neighborhood = neighborhood, directed = FALSE,
#'                    type_x = "binomial", type_y = "binomial")
#' # Define ground truth coefficients
#' gt_coef <- c("edges_local" = 3, "attribute_y" = -1, "attribute_x" = -1)
#' gt_coef_pop <- rnorm(n = n_actors, -2, 1)
#' # Define MCMC sampler
#' sampler_new <- sampler.iglm(n_burn_in = 100, n_simulation = 10,
#'                                sampler.x = sampler.net_attr(n_proposals = n_actors * 10, seed = 13),
#'                                sampler.y = sampler.net_attr(n_proposals = n_actors * 10, seed = 32),
#'                                sampler.z = sampler.net_attr(n_proposals = sum(neighborhood > 0
#'                                ) * 10, seed = 134),
#'                                init_empty = FALSE)
#' # Create iglm model specification
#' model_tmp_new <- iglm(formula = xyz_obj ~ edges(mode = "local") +
#'                           attribute_y + attribute_x + popularity,
#'                           coef = gt_coef,
#'                           coef_popularity = gt_coef_pop,
#'                           sampler = sampler_new,
#'                           control = control.iglm(accelerated = FALSE,
#'                           max_it = 200, display_progress = FALSE, var = TRUE))
#' # Simulate from the model
#' model_tmp_new$simulate()
#' model_tmp_new$set_target(model_tmp_new$get_samples()[[1]])
#' 
#' # Estimate model parameters
#' model_tmp_new$estimate()
#' 
#' # Model Assessment
#' model_tmp_new$model_assessment(formula = ~  degree_distribution )
#' # model_tmp_new$results$plot(model_assessment = TRUE)
#'                                                    
#' @references 
#' 
#' Fritz, C., Schweinberger, M., Bhadra, S., and D.R. Hunter (2025). A Regression Framework for Studying Relationships among Attributes under Network Interference. Journal of the American Statistical Association, to appear.
#' 
#' Schweinberger, M. and M.S. Handcock (2015). Local Dependence in Random Graph Models: Characterization, Properties, and Statistical Inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.
#' 
#' Schweinberger, M. and J.R. Stewart (2020). Concentration and Consistency Results for Canonical and Curved Exponential-Family Models of Random Graphs. The Annals of Statistics, 48, 374-396.
#' 
#' Stewart, J.R. and M. Schweinberger (2025). Pseudo-Likelihood-Based M-Estimation of Random Graphs with Dependent Edges and Parameter Vectors of Increasing Dimension. The Annals of Statistics, to appear. 
#' 
#' @export
iglm <- function(formula = NULL, coef= NULL, coef_popularity = NULL, sampler = NULL, 
                 control = NULL, file = NULL) {
  # browser()
  iglm_object_generator$new(
    formula = formula, 
    coef = coef, 
    coef_popularity = coef_popularity, 
    sampler = sampler, 
    control = control, 
    file = file
  )
}

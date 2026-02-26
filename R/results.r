#' @docType class
#' @title R6 Class for Storing iglm Estimation and Simulation Results
#' @description
#' The `results` class stores estimation (`$estimate()`) and simulation (`$simulate()`) results.
#'
#' This class is primarily intended for internal use within the `iglm`
#' framework but provides structured access to the results via the active
#' bindings of the main `iglm_object`.
#' @import R6
#' @import RcppProgress
#' @importFrom graphics plot lines abline layout title par axis boxplot
results.generator <- R6::R6Class("results",
                                 private = list(
                                   .coefficients_path = NULL,
                                   .samples = NULL,
                                   .stats= NULL,
                                   .var= NULL,
                                   .fisher_degrees= NULL,
                                   .fisher_nondegrees= NULL,
                                   .score_degrees= NULL,
                                   .score_nondegrees= NULL,
                                   .llh= NULL, 
                                   .model_assessment= NULL, 
                                   .prediction= NULL, 
                                   .estimated = FALSE
                                 ),public = list(
                                   #' @description
                                   #' Creates a new `results` object. Initializes internal fields, primarily
                                   #' setting up an empty matrix for the `coefficients_path` based on the
                                   #' expected number of coefficients.
                                   #' @param size_coef (integer) The number of non-degrees (structural)
                                   #'   coefficients in the model.
                                   #' @param size_coef_degrees (integer) The number of degrees coefficients
                                   #'   in the model (0 if none).
                                   #' @param file (character or `NULL`) If provided, loads the sampler state from
                                   #'  the specified .rds file instead of initializing from parameters.
                                   #' @return A new `results` object, initialized to hold results for a model
                                   #'   with the specified dimensions.
                                   initialize = function(size_coef, size_coef_degrees, file) {
                                     if(is.null(file)){
                                       private$.coefficients_path <- matrix(NA, 
                                                                            nrow = 0, 
                                                                            ncol = size_coef + size_coef_degrees)
                                       private$.llh <- numeric(0)
                                       private$.samples <- list() 
                                       private$.prediction <- list() 
                                     } else {
                                       data_loaded <- readRDS(file)
                                       required_fields <- c("coefficients_path", "samples", "stats", "var",
                                                            "fisher_degrees", "fisher_nondegrees",
                                                            "score_degrees", "score_nondegrees",
                                                            "llh", "model_assessment", "estimated","prediction")
                                       if (!is.list(data_loaded) || !all(required_fields %in% names(data_loaded))) {
                                         stop("File does not contain a valid results state.", call. = FALSE)
                                       }
                                       private$.coefficients_path <- data_loaded$coefficients_path
                                       private$.samples <- data_loaded$samples
                                       private$.stats <- data_loaded$stats
                                       private$.var <- data_loaded$var
                                       private$.fisher_degrees <- data_loaded$fisher_degrees
                                       private$.fisher_nondegrees <- data_loaded$fisher_nondegrees
                                       private$.score_degrees <- data_loaded$score_degrees
                                       private$.score_nondegrees <- data_loaded$score_nondegrees
                                       private$.llh <- data_loaded$llh
                                       private$.model_assessment <- data_loaded$model_assessment
                                       private$.estimated <- data_loaded$estimated
                                       private$.prediction <- data_loaded$prediction
                                     }
                                     
                                     invisible(self)
                                   }, 
                                   #' @description
                                   #' Stores the results object generated by a model assessment (goodness-of-fit)
                                   #' procedure within this `results` container.
                                   #' @param res An object containing the model assessment results, expected to
                                   #'   have the class `iglm_model_assessment`.
                                   #' @return The `results` object itself (`self`), invisibly. Called for its
                                   #'   side effect of storing the assessment results.
                                   set_model_assessment = function(res){
                                     if(!inherits(res, "iglm_model_assessment")){
                                       stop("`res` must be of class `iglm_model_assessment`.", call. = FALSE)
                                     }
                                     private$.model_assessment <- res
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Stores prediction results.
                                   #' @param prediction An object containing the prediction results (is a list of class `iglm.prediction`.
                                   set_prediction = function(prediction){
                                     if(!inherits(prediction, "iglm.prediction")){
                                       stop("`prediction` must be of class `iglm.prediction`.", call. = FALSE)
                                     }
                                     private$.prediction <- prediction
                                   },
                                   #' @description
                                   #' Gathers the current state of the `results` object into a list for saving
                                   #' or inspection. This includes all internal fields such as coefficient paths,
                                   #' samples, statistics, variance-covariance matrix, Fisher information,
                                   #' score vectors, log-likelihood values, model assessment results, and
                                   #' estimation status.
                                   #' @return A list containing all the internal fields of the `results` object.
                                   gather = function(){
                                     data_to_save <- list(
                                       coefficients_path = private$.coefficients_path,
                                       samples = private$.samples,
                                       stats = private$.stats,
                                       var = private$.var,
                                       fisher_degrees = private$.fisher_degrees,
                                       fisher_nondegrees = private$.fisher_nondegrees,
                                       score_degrees = private$.score_degrees,
                                       score_nondegrees = private$.score_nondegrees,
                                       llh = private$.llh,
                                       prediction = private$.prediction,
                                       model_assessment = private$.model_assessment,
                                       estimated = private$.estimated
                                     )
                                     return(data_to_save)
                                   }, 
                                   #' @description
                                   #' Saves the current state of the `results` object to a specified file path
                                   #' in RDS format. This allows for persisting the results for later
                                   #' retrieval and analysis.
                                   #' @param file (character) The file path where the results state should be
                                   #'   saved. Must be a valid character string.
                                   #' @return The `results` object itself (`self`), invisibly. 
                                   save = function(file){
                                     if (missing(file) || !is.character(file) || length(file) != 1) {
                                       stop("A valid 'file' (character string) must be provided.", call. = FALSE)
                                     }
                                     data_to_save <- self$gather()
                                     saveRDS(data_to_save, file = file)
                                     message(paste("Object state saved to:", file))
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Resizes the internal storage for the coefficient paths to accommodate a
                                   #' different number of coefficients. This is useful if the model structure
                                   #' changes and the results object needs to be reset.
                                   #' @param size_coef (integer) The new number of non-degrees coefficients.
                                   #' @param size_coef_degrees (integer) The new number of degrees
                                   #'  coefficients.
                                   #'  @return The `results` object itself (`self`), invisibly. 
                                   resize = function(size_coef, size_coef_degrees){
                                     private$.coefficients_path <- matrix(NA, 
                                                                          nrow = 0, 
                                                                          ncol = size_coef + size_coef_degrees)
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Updates the internal fields of the `results` object with new outputs,
                                   #' typically after an estimation run (`$estimate()`) or simulation run
                                   #' (`$simulate()`). Allows selectively updating components. Appends to
                                   #' `coefficients_path` and `llh` if called multiple times after estimation.
                                   #' Replaces `samples` and `stats`.
                                   #' @param coefficients_path (matrix) A matrix where rows represent iterations
                                   #'   and columns represent all coefficients (non-degrees then degrees),
                                   #'   showing their values during estimation. If provided, appends to any
                                   #'   existing path.
                                   #' @param samples (list) A list of simulated `iglm.data` objects (class `iglm.data.list`).
                                   #'   If provided, replaces any existing samples.
                                   #' @param var (matrix) The estimated variance-covariance matrix for the
                                   #'   non-degrees coefficients. Replaces existing matrix.
                                   #' @param fisher_degrees (matrix) The Fisher information matrix for
                                   #'   degrees coefficients. Replaces existing matrix.
                                   #' @param fisher_nondegrees (matrix) The Fisher information matrix for
                                   #'   non-degrees coefficients. Replaces existing matrix.
                                   #' @param score_degrees (numeric) The score vector for degrees coefficients.
                                   #'   Replaces existing vector.
                                   #' @param score_nondegrees (numeric) The score vector for non-degrees
                                   #'   coefficients. Replaces existing vector.
                                   #' @param llh (numeric) Log-likelihood value(s). If provided, appends to the
                                   #'   existing vector of log-likelihoods.
                                   #' @param stats (matrix) A matrix of summary statistics from simulations,
                                   #'   where rows correspond to simulations and columns to statistics. Replaces
                                   #'   or extends the existing matrix and will be turned into a mcmc object from the `coda` package.
                                   #' @param estimated (logical) A flag indicating whether these results come
                                   #'   from a completed estimation run. Updates the internal status.
                                   #' @return The `results` object itself (`self`), invisibly. Called for its
                                   #'   side effects.
                                   #' @importFrom coda mcmc thin
                                   update = function(coefficients_path= NULL, 
                                                     samples= NULL, 
                                                     var= NULL,
                                                     fisher_degrees= NULL,
                                                     fisher_nondegrees= NULL,
                                                     score_degrees= NULL,
                                                     score_nondegrees= NULL,
                                                     llh= NULL,
                                                     stats= NULL, 
                                                     estimated = FALSE) {
                                     
                                     private$.estimated <- estimated
                                     if(!is.null(var)){
                                       private$.var <- var
                                     }
                                     if(!is.null(fisher_degrees)){
                                       private$.fisher_degrees <- fisher_degrees
                                     }
                                     if(!is.null(fisher_nondegrees)){
                                       private$.fisher_nondegrees <- fisher_nondegrees
                                     }
                                     if(!is.null(score_degrees)){
                                       private$.score_degrees <- score_degrees
                                     }
                                     if(!is.null(score_nondegrees)){
                                       private$.score_nondegrees <- score_nondegrees
                                     }
                                     if(!is.null(llh)){
                                       private$.llh <- c(private$.llh, llh)
                                     }
                                     if(!is.null(stats)){
                                       
                                       if(length(private$.stats) ==0){
                                         private$.stats <- mcmc(stats)
                                       } else {
                                         combined_data <- rbind(private$.stats, mcmc(stats))
                                         private$.stats <-  combined_chain <- mcmc(
                                           data = combined_data, 
                                           start = start(private$.stats), 
                                           thin = thin(private$.stats)
                                         )
                                         
                                         
                                       }
                                       
                                     }
                                     if(!is.null(coefficients_path)){
                                       private$.coefficients_path <- 
                                         rbind(private$.coefficients_path, coefficients_path)
                                     }
                                     if(!is.null(samples)){
                                       if(! "iglm.data.list" %in% class(samples)){
                                         stop("`samples` must be of class `iglm.data.list`.", call. = FALSE)
                                       } else {
                                         if(length(private$.samples) ==0){
                                           private$.samples <- samples
                                         } else {
                                           private$.samples <- append_iglm.data(private$.samples, samples)
                                         }
                                         
                                         
                                       }
                                     }
                                     invisible(self)
                                   }, 
                                   #' @description
                                   #' Clears the stored simulation samples (`.samples`) and statistics (`.stats`) from the object,
                                   #' resetting it to an empty list. This might be used to save memory or
                                   #' before running new simulations. 
                                   #' @return The `results` object itself (`self`), invisibly.
                                   remove_samples = function(){
                                     private$.samples <- NULL
                                     private$.stats <- NULL
                                   },
                                   #' @description
                                   #' Generates diagnostic plots for the estimation results. Currently plots:
                                   #' \itemize{
                                   #'   \item The log-likelihood path across iterations.
                                   #'   \item The convergence paths for degrees coefficients (if present).
                                   #'   \item The convergence paths for non-degrees coefficients.
                                   #' }
                                   #' Optionally, can also trigger plotting of model assessment results if available.
                                   #' @param trace (logical) If `TRUE` (default), plot the trace plots of the estimation
                                   #'   (log-likelihood and coefficient paths). Requires
                                   #'   model to be estimated.
                                   #' @param model_assessment (logical) If `TRUE`, attempts to plot the results
                                   #'   stored in the `.model_assessment` field. Requires model assessment to
                                   #'   have been run and a suitable `plot` method for `iglm_model_assessment`
                                   #'   objects to exist. Default is `FALSE`.
                                   #' @param stats (logical) If `TRUE`, plots the normalized statistics from simulations. 
                                   #'    Default is `FALSE`.
                                   #' @param ... Additional fits with identical model_assessment terms are currently identified from this argument. 
                                   #'    The names of the arguments are shown as the legend in the model assessment plots. 
                                   #' 
                                   #' @details Requires estimation results (`private$.estimated == TRUE`) to plot
                                   #'   convergence diagnostics. Requires model assessment results for the
                                   #'   model assessment plots. 
                                   plot = function(trace = FALSE, stats = FALSE, model_assessment = FALSE, ...){
                                     
                                     
                                     if(stats + trace + model_assessment == 0){
                                       stop("At least one of `stats`, `trace`, or `model_assessment` must be TRUE.", call. = FALSE)
                                     }
                                     if(stats){
                                       if(length(private$.stats) == 0){
                                         stop("No samples available to plot.", call. = FALSE)
                                       } else {
                                         normalized <- scale(private$.stats)
                                         plot(NA, xlim=c(1,nrow(normalized)), ylim=range(normalized), xlab="Sample", ylab="Normalized Statistic", bty ="l")
                                         for (tmp in 1:ncol(normalized)){
                                           lines(y = normalized[,tmp], x = 1:nrow(normalized), col = tmp)
                                         }
                                         
                                       }
                                     }
                                     if(trace){
                                       if(!private$.estimated){
                                         stop("Model has not been estimated yet. Cannot plot results.", call. = FALSE)
                                       }
                                       
                                       plot(private$.llh, type = "l", xlab = "Iteration", ylab = "Log-likelihood", bty ="l")
                                       
                                       if(!is.null(private$.score_degrees)){
                                         coefficients_path_np <- private$.coefficients_path[,1:nrow(private$.var)]
                                         coefficients_path_p <- private$.coefficients_path[,(nrow(private$.var)+1):ncol(private$.coefficients_path)]
                                         
                                         plot(NA, xlim=c(1,nrow(coefficients_path_p)), ylim=range(coefficients_path_p),
                                              xlab="Iteration", ylab="Degree Coefficients", bty ="l")
                                         for (tmp in 1:ncol(coefficients_path_p)){
                                           lines(y = coefficients_path_p[,tmp], x = 1:nrow(coefficients_path_p), col = tmp)
                                         }
                                         
                                         plot(NA, xlim=c(1,nrow(coefficients_path_np)), 
                                              ylim=range(coefficients_path_np), 
                                              xlab="Iteration", ylab="Coefficients", bty ="l")
                                         for (tmp in 1:ncol(coefficients_path_np)){
                                           lines(y = coefficients_path_np[,tmp], x = 1:nrow(coefficients_path_np), col = tmp)
                                         }
                                       } else {
                                         plot(NA, xlim=c(1,nrow(private$.coefficients_path)), ylim=range(private$.coefficients_path), xlab="Iteration", ylab="Coefficients")
                                         for (tmp in 1:ncol(private$.coefficients_path)){
                                           lines(y = private$.coefficients_path[,tmp], x = 1:nrow(private$.coefficients_path), col = tmp)
                                         }
                                       }
                                     }
                                     if(model_assessment){
                                       if(is.null(private$.model_assessment)){
                                         stop("No model assessment available to plot.", call. = FALSE)
                                       }
                                       
                                       dot_list <- list(...)
                                       # Check what parts of the dot_list are of class "iglm_model_assessment" and whose names conincide with the planned "model_assessment"
                                       good_ind <- unlist(lapply(dot_list, function(x){
                                         if(inherits(x, "iglm_model_assessment")){
                                           return(TRUE)
                                         } else {
                                           return(FALSE)
                                         }
                                       }))
                                       dot_list <- dot_list[good_ind]
                            
                                       good_ind <- unlist(lapply(dot_list, function(x){
                                         if(identical(sort(x$names), sort(private$.model_assessment$names))){
                                           return(TRUE)
                                         } else {
                                           return(FALSE)
                                         }
                                       }))
                                       dot_list <- dot_list[good_ind]
                                       if(length(dot_list) >0){
                                         add <- TRUE
                                         colors_tmp <- generate_colors(length(dot_list)+1)
                                         if(is.null(names(dot_list))){
                                           names_tmp <- c(private$.model_assessment$name, paste0("Model ", 1:length(dot_list)))
                                         } else {
                                           names_tmp <- c(private$.model_assessment$name, names(dot_list))
                                         }
                                       } else {
                                         add <- FALSE
                                       }
                                       
                                       if(private$.model_assessment$include_mcmc){
                                         if(add){
                                           stop("Adding model assessment plots when MCMC diagnostics are included is not supported.", call. = FALSE)
                                         } 
                                         normalized <- private$.stats
                                         normalized <- sweep(normalized, 2, private$.model_assessment$sufficient_statistics, "/")
                                         for(i in 1:ncol(normalized)){
                                           plot(density(normalized[,i]), main = paste0(names(private$.model_assessment$sufficient_statistics)[i]), 
                                                bty ="l", xlab = "Ratio between Simulated and Observed Sufficient Statistics")
                                           rug(normalized[,i], lwd = 1)
                                         }
                                       }
                                       tmp_names <- names(private$.model_assessment$observed)
                                       base_names <- private$.model_assessment$base_name
                                       k <- 0
                                       for(i in base_names){
                                         k <- k+1
                                         if(i == "degree_distribution"){
                                           # Degree -----
                                           # Directed 
                                           # In degree
                                           if(is.list(private$.model_assessment$observed$degree_distribution)){
                                             if(add){
                                               x_positions <- private$.model_assessment$observed$degree_distribution$in_degree
                                               
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution$in_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               simulated_list <- list()
                                               for(j in 1:length(dot_list)){
                                                 simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution$in_degree
                                                 })
                                                 simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                               }
                                               simulated_list <- do.call("rbind",simulated_list)
                                               ylim <- range(c(simulated_list, simulated,
                                                               private$.model_assessment$observed$degree_distribution$in_degree))
                                               
                                               plot(private$.model_assessment$observed$degree_distribution$in_degree, 
                                                    xlab = "Indegree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution$in_degree))),
                                                                n = 10))
                                               axis(side = 2)       
                                               x <- as.numeric(names(private$.model_assessment$observed$degree_distribution$in_degree))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 2)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               for(j in 1:length(dot_list)){
                                                 x_positions <- dot_list[[j]]$observed$degree_distribution$in_degree
                                                 simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution$in_degree
                                                 })
                                                 simulated <- do.call("rbind",simulated)
                                                 
                                                 x <- as.numeric(names(dot_list[[j]]$observed$degree_distribution$in_degree))
                                                 x_polygon <- c(x, rev(x)) 
                                                 y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                                 polygon(x_polygon, y_polygon, 
                                                         col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                         border = NA)    
                                                 lines(x, colMeans(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 2)
                                                 lines(x, colMins(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                                 lines(x, colMaxs(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                               }
                                               lines(private$.model_assessment$observed$degree_distribution$in_degree, type = "l", 
                                                     col = "black",lwd = 2)
                                               legend("topright", legend = c("Observed",names_tmp),
                                                      col = c("black",colors_tmp),
                                                      lwd = 2, bty = "n")
                                               
                                             } else {
                                               
                                               x_positions <- private$.model_assessment$observed$degree_distribution$in_degree
                                               
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution$in_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$degree_distribution$in_degree))
                                               
                                               plot(private$.model_assessment$observed$degree_distribution$in_degree, 
                                                    xlab = "Indegree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution$in_degree))),
                                                                n = 10))
                                               axis(side = 2)       
                                               boxplot(simulated, at = as.numeric(names(private$.model_assessment$observed$degree_distribution$in_degree)),
                                                       add = TRUE, col = "#87CEEB80", boxwex = 0.5,axes = FALSE)
                                               lines(private$.model_assessment$observed$degree_distribution$in_degree, type = "l", 
                                                     col = "#D55E00",lwd = 2)
                                             }
                                             # Out degree 
                                             if(add){
                                               x_positions <- private$.model_assessment$observed$degree_distribution$out_degree
                                               
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution$out_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               simulated_list <- list()
                                               for(j in 1:length(dot_list)){
                                                 simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution$out_degree
                                                 })
                                                 simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                               }
                                               simulated_list <- do.call("rbind",simulated_list)
                                               ylim <- range(c(simulated_list, simulated,
                                                               private$.model_assessment$observed$degree_distribution$out_degree))
                                            
                                               plot(private$.model_assessment$observed$degree_distribution$out_degree, 
                                                    xlab = "Outdegree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution$out_degree))),
                                                                n = 10))
                                               axis(side = 2)       
                                               x <- as.numeric(names(private$.model_assessment$observed$degree_distribution$out_degree))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 2)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               for(j in 1:length(dot_list)){
                                                 x_positions <- dot_list[[j]]$observed$degree_distribution$out_degree
                                                 simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution$out_degree
                                                 })
                                                 simulated <- do.call("rbind",simulated)
                                                 
                                                 x <- as.numeric(names(dot_list[[j]]$observed$degree_distribution$out_degree))
                                                 x_polygon <- c(x, rev(x)) 
                                                 y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                                 polygon(x_polygon, y_polygon, 
                                                         col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                         border = NA)    
                                                 lines(x, colMeans(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 2)
                                                 lines(x, colMins(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                                 lines(x, colMaxs(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                               }
                                               lines(private$.model_assessment$observed$degree_distribution$out_degree, type = "l", 
                                                     col = "black",lwd = 2)
                                               legend("topright", legend = c("Observed",names_tmp),
                                                      col = c("black",colors_tmp),
                                                      lwd = 2, bty = "n")
                                               
                                             } else {
                                               x_positions <- private$.model_assessment$observed$degree_distribution$out_degree
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution$out_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$degree_distribution$out_degree))
                                               
                                               plot(private$.model_assessment$observed$degree_distribution$out_degree, 
                                                    xlab = "Outdegree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution$out_degree))),
                                                                n = 10))
                                               axis(side = 2)     
                                               boxplot(simulated, at = as.numeric(names(private$.model_assessment$observed$degree_distribution$out_degree)),
                                                       add = TRUE, col = "#87CEEB80",boxwex = 0.5, axes = FALSE)
                                               lines(private$.model_assessment$observed$degree_distribution$out_degree, type = "l", 
                                                     col = "#D55E00",lwd = 2)
                                             }
                                           } else {
                                             # Undirected 
                                             if(add){
                                               x_positions <- private$.model_assessment$observed$degree_distribution
                                               
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               simulated_list <- list()
                                               for(j in 1:length(dot_list)){
                                                 simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution
                                                 })
                                                 simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                               }
                                               simulated_list <- do.call("rbind",simulated_list)
                                               ylim <- range(c(simulated_list, simulated,
                                                               private$.model_assessment$observed$degree_distribution))
                                               
                                               
                                               plot(private$.model_assessment$observed$degree_distribution, 
                                                    xlab = "Degree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution))),
                                                                n = 10))
                                               axis(side = 2)       
                                               x <- as.numeric(names(private$.model_assessment$observed$degree_distribution))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[1],lwd = 1)
                                               for(j in 1:length(dot_list)){
                                                 x_positions <- dot_list[[j]]$observed$degree_distribution
                                                 simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                   x$degree_distribution
                                                 })
                                                 simulated <- do.call("rbind",simulated)
                                                 
                                                 x <- as.numeric(names(dot_list[[j]]$observed$degree_distribution))
                                                 x_polygon <- c(x, rev(x)) 
                                                 y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                                 polygon(x_polygon, y_polygon, 
                                                         col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                         border = NA)    
                                                 lines(x, colMeans(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                                 lines(x, colMins(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                                 lines(x, colMaxs(simulated), type = "l", 
                                                       col = colors_tmp[j+1],lwd = 1)
                                                 
                                               }
                                               lines(private$.model_assessment$observed$degree_distribution, type = "l", 
                                                     col = "black",lwd = 2)
                                               legend("topright", legend = c("Observed",names_tmp),
                                                      col = c("black",colors_tmp),
                                                      lwd = 2, bty = "n")
                                               
                                             } else {
                                               x_positions <- private$.model_assessment$observed$degree_distribution
                                               
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$degree_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$degree_distribution))
                                               
                                               plot(private$.model_assessment$observed$degree_distribution, 
                                                    xlab = "Degree",ylim=ylim, 
                                                    xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3), 
                                                    ylab = "Percentage",type = "n", bty ="l", axes = FALSE)
                                               axis(side = 1,              
                                                    at = pretty(range(as.numeric(names(private$.model_assessment$observed$degree_distribution))),
                                                                n = 10))
                                               axis(side = 2)     
                                               
                                               
                                               boxplot(simulated, at = as.numeric(names(private$.model_assessment$observed$degree_distribution)),
                                                       add = TRUE, col = "#87CEEB80", boxwex = 0.5,axes = FALSE)
                                               lines(private$.model_assessment$observed$degree_distribution, type = "l", 
                                                     col = "#D55E00",lwd = 2)
                                             }
                                             
                                           }
                                           
                                           
                                           
                                           
                                         } else if(i  %in% c("dyadwise_shared_partner_distribution", 
                                                             "edgewise_shared_partner_distribution")){
                                           
                                           if(i  == "dyadwise_shared_partner_distribution"){
                                             xlab_tmp <- "Dyadwise Shared Partner"
                                           } else{
                                             xlab_tmp <- "Edgewise Shared Partner"
                                           }
                                           # ESP/DSP -----
                                           if(add){
                                             x_positions <- eval(parse(text = paste0("private$.model_assessment$observed$",tmp_names[k])))
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               eval(parse(text = paste0("x$",tmp_names[k])))
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             
                                             simulated_list <- list()
                                             for(j in 1:length(dot_list)){
                                               simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                 eval(parse(text = paste0("x$",tmp_names[k])))
                                               })
                                               simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                             }
                                             simulated_list <- do.call("rbind",simulated_list)
                                             ylim <- range(c(simulated_list, simulated,
                                                             x_positions))
                                             
                                             
                                             ylim <- range(c(simulated,
                                                             x_positions))
                                             plot(x_positions, 
                                                  xlab = xlab_tmp,ylab = "Percentage",type = "n",
                                                  ylim = ylim, axes = F,
                                                  bty ="l",
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3))
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(x_positions))),
                                                              n = 10))
                                             axis(side = 2)        
                                             
                                             x <- as.numeric(names(x_positions))
                                             x_polygon <- c(x, rev(x)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             lines(x, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               
                                               x_positions <- eval(parse(text = paste0("dot_list[[j]]$observed$",tmp_names[k])))
                                               
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 eval(parse(text = paste0("x$",tmp_names[k])))
                                               })
                                               
                                               simulated <- do.call("rbind",simulated)
                                               
                                               x <-  as.numeric(names(x_positions))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             
                                             lines(x_positions, type = "l", 
                                                   col = "black",lwd = 2)
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                           } else {
                                             x_positions <- eval(parse(text = paste0("private$.model_assessment$observed$",tmp_names[k])))
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               # x$edgewise_shared_partner_distribution
                                               eval(parse(text = paste0("x$",tmp_names[k])))
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             ylim <- range(c(simulated,
                                                             x_positions))
                                             plot(x_positions, 
                                                  xlab = xlab_tmp,ylab = "Percentage",type = "n",
                                                  ylim = ylim, axes = F,
                                                  bty ="l",
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.3, max(as.numeric(names(x_positions)))+0.3))
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(x_positions))),
                                                              n = 10))
                                             axis(side = 2)          
                                             
                                             boxplot(simulated, at = as.numeric(names(x_positions)),
                                                     add = TRUE, col = "#87CEEB80",boxwex = 0.5,axes = FALSE)
                                             lines(x_positions, type = "l", 
                                                   col = "#D55E00",lwd = 2)
                                           }
                                           
                                         } else if(i == "spillover_degree_distribution"){
                                           # Spillover degree -----
                                           if(add){
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$spillover_degree_distribution$in_spillover_degree
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             
                                             
                                             simulated_list <- list()
                                             for(j in 1:length(dot_list)){
                                               simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$spillover_degree_distribution$in_spillover_degree
                                               })
                                               simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                             }
                                             simulated_list <- do.call("rbind",simulated_list)
                                             ylim <- range(c(simulated_list,simulated,
                                                             private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree), na.rm = TRUE)
                                             x_positions <- private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree
                                             
                                             plot(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree, 
                                                  xlab = "Spillover Indegree",ylab = "Percentage",type = "n", 
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.5, max(as.numeric(names(x_positions)))+0.5),
                                                  ylim = ylim, bty ="l", axes = FALSE)
                                             
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree))),
                                                              n = 10))
                                             axis(side = 2)    
                                             x <- as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree))
                                             x_polygon <- c(x, rev(x)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             lines(x, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               x_positions <- dot_list[[j]]$observed$spillover_degree_distribution$in_spillover_degree
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$spillover_degree_distribution$in_spillover_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               x <- as.numeric(names(dot_list[[j]]$observed$spillover_degree_distribution$in_spillover_degree))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             lines(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree, type = "l", 
                                                   col = "black",lwd = 2)
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$spillover_degree_distribution$out_spillover_degree
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             
                                             ylim <- range(c(simulated,
                                                             private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree), na.rm = TRUE)
                                             x_positions <- private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree
                                             
                                             
                                             plot(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree, 
                                                  xlab = "Spillover Outdegree",ylab = "Percentage",type = "n", 
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.5, max(as.numeric(names(x_positions)))+0.5),
                                                  ylim = ylim, bty ="l", axes = FALSE)
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree))),
                                                              n = 10))
                                             axis(side = 2)    
                                             x <- as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree))
                                             x_polygon <- c(x, rev(x)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             lines(x, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               x_positions <- dot_list[[j]]$observed$spillover_degree_distribution$out_spillover_degree
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$spillover_degree_distribution$out_spillover_degree
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               x <- as.numeric(names(dot_list[[j]]$observed$spillover_degree_distribution$out_spillover_degree))
                                               x_polygon <- c(x, rev(x)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             
                                             lines(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree, type = "l", 
                                                   col = "black",lwd = 2)
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                           } else {
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$spillover_degree_distribution$in_spillover_degree
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             
                                             ylim <- range(c(simulated,
                                                             private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree), na.rm = TRUE)
                                             x_positions <- private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree
                                             
                                             plot(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree, 
                                                  xlab = "Spillover Indegree",ylab = "Percentage",type = "n", 
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.5, max(as.numeric(names(x_positions)))+0.5),
                                                  ylim = ylim, bty ="l", axes = FALSE)
                                             
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree))),
                                                              n = 10))
                                             axis(side = 2)    
                                             
                                             boxplot(simulated, 
                                                     at = as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree)),
                                                     add = TRUE, col = "#87CEEB80", boxwex = 0.5,axes = FALSE)
                                             lines(private$.model_assessment$observed$spillover_degree_distribution$in_spillover_degree, type = "l", 
                                                   col = "#D55E00",lwd = 2)
                                             
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$spillover_degree_distribution$out_spillover_degree
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             ylim <- range(c(simulated,
                                                             private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree))
                                             x_positions <- private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree
                                             
                                             plot(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree, 
                                                  xlab = "Spillover Outdegree",ylab = "Percentage",type = "n", 
                                                  xlim = c(min(as.numeric(names(x_positions)))-0.5, max(as.numeric(names(x_positions)))+0.5),
                                                  ylim = ylim, bty ="l", axes = FALSE)
                                             
                                             axis(side = 1,              
                                                  at = pretty(range(as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree))),
                                                              n = 10))
                                             axis(side = 2)    
                                             
                                             boxplot(simulated, 
                                                     at = as.numeric(names(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree)),
                                                     add = TRUE, col = "#87CEEB80", boxwex = 0.5,axes = FALSE)
                                             lines(private$.model_assessment$observed$spillover_degree_distribution$out_spillover_degree, type = "l", 
                                                   col = "#D55E00",lwd = 2)
                                             
                                           }
                                           
                                           
                                           
                                         } else if (i == "geodesic_distances_distribution"){
                                           # Geodesic distances -----
                                           if(add){
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$geodesic_distances_distribution
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             
                                             simulated_list <- list()
                                             for(j in 1:length(dot_list)){
                                               simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$geodesic_distances_distribution
                                               })
                                               simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                             }
                                             simulated_list <- do.call("rbind",simulated_list)
                                             ylim <- range(c(simulated_list,simulated,
                                                             private$.model_assessment$observed$geodesic_distances_distribution))
                                             
                                             
                                             x_positions <- 1:length(private$.model_assessment$observed$geodesic_distances_distribution)
                                             plot(x_positions, as.vector(private$.model_assessment$observed$geodesic_distances_distribution), 
                                                  xlab = "Geodesic Distance",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                  xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                             
                                             
                                             colnames(simulated) <- x_positions
                                             
                                             lines(x_positions, as.vector(private$.model_assessment$observed$geodesic_distances_distribution),
                                                   type = "l", col = "black",lwd = 2)
                                             axis(side = 1,
                                                  at = x_positions,             
                                                  labels = names(private$.model_assessment$observed$geodesic_distances_distribution)) 
                                             x_polygon <- c(x_positions, rev(x_positions)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             
                                             lines(x_positions, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x_positions, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x_positions, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               x_positions <- 1:length(dot_list[[j]]$observed$geodesic_distances_distribution)
                                               
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$geodesic_distances_distribution
                                               })
                                               
                                               simulated <- do.call("rbind",simulated)
                                               x_polygon <- c(x_positions, rev(x_positions)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x_positions, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x_positions, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x_positions, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                             
                                           } else {
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$geodesic_distances_distribution
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             ylim <- range(c(simulated,
                                                             private$.model_assessment$observed$geodesic_distances_distribution))
                                             
                                             
                                             x_positions <- 1:length(private$.model_assessment$observed$geodesic_distances_distribution)
                                             plot(x_positions, as.vector(private$.model_assessment$observed$geodesic_distances_distribution), 
                                                  xlab = "Geodesic Distance",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                  xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                             
                                             
                                             colnames(simulated) <- x_positions
                                             boxplot(simulated, at = x_positions,boxwex = 0.5,
                                                     add = TRUE, col = "#87CEEB80", axes = FALSE)
                                             
                                             lines(x_positions, as.vector(private$.model_assessment$observed$geodesic_distances_distribution),
                                                   type = "l", col = "#D55E00",lwd = 2)
                                             axis(side = 1,
                                                  at = x_positions,             
                                                  labels = names(private$.model_assessment$observed$geodesic_distances_distribution)) 
                                           }
                                           
                                         } else if (i == "y_distribution"){
                                           # y_distribution -----
                                           if(add){
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$y_distribution
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             simulated_list <- list()
                                             for(j in 1:length(dot_list)){
                                               simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$y_distribution
                                               })
                                               simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                             }
                                             simulated_list <- do.call("rbind",simulated_list)
                                             
                                             ylim <- range(c(simulated_list, simulated,
                                                             private$.model_assessment$observed$y_distribution))
                                             
                                             
                                             x_positions <- 1:length(private$.model_assessment$observed$y_distribution)
                                             plot(x_positions, as.vector(private$.model_assessment$observed$y_distribution), 
                                                  xlab = "Distribution of Y",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                  xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                             
                                             
                                             colnames(simulated) <- x_positions
                                             
                                             lines(x_positions, as.vector(private$.model_assessment$observed$y_distribution),
                                                   type = "l", col = "black",lwd = 2)
                                             axis(side = 1,
                                                  at = pretty(x_positions),             
                                                  labels = pretty(x_positions)) 
                                             x_polygon <- c(x_positions, rev(x_positions)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             
                                             lines(x_positions, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x_positions, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x_positions, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               x_positions <- 1:length(dot_list[[j]]$observed$y_distribution)
                                               
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$y_distribution
                                               })
                                               
                                               simulated <- do.call("rbind",simulated)
                                               x_polygon <- c(x_positions, rev(x_positions)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x_positions, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x_positions, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x_positions, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                             
                                           } else {
                                             if(private$.samples[[1]]$type_y == "normal"){
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$y_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$y_distribution))
                                               
                                               
                                               x_positions <- as.numeric(names(private$.model_assessment$observed$y_distribution))
                                               plot(x_positions, as.vector(private$.model_assessment$observed$y_distribution), 
                                                    xlab = "Distribution of Y",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                    xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                               
                                               
                                               colnames(simulated) <- x_positions
                                               
                               
                                               axis(side = 1,
                                                    at = pretty(x_positions),             
                                                    labels = pretty(x_positions)) 
                                               x_polygon <- c(x_positions, rev(x_positions)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha("#87CEEB80",alpha_level = 0.4),
                                                       border = NA)    
                                               
                                               lines(x_positions, colMeans(simulated), type = "l", 
                                                     col ="#87CEEB80",lwd = 2)
                                               lines(x_positions, colMins(simulated), type = "l", 
                                                     col = "#87CEEB80",lwd = 1)
                                               lines(x_positions, colMaxs(simulated), type = "l", 
                                                     col = "#87CEEB80",lwd = 1)

                                               lines(x_positions, as.vector(private$.model_assessment$observed$y_distribution),
                                                     type = "l", col = "#D55E00",lwd = 2)
                                               # legend("topright", legend = c("Observed","Simulated"),
                                               #        col = c("#D55E00","#87CEEB80"),
                                               #        lwd = 2, bty = "n")
                                               
                                             } else{
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$y_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$y_distribution))
                                               
                                               
                                               x_positions <- 1:length(private$.model_assessment$observed$y_distribution)
                                               plot(x_positions, as.vector(private$.model_assessment$observed$y_distribution), 
                                                    xlab = "Distribution of Y",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                    xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                               
                                               
                                               colnames(simulated) <- x_positions
                                               boxplot(simulated, at = x_positions,boxwex = 0.5,
                                                       add = TRUE, col = "#87CEEB80", axes = FALSE)
                                               
                                               lines(x_positions, as.vector(private$.model_assessment$observed$y_distribution),
                                                     type = "l", col = "#D55E00",lwd = 2)
                                               axis(side = 1,
                                                    at = pretty(x_positions),             
                                                    labels = pretty(x_positions)) 
                                             }
                                           
                                           }
                                           
                                         }else if (i == "x_distribution"){
                                           # x_distribution -----
                                           if(add){
                                             simulated <- lapply(private$.model_assessment$simulated, function(x){
                                               x$x_distribution
                                             })
                                             simulated <- do.call("rbind",simulated)
                                             simulated_list <- list()
                                             for(j in 1:length(dot_list)){
                                               simulated_list[[j]] <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$x_distribution
                                               })
                                               simulated_list[[j]] <- do.call("rbind",simulated_list[[j]])
                                             }
                                             simulated_list <- do.call("rbind",simulated_list)
                                             
                                             ylim <- range(c(simulated_list, simulated,
                                                             private$.model_assessment$observed$x_distribution))
                                             
                                             x_positions <- 1:length(private$.model_assessment$observed$x_distribution)
                                             plot(x_positions, as.vector(private$.model_assessment$observed$x_distribution), 
                                                  xlab = "Distribution of Y",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                  xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                             
                                             
                                             colnames(simulated) <- x_positions
                                             
                                             lines(x_positions, as.vector(private$.model_assessment$observed$x_distribution),
                                                   type = "l", col = "black",lwd = 2)
                                             axis(side = 1,
                                                  at = pretty(x_positions),             
                                                  labels = pretty(x_positions)) 
                                             x_polygon <- c(x_positions, rev(x_positions)) 
                                             y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                             polygon(x_polygon, y_polygon, 
                                                     col = add_alpha(colors_tmp[1],alpha_level = 0.4),
                                                     border = NA)    
                                             
                                             lines(x_positions, colMeans(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 2)
                                             lines(x_positions, colMins(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             lines(x_positions, colMaxs(simulated), type = "l", 
                                                   col = colors_tmp[1],lwd = 1)
                                             
                                             for(j in 1:length(dot_list)){
                                               x_positions <- 1:length(dot_list[[j]]$observed$x_distribution)
                                               
                                               simulated <- lapply(dot_list[[j]]$simulated, function(x){
                                                 x$x_distribution
                                               })
                                               
                                               simulated <- do.call("rbind",simulated)
                                               x_polygon <- c(x_positions, rev(x_positions)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha(colors_tmp[j+1],alpha_level = 0.4),
                                                       border = NA)    
                                               lines(x_positions, colMeans(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 2)
                                               lines(x_positions, colMins(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                               lines(x_positions, colMaxs(simulated), type = "l", 
                                                     col = colors_tmp[j+1],lwd = 1)
                                             }
                                             
                                             legend("topright", legend = c("Observed",names_tmp),
                                                    col = c("black",colors_tmp),
                                                    lwd = 2, bty = "n")
                                             
                                           } else {
                                             if(private$.samples[[1]]$type_y == "normal"){
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$x_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$x_distribution))
                                               
                                               x_positions <- as.numeric(names(private$.model_assessment$observed$x_distribution))
                                               plot(x_positions, as.vector(private$.model_assessment$observed$x_distribution), 
                                                    xlab = "Distribution of Y",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                    xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                               
                                               
                                               colnames(simulated) <- x_positions
                                               
                                            
                                               axis(side = 1,
                                                    at = pretty(x_positions),             
                                                    labels = pretty(x_positions)) 
                                               x_polygon <- c(x_positions, rev(x_positions)) 
                                               y_polygon <- c(colMins(simulated), rev(colMaxs(simulated))) 
                                               polygon(x_polygon, y_polygon, 
                                                       col = add_alpha("#87CEEB80",alpha_level = 0.4),
                                                       border = NA)    
                                               
                                               lines(x_positions, colMeans(simulated), type = "l", 
                                                     col = "#87CEEB80",lwd = 2)
                                               lines(x_positions, colMins(simulated), type = "l", 
                                                     col = "#87CEEB80",lwd = 1)
                                               lines(x_positions, colMaxs(simulated), type = "l", 
                                                     col = "#87CEEB80",lwd = 1) 
                                               lines(x_positions, as.vector(private$.model_assessment$observed$x_distribution),
                                                     type = "l", col = "#D55E00",lwd = 2)
                                               # legend("topright", legend = c("Observed","Simulated"),
                                               #        col = c("black","#87CEEB80"),
                                               #        lwd = 2, bty = "n")
                                             } else{
                                               simulated <- lapply(private$.model_assessment$simulated, function(x){
                                                 x$x_distribution
                                               })
                                               simulated <- do.call("rbind",simulated)
                                               ylim <- range(c(simulated,
                                                               private$.model_assessment$observed$x_distribution))
                                               
                                               
                                               x_positions <- 1:length(private$.model_assessment$observed$x_distribution)
                                               plot(x_positions, as.vector(private$.model_assessment$observed$x_distribution), 
                                                    xlab = "Distribution of X",ylab = "Percentage",type = "n", xaxt = "n", ylim=ylim, 
                                                    xlim = c(min(x_positions)-0.3, max(x_positions)+0.3), bty ="l")
                                               
                                               
                                               colnames(simulated) <- x_positions
                                               boxplot(simulated, at = x_positions,boxwex = 0.5,
                                                       add = TRUE, col = "#87CEEB80", axes = FALSE)
                                               
                                               lines(x_positions, as.vector(private$.model_assessment$observed$x_distribution),
                                                     type = "l", col = "#D55E00",lwd = 2)
                                               axis(side = 1,
                                                    at = pretty(x_positions),             
                                                    labels = pretty(x_positions)) 
                                             }
                                             
                                           }
                                           
                                         } 
                                         
                                       }
                                     }
                                   },
                                   #' @description
                                   #' Prints a concise summary of the contents of the `results` object,
                                   #' indicating whether various components (coefficients path, variance matrix,
                                   #' Fisher info, score, samples, stats, etc.) are available.
                                   #' @param ... Additional arguments (currently ignored).
                                   #' @return The `results` object itself (`self`), invisibly.
                                   print = function(...) {
                                     cat("Results Summary:\n")
                                     cat("----------------\n")
                                     cat("Number of Coefficient Paths Recorded:", nrow(private$.coefficients_path), "\n")
                                     cat("Log-Likelihoods Recorded:", length(private$.llh), "\n")
                                     if(!is.null(private$.var)){
                                       cat("Variance-Covariance Matrix Available\n")
                                     } else {
                                       cat("Variance-Covariance Matrix Not Available\n")
                                     }
                                     if(!is.null(private$.fisher_degrees)){
                                       cat("Fisher Information for degrees Available\n")
                                     } else {
                                       cat("Fisher Information for degrees Not Available\n")
                                     }
                                     if(!is.null(private$.fisher_nondegrees)){
                                       cat("Fisher Information for Non-degrees Available\n")
                                     } else {
                                       cat("Fisher Information for Non-degrees Not Available\n")
                                     }
                                     if(!is.null(private$.score_degrees)){
                                       cat("Score for degrees Available\n")
                                     } else {
                                       cat("Score for degrees Not Available\n")
                                     }
                                     if(!is.null(private$.score_nondegrees)){
                                       cat("Score for Non-degrees Available\n")
                                     } else {
                                       cat("Score for Non-degrees Not Available\n")
                                     }
                                     if(!is.null(private$.stats)){
                                       cat("Statistics from Simulations Available\n")
                                     } else {
                                       cat("Statistics from Simulations Not Available\n")
                                     }
                                     if(!is.null(private$.samples)){
                                       cat("Samples Available of class:", class(private$.samples), "\n")
                                     } else {
                                       cat("Samples Not Available\n")
                                     }
                                     invisible(self)
                                   }
                                 ),
                                 active = list(
                                   #' @field coefficients_path (`matrix` or `NULL`) Read-only. The path of all estimated coefficients across iterations.
                                   coefficients_path = function(value) { if(missing(value)) private$.coefficients_path else stop("`coefficients_path` is read-only.", call. = FALSE) },
                                   #' @field samples (`list` or `NULL`) Read-only. A list of simulated `iglm.data` objects (class `iglm.data.list`).
                                   samples = function(value) { if(missing(value)) private$.samples else stop("`samples` is read-only.", call. = FALSE)},
                                   #' @field stats (`matrix` or `NULL`) Read-only. Matrix of summary statistics for simulated samples, which are an `mcmc` obect from `coda`.
                                   stats = function(value) { if(missing(value)) private$.stats else stop("`stats` is read-only.", call. = FALSE)},
                                   #' @field var (`matrix` or `NULL`) Read-only. Estimated variance-covariance matrix for non-degrees coefficients.
                                   var = function(value) { if(missing(value)) private$.var else stop("`var` is read-only.", call. = FALSE) },
                                   #' @field fisher_degrees (`matrix` or `NULL`) Read-only. Fisher information matrix for degrees coefficients.
                                   fisher_degrees = function(value) { if(missing(value)) private$.fisher_degrees else stop("`fisher_degrees` is read-only.", call. = FALSE) },
                                   #' @field fisher_nondegrees (`matrix` or `NULL`) Read-only. Fisher information matrix for non-degrees coefficients.
                                   fisher_nondegrees = function(value) { if(missing(value)) private$.fisher_nondegrees else stop("`fisher_nondegrees` is read-only.", call. = FALSE) },
                                   #' @field score_degrees (`numeric` or `NULL`) Read-only. Score vector for degrees coefficients.
                                   score_degrees = function(value) { if(missing(value)) private$.score_degrees else stop("`score_degrees` is read-only.", call. = FALSE) },
                                   #' @field score_nondegrees (`numeric` or `NULL`) Read-only. Score vector for non-degrees coefficients.
                                   score_nondegrees = function(value) { if(missing(value)) private$.score_nondegrees else stop("`score_nondegrees` is read-only.", call. = FALSE) },
                                   #' @field llh (`numeric` or `NULL`) Read-only. Vector of log-likelihood values recorded during estimation.
                                   llh = function(value) { if(missing(value)) private$.llh else stop("`llh` is read-only.", call. = FALSE) },
                                   #' @field model_assessment (`list` or `NULL`) Read-only. Results from model assessment (goodness-of-fit).
                                   model_assessment = function(value) { if(missing(value)) private$.model_assessment else stop("`model_assessment` is read-only.", call. = FALSE)},
                                   #' @field estimated (`logical`) Read-only. Flag indicating if estimation has been completed.
                                   estimated = function(value) { if(missing(value)) private$.estimated else stop("`estimated` is read-only.", call. = FALSE) }
                                 )
)

#' Constructor for the results R6 Object
#'
#' @description
#' Creates a new instance of the `results` R6 class. This class is designed to
#' store various outputs from `iglm` model estimation and simulation. Users
#' typically do not need to call this constructor directly; it is used internally
#' by the `iglm_object`.
#'
#' @param size_coef (integer) The number of non-degrees coefficients the object
#'   should be initialized to accommodate.
#' @param size_coef_degrees (integer) The number of degrees coefficients
#'   the object should be initialized to accommodate.
#' @param file (character or NULL) Optional file path to load a previously saved
#'  `results` object. If provided, the object will be initialized by loading
#'  from this file.
#' @return An object of class `results` (and `R6`), initialized with empty or
#'   NA structures appropriately sized based on the input dimensions.
#' @export 
results <- function(size_coef,size_coef_degrees, file = NULL) {
  results.generator$new(size_coef = size_coef, 
                        size_coef_degrees =size_coef_degrees, 
                        file = file)
}

#' @importFrom grDevices hcl
generate_colors <- function(n, chroma = 35, luminance = 85, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  # Construct hex codes using the HCL space
  colors <- hcl(h = hues, c = chroma, l = luminance)
  return(colors)
}


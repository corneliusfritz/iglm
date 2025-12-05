#' @docType class
#' @title A R6 class to represent networks with unit-level attributes
#' @description
#' The `iglm.data` class is a container for storing, validating, and analyzing
#' unit-level attributes (x_attribute, y_attribute) and connections (z_network).
#' @import R6
#' @importFrom igraph graph_from_edgelist layout_with_fr vcount add_vertices V plot.igraph
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom graphics legend plot
#' @importFrom Matrix sparseMatrix spMatrix bdiag diag
#' @export
iglm.data_generator <- R6::R6Class("iglm.data",
                                 private = list(
                                   .x_attribute = NULL,
                                   .y_attribute = NULL,
                                   .z_network = NULL,
                                   .neighborhood = NULL,
                                   .overlap = NULL,
                                   .directed = NULL,
                                   .n_actor = NULL,
                                   .type_x = NULL,
                                   .type_y = NULL,
                                   .scale_x = NULL,
                                   .scale_y = NULL,
                                   .descriptives = NULL,
                                   .validate = function() {
                                     errors <- character()
                                     # Check scales
                                     if (length(private$.scale_x) != 1 || private$.scale_x <= 0) {
                                       errors <- c(errors, "'scale_x' must be a single positive number.")
                                     }
                                     if (length(private$.scale_y) != 1 || private$.scale_y <= 0) {
                                       errors <- c(errors, "'scale_y' must be a single positive number.")
                                     }
                                     
                                     # Check types
                                     valid_types <- c("binomial", "poisson", "normal")
                                     if (!private$.type_x %in% valid_types) {
                                       errors <- c(errors, "type_x must be one of 'binomial', 'poisson', or 'normal'.")
                                     }
                                     if (!private$.type_y %in% valid_types) {
                                       errors <- c(errors, "type_y must be one of 'binomial', 'poisson', or 'normal'.")
                                     }
                                     
                                     # Check attribute lengths
                                     if (length(private$.x_attribute) != private$.n_actor) {
                                       errors <- c(errors, "Length of 'x_attribute' must be equal to 'n_actor'.")
                                     }
                                     if (length(private$.y_attribute) != private$.n_actor) {
                                       errors <- c(errors, "Length of 'y_attribute' must be equal to 'n_actor'.")
                                     }
                                     if(!inherits(private$.descriptives,"list")){
                                       errors <- c(errors, "'descriptives' must be a list.")
                                     }
                                     
                                     # Check attribute value constraints
                                     if (private$.type_x == "binomial" && !all(private$.x_attribute %in% c(0, 1))) {
                                       errors <- c(errors, "For 'binomial' type, 'x_attribute' must be a binary vector.")
                                     }
                                     # browser()
                                     if (private$.type_x == "poisson" && !all(floor(private$.x_attribute) == private$.x_attribute & private$.x_attribute >= 0)) {
                                       errors <- c(errors, "For 'poisson' type, 'x_attribute' must be a vector of non-negative integers.")
                                     }
                                     if (private$.type_y == "binomial" && !all(private$.y_attribute %in% c(0, 1))) {
                                       errors <- c(errors, "For 'binomial' type, 'y_attribute' must be a binary vector.")
                                     }
                                     # browser()
                                     if (private$.type_y == "poisson" && !all(floor(private$.y_attribute) == private$.y_attribute & private$.y_attribute >= 0)) {
                                       errors <- c(errors, "For 'poisson' type, 'y_attribute' must be a vector of non-negative integers.")
                                     }
                                     # Check z_network format
                                     if (!is.matrix(private$.z_network) && !inherits(private$.z_network, "Matrix")) {
                                       errors <- c(errors, "'z_network' must be a matrix or a sparse Matrix object.")
                                     } else {
                                       if (ncol(private$.z_network) == 2) {
                                         if (any(private$.z_network < 1) || any(private$.z_network > private$.n_actor)) {
                                           errors <- c(errors, "'z_network' edge list contains invalid actor indices.")
                                         }
                                       } else if (ncol(private$.z_network) != private$.n_actor) {
                                         errors <- c(errors, "'z_network' must be either an edge list with 2 columns or an adjacency matrix of size n_actor x n_actor.")
                                       }
                                     }
                                     # Check neighborhood format
                                     if (!is.null(private$.neighborhood)) {
                                       if (!is.matrix(private$.neighborhood) && !inherits(private$.neighborhood, "Matrix")) {
                                         errors <- c(errors, "'neighborhood' must be a matrix or a sparse Matrix object.")
                                       } else {
                                         if (ncol(private$.neighborhood) == 2) {
                                           if (any(private$.neighborhood < 1) || any(private$.neighborhood > private$.n_actor)) {
                                             errors <- c(errors, "'neighborhood' edge list contains invalid actor indices.")
                                           }
                                         } else if (ncol(private$.neighborhood) != private$.n_actor) {
                                           errors <- c(errors, "'neighborhood' must be either an edge list with 2 columns or an adjacency matrix of size n_actor x n_actor.")
                                         }
                                       }
                                     }
                                     # Check directed flag
                                     if (!is.logical(private$.directed) || length(private$.directed) != 1) {
                                       errors <- c(errors, "'directed' must be a single logical value (TRUE or FALSE).")
                                     }
                                     if (length(errors) > 0) stop(paste(errors, collapse = "\n"))
                                     
                                     invisible(self)
                                   }
                                 ),
                                 
                                 public = list(
                                   #' @description
                                   #' Create a new `iglm.data` object, that includes data on two attributes and one network. 
                                   #'
                                   #' @param x_attribute A numeric vector for the first unit-level attribute.
                                   #' @param y_attribute A numeric vector for the second unit-level attribute.
                                   #' @param z_network A matrix representing the network. Can be a 2-column
                                   #'   edgelist or a square adjacency matrix.
                                   #' @param neighborhood An optional matrix for the neighborhood representing local dependence.
                                   #'   Can be a 2-column edgelist or a square adjacency matrix. 
                                   #'   A tie in `neighborhood` between actor i and j indicates that j is in the neighborhood of i, 
                                   #'   implying dependence between the respective actors. 
                                   #' @param directed A logical value indicating if `z_network` is directed.
                                   #'   If `NA` (default), directedness is inferred from the symmetry of
                                   #'   `z_network`.
                                   #' @param n_actor An integer for the number of actors in the system.
                                   #'   If `NA` (default), `n_actor` is inferred from the attributes or
                                   #'   network matrices.
                                   #' @param type_x Character string for the type of `x_attribute`.
                                   #'   Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
                                   #'   Default is `"binomial"`.
                                   #' @param type_y Character string for the type of `y_attribute`.
                                   #'   Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
                                   #'   Default is `"binomial"`.
                                   #' @param scale_x A positive numeric value for scaling (e.g., variance
                                   #'   for "normal" type). Default is 1.
                                   #' @param scale_y A positive numeric value for scaling (e.g., variance
                                   #'   for "normal" type). Default is 1.
                                   #' @param return_neighborhood Logical. If `TRUE` (default) and
                                   #'   `neighborhood` is `NULL`, a full neighborhood (all dyads) is
                                   #'   generated implying global dependence. If `FALSE`, no neighborhood is set.
                                   #' @param file (character) Optional file path to load a saved `iglm.data` object state.  
                                   #' @return A new `iglm.data` object.
                                   initialize = function(x_attribute = NULL, y_attribute = NULL, z_network = NULL,
                                                         neighborhood = NULL, directed = NA, n_actor = NA, 
                                                         type_x = "binomial", type_y = "binomial",
                                                         scale_x = 1, scale_y = 1, return_neighborhood = TRUE, 
                                                         file = NULL) {
                                     if(!is.null(file)){
                                       if(!file.exists(file)){
                                         stop(paste("File",file,"does not exist."))
                                       }
                                       data_loaded <- readRDS(file)
                                       required_fields <- c("x_attribute", "y_attribute", "z_network",
                                                            "neighborhood", "directed", "n_actor",
                                                            "type_x", "type_y", "scale_x", "scale_y")
                                       if (!is.list(data_loaded) || !all(required_fields %in% names(data_loaded))) {
                                         stop("File does not contain a valid iglm.data state.", call. = FALSE)
                                       }
                                       x_attribute = data_loaded$x_attribute
                                       y_attribute = data_loaded$y_attribute
                                       z_network = data_loaded$z_network
                                       neighborhood = data_loaded$neighborhood
                                       directed = data_loaded$directed
                                       n_actor = data_loaded$n_actor
                                       type_x = data_loaded$type_x
                                       type_y = data_loaded$type_y
                                       scale_x = data_loaded$scale_x
                                       scale_y = data_loaded$scale_y
                                     } 
                                     private$.type_x <- type_x
                                     private$.type_y <- type_y
                                     private$.scale_x <- scale_x
                                     private$.scale_y <- scale_y
                                     private$.descriptives <- list()
                                   
                                     if(return_neighborhood){
                                       if(is.null(neighborhood)){
                                         if(is.na(n_actor)){
                                           stop("n_actor must be provided if neighborhood is not provided.")
                                         }
                                         neighborhood <- expand.grid(1:n_actor,1:n_actor)
                                         neighborhood <- neighborhood[neighborhood$Var1 != neighborhood$Var2, ]
                                         neighborhood <- as.matrix(neighborhood)
                                       }
                                     }
                                     if(is.null(z_network)){
                                       private$.z_network <- matrix(0,nrow = 0,ncol = 2)
                                     } else {
                                       private$.z_network<- z_network
                                     }
                                     if(ncol(private$.z_network) > 2){
                                       if(!directed){
                                         private$.z_network[lower.tri(private$.z_network)] <- 0
                                       }
                                       private$.z_network <- which(private$.z_network == 1, arr.ind = T)
                                     }
                                     
                                     if(!directed){
                                       wrong_tmp <- private$.z_network[,1] > private$.z_network[,2]
                                       correct_tmp <- private$.z_network[,1] < private$.z_network[,2]
                                       private$.z_network <- rbind(private$.z_network[correct_tmp, c(1,2)], 
                                                                   private$.z_network[wrong_tmp, c(2,1)])
                                       private$.z_network <- private$.z_network[!duplicated(private$.z_network), ]
                                     }
                                     
                                     if(is.na(n_actor)){
                                       if(return_neighborhood){
                                         if(ncol(neighborhood)>2){
                                           private$.n_actor <- nrow(neighborhood)  
                                         } else {
                                           private$.n_actor <- max(neighborhood)
                                         }
                                       } else if(!is.null(x_attribute)) {
                                         private$.n_actor = length(x_attribute)
                                       } else if(!is.null(y_attribute)) {
                                         private$.n_actor = length(y_attribute)
                                       } else if(ncol(z_network)>2){
                                         private$.n_actor = nrow(z_network)  
                                       } else {
                                         private$.n_actor = max(z_network)
                                       }
                                     } else {
                                       private$.n_actor = n_actor
                                     }
                                     if(is.na(private$.n_actor)){
                                       stop("n_actor could not be inferred. Please provide n_actor.")
                                     }
                                     if(is.null(x_attribute)){
                                       private$.x_attribute = numeric(length = private$.n_actor)
                                     } else {
                                       private$.x_attribute = x_attribute
                                     }
                                     if(is.null(y_attribute)){
                                       private$.y_attribute = numeric(length = private$.n_actor)
                                     } else {
                                       private$.y_attribute = y_attribute
                                     }
                                     if(is.na(directed)){
                                       if(ncol(private$.z_network) == 2){
                                         z_network_tmp = matrix(0,nrow = private$.n_actor,ncol = private$.n_actor)
                                         z_network_tmp[private$.z_network] = 1
                                         private$.directed = !isSymmetric(z_network_tmp)
                                       } else {
                                         private$.directed = !isSymmetric(private$.z_network)  
                                       }
                                     } else {
                                       private$.directed = directed
                                     }
                                     if(return_neighborhood){
                                       if(ncol(neighborhood) == 2){
                                         sp_nb = spMatrix(nrow = private$.n_actor,ncol = private$.n_actor, 
                                                          i = neighborhood[,1], j = neighborhood[,2],
                                                          x = rep(1,length(neighborhood[,2])))
                                         sp_nb_trans = sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)
                                         
                                         overlap =  sp_nb %*%sp_nb_trans
                                         overlap = as(overlap, "TsparseMatrix")
                                         overlap = cbind(overlap@i + 1,
                                                         overlap@j+1)
                                         private$.overlap = overlap[overlap[,1] != overlap[,2],]
                                         private$.neighborhood = neighborhood
                                         
                                       } else {
                                         positions = which(neighborhood == 1, arr.ind = T)
                                         sp_nb = spMatrix(nrow = ncol(neighborhood),ncol = ncol(neighborhood), 
                                                          i = positions[,1], j = positions[,2],x = rep(1,length(positions[,2])))
                                         sp_nb_trans = sparseMatrix(i = sp_nb@j + 1, j = sp_nb@i + 1, dims = sp_nb@Dim)
                                         
                                         overlap =  as.matrix(sp_nb %*%sp_nb_trans>0)
                                         diag(overlap) = 0
                                         private$.overlap = which(overlap == 1, arr.ind = T)
                                         private$.neighborhood = which(neighborhood == 1, arr.ind = T)
                                       }
                                     } else{
                                       private$.overlap = matrix(0,nrow = 0,ncol = 2)
                                     }
                                     private$.validate()
                                     
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `z_network` of the `iglm.data` object.
                                   #' @param z_network A matrix representing the network. Can be a 2-column
                                   #'  edgelist or a square adjacency matrix.
                                   #'  @return The `iglm.data` object itself (`self`), invisibly.
                                   set_z_network = function(z_network){
                                     if(ncol(z_network) > 2){
                                       if(!private$.directed){
                                         z_network[lower.tri(z_network)] <- 0
                                       }
                                       z_network <- which(z_network == 1, arr.ind = T)
                                     }
                                     private$.z_network <- z_network
                                     private$.validate()
                                     invisible(self)
                                   },
                                   
                                   #' @description
                                   #' Sets the `type_x` of the `iglm.data` object.
                                   #' @param type_x A character string for the type of `x_attribute`.
                                   #'  Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
                                   #'  @return The `iglm.data` object itself (`self`), invisibly.
                                   set_type_x = function(type_x){
                                     private$.type_x <- type_x
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `type_y` of the `iglm.data` object.
                                   #' @param type_y A character string for the type of `y_attribute`.
                                   #' Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
                                   #' @return The `iglm.data` object itself (`self`), invisibly.
                                   set_type_y = function(type_y){
                                     private$.type_y <- type_y
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `scale_x` of the `iglm.data` object.
                                   #' @param scale_x A positive numeric value for scaling (e.g., variance
                                   #' for "normal" type).
                                   #' @return The `iglm.data` object itself (`self`), invisibly.
                                   set_scale_x = function(scale_x){
                                     private$.scale_x <- scale_x
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `scale_y` of the `iglm.data` object.
                                   #' @param scale_y A positive numeric value for scaling (e.g., variance
                                   #' for "normal" type).
                                   #' @return The `iglm.data` object itself (`self`), invisibly.
                                   set_scale_y = function(scale_y){
                                     private$.scale_y <- scale_y
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `x_attribute` of the `iglm.data` object.
                                   #' @param x_attribute A numeric vector for the first unit-level attribute.
                                   #' @return The `iglm.data` object itself (`self`), invisibly.
                                   set_x_attribute = function(x_attribute){
                                     private$.x_attribute <- x_attribute
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Sets the `y_attribute` of the `iglm.data` object.
                                   #' @param y_attribute A numeric vector for the first unit-level attribute.
                                   #' @return The `iglm.data` object itself (`self`), invisibly.
                                   set_y_attribute = function(y_attribute){
                                     private$.y_attribute <- y_attribute
                                     private$.validate()
                                     invisible(self)
                                   },
                                   #' @description
                                   #' Gathers the current state of the `iglm.data` object into a list.
                                   #' This includes all attributes, network, and configuration
                                   #' details necessary to reconstruct the object later.
                                   #' @return A list containing the current state of the `iglm.data` object.
                                   gather = function(){
                                     data_to_save <- list(
                                       x_attribute = private$.x_attribute,
                                       y_attribute = private$.y_attribute,
                                       z_network = private$.z_network,
                                       neighborhood = private$.neighborhood,
                                       directed = private$.directed,
                                       n_actor = private$.n_actor,
                                       type_x = private$.type_x,
                                       type_y = private$.type_y,
                                       scale_x = private$.scale_x,
                                       scale_y = private$.scale_y
                                     )
                                     return(data_to_save)
                                   },
                                   #' @description
                                   #' Saves the current state of the `iglm.data` object to a specified file path
                                   #' in RDS format. This includes all attributes, network, and configuration
                                   #' details necessary to reconstruct the object later.
                                   #' @param file (character) The file where the object state should be saved.
                                   #' @return The `iglm.data` object itself (`self`), invisibly. 
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
                                   #' Calculates the density of the `z_network`.
                                   #' @return A numeric value for the network density.
                                   mean_z = function(){
                                     m = nrow(private$.z_network)/(private$.n_actor*(private$.n_actor-1)/(2 - private$.directed))
                                     private$.descriptives$density_z <- m
                                     invisible(m)
                                   },
                                   #' @description
                                   #' Calculates the mean of the `x_attribute`.
                                   #' @return A numeric value for the mean of `x_attribute`.
                                   mean_x = function(){
                                     m = mean(private$.x_attribute)
                                     private$.descriptives$density_x <- m
                                     invisible(m)
                                   },
                                   #' @description
                                   #' Calculates the mean of the `y_attribute`.
                                   #' @return A numeric value for the mean of `y_attribute`.
                                   mean_y = function(){
                                     m = mean(private$.y_attribute)
                                     private$.descriptives$density_y <- m
                                     invisible(m)
                                   },
                                   #' @description
                                   #' Calculates the distribution of the `x_attribute`.
                                   #' @param value_range (numeric vector) Optional range of values to consider for the distribution. If `NULL` (default), the range is inferred from the data.
                                   #' @param prob (logical) If `TRUE` (default), returns probabilities; if `FALSE`, returns frequencies.
                                   #' @param plot (logical) If `TRUE` (default), plots the distribution using a density plot for continuous data or a bar plot for discrete data.
                                   #' @return A numeric vector representing the distribution of `x_attribute` (invisible).
                                   x_distribution = function(value_range = NULL, prob = TRUE, plot = FALSE){
                                     if(private$.type_x == "gaussian"){
                                       tmp_density <- density(private$.x_attribute, from = value_range[1], to = value_range[2])
                                       names(tmp_density$y) <- tmp_density$x
                                       
                                       if(plot){
                                         plot(tmp_density, main = "Density of x_attribute", 
                                              xlab = "x_attribute values", ylab = "Density")
                                       }
                                       private$.descriptives$x_distribution <- tmp_density$y
                                     } else {
                                       if(is.null(value_range)){
                                         value_range <- range(private$.x_attribute)
                                       }
                                       info <- factor(as.numeric(private$.x_attribute), 
                                                      levels = seq(from = value_range[1], to = value_range[2]))
                                       info <- table(info)
                                       if(sum(info) > 0){
                                         info <- info/(sum(info)*prob + (!prob))
                                       }
                                       private$.descriptives$x_distribution <- info
                                       if(plot){
                                         barplot(info, main = "Distribution of x_attribute", 
                                                 xlab = "x_attribute values", ylab = ifelse(prob, "Probability", "Frequency"))
                                       }
                                     }
                                     invisible(private$.descriptives$x_distribution)
                                   },
                                   #' @description
                                   #' Calculates the distribution of the `y_attribute`.
                                   #' @param value_range (numeric vector) Optional range of values to consider for the distribution. If `NULL` (default), the range is inferred from the data.
                                   #' @param prob (logical) If `TRUE` (default), returns probabilities; if `FALSE`, returns frequencies.
                                   #' @param plot (logical) If `TRUE` (default), plots the distribution using a density plot for continuous data or a bar plot for discrete data.
                                   #' @return A numeric vector representing the distribution of `y_attribute` (invisible).
                                   y_distribution = function(value_range = NULL, prob = TRUE, plot = FALSE){
                                     if(private$.type_y == "gaussian"){
                                       tmp_density <- density(private$.y_attribute, from = value_range[1], to = value_range[2])
                                       names(tmp_density$y) <- tmp_density$x
                                       
                                       if(plot){
                                         plot(tmp_density, main = "Density of y_attribute", 
                                              xlab = "y_attribute values", ylab = "Density")
                                       }
                                       private$.descriptives$y_distribution <- tmp_density$y
                                      } else {
                                        if(is.null(value_range)){
                                          value_range <- range(private$.y_attribute)
                                        }
                                        info <- factor(as.numeric(private$.y_attribute), 
                                                      levels = seq(from = value_range[1], to = value_range[2]))
                                        info <- table(info)
                                        if(sum(info) > 0){
                                          info <- info/(sum(info)*prob + (!prob))
                                        }
                                        private$.descriptives$y_distribution <- info
                                        if(plot){
                                          barplot(info, main = "Distribution of y_attribute", 
                                                  xlab = "y_attribute values", ylab = ifelse(prob, "Probability", "Frequency"))
                                        }
                                      }
                                     invisible(private$.descriptives$y_distribution)
                                   },
                                   #' @description
                                   #' Calculates the matrix of edgewise shared partners.
                                   #' This is a two-path matrix (e.g., $A A^T$ or $A^T A$).
                                   #'
                                   #' @param type (character) The type of two-path to calculate for directed
                                   #'   networks. Ignored if network is undirected.
                                   #'   Must be one of:
                                   #'   `"OTP"` (Outgoing Two-Path),
                                   #'   `"ISP"` (In-Star),
                                   #'   `"OSP"` (Out-Star),
                                   #'   `"ITP"` (Incoming Two-Path),
                                   #'   `"ALL"` (Symmetric all-partner).
                                   #'   Default is `"ALL"`.
                                   #' @return A sparse matrix (`dgCMatrix`) of shared partner counts.
                                   edgewise_shared_partner = function(type = "ALL"){
                                     if(is.null(private$.descriptives$edgewise_shared_partner)){
                                       private$.descriptives$edgewise_shared_partner <- list()
                                     }
                                     if(!private$.directed){
                                       res <- self$dyadwise_shared_partner(type = "ALL")
                                       res <- res[private$.z_network]
                                       private$.descriptives$edgewise_shared_partner$ALL<- res
                                     } else {
                                       if(type == "OTP"){
                                         res <- self$dyadwise_shared_partner(type = "OTP")
                                         res <- res[private$.z_network]
                                         private$.descriptives$edgewise_shared_partner$OTP<- res
                                         return(res)
                                       }  else if(type == "ISP"){
                                         res <- self$dyadwise_shared_partner(type = "ISP")
                                         res <- res[private$.z_network]
                                         private$.descriptives$edgewise_shared_partner$ISP<- res
                                       }else if(type == "OSP"){
                                         res <- self$dyadwise_shared_partner(type = "OSP")
                                         res <- res[private$.z_network]
                                         private$.descriptives$edgewise_shared_partner$OSP<- res
                                       } else if(type == "ITP"){
                                         res <- self$dyadwise_shared_partner(type = "ITP")
                                         res <- res[private$.z_network]
                                         private$.descriptives$edgewise_shared_partner$ITP<- res
                                       } else if(type == "ALL"){
                                         res <- self$dyadwise_shared_partner(type = "ALL")
                                         res <- res[private$.z_network]
                                         private$.descriptives$edgewise_shared_partner$ALL<- res
                                       }  else{
                                         stop("type must be one of 'OTP', 'ISP', 'OSP', 'ITP', or 'ALL'.")
                                       }
                                     }
                                     invisible(res)
                                   }
                                   ,
                                   #' @description
                                   #' Sets the neighborhood and overlap matrices.
                                   #' @param neighborhood A matrix for a secondary neighborhood.
                                   #'  Can be a 2-column edgelist or a square adjacency matrix.
                                   #' @param overlap A matrix for the overlap network.
                                   #'  Can be a 2-column edgelist or a square adjacency matrix.
                                   #' @return None. Updates the internal neighborhood and overlap matrices.
                                   set_neighborhood_overlap = function(neighborhood, overlap){
                                     if(!is.matrix(neighborhood)){
                                       stop("'neighborhood' must be a matrix or a sparse Matrix object.")
                                     }
                                     if(ncol(neighborhood) > 2){
                                       private$.neighborhood = which(neighborhood == 1, arr.ind = T)
                                     } else {
                                       private$.neighborhood = neighborhood
                                     }
                                     
                                     if(!is.matrix(overlap)){
                                       stop("'overlap' must be a matrix or a sparse Matrix object.")
                                     }
                                     if(ncol(overlap) > 2){
                                       private$.overlap = which(overlap == 1, arr.ind = T)
                                     } else {
                                       private$.overlap = overlap
                                     }
                                     
                                   },
                                   #' @description
                                   #' Calculates the matrix of edgewise shared partners.
                                   #' This is a two-path matrix (e.g., $A A^T$ or $A^T A$).
                                   #'
                                   #' @param type (character) The type of two-path to calculate for directed
                                   #'   networks. Ignored if network is undirected.
                                   #'   Must be one of:
                                   #'   `"OTP"` (Outgoing Two-Path, z_i,j*z_j,h ),
                                   #'   `"ISP"` (In-Star),
                                   #'   `"OSP"` (Out-Star),
                                   #'   `"ITP"` (Incoming Two-Path),
                                   #'   `"ALL"` (Symmetric all-partner).
                                   #'   Default is `"ALL"`.
                                   #' @return A sparse matrix (`dgCMatrix`) of shared partner counts.
                                   dyadwise_shared_partner = function(type = "ALL"){
                                     
                                     adj_mat <- sparseMatrix(i = private$.z_network[,1],
                                                             j = private$.z_network[,2],symmetric = !private$.directed,
                                                             dims = c(private$.n_actor, private$.n_actor))
                                     if(is.null(private$.descriptives$dyadwise_shared_partner)){
                                       private$.descriptives$dyadwise_shared_partner <- list()
                                     }
                                     if(!private$.directed){
                                       res <- Matrix::t(adj_mat) %*% t(adj_mat)
                                       private$.descriptives$dyadwise_shared_partner$ALL<- res
                                     } else {
                                       if(type == "OTP"){
                                         res <- adj_mat %*% Matrix::t(adj_mat)
                                         private$.descriptives$dyadwise_shared_partner$OTP<- res
                                         return(res)
                                       }  else if(type == "ISP"){
                                         res <- adj_mat %*% adj_mat
                                         private$.descriptives$dyadwise_shared_partner$ISP<- res
                                       }else if(type == "OSP"){
                                         res <- Matrix::t(adj_mat) %*% Matrix::t(adj_mat)
                                         private$.descriptives$dyadwise_shared_partner$OSP<- res
                                       } else if(type == "ITP"){
                                         res <- Matrix::t(adj_mat) %*% adj_mat
                                         private$.descriptives$dyadwise_shared_partner$ITP<- res
                                       } else if(type == "ALL"){
                                         adj_mat_symm = pmax(adj_mat, t(adj_mat))
                                         res <- Matrix::t(adj_mat_symm) %*% adj_mat_symm
                                         private$.descriptives$dyadwise_shared_partner$ALL<- res
                                       }  else{
                                         stop("type must be one of 'OTP', 'ISP', 'OSP', 'ITP', or 'ALL'.")
                                       }
                                     }
                                     invisible(res)
                                   }
                                   ,
                                   #' @description
                                   #' Calculates the geodesic distance distribution of the symmetrized
                                   #' `z_network`.
                                   #'
                                   #' @param value_range (numeric vector) A vector `c(min, max)` specifying
                                   #'   the range of distances to tabulate. If `NULL` (default), the range
                                   #'   is inferred from the data.
                                   #' @param plot (logical) If `TRUE`, plots the distribution.
                                   #' @param prob (logical) If `TRUE` (default), returns a probability
                                   #'   distribution (proportions). If `FALSE`, returns raw counts.
                                   #' @return A named vector (a `table` object) with the distribution of
                                   #'   geodesic distances. Includes `Inf` for unreachable pairs.
                                   geodesic_distances_distribution = function(value_range = NULL, prob = TRUE, plot= FALSE){
                                     if(is.null(private$.descriptives$geodesic_distances)){
                                       self$geodesic_distances()
                                     }
                                     D <- private$.descriptives$geodesic_distances
                                     D_vec <- as.vector(D)
                                     if(is.null(value_range)){
                                       if(length(D_vec[is.finite(D_vec) & D_vec>0])>0){
                                         value_range = range(D_vec[is.finite(D_vec) & D_vec>0])  
                                       } else {
                                         value_range = c(1,1)
                                       }
                                       
                                     }
                                     info_factor <- factor(D_vec, 
                                                           levels = c(seq(from = value_range[1], to = value_range[2]), Inf))
                                     info <- table(info_factor)
                                     
                                   
                                     if(sum(info) > 0){
                                       info <- info/(sum(info)*prob + (!prob))
                                     }
                                     private$.descriptives$geodesic_distances_distribution <- info
                                     if(plot){
                                       barplot(info,
                                               xlab = "Geodesic Distance", ylab = ifelse(prob,"Proportion","Count"), 
                                               las = 1 )
                                     }
                                     invisible(info)
                                   },
                                   #' @description
                                   #' Calculates the all-pairs geodesic distance matrix for the
                                   #' symmetrized `z_network` using a matrix-based BFS algorithm.
                                   #' @return A sparse matrix (`dgCMatrix`) where `D[i, j]` is the
                                   #'   shortest path distance from i to j. `Inf` indicates no path.
                                   #' @importFrom Matrix sparseMatrix t Matrix diag nnzero Diagonal
                                   geodesic_distances = function() {
                                     adj_mat <- Matrix::sparseMatrix(i = private$.z_network[,1],
                                                             j = private$.z_network[,2],
                                                             dims = c(private$.n_actor, private$.n_actor))
                                     # Make symmetric
                                     adj_mat = pmax(adj_mat, Matrix::t(adj_mat))
                                     D <- Matrix::Matrix(Inf, nrow(adj_mat), nrow(adj_mat), dimnames = dimnames(adj_mat))
                                     # Distance to oneself is 0
                                     Matrix::diag(D) <- 0
                                     F_tmp <- Matrix::Diagonal( nrow(adj_mat))
                                     # k = current path length
                                     k <- 0
                                     # Loop while any search has a non-empty frontier
                                     while(Matrix::nnzero(F_tmp) > 0) {
                                       k <- k + 1
                                       F_next <- (adj_mat %*% F_tmp) > 0
                                       F_new <- F_next & (D == Inf)
                                       D[F_new] <- k
                                       F_tmp <- F_new
                                     }
                                     private$.descriptives$geodesic_distances <- D
                                     invisible(D)
                                   },
                                   #' @description
                                   #' Calculates the distribution of edgewise shared partners.
                                   #'
                                   #' @param type (character) The type of shared partner matrix to use.
                                   #'   See `edgewise_shared_partner` for details. Default is `"ALL"`.
                                   #' @param value_range (numeric vector) A vector `c(min, max)` specifying
                                   #'   the range of counts to tabulate. If `NULL` (default), the range
                                   #'   is inferred from the data.
                                   #' @param prob (logical) If `TRUE` (default), returns a probability
                                   #'   distribution (proportions). If `FALSE`, returns raw counts.
                                   #' @param plot (logical) If `TRUE`, plots the distribution.
                                   #' @return A named vector (a `table` object) with the distribution of
                                   #'   shared partner counts.
                                   edgewise_shared_partner_distribution = function(type  ="ALL",
                                                                                   value_range = NULL, 
                                                                                   prob = TRUE, 
                                                                                   plot = FALSE){
                                     
                                     if(!type %in% c("OTP","ITP", "ISP", "OSP", "ALL")){
                                       stop("type must be one of 'OTP', 'ISP', 'OSP', 'ITP', or 'ALL'.")
                                     }
                                     if(length(value_range) != 2 & !is.null(value_range)){
                                       stop("'value_range' must be a numeric vector of length 2.")
                                     }
                                     if(sum(value_range<0)>0){
                                       stop("'value_range' values must be non-negative.")
                                     }
                                     if(is.null(private$.descriptives$edgewise_shared_partner_dist)){
                                       private$.descriptives$edgewise_shared_partner_dist <- list()
                                     }
                                     
                                     if(type  == "OTP"){
                                       # Calculate the ESP if not already done so
                                       if(is.null(private$.descriptives$edgewise_shared_partner$OTP)){
                                         self$edgewise_shared_partner(type = "OTP")
                                       }
                                       info <- private$.descriptives$edgewise_shared_partner$OTP
                                     } else if(type  == "ITP"){
                                       if(is.null(private$.descriptives$edgewise_shared_partner$ITP)){
                                         self$edgewise_shared_partner(type = "ITP")
                                       }
                                       info <- private$.descriptives$edgewise_shared_partner$ITP
                                     } else if(type  == "ISP"){
                                       if(is.null(private$.descriptives$edgewise_shared_partner$ISP)){
                                         self$edgewise_shared_partner(type = "ISP")
                                       }
                                       info <- private$.descriptives$edgewise_shared_partner$ISP
                                     } else if(type  == "OSP"){
                                       if(is.null(private$.descriptives$edgewise_shared_partner$OSP)){
                                         self$edgewise_shared_partner(type = "OSP")
                                       }
                                       info <- private$.descriptives$edgewise_shared_partner$OSP
                                     } else if(type  == "ALL"){
                                       if(is.null(private$.descriptives$edgewise_shared_partner$ALL)){
                                         self$edgewise_shared_partner(type = "ALL")
                                       }
                                       info <- private$.descriptives$edgewise_shared_partner$ALL
                                     }
                                     if(is.null(value_range)){
                                       value_range = range(as.numeric(info))
                                     }
                                     info <- factor(as.numeric(info), 
                                                    levels = seq(from = value_range[1], to = value_range[2]))
                                     info <- table(info)
                                     # Transform to probability from frequency 
                                     if(sum(info) > 0){
                                       info <- info/(sum(info)*prob + (!prob))
                                     }
                                     if(type  == "OTP"){
                                       private$.descriptives$edgewise_shared_partner_distribution$OTP <- info
                                     } else if(type  == "ITP"){
                                       private$.descriptives$edgewise_shared_partner_distribution$ITP<- info
                                     } else if(type  == "ISP"){
                                       private$.descriptives$edgewise_shared_partner_distribution$ISP<- info
                                     } else if(type  == "OSP"){
                                       private$.descriptives$edgewise_shared_partner_distribution$OSP<- info
                                     } else if(type  == "ALL"){
                                       private$.descriptives$edgewise_shared_partner_distribution$ALL<- info
                                     }
                                     if(plot){
                                       barplot(info, 
                                               xlab = paste0("Number of ", type, "- Edgewise Shared Partners"),
                                               ylab = ifelse(prob, "Proportion", "Count"), 
                                               las = 1 )
                                     }
                                     invisible(info)
                                   },
                                   #' @description
                                   #' Calculates the distribution of edgewise shared partners.
                                   #'
                                   #' @param type (character) The type of shared partner matrix to use.
                                   #'   See `edgewise_shared_partner` for details. Default is `"ALL"`.
                                   #' @param value_range (numeric vector) A vector `c(min, max)` specifying
                                   #'   the range of counts to tabulate. If `NULL` (default), the range
                                   #'   is inferred from the data.
                                   #' @param plot (logical) If `TRUE`, plots the distribution.
                                   #' @param prob (logical) If `TRUE` (default), returns a probability
                                   #'   distribution (proportions). If `FALSE`, returns raw counts.
                                   #' @return A named vector (a `table` object) with the distribution of
                                   #'   shared partner counts.
                                   dyadwise_shared_partner_distribution = function(type  ="ALL", 
                                                                                   value_range = NULL,
                                                                                   prob = TRUE, plot = FALSE){
                                     
                                     if(!type %in% c("OTP","ITP", "ISP", "OSP", "ALL")){
                                       stop("type must be one of 'OTP', 'ISP', 'OSP', 'ITP', or 'ALL'.")
                                     }
                                     if(length(value_range) != 2 & !is.null(value_range)){
                                       stop("'value_range' must be a numeric vector of length 2.")
                                     }
                                     if(sum(value_range<0)>0){
                                       stop("'value_range' values must be non-negative.")
                                     }
                                     if(is.null(private$.descriptives$dyadwise_shared_partner_dist)){
                                       private$.descriptives$dyadwise_shared_partner_dist <- list()
                                     }
                                     
                                     if(type  == "OTP"){
                                       # Calculate the ESP if not already done so
                                       if(is.null(private$.descriptives$dyadwise_shared_partner$OTP)){
                                         self$dyadwise_shared_partner(type = "OTP")
                                       }
                                       info <- private$.descriptives$dyadwise_shared_partner$OTP
                                     } else if(type  == "ITP"){
                                       if(is.null(private$.descriptives$dyadwise_shared_partner$ITP)){
                                         self$dyadwise_shared_partner(type = "ITP")
                                       }
                                       info <- private$.descriptives$dyadwise_shared_partner$ITP
                                     } else if(type  == "ISP"){
                                       if(is.null(private$.descriptives$dyadwise_shared_partner$ISP)){
                                         self$dyadwise_shared_partner(type = "ISP")
                                       }
                                       info <- private$.descriptives$dyadwise_shared_partner$ISP
                                     } else if(type  == "OSP"){
                                       if(is.null(private$.descriptives$dyadwise_shared_partner$OSP)){
                                         self$dyadwise_shared_partner(type = "OSP")
                                       }
                                       info <- private$.descriptives$dyadwise_shared_partner$OSP
                                     } else if(type  == "ALL"){
                                       if(is.null(private$.descriptives$dyadwise_shared_partner$ALL)){
                                         self$dyadwise_shared_partner(type = "ALL")
                                       }
                                       info <- private$.descriptives$dyadwise_shared_partner$ALL
                                     }
                                     if(is.null(value_range)){
                                       value_range = range(as.numeric(info))
                                     }
                                     info <- factor(as.numeric(info), 
                                                    levels = seq(from = value_range[1], to = value_range[2]))
                                     info <- table(info)
                                     # Transform to probability from frequency
                                     if(sum(info) > 0){
                                       info <- info/(sum(info)*prob + (!prob))
                                     }
                                     
                                     if(type  == "OTP"){
                                       private$.descriptives$dyadwise_shared_partner_distribution$OTP <- info
                                     } else if(type  == "ITP"){
                                       private$.descriptives$dyadwise_shared_partner_distribution$ITP<- info
                                     } else if(type  == "ISP"){
                                       private$.descriptives$dyadwise_shared_partner_distribution$ISP<- info
                                     } else if(type  == "OSP"){
                                       private$.descriptives$dyadwise_shared_partner_distribution$OSP<- info
                                     } else if(type  == "ALL"){
                                       private$.descriptives$dyadwise_shared_partner_distribution$ALL<- info
                                     }
                                     if(plot){
                                       barplot(info, 
                                               xlab = paste0("Number of ",type, "- Dyadwise Shared Partners"), 
                                               ylab = ifelse(prob, "Proportion", "Count"), 
                                               las = 1 )
                                     }
                                     invisible(info)
                                   },
                                   #' @description
                                   #' Calculates the degree distribution of the `z_network`.
                                   #'
                                   #' @param value_range (numeric vector) A vector `c(min, max)` specifying
                                   #'   the range of degrees to tabulate. If `NULL` (default), the range
                                   #'   is inferred from the data.
                                   #' @param prob (logical) If `TRUE` (default), returns a probability
                                   #'   distribution (proportions). If `FALSE`, returns raw counts.
                                   #' @param plot (logical) If `TRUE`, plots the degree distribution.
                                   #' @return If the network is directed, a list containing two `table`
                                   #'   objects: `in_degree` and `out_degree`. If undirected, a single
                                   #'   `table` object with the degree distribution.
                                   degree_distribution = function(value_range = NULL, prob = TRUE, plot = FALSE){
                                     if(is.null(private$.descriptives$degree)){
                                       self$degree()
                                     }
                                     if(is.null(value_range)){
                                       value_range = range(unlist(private$.descriptives$degree))
                                     }
                                     if(private$.directed){
                                      
                                       info_in <- factor(private$.descriptives$degree$in_degree_seq, 
                                                         levels = seq(from = value_range[1], to = value_range[2]))
                                       
                                       info_out <- factor(private$.descriptives$degree$out_degree_seq, 
                                                          levels = seq(from = value_range[1], to = value_range[2]))
                                       info_in <- table(info_in)
                                       info_out <- table(info_out)
                                       
                                       if(sum(info_in) > 0){
                                         info_in <- info_in/(sum(info_in)*prob + (!prob))
                                       } 
                                       if(sum(info_out) > 0){
                                         info_out <- info_out/(sum(info_out)*prob + (!prob))
                                       }
                                       info <- list(in_degree = info_in,
                                                    out_degree = info_out)
                                       private$.descriptives$degree_distribution <- info
                                     } else {
                                       info <- factor(private$.descriptives$degree$degree_seq, 
                                                      levels = seq(from = value_range[1], to = value_range[2]))
                                       if(is.null(value_range)){
                                         value_range = range(as.numeric(info))
                                       }
                                       info <-table(info)
                                       if(sum(info) > 0){
                                         info <- info/(sum(info)*prob + (!prob))  
                                       } 
                                       private$.descriptives$degree_distribution <- info
                                     }
                                     if(plot){
                                       if(private$.directed){
                                         barplot(info$in_degree, 
                                                 xlab = "In-Degree", 
                                                 ylab = ifelse(prob, "Proportion", "Count"), 
                                                 las = 1 )
                                         barplot(info$out_degree,
                                                 xlab = "Out-Degree", 
                                                 ylab = ifelse(prob, "Proportion", "Count"), 
                                                 las = 1 )
                                       } else {
                                         barplot(info, 
                                                 xlab = "Degree",
                                                 ylab = ifelse(prob, "Proportion", "Count"), 
                                                 las = 1 )
                                       }
                                     }
                                     invisible(info)
                                   },
                                   #' @description
                                   #' Calculates the degree sequence(s) of the `z_network`.
                                   #'
                                   #' @return If the network is directed, a list containing two vectors:
                                   #'   `in_degree_seq` and `out_degree_seq`. If undirected, a single
                                   #'   list containing the vector `degree_seq`.
                                   degree = function(){
                                     res = list()
                                     if(private$.directed){
                                       if(ncol(private$.z_network) == 2){
                                         in_degree_seq_res = table(private$.z_network[,2])
                                         out_degree_seq_res = table(private$.z_network[,1])
                                         in_degree_seq = numeric(private$.n_actor)
                                         in_degree_seq[as.numeric(names(in_degree_seq_res))] <- in_degree_seq_res
                                         
                                         out_degree_seq = numeric(private$.n_actor)
                                         out_degree_seq[as.numeric(names(out_degree_seq_res))] <- out_degree_seq_res
                                         # in_zero_values = which(! (1:private$.n_actor %in% names(in_degree_seq)))
                                         # in_degree_seq = c(in_degree_seq, rep(0, length(in_zero_values)))
                                         # out_zero_values = which(! (1:private$.n_actor %in% names(out_degree_seq)))
                                         # out_degree_seq = c(in_degree_seq, rep(0, length(out_zero_values)))
                                       } else {
                                         in_degree_seq = colSums(private$.z_network)
                                         out_degree_seq = rowSums(private$.z_network)
                                         in_zero_values = which(in_degree_seq == 0)
                                         out_zero_values = which(out_degree_seq == 0)
                                       }
                                       res$in_degree_seq = in_degree_seq
                                       res$out_degree_seq = out_degree_seq
                                       
                                     } else {
                                       if(ncol(private$.z_network) == 2){
                                         tmp <- table(private$.z_network)
                                         degree_seq <- numeric(private$.n_actor)
                                         degree_seq[as.numeric(names(tmp))] <- tmp
                                         
                                       }  else {
                                         degree_seq = colSums(private$.z_network)
                                       }
                                       res$degree_seq = degree_seq
                                     }
                                     private$.descriptives$degree = res
                                     invisible(res)
                                   },
                                   #' @description
                                   #' Calculates the spillover degree distribution between actors with
                                   #' `x_attribute == 1` and actors with `y_attribute == 1`.
                                   #'
                                   #' @param prob (logical) If `TRUE` (default), returns a probability
                                   #'   distribution (proportions). If `FALSE`, returns raw counts.
                                   #' @param value_range (numeric vector) A vector `c(min, max)` specifying
                                   #'   the range of degrees to tabulate. If `NULL` (default), the range
                                   #'   is inferred from the data.
                                   #' @param plot (logical) If `TRUE`, plots the distributions.
                                   #' @return A list containing two `table` objects:
                                   #'   `out_spillover_degree` (from x_i=1 to y_j=1) and
                                   #'   `in_spillover_degree` (from y_i=1 to x_j=1).
                                   spillover_degree_distribution = function( prob = TRUE, value_range = NULL, plot = FALSE){
                                     
                                     # actors_x = which(x$x_attribute == 1)
                                     # actors_y = which(x$y_attribute == 1)
                                     # 
                                     # adj_mat_x_y = matrix(data = 0, nrow = length(actors_x), ncol = length(actors_y), dimnames = list(actors_x, actors_y))
                                     # edges_x_y = x$z_network[(x$z_network[,1] %in% actors_x) & (x$z_network[,2] %in% actors_y),]
                                     # edges_x_y[,1] = match(edges_x_y[,1], rownames(adj_mat_x_y))
                                     # edges_x_y[,2] = match(edges_x_y[,2], colnames(adj_mat_x_y))
                                     # adj_mat_x_y[edges_x_y] = 1
                                     # adj_mat_x_y[which(neighborhood[actors_x,as.numeric(colnames(adj_mat_x_y))] == 0, arr.ind = T)] = NA
                                     # adj_mat_x_y = adj_mat_x_y[rowSums(neighborhood[actors_x,as.numeric(colnames(adj_mat_x_y))],na.rm = T) != 0,]
                                     # adj_mat_x_y = adj_mat_x_y[,colSums(neighborhood[actors_x,as.numeric(colnames(adj_mat_x_y))],na.rm = T) != 0]
                                     # 
                                     # out_degree_x_y = table(factor(rowSums(adj_mat_x_y,na.rm = T), levels = 0:max_outdegree_x_y))/nrow(adj_mat_x_y)
                                     # in_degree_x_y = table(factor(colSums(adj_mat_x_y,na.rm = T), levels = 0:max_indegree_x_y))/ncol(adj_mat_x_y)
                                     # 
                                     # browser()
                                     actors_x = which(private$.x_attribute == 1)
                                     actors_y = which(private$.y_attribute == 1)
                                     
                                     adj_mat_x_y = matrix(data = 0, nrow = length(actors_x),
                                                          ncol = length(actors_y),
                                                          dimnames = list(actors_x, actors_y))
                                     edges_x_y = matrix(private$.z_network[(private$.z_network[,1] %in% actors_x) & (private$.z_network[,2] %in% actors_y),],
                                                        ncol = 2)
                                     if(nrow(edges_x_y) > 0){
                                       edges_x_y[,1] = match(edges_x_y[,1], rownames(adj_mat_x_y))
                                       edges_x_y[,2] = match(edges_x_y[,2], colnames(adj_mat_x_y))
                                       which_overlap <- check_overlap(edges_x_y,private$.neighborhood)
                                       edges_x_y_overlap <- edges_x_y[which_overlap,]
                                       edges_x_y_nonoverlap <- edges_x_y[!which_overlap,]
                                       
                                       adj_mat_x_y[edges_x_y_overlap] = 1
                                       adj_mat_x_y[edges_x_y_nonoverlap] = NA
                                       
                                       tmp1 <- rowSums(adj_mat_x_y,na.rm = T)
                                       tmp2 <- colSums(adj_mat_x_y,na.rm = T)
                                       if(is.null(value_range)){
                                         value_range = range(unique(c(tmp1,tmp2)))
                                       }
                                       out_degree_x_y = table(factor(tmp1, levels = seq(from = value_range[1], to = value_range[2])))/
                                         (nrow(adj_mat_x_y)*prob + (!prob))
                                       in_degree_x_y = table(factor(tmp2, levels = seq(from = value_range[1], to = value_range[2])))/
                                         (ncol(adj_mat_x_y)*prob + (!prob))
                                       res <- list(out_spillover_degree = out_degree_x_y,
                                                   in_spillover_degree = in_degree_x_y)
                                     } else {
                                       if(is.null(value_range)){
                                         value_range = c(0,1)
                                       }
                                       out_degree_x_y = table(factor(0, levels = seq(from = value_range[1], to = value_range[2])))
                                       in_degree_x_y = table(factor(0, levels = seq(from = value_range[1], to = value_range[2])))
                                       res <- list(out_spillover_degree = out_degree_x_y,
                                                   in_spillover_degree = in_degree_x_y)
                                     }
                                     
                                     private$.descriptives$spillover_degree_distribution <- res
                                     if(plot){
                                       barplot(out_degree_x_y, 
                                               xlab = "Out-Spillover Degree", 
                                               ylab = ifelse(prob, "Proportion", "Count"), 
                                               las = 1 )
                                       barplot(in_degree_x_y, 
                                               xlab = "In-Spillover Degree", 
                                               ylab = ifelse(prob, "Proportion", "Count"), 
                                               las = 1)
                                     }
                                     invisible(res)
                                    },
                                   #' @description
                                   #' Plot the network using `igraph`.
                                   #'
                                   #' Visualizes the `z_network` using the `igraph` package.
                                   #' Nodes can be colored by `x_attribute` and sized by `y_attribute`.
                                   #' `neighborhood` edges can be plotted as a background layer.
                                   #'
                                   #' @param node_color (character) Attribute to map to node color.
                                   #'   One of `"x"` (default), `"y"`, or `"none"`.
                                   #' @param node_size (character) Attribute to map to node size.
                                   #'   One of `"y"` (default), `"x"`, or `"constant"`.
                                   #' @param show_overlap (logical) If `TRUE` (default), plot the
                                   #'   `neighborhood` edges as a background layer.
                                   #' @param layout An `igraph` layout function (e.g., `igraph::layout_with_fr`).
                                   #' @param network_edges_col (character) Color for the `z_network` edges.
                                   #' @param neighborhood_edges_col (character) Color for the `neighborhood` edges.
                                   #' @param main (character) The main title for the plot.
                                   #' @param legend_col_n_levels (integer) Number of levels for the color legend.
                                   #' @param legend_size_n_levels (integer) Number of levels for the size legend.
                                   #' @param legend_pos (character) Position of the legend (e.g., `"right"`).
                                   #' @param alpha_neighborhood (numeric) Alpha transparency for neighborhood edges.
                                   #' @param edge.width (numeric) Width of the network edges.
                                   #' @param vertex.frame.width (numeric) Width of the vertex frame.
                                   #' @param edge.arrow.size (numeric) Size of the arrowheads for directed edges.
                                   #' @param coords (matrix) Optional matrix of x-y coordinates for node layout.
                                   #' @param legend_size (numeric) Scaling factor for the size legend.
                                   #' @param ... Additional arguments passed to `plot.igraph`.
                                   #' @return A list containing the `igraph` object (`graph`) and the
                                   #'   layout coordinates (`coords`), invisibly.
                                   plot = function(node_color = "x",
                                                   node_size = "y",
                                                   show_overlap = TRUE,
                                                   layout = igraph::layout_with_fr,
                                                   network_edges_col = "grey60",
                                                   neighborhood_edges_col = "orange",
                                                   main = "",
                                                   legend_col_n_levels = NULL,
                                                   legend_size_n_levels = NULL,
                                                   legend_pos = "right",
                                                   alpha_neighborhood = 0.2,
                                                   edge.width = 1, 
                                                   edge.arrow.size = 1, 
                                                   vertex.frame.width = 0.5,
                                                   coords = NULL, legend_size = 0.5,
                                                   ...) {
                                     # browser()
                                     if(is.null(legend_col_n_levels)){
                                       if(node_color == "x"){
                                         if(private$.type_x == "binomial"){
                                           legend_col_n_levels <- 2
                                         } else {
                                           legend_col_n_levels <- 4
                                         }
                                       } else {
                                         if(private$.type_y == "binomial"){
                                           legend_col_n_levels <- 2
                                         } else {
                                           legend_col_n_levels <- 4
                                         }
                                       }
                                     } 
                                     
                                     if(is.null(legend_size_n_levels)){
                                       if(node_size == "x"){
                                         if(private$.type_x == "binomial"){
                                           legend_size_n_levels <- 2
                                         } else {
                                           legend_size_n_levels <- 3
                                         }
                                       } else {
                                         if(private$.type_y == "binomial"){
                                           legend_size_n_levels <- 2
                                         } else {
                                           legend_size_n_levels <- 3
                                         }
                                       }
                                     } 
                                     
                                     # --- build igraph object from edge list
                                     
                                     if(!private$.directed){
                                       g <- igraph::graph_from_edgelist(as.matrix(private$.z_network[private$.z_network[,1] > private$.z_network[,2],]),
                                                                        directed = private$.directed)
                                     }else {
                                       g <- igraph::graph_from_edgelist(as.matrix(private$.z_network), directed = private$.directed)
                                     }
                                     n <- private$.n_actor
                                     if (igraph::vcount(g) < n) {
                                       g <- igraph::add_vertices(g, n - igraph::vcount(g))
                                     }
                                     # set canonical vertex names = 1:n so mapping is stable
                                     igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
                                     
                                     # --- map attributes to nodes
                                     xa <- private$.x_attribute
                                     ya <- private$.y_attribute
                                     
                                     # normalize attributes for visual mapping
                                     # norm <- function(v) (v - min(v, na.rm = TRUE)) / (max(v, na.rm = TRUE) - min(v, na.rm = TRUE) + 1e-9)
                                     xa_n <- find_ranks(xa)
                                     xa_n <- xa_n/max(xa_n)
                                     ya_n <- find_ranks(ya)
                                     ya_n <- ya_n/max(ya_n)
                                     
                                     
                                     # node color
                                     if (node_color == "x") {
                                       pal <-  grDevices::colorRampPalette(c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"))(100)
                                       node_cols <- pal[as.numeric(cut(xa_n, 100))]
                                       color_attr <- xa
                                       color_lab  <-  "x"
                                     } else if (node_color == "y") {
                                       pal <-  grDevices::colorRampPalette(c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"))(100)
                                       node_cols <- pal[as.numeric(cut(ya_n, 100))]
                                       color_attr <- ya
                                       color_lab  <- "y"
                                     } else {
                                       node_cols <- rep("grey70", n)
                                       color_attr <- NULL
                                       color_lab  <- NULL
                                     }
                                     
                                     if (node_size == "x") {
                                       node_sizes <- xa_n
                                       size_attr <- xa
                                       
                                       size_lab  <- "x"
                                     } else if (node_size == "y") {
                                       node_sizes <- ya_n
                                       size_attr <- ya
                                       
                                       size_lab  <-"y"
                                     } else {
                                       node_sizes <- rep(6, n)
                                       size_attr <- NULL
                                       size_lab  <- NULL
                                     }
                                     
                                     if (is.null(coords)) {
                                       coords <- layout(g)
                                     } else {
                                       # sanity: make sure dimensions match
                                       if (!is.matrix(coords) || nrow(coords) != igraph::vcount(g) || ncol(coords) < 2) {
                                         stop("`coords` must be a matrix with nrow = vcount(g) and at least 2 columns.")
                                       }
                                     }
                                     if (show_overlap && !is.null(private$.neighborhood) && nrow(private$.neighborhood) > 0) {
                                       overlap_edges <- as.matrix(private$.neighborhood)
                                       overlap_edges <- overlap_edges[overlap_edges[, 1] != overlap_edges[, 2], , drop = FALSE]
                                       g2 <- igraph::graph_from_edgelist(overlap_edges, directed = FALSE)
                                       igraph::V(g2)$name <- as.character(seq_len(igraph::vcount(g2)))
                                       igraph::plot.igraph(
                                         g2,
                                         edge.color = grDevices::adjustcolor(neighborhood_edges_col, alpha.f = alpha_neighborhood ),
                                         vertex.color = "white",
                                         edge.width = 1, 
                                         vertex.size = 0,
                                         vertex.label = NA,
                                         layout = coords
                                       )
                                       add_mode <- TRUE
                                     } else {
                                       add_mode <- FALSE
                                     }
                                     # browser()
                                     plot(
                                       g,
                                       vertex.frame.width =vertex.frame.width,
                                       vertex.color = node_cols,
                                       vertex.size  = log(node_sizes+1)*10+2,
                                       vertex.label = NA,
                                       edge.color   = network_edges_col,
                                       edge.width = edge.width, 
                                       edge.arrow.size = edge.arrow.size, 
                                       add = add_mode, 
                                       layout = coords, 
                                     )
                                     
                                     if (!is.null(color_attr) || !is.null(size_attr)) {
                                       
                                       if (!is.null(color_attr)) {
                                         # browser()
                                         rng <- quantile(color_attr, na.rm = TRUE, probs = seq(from = 0, to = 1, length.out = legend_col_n_levels))
                                         color_vals <- round(rng, 2)
                                         color_cols <- c(pal[1],pal[round(seq(from = 0, to = 1, length.out = legend_col_n_levels)*100)])
                                       } else{
                                         color_cols <- NULL
                                         color_vals <- NULL
                                       }
                                       if (!is.null(size_attr)) {
                                         
                                         vals <- sort(unique(size_attr))[findInterval(vec = sort(unique(node_sizes)),
                                                                              x = seq(from = 0, to = 0.99, length.out = legend_size_n_levels))+1]
                                         sizes <- sort(unique(node_sizes))[findInterval(vec = sort(unique(node_sizes)),
                                                                                      x = seq(from = 0, to = 0.99, length.out = legend_size_n_levels))+1]
                                         sizes <- log(sizes+1)*10+2
                                       } else {
                                         sizes <- NULL
                                         vals <- NULL
                                       }
                                       length_a <- length(color_vals)
                                       length_b <- length(vals)  
                                       if(length_a > 0 ){
                                         labels_a <- paste(color_lab, "=", color_vals)
                                       } else {
                                         labels_a <- c()
                                       }
                                       if(length_b > 0 ){
                                         labels_b <- paste(size_lab, "=", round(vals, 2))
                                         
                                         color_cols_b <- rep("grey70", times = length_b)
                                       } else {
                                         labels_b <- c()
                                         color_cols_b <- c()
                                       }
                                       
                                       labels <- c(labels_a, labels_b)
                                       colors_tmp <- c(color_cols, color_cols_b)
                                       sizes <- c(rep(1, times = length_a), sizes/8)
                                       legend(title = "Legend", 
                                              legend_pos,
                                              legend = labels,
                                              pt.bg   = colors_tmp,
                                              border = "black",
                                              pch = 21,
                                              pt.cex = sizes,
                                              cex    = legend_size )
                                     }
                                     
                                     invisible(list(graph = g, coords = coords))
                                   },
                                   #' @description
                                   #' Print a summary of the `iglm.data` object to the console.
                                   #'
                                   #' @param digits (integer) Number of digits to round numeric output to.
                                   #' @param ... Additional arguments (not used).
                                   #' @return The object's private environment, invisibly.
                                   print = function(digits = 3, ...) {
                                     n <- as.integer(private$.n_actor)
                                     dir_flag <- isTRUE(private$.directed)
                                     
                                     # --- helper to coerce stored edge list to a clean 2-col integer matrix
                                     edge_mat <- function(E) {
                                       if (is.null(E)) return(matrix(integer(0), ncol = 2))
                                       E <- as.matrix(E)
                                       if (!is.numeric(E) || ncol(E) != 2) {
                                         stop("Edge list must be numeric matrix with 2 columns.")
                                       }
                                       storage.mode(E) <- "integer"
                                       # drop self-loops for counting/printing
                                       E <- E[ E[,1] != E[,2], , drop = FALSE ]
                                       E
                                     }
                                     
                                     zE   <- edge_mat(private$.z_network)
                                     nbE  <- edge_mat(private$.neighborhood)
                                     ovE  <- edge_mat(private$.overlap)
                                     
                                     # --- edge counting: for undirected, count unique unordered pairs
                                     count_edges <- function(E, directed) {
                                       if (nrow(E) == 0) return(0L)
                                       if (directed) {
                                         return(nrow(E))
                                       } else {
                                         U <- cbind(pmin(E[,1], E[,2]), pmax(E[,1], E[,2]))
                                         return(nrow(unique(U)))
                                       }
                                     }
                                     
                                     m_z  <- count_edges(zE,  dir_flag)
                                     m_nb <- count_edges(nbE, dir_flag)  # "neighborhood" may be a relation; still count pairs
                                     
                                     # --- density computed on z_network only (primary network)
                                     denom <- if (dir_flag) n * (n - 1) else n * (n - 1) / 2
                                     dens  <- if (denom > 0) m_z / denom else NA_real_
                                     
                                     # --- quick attribute summaries based on declared types
                                     numfmt <- function(v) format(round(v, digits), nsmall = digits, trim = TRUE)
                                     
                                     summarize_attr <- function(v, type) {
                                       v <- as.vector(v)
                                       if (type == "binomial") {
                                         n1 <- sum(v == 1, na.rm = TRUE)
                                         n0 <- sum(v == 0, na.rm = TRUE)
                                         p1 <- mean(v == 1, na.rm = TRUE)
                                         paste0("(binomial): 1s=", n1, ", 0s=", n0, ", P(1)=", numfmt(p1))
                                       } else if (type == "poisson") {
                                         qs <- stats::quantile(v, c(.00,.25,.50,.75,1), na.rm = TRUE, names = FALSE)
                                         paste0(
                                           "(poisson): min=", numfmt(qs[1]),
                                           ", q1=", numfmt(qs[2]),
                                           ", med=", numfmt(qs[3]),
                                           ", q3=", numfmt(qs[4]),
                                           ", max=", numfmt(qs[5]),
                                           ", mean=", numfmt(mean(v, na.rm = TRUE))
                                         )
                                       } else if (type == "normal") {
                                         paste0(
                                           "(normal): mean=", numfmt(mean(v, na.rm = TRUE)),
                                           ", sd=", numfmt(stats::sd(v, na.rm = TRUE)),
                                           ", min=", numfmt(min(v, na.rm = TRUE)),
                                           ", max=", numfmt(max(v, na.rm = TRUE))
                                         )
                                       } else {
                                         paste0("unknown type; length=", length(v))
                                       }
                                     }
                                     
                                     x_sum <- summarize_attr(private$.x_attribute, private$.type_x)
                                     y_sum <- summarize_attr(private$.y_attribute, private$.type_y)
                                     
                                     
                                     cat("iglm.data object\n")
                                     cat("  vertices   :", n, "\n")
                                     cat("  directed   :", if (dir_flag) "TRUE" else "FALSE", "\n")
                                     cat("  scales     : scale_x=", numfmt(private$.scale_x),
                                         ", scale_y=", numfmt(private$.scale_y), "\n", sep = "")
                                     cat("  edges (z)  :", m_z, "\n")
                                     cat("  density (z):", if (is.na(dens)) "NA" else numfmt(dens), "\n")
                                     cat("  neighborhood edges :", m_nb, "\n")
                                     cat("\n")
                                     cat("Attribute summaries\n")
                                     cat("  x_attribute ", x_sum, "\n", sep = "")
                                     cat("  y_attribute ", y_sum, "\n", sep = "")
                                     cat("\n")
                                     invisible(private)
                                   }
                                 ),
                                 
                                 # --- Active Bindings for Read-Only Access ---
                                 active = list(
                                   #' @field x_attribute (`numeric`) Read-only. The vector for the first unit-level attribute.
                                   x_attribute = function(value) { if(missing(value)) private$.x_attribute else stop("`x_attribute` is read-only.", call. = FALSE) },
                                   
                                   #' @field y_attribute (`numeric`) Read-only. The vector for the second unit-level attribute.
                                   y_attribute = function(value) { if(missing(value)) private$.y_attribute else stop("`y_attribute` is read-only.", call. = FALSE) },
                                   
                                   #' @field z_network (`matrix`) Read-only. The primary network structure as a 2-column integer edgelist.
                                   z_network = function(value) { if(missing(value)) private$.z_network else stop("`z_network` is read-only.", call. = FALSE) },
                                   
                                   #' @field neighborhood (`matrix` or `NULL`) Read-only. The secondary/neighborhood structure as a 2-column integer edgelist. `NULL` if not provided.
                                   neighborhood = function(value) { if(missing(value)) private$.neighborhood else stop("`neighborhood` is read-only.", call. = FALSE) },
                                   
                                   #' @field overlap (`matrix`) Read-only. The calculated overlap relation (dyads with shared neighbors in `neighborhood`) as a 2-column integer edgelist.
                                   overlap = function(value) { if(missing(value)) private$.overlap else stop("`overlap` is read-only.", call. = FALSE) },
                                   
                                   #' @field directed (`logical`) Read-only. Indicates if the `z_network` is treated as directed.
                                   directed = function(value) { if(missing(value)) private$.directed else stop("`directed` is read-only.", call. = FALSE) },
                                   
                                   #' @field n_actor (`integer`) Read-only. The total number of actors (nodes) in the network.
                                   n_actor = function(value) { if(missing(value)) private$.n_actor else stop("`n_actor` is read-only.", call. = FALSE) },
                                   
                                   #' @field type_x (`character`) Read-only. The specified distribution type for the `x_attribute`.
                                   type_x = function(value) { if(missing(value)) private$.type_x else stop("`type_x` is read-only.", call. = FALSE) },
                                   
                                   #' @field type_y (`character`) Read-only. The specified distribution type for the `y_attribute`.
                                   type_y = function(value) { if(missing(value)) private$.type_y else stop("`type_y` is read-only.", call. = FALSE) },
                                   
                                   #' @field scale_x (`numeric`) Read-only. The scale parameter associated with the `x_attribute`.
                                   scale_x = function(value) { if(missing(value)) private$.scale_x else stop("`scale_x` is read-only.", call. = FALSE) },
                                   
                                   #' @field scale_y (`numeric`) Read-only. The scale parameter associated with the `y_attribute`.
                                   scale_y = function(value) { if(missing(value)) private$.scale_y else stop("`scale_y` is read-only.", call. = FALSE) }
                                 )
)

#' Constructor for the iglm.data R6 object
#'
#' @description
#' Creates a `iglm.data` object, which stores network and attribute data.
#' This function acts as a user-friendly interface to the `iglm.data` R6 class generator.
#' It handles data input, infers parameters like the number of actors (`n_actor`)
#' and network directedness (`directed`) if not explicitly provided, processes
#' network data into a consistent edgelist format, calculates the overlap
#' relation based on an optional neighborhood definition, and performs
#' extensive validation of all inputs.
#'
#' @param x_attribute A numeric vector for the first unit-level attribute.
#' @param y_attribute A numeric vector for the second unit-level attribute.
#' @param z_network A matrix representing the network. Can be a 2-column
#'   edgelist or a square adjacency matrix.
#' @param neighborhood An optional matrix for the neighborhood representing local dependence.
#'   Can be a 2-column edgelist or a square adjacency matrix. 
#'   A tie in `neighborhood` between actor i and j indicates that j is in the neighborhood of i, 
#'   implying dependence between the respective actors. 
#' @param directed A logical value indicating if `z_network` is directed.
#'   If `NA` (default), directedness is inferred from the symmetry of
#'   `z_network`.
#' @param n_actor An integer for the number of actors in the system.
#'   If `NA` (default), `n_actor` is inferred from the attributes or
#'   network matrices.
#' @param type_x Character string for the type of `x_attribute`.
#'   Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
#'   Default is `"binomial"`.
#' @param type_y Character string for the type of `y_attribute`.
#'   Must be one of `"binomial"`, `"poisson"`, or `"normal"`.
#'   Default is `"binomial"`.
#' @param scale_x A positive numeric value for scaling (e.g., variance
#'   for "normal" type). Default is 1.
#' @param scale_y A positive numeric value for scaling (e.g., variance
#'   for "normal" type). Default is 1.
#' @param return_neighborhood Logical. If `TRUE` (default) and
#'   `neighborhood` is `NULL`, a full neighborhood (all dyads) is
#'   generated implying global dependence. If `FALSE`, no neighborhood is set.
#' @param file (character) Optional file path to load a saved `iglm.data` object state.   
#' @return An object of class `iglm.data` (and `R6`).
#' @examples
#' data(state_twitter)
#' state_twitter$iglm.data$degree_distribution(prob = FALSE, plot = TRUE)
#' state_twitter$iglm.data$geodesic_distances_distribution(prob = FALSE, plot = TRUE)
#' state_twitter$iglm.data$density_x()
#' state_twitter$iglm.data$density_y()
#' 
#' # Generate a small iglm data object either via adjacency matrix or edgelist
#' tmp_adjacency <- iglm.data(z_network = matrix(c(0,1,1,0,
#'                                                 1,0,0,1,
#'                                                 1,0,0,1,
#'                                                 0,1,1,0), nrow=4, byrow=TRUE),
#'                            directed = FALSE,
#'                            n_actor = 4,
#'                            type_x = "binomial",
#'                            type_y = "binomial")
#' 
#' 
#' tmp_edgelist <- iglm.data(z_network = tmp_adjacency$z_network, 
#'                           directed = FALSE,
#'                        n_actor = 4,
#'                        type_x = "binomial",
#'                        type_y = "binomial")
#' 
#' tmp_edgelist$mean_z()
#' tmp_adjacency$mean_z()
#' @export
iglm.data <- function(x_attribute = NULL, y_attribute = NULL, z_network = NULL,
                    neighborhood = NULL, directed = TRUE, n_actor = NA, 
                    type_x = "binomial", type_y = "binomial",
                    scale_x = 1, scale_y = 1, 
                    return_neighborhood = TRUE, file = NULL) {
  
  iglm.data_generator$new(x_attribute = x_attribute,
                        y_attribute = y_attribute,
                        z_network = z_network,
                        neighborhood = neighborhood,
                        directed = directed,
                        n_actor = n_actor,
                        type_x = type_x,
                        type_y = type_y,
                        scale_x = scale_x,
                        scale_y = scale_y, 
                        return_neighborhood = return_neighborhood, 
                        file = file)
}

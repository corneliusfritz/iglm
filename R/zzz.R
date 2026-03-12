#' @useDynLib iglm, .registration=TRUE

# Rcpp::loadModule("double_cpp", TRUE)
.onUnload <- function(libpath) {
  library.dynam.unload("iglm", libpath)
}
# Small helper functions
rowMins <- function(x) {
  apply(x, 1, min)
}
rowMaxs <- function(x) {
  apply(x, 1, max)
}
colMins <- function(x) {
  apply(x, 2, min)
}
colMaxs <- function(x) {
  apply(x, 2, max)
}
add_alpha <- function(color_code, alpha_level) {
  if (alpha_level < 0 || alpha_level > 1) {
    stop("Alpha level must be between 0 and 1.")
  }
  alpha_int <- round(alpha_level * 255)
  rgb_matrix <- col2rgb(color_code)
  new_color <- rgb(
    red = rgb_matrix[1, ],
    green = rgb_matrix[2, ],
    blue = rgb_matrix[3, ],
    alpha = alpha_int,
    maxColorValue = 255
  )

  return(new_color)
}

#' @title Model specification for a `iglm' object
#'
#' @description
#'
#' The help pages of \code{\link{iglm}} describe the model with details on model fitting
#' and estimation.
#' Generally, a model is specified via it's sufficient statistics,
#' that can be further decomposed into two parts:
#' \itemize{
#'   \item \strong{\eqn{\mathbf{g}_i(x_i^*,y_i^*) = \mathbf{g}_i(x_i,y_i)= (g_i(x_i,y_i))}}: A vector of unit-level functions (or "g-terms")
#'     that describe the relationship between an individual actor \eqn{i}'s
#'     predictors (\eqn{x_i}) and their own response (\eqn{y_i}).
#'   \item \strong{\eqn{\mathbf{h}_{i,j}(x_i^*,x_j^*, y_i^*, y_j^*, z) = \mathbf{h}_{i,j}(x,y,z)= (h_{i,j}(x,y,z))}}: A vector of pair-level functions (or "h-terms")
#'     that specify how the connections (\eqn{z}) and responses (\eqn{y_i, y_j})
#'     of a pair of units \eqn{\{i,j\}} depend on each other and the wider
#'     network structure.
#' }
#' Each term defines a component for the model's features, which
#' are a sum of unit-level components, \eqn{\sum_i g_i(x_i,y_i)}, and/or
#' pair-level components, \eqn{\sum_{i \ne j} h_{i,j}(x,y,z)}.
#' The implemented terms are grouped into three categories:
#' \enumerate{
#' \item  \eqn{{g}_i} terms for attribute dependence,
#' \item \eqn{{h}_{i,j}} terms for network dependence,
#' \item \eqn{{h}_{i,j}} Terms for joint attribute/network dependence.
#' }
#' Below is a detailed description of each term that can be specified in the model formula (see, \code{\link{iglm}}):
#' \code{iglm.data ~ <term_1> + <term_2> + ... }, where the left-hand side of the formula has to be a \code{\link{iglm.data}} object,
#' and the right-hand side lists some of the terms from the list below, which should be included in the model.
#' Setting, e.g., \code{<term_1>} to  \code{attribute_x} includes the term \eqn{g_i(x,y,z) = x_i} in \eqn{g_i}.
#' Note that the function \code{\link{create_userterms_skeleton}} can be used to create custom terms.
#'
#' @section Notation:
#'
#' Here, \eqn{x_i} and \eqn{y_i} are the attributes for actor \eqn{i},
#' and \eqn{z_{i,j}} indicates the presence (1) or absence (0) of a tie
#' from actor \eqn{i} to actor \eqn{j}.
#' The local neighborhood of actor \eqn{i} is denoted \eqn{\mathcal{N}_i},
#' and the indicator for whether actors \eqn{i} and \eqn{j} share a local
#' neighborhood is given by
#' \eqn{c_{i,j} = \mathbb{I}(\mathcal{N}_i \cap \mathcal{N}_j \neq \emptyset)}.
#' The functions below specify the forms of \eqn{g_i(x_i,y_i)} and \eqn{h_{i,j}(x,y,z)}
#' for each term.
#' Some terms also depend on other covariates, which are denoted by
#' \eqn{v = (v_1, ..., v_N)} (unit-level) and \eqn{w = (w_{i,j}) \in \mathbb{R}^{N \times N}}(dyadic).
#' These covariates must be provided by the user via the \code{data} argument of the terms.
#' Assuming that the matriv \code{x} exists in the environment associated with the
#' used formula, \code{cov_z(data = v)} includes the dyadic covariable \eqn{v = (v_{i,j})} in the model.
#' Some terms also have a \code{mode} argument, which can take values
#' \code{"global"}, \code{"local"}, or \code{"alocal"}.
#' The \code{"global"} mode indicates that the statistic is computed
#' over the entire network, while \code{"local"} restricts the statistic
#' to local neighborhoods only (i.e., edges where \eqn{c_{i,j} = 1}).
#' The \code{"alocal"} mode restricts the statistic to non-local edges
#' only (i.e., edges where \eqn{c_{i,j} = 0}).
#' For instance, \code{edges(mode = "local")} counts the number of edges
#' that connect actors with overlapping neighborhoods.
#' See underneath for which options are implemented for each term.
#' See the documentation for \code{\link{iglm}} for details on model fitting
#' and estimation.
#'
#' @section List of Terms:
#'
#' \describe{
#'    \item{\strong{1. \eqn{g_i} Terms for Attribute Dependence}}{}
#'    \item{\code{attribute_x}}{
#'      \strong{Attribute (X) [g-term]:} Intercept for attribute \eqn{x}.
#'      \eqn{g_i(x_i,y_i) = x_i}
#'    }
#'    \item{\code{attribute_y}}{
#'      \strong{Attribute (Y) [g-term]:} Intercept for attribute \eqn{y}.
#'      \eqn{g_i(x_i,y_i) = y_i}
#'    }
#'    \item{\code{cov_x}}{
#'      \strong{Nodal Covariate (X) [g-term]:} Effect of a unit-level exogenous covariate \eqn{v_i} on endogenous attribute \eqn{x_i}.
#'      \eqn{g_i(x_i,y_i) = v_i x_i}
#'    }
#'    \item{\code{cov_y(data = v)}}{
#'      \strong{Nodal Covariate (Y) [g-term]:} Effect of a unit-level exogenous covariate \eqn{v_i} on endogenous attribute \eqn{y_i}.
#'      \eqn{g_i(x_i,y_i) = v_i y_i}
#'    }
#'    \item{\code{attribute_xy(mode = "global")}}{
#'      \strong{Nodal Attribute Interaction (X-Y) [g-term]:} Interaction of attributes \eqn{x_i} and \eqn{y_i}. Spatial or structural boundaries determine the interaction scale.
#'      \itemize{
#'        \item \code{global}: \eqn{g_i(x_i,y_i) = x_i y_i}
#'        \item \code{local}: \eqn{g_i(x_i,y_i) = x_i \sum_{j \in \mathcal{N}_i} y_j + y_i \sum_{j \in \mathcal{N}_i} x_j}
#'        \item \code{alocal}: \eqn{g_i(x_i,y_i) = x_i \sum_{j \notin \mathcal{N}_i} y_j + y_i \sum_{j \notin \mathcal{N}_i} x_j}
#'      }
#'    }
#'    \item{\strong{2. \eqn{h_{i,j}} Terms for Network Dependence}}{}
#'    \item{\code{degrees}}{
#'      \strong{Degrees [h-term]:} Specifies node-level fixed effects. Estimation requires an MM algorithm constraint.
#'      Directed networks assign distinct sender and receiver effects (constraining node N's out-effect to 0 for identifiability).
#'      Undirected networks assign a singular degree effect per node.
#'    }
#'    \item{\code{edges(mode = "global")}}{
#'      \strong{Edges [h-term]:} Captures the baseline propensity of tie formation \eqn{z_{i,j}}, partitioned by structural boundary \eqn{c_{i,j}}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) z_{i,j}}
#'      }
#'    }
#'    \item{\code{mutual(mode = "global")}}{
#'      \strong{Mutual Reciprocity [h-term]:} Evaluates reciprocal tie formation in directed networks.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = z_{i,j} z_{j,i}} (for \eqn{i < j})
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} z_{i,j} z_{j,i}} (for \eqn{i < j})
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) z_{i,j} z_{j,i}} (for \eqn{i < j})
#'      }
#'    }
#'    \item{\code{cov_z(data, mode = "global")}}{
#'      \strong{Dyadic Covariate [h-term]:} Exogenous dyadic covariate \eqn{w_{i,j}} influence on edge formation.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = w_{i,j} z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} w_{i,j} z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) w_{i,j} z_{i,j}}
#'      }
#'    }
#'    \item{\code{isolates}}{
#'      \strong{Isolates [z-term]:} Captures the frequency of nodes with degree zero \eqn{\mathbb{I} (\sum_{j\neq i} z_{i,j} = 0)}.
#'    }
#'    \item{\code{nonisolates}}{
#'      \strong{Non-Isolates [z-term]:} Captures the frequency of nodes with degree strictly greater than zero \eqn{\mathbb{I} (\sum_{j\neq i} z_{i,j} \neq 0)}.
#'    }
#'    \item{\code{gwodegree(mode = "global", decay = 2.0)}}{
#'      \strong{Geometrically Weighted Out-Degree [z-term]:} Captures the out-degree distribution utilizing an exponential decay parameter \eqn{\alpha}.
#'      \itemize{
#'        \item \code{global}: Evaluated across the entire adjacency matrix.
#'        \item \code{local}: Evaluated strictly within local neighborhood partitions \eqn{c_{i,j}}.
#'      }
#'    }
#'    \item{\code{gwidegree(mode = "global", decay = 2.0)}}{
#'      \strong{Geometrically Weighted In-Degree [z-term]:} Captures the in-degree distribution utilizing an exponential decay parameter \eqn{\alpha}.
#'      \itemize{
#'        \item \code{global}: Evaluated across the entire adjacency matrix.
#'        \item \code{local}: Evaluated strictly within local neighborhood partitions \eqn{c_{i,j}}.
#'      }
#'    }
#'    \item{\code{gwesp(data, mode = "global", type = "OTP", decay = 2.0)}}{
#'      \strong{Geometrically Weighted Edgewise-Shared Partners [h-term]:} Models triadic closure propensity conditioning on existing edges, scaled by decay parameter \eqn{\alpha}.
#'      Types dictate the structural path constraint: OTP (outgoing two-paths, \eqn{z_{i,h}\,z_{h,j}\, z_{i,j}}), ITP (incoming two-paths, \eqn{z_{h,i}\,z_{j,h}\, z_{i,j}}),
#'      OSP (outgoing shared partners, \eqn{z_{i,h}\,z_{j,h}\, z_{i,j}}), ISP (incoming shared partners, \eqn{z_{h,i}\,z_{h,j}\, z_{i,j}}).
#'      \itemize{
#'        \item \code{global}: Path evaluations span the entire network topology.
#'        \item \code{local}: Path evaluations are constrained exclusively to \eqn{c_{i,j}} overlapping neighborhoods.
#'      }
#'    }
#'    \item{\code{gwdsp(data, mode = "global", type = "OTP", decay = 2.0)}}{
#'      \strong{Geometrically Weighted Dyadwise-Shared Partners [h-term]:} Models the baseline triadic potential irrespective of the closing \eqn{z_{i,j}} edge, scaled by decay parameter \eqn{\alpha}.
#'      Types dictate the structural path constraint: OTP (outgoing two-paths, \eqn{z_{i,h}\,z_{h,j}}), ITP (incoming two-paths, \eqn{z_{h,i}\,z_{j,h}}),
#'      OSP (outgoing shared partners, \eqn{z_{i,h}\,z_{j,h}}), ISP (incoming shared partners, \eqn{z_{h,i}\,z_{h,j}}).
#'      \itemize{
#'        \item \code{global}: Path evaluations span the entire network topology.
#'        \item \code{local}: Path evaluations are constrained exclusively to \eqn{c_{i,j}} overlapping neighborhoods.
#'      }
#'    }
#'    \item{\code{cov_z_out(data, mode = "global")}}{
#'      \strong{Covariate Sender [h-term]:} Exogenous monadic covariate \eqn{v_{i}} influence on generating an outgoing tie.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = v_i z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} v_i z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) v_i z_{i,j}}
#'      }
#'    }
#'    \item{\code{cov_z_in(data, mode = "global")}}{
#'      \strong{Covariate Receiver [h-term]:} Exogenous monadic covariate \eqn{v_{j}} influence on receiving an incoming tie.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = v_j z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} v_j z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) v_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{transitive}}{
#'      \strong{Transitivity (Local) [Joint]:} Indicator evaluating the presence of a local transitive triad configuration. An actor \eqn{h \neq i,j} must exist fulfilling
#'      \eqn{h\in \mathcal{N}_i, h\in \mathcal{N}_j} mapping to \eqn{z_{i,j} = z_{i,h} = z_{h,j}}: \eqn{h_{i,j} = c_{i,j} z_{i,j} \mathbb{I}(\sum_{k} c_{i,k} c_{j,k} z_{i,k} z_{k,j}>1)}.
#'    }
#'    \item{\strong{3. \eqn{h_{i,j}} Terms for Joint Attribute/Network Dependence}}{}
#'    \item{\code{attribute_xz(mode = "local")}}{
#'      \strong{Attribute Sum (X-Z) [Joint h-term]:} Models the additive effect of \eqn{x_i} and \eqn{x_j} on edge formation within local neighborhoods.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} (x_i + x_j) z_{i,j}}
#'      }
#'    }
#'    \item{\code{attribute_yz(mode = "local")}}{
#'      \strong{Attribute Sum (Y-Z) [Joint h-term]:} Models the additive effect of \eqn{y_i} and \eqn{y_j} on edge formation within local neighborhoods.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} (y_i + y_j) z_{i,j}}
#'      }
#'    }
#'    \item{\code{edges_x_match(mode = "global")}}{
#'      \strong{Attribute Match (X-Z) [Joint h-term]:} Models homophily/matching on the binary attribute \eqn{x}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = \mathbb{I}(x_i = x_j) z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} \mathbb{I}(x_i = x_j) z_{i,j}}
#'      }
#'    }
#'    \item{\code{edges_y_match(mode = "global")}}{
#'      \strong{Attribute Match (Y-Z) [Joint h-term]:} Models homophily/matching on the binary attribute \eqn{y}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = \mathbb{I}(y_i = y_j) z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} \mathbb{I}(y_i = y_j) z_{i,j}}
#'      }
#'    }
#'    \item{\code{outedges_x(mode = "global")}}{
#'      \strong{Attribute Out-Degree (X-Z) [Joint h-term]:} Influence of endogenous \eqn{x_i} on out-degree formation.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = x_i z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_i z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) x_i z_{i,j}}
#'      }
#'    }
#'    \item{\code{inedges_x(mode = "global")}}{
#'      \strong{Attribute In-Degree (X-Z) [Joint h-term]:} Influence of endogenous \eqn{x_j} on in-degree reception.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = x_j z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_j z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) x_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{outedges_y(mode = "global")}}{
#'      \strong{Attribute Out-Degree (Y-Z) [Joint h-term]:} Influence of endogenous \eqn{y_i} on out-degree formation.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = y_i z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) y_i z_{i,j}}
#'      }
#'    }
#'    \item{\code{inedges_y(mode = "global")}}{
#'      \strong{Attribute In-Degree (Y-Z) [Joint h-term]:} Influence of endogenous \eqn{y_j} on in-degree reception.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = y_j z_{i,j}}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_j z_{i,j}}
#'        \item \code{alocal}: \eqn{h_{i,j}(x,y,z) = (1 - c_{i,j}) y_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_xx(mode = "local")}}{
#'      \strong{Symmetric X-X-Z Outcome Spillover [h-term]:} Propagates \eqn{x}-outcome spillover effects through explicit network ties.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_i x_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_xx_scaled(mode = "local")}}{
#'      \strong{Scaled X-X-Z Outcome Spillover [h-term]:} Normalizes the \eqn{x}-outcome spillover influence by the relevant out-degree topology of actor \eqn{i}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = x_i x_j z_{i,j} / \text{deg}(i)}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_i x_j z_{i,j} / \text{deg}(i, \text{local})}
#'      }
#'    }
#'    \item{\code{spillover_yy(mode = "local")}}{
#'      \strong{Symmetric Y-Y-Z Outcome Spillover [h-term]:} Propagates \eqn{y}-outcome spillover effects through explicit network ties.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i y_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_yy_scaled(mode = "local")}}{
#'      \strong{Scaled Y-Y-Z Outcome Spillover [h-term]:} Normalizes the \eqn{y}-outcome spillover influence by the relevant out-degree topology of actor \eqn{i}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = y_i y_j z_{i,j} / \text{deg}(i)}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i y_j z_{i,j} / \text{deg}(i, \text{local})}
#'      }
#'    }
#'    \item{\code{spillover_xy(mode = "local")}}{
#'      \strong{Directed X-Y-Z Treatment Spillover [h-term]:} Maps the cross-attribute \eqn{x_i \to y_j} treatment assignment across the network graph.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_i y_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_xy_scaled(mode = "local")}}{
#'      \strong{Scaled X-Y-Z Treatment Spillover [h-term]:} Normalizes the cross-attribute \eqn{x_i \to y_j} spillover influence by the relevant out-degree topology of actor \eqn{i}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = x_i y_j z_{i,j} / \text{deg}(i)}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} x_i y_j z_{i,j} / \text{deg}(i, \text{local})}
#'      }
#'    }
#'    \item{\code{spillover_yx(mode = "local")}}{
#'      \strong{Directed Y-X-Z Treatment Spillover [h-term]:} Maps the cross-attribute \eqn{y_i \to x_j} treatment assignment across the network graph.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i x_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_yx_scaled(mode = "local")}}{
#'      \strong{Scaled Y-X-Z Treatment Spillover [h-term]:} Normalizes the cross-attribute \eqn{y_i \to x_j} spillover influence by the relevant out-degree topology of actor \eqn{i}.
#'      \itemize{
#'        \item \code{global}: \eqn{h_{i,j}(x,y,z) = y_i x_j z_{i,j} / \text{deg}(i)}
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i x_j z_{i,j} / \text{deg}(i, \text{local})}
#'      }
#'    }
#'    \item{\code{spillover_yc(mode = "local")}}{
#'      \strong{Directed Y-C-Z Treatment Spillover [h-term]:} Specifies exogenous covariate \eqn{v} interacting with endogenous trait \eqn{y} over the network.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} y_i v_j z_{i,j}}
#'      }
#'    }
#'    \item{\code{spillover_yc_symm(data = v, mode = "local")}}{
#'      \strong{Symmetric Y-C-Z Treatment Spillover [h-term]:} Bidirectional mapping of the exogenous covariate \eqn{v} and endogenous trait \eqn{y} interaction.
#'      \itemize{
#'        \item \code{local}: \eqn{h_{i,j}(x,y,z) = c_{i,j} (v_i y_j + v_j y_i) z_{i,j}}
#'      }
#'    }
#' }
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
#' @aliases terms
#' @name model.terms
NULL

#' Generate the Skeleton for an R package to implement additional iglm terms
#'
#' @description
#' This function generates the directory structure and source files for a new R package
#' named \code{iglm.userterms} (or whatever name is provided in the parameter \code{pkg_name}).
#' This auxiliary package serves as a template for extending the
#' \code{iglm} framework to user-defined sufficient statistics.
#' By compiling this package, users can link custom C++ implementations of change statistics
#' directly with the \code{iglm} package, enabling seamless integration of new model terms.
#'
#' @param path A character string specifying the path where the package directory
#'   should be created. Defaults to the current working directory (\code{"."}).
#' @param pkg_name A character string specifying the name of the package to be created.
#'
#' @details
#' The function creates a directory with the name specified in \code{pkg_name}
#' at the specified location.
#' As an example for a possible statistic, the statistic counting mutual
#' connections in the network is implemented.
#' After defining all possible change-statistics in the c++ function (this has to include a change for
#' \code{z_ij} (network), \code{x_i} (attribute x), and \code{y_i} (attribute y) all toggling from 0 to 1),
#' the function has to be registered using the \code{EFFECT_REGISTER} macro.
#' After compiling the package,
#' users have to load the package using \code{library(pkg_name)} before using it in \code{iglm}.
#'
#' @export
create_userterms_skeleton <- function(path = ".", pkg_name = "iglm.userterms") {
  pkg_path <- file.path(path, pkg_name)

  if (dir.exists(pkg_path)) {
    stop(paste(
      "Directory", pkg_path,
      "already exists. Please remove it or choose a different location."
    ))
  }

  # 2. Create Directory Structure
  dir.create(pkg_path, recursive = TRUE)
  dir.create(file.path(pkg_path, "R"))
  dir.create(file.path(pkg_path, "src"))

  # 3. Define File Contents based on your upload

  # --- DESCRIPTION  ---
  desc_content <- c(
    paste("Package: ", pkg_name),
    "Type: Package",
    "Title: Userterms for Regression under Network Interference",
    "Version: 1.0",
    "Date: 2025-11-09",
    "Authors@R: c(person(given = \"Your\", family = \"Name\", role = c(\"aut\", \"cre\"), email = \"corneliusfritz2010@gmail.com\"))",
    "Description: Userterms for generalized linear models (GLMs) for studying relationships among attributes in connected populations.",
    "License: GPL-3",
    "Encoding: UTF-8",
    "Imports: iglm, Rcpp",
    "LinkingTo: Rcpp, RcppArmadillo, iglm",
    "NeedsCompilation: yes",
    "RoxygenNote: 7.3.3"
  )

  # --- NAMESPACE ---
  ns_content <- c(
    "# Generated by roxygen2: do not edit by hand",
    "",
    "import(iglm)",
    "importFrom(Rcpp,evalCpp)",
    paste0("useDynLib(", pkg_name, ", .registration=TRUE)")
  )

  zzz_content <- c(
    paste("#' @useDynLib", pkg_name, ", .registration=TRUE"),
    "#' @import iglm",
    "#' @importFrom Rcpp evalCpp",
    "NULL"
  )


  # --- src/Makevars  ---
  makevars_content <- c(
    "IGLM_LIB_PATH = $(shell \"${R_HOME}/bin/Rscript\" -e \"cat(file.path(system.file(package='iglm'), 'libs'))\")",
    "IGLM_SO_PATH = $(shell \"${R_HOME}/bin/Rscript\" -e \"cat(file.path(system.file(package = 'iglm'), 'libs',paste0('iglm.so')))\")",
    "PKG_LIBS = $(IGLM_SO_PATH) -Wl,-rpath,$(IGLM_LIB_PATH) -lR $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"
  )

  # --- src/additional_cc.cpp ---
  cpp_content <- c(
    "#include <RcppArmadillo.h>",
    "#include \"iglm/extension_api.hpp\"",
    "#include \"iglm/xyz_class.h\"",
    "",
    "double xyz_stat_my_mutual(const XYZ_class &object,",
    "                          const int &actor_i,",
    "                          const int &actor_j,",
    "                          const arma::mat &data,",
    "                          const double &type,",
    "                          const std::string &mode,",
    "                          const bool &is_full_neighborhood){",
    "  // Define change statistic if Z_i,j: 0 to 1",
    "  if(mode == \"z\"){",
    "    return(object.z_network.get_val(actor_j, actor_i));",
    "  } else if(mode == \"x\"){",
    "    // Define change statistic if X_i: 0 to 1",
    "    return(0);",
    "  }else {",
    "    // Define change statistic if Y_i: 0 to 1",
    "    return(0);",
    "  }",
    "};",
    "",
    "// Register function",
    "EFFECT_REGISTER(\"my_mutual\", ::xyz_stat_my_mutual, \"my_mutual\", 0);"
  )

  # 4. Write Files
  writeLines(desc_content, file.path(pkg_path, "DESCRIPTION"))
  writeLines(ns_content, file.path(pkg_path, "NAMESPACE"))
  writeLines(zzz_content, file.path(pkg_path, "R", "zzz.R"))
  writeLines(makevars_content, file.path(pkg_path, "src", "Makevars"))
  writeLines(cpp_content, file.path(pkg_path, "src", "additional_cc.cpp"))

  message(paste0("Package skeleton ", pkg_name, " created at:", normalizePath(pkg_path)))
}

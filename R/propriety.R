#' @title Symmetrise a two-column edgelist
#' @param edges A two-column integer edgelist.
#' @return The edgelist with both orientations of every pair, duplicates removed.
#' @keywords internal
#' @noRd
.symmetrise_edges <- function(edges) {
  if (is.null(edges) || length(edges) == 0 || nrow(edges) == 0) {
    return(matrix(integer(0), nrow = 0, ncol = 2))
  }
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  unique(rbind(edges, edges[, 2:1, drop = FALSE]))
}

#' @title Largest eigenvalue of a symmetric non-negative adjacency structure
#'
#' @description Power iteration on a sparse 0/1 matrix. For a non-negative
#'   symmetric matrix the Perron root equals the spectral radius, so a
#'   non-negative starting vector converges to the largest eigenvalue. Used
#'   instead of \code{base::eigen}, which is O(N^3) and unusable at the
#'   population sizes this package targets; each iteration here is O(|E|).
#'
#' @param edges A two-column edgelist, assumed already symmetrised.
#' @param n_actor Number of units.
#' @param tol Relative convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @return The largest eigenvalue, or 0 if the structure is empty.
#' @keywords internal
#' @noRd
.spectral_radius_edges <- function(edges, n_actor, tol = 1e-9, max_iter = 1000L) {
  if (nrow(edges) == 0 || n_actor < 1) {
    return(0)
  }
  A <- Matrix::sparseMatrix(
    i = as.integer(edges[, 1]), j = as.integer(edges[, 2]),
    x = rep(1, nrow(edges)), dims = c(n_actor, n_actor)
  )
  v <- rep(1 / sqrt(n_actor), n_actor)
  lambda <- 0
  for (iter in seq_len(max_iter)) {
    w <- as.numeric(A %*% v)
    nw <- sqrt(sum(w^2))
    if (!is.finite(nw) || nw == 0) {
      return(0)
    }
    v <- w / nw
    if (abs(nw - lambda) <= tol * max(1, nw)) {
      return(nw)
    }
    lambda <- nw
  }
  lambda
}

#' @title Guard against an improper Gaussian joint distribution
#'
#' @description With \code{type_y = "normal"} the reference measure for
#'   \eqn{Y} is Gaussian, so a positive outcome-outcome coupling makes the
#'   joint density integrable only up to a threshold. Writing \eqn{U} for the
#'   spillover structure \eqn{u_{i,j} = c_{i,j} z_{i,j}}, the conditional
#'   specification is compatible with a joint distribution only while
#'
#'   \deqn{\theta_{y,y,z} \lambda_{\max}(U) / \psi_y < 1.}
#'
#'   Past that point the Gibbs sampler diverges rather than failing: it returns
#'   finite but astronomically large draws before eventually producing
#'   \code{Inf}, with no diagnostic. This check makes the failure explicit.
#'
#'   The test is applied only when it is exact: a single \code{spillover_yy}
#'   term with a positive coefficient and no \code{spillover_yy_scaled}. For a
#'   negative coefficient the binding eigenvalue is the most negative one
#'   rather than the largest, and no other built-in term contributes an
#'   unbounded quadratic form in \eqn{y}, so every other specification is left
#'   to the post-hoc check in \code{simulate_iglm}.
#'
#' @param preprocessed The result of \code{formula_preprocess}.
#' @param coef The coefficient vector, aligned with \code{preprocessed$term_names}.
#' @return Invisibly \code{NULL}. Raises an error if the joint is improper.
#' @keywords internal
#' @noRd
.check_gaussian_propriety <- function(preprocessed, coef) {
  data_object <- preprocessed$data_object
  if (!identical(data_object$type_y, "normal")) {
    return(invisible(NULL))
  }

  terms <- preprocessed$term_names
  if (any(grepl("^spillover_yy_scaled", terms))) {
    return(invisible(NULL))
  }
  idx <- which(terms == "spillover_yy")
  if (length(idx) != 1L) {
    return(invisible(NULL))
  }

  theta <- coef[idx]
  if (!is.finite(theta) || theta <= 0) {
    return(invisible(NULL))
  }

  psi <- data_object$scale_y
  if (!is.finite(psi) || psi <= 0) {
    return(invisible(NULL))
  }

  z_edges <- .symmetrise_edges(data_object$z_network)
  if (nrow(z_edges) == 0) {
    return(invisible(NULL))
  }
  overlap_edges <- .symmetrise_edges(data_object$overlap)
  if (nrow(overlap_edges) == 0) {
    return(invisible(NULL))
  }

  # u_ij = c_ij z_ij: the edges that actually carry spillover.
  z_key <- paste(z_edges[, 1], z_edges[, 2], sep = "-")
  o_key <- paste(overlap_edges[, 1], overlap_edges[, 2], sep = "-")
  u_edges <- z_edges[z_key %in% o_key, , drop = FALSE]
  if (nrow(u_edges) == 0) {
    return(invisible(NULL))
  }

  lambda <- .spectral_radius_edges(u_edges, data_object$n_actor)
  rho <- theta * lambda / psi

  if (is.finite(rho) && rho >= 1) {
    stop(sprintf(
      paste0(
        "The joint distribution implied by this model is improper: ",
        "spillover_yy coefficient (%.4g) x largest eigenvalue of the spillover ",
        "structure (%.4g) / scale_y (%.4g) = %.4g, which must be < 1 when ",
        "type_y = \"normal\". Simulation would diverge rather than converge. ",
        "Reduce the spillover_yy coefficient below %.4g, or increase scale_y."
      ),
      theta, lambda, psi, rho, psi / lambda
    ), call. = FALSE)
  }

  invisible(NULL)
}

#' @title Detect a diverged Gaussian simulation
#'
#' @description Backstop for specifications the a-priori check does not cover:
#'   any non-finite simulated outcome means the chain diverged, which for
#'   \code{type_y = "normal"} means the joint density was not integrable.
#'
#' @param values Simulated outcomes, or any numeric summary of them.
#' @param type_y The outcome type of the data object.
#' @return Invisibly \code{NULL}. Raises an error if the simulation diverged.
#' @keywords internal
#' @noRd
.check_simulation_finite <- function(values, type_y) {
  if (!identical(type_y, "normal")) {
    return(invisible(NULL))
  }
  values <- suppressWarnings(as.numeric(unlist(values, use.names = FALSE)))
  if (length(values) == 0 || !any(!is.finite(values))) {
    return(invisible(NULL))
  }
  stop(paste0(
    "Simulation diverged: non-finite outcomes were produced. With ",
    "type_y = \"normal\" this means the outcome-outcome coupling is too ",
    "strong for the joint distribution to be proper. Reduce the spillover ",
    "coefficient, or increase scale_y."
  ), call. = FALSE)
}

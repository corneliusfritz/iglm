make_test_data <- function(n_actor = 30, degree = 6, seed = 11, scale_y = 1) {
  set.seed(seed)
  p <- degree / (n_actor - 1)
  z <- matrix(0, n_actor, n_actor)
  upper <- which(upper.tri(z))
  z[upper] <- stats::rbinom(length(upper), 1, p)
  z <- z + t(z)
  diag(z) <- 0
  neighborhood <- pmin(z + diag(1, n_actor), 1)
  iglm.data(
    x_attribute = stats::rbinom(n_actor, 1, 0.4),
    y_attribute = stats::rnorm(n_actor),
    z_network = z, neighborhood = neighborhood, directed = FALSE,
    type_x = "binomial", type_y = "normal",
    fix_x = TRUE, fix_z = TRUE, scale_y = scale_y
  )
}

simulate_at <- function(basis, theta_yy, n_actor = 30) {
  sampler <- sampler.iglm(
    n_burn_in = 100, n_simulation = 2, seed = 1,
    sampler_x = sampler.net.attr(n_proposals = 1),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 10),
    sampler_z = sampler.net.attr(n_proposals = 1),
    init_empty = FALSE
  )
  simulate_iglm(
    formula = basis ~ attribute_y + attribute_xy(mode = "global") +
      spillover_yy(mode = "local"),
    basis = basis, coef = c(0, 1, theta_yy), sampler = sampler,
    only_stats = FALSE, fix_x = TRUE, fix_z = TRUE, display_progress = FALSE
  )
}

test_that("power iteration recovers the largest eigenvalue", {
  set.seed(4)
  n <- 25
  a <- matrix(stats::rbinom(n * n, 1, 0.3), n, n)
  a[lower.tri(a)] <- t(a)[lower.tri(a)]
  diag(a) <- 0
  edges <- which(a == 1, arr.ind = TRUE)
  expected <- max(eigen(a, symmetric = TRUE, only.values = TRUE)$values)
  got <- iglm:::.spectral_radius_edges(iglm:::.symmetrise_edges(edges), n)
  expect_equal(got, expected, tolerance = 1e-6)

  # Empty structures are handled rather than erroring.
  expect_equal(iglm:::.spectral_radius_edges(matrix(integer(0), 0, 2), n), 0)
})

test_that("an improper Gaussian specification is rejected before simulating", {
  basis <- make_test_data()
  z <- matrix(0, basis$n_actor, basis$n_actor)
  z[basis$z_network] <- 1
  z[basis$z_network[, 2:1, drop = FALSE]] <- 1
  overlap <- matrix(0, basis$n_actor, basis$n_actor)
  overlap[basis$overlap] <- 1
  lambda <- max(eigen(z * overlap, symmetric = TRUE, only.values = TRUE)$values)
  threshold <- basis$scale_y / lambda

  # Comfortably below the threshold: unchanged behaviour, finite draws.
  below <- simulate_at(basis, threshold * 0.6)
  expect_true(all(is.finite(unlist(lapply(below$samples, function(s) s$y_attribute)))))

  # Above the threshold: an explicit error instead of draws of order 1e96.
  expect_error(
    simulate_at(basis, threshold * 1.2),
    "joint distribution implied by this model is improper"
  )

  # The message names the largest admissible coefficient.
  expect_error(simulate_at(basis, threshold * 1.2), sprintf("%.4g", threshold))
})

test_that("the propriety check does not fire where it does not apply", {
  # Negative coupling: the binding eigenvalue is not the largest one, so the
  # a-priori test is skipped rather than guessed at.
  basis <- make_test_data()
  expect_silent(.check <- iglm:::.check_gaussian_propriety(
    list(
      data_object = basis,
      term_names = c("attribute_y", "attribute_xy_global", "spillover_yy")
    ),
    c(0, 1, -50)
  ))

  # Binary outcomes are bounded, so no threshold exists.
  set.seed(5)
  n <- 20
  z <- matrix(stats::rbinom(n * n, 1, 0.3), n, n)
  z[lower.tri(z)] <- t(z)[lower.tri(z)]
  diag(z) <- 0
  bin <- iglm.data(
    x_attribute = stats::rbinom(n, 1, 0.4), y_attribute = stats::rbinom(n, 1, 0.5),
    z_network = z, neighborhood = pmin(z + diag(1, n), 1), directed = FALSE,
    type_x = "binomial", type_y = "binomial", fix_x = TRUE, fix_z = TRUE
  )
  expect_silent(.check <- iglm:::.check_gaussian_propriety(
    list(data_object = bin, term_names = c("spillover_yy")), c(500)
  ))
})

test_that("scale_y moves the threshold proportionally", {
  # The condition is theta * lambda / psi < 1, so doubling psi doubles the
  # admissible coefficient.
  for (psi in c(1, 4)) {
    basis <- make_test_data(scale_y = psi)
    z <- matrix(0, basis$n_actor, basis$n_actor)
    z[basis$z_network] <- 1
    z[basis$z_network[, 2:1, drop = FALSE]] <- 1
    overlap <- matrix(0, basis$n_actor, basis$n_actor)
    overlap[basis$overlap] <- 1
    lambda <- max(eigen(z * overlap, symmetric = TRUE, only.values = TRUE)$values)
    threshold <- psi / lambda

    expect_error(simulate_at(basis, threshold * 1.05), "improper")
    expect_silent(simulate_at(basis, threshold * 0.5))
  }
})

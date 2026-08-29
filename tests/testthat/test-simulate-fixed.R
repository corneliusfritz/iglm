test_that("Simulation with fix_x preserves observed x", {
  n_actor <- 20
  neighborhood <- matrix(1, nrow = n_actor, ncol = n_actor)
  diag(neighborhood) <- 0

  obs_x <- rep(c(1, 0), length.out = n_actor)
  obs_y <- rep(c(0, 1), length.out = n_actor)
  obs_z <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    x_attribute = obs_x,
    y_attribute = obs_y,
    z_network = obs_z,
    neighborhood = neighborhood,
    directed = FALSE,
    fix_x = TRUE
  )

  # Test with init_empty = TRUE
  sampler_empty <- sampler.iglm(
    init_empty = TRUE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 10)
  )

  mod <- iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(-1, 0.5, 0.5),
    sampler = sampler_empty
  )

  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()
  expect_equal(length(samples), 3)
  for (s in samples) {
    expect_equal(s$x_attribute, obs_x)
  }

  # Test with init_empty = FALSE
  sampler_obs <- sampler.iglm(
    init_empty = FALSE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 10)
  )
  mod$set_sampler(sampler_obs)
  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()
  for (s in samples) {
    expect_equal(s$x_attribute, obs_x)
  }
})

test_that("Simulation with fix_z preserves observed network", {
  n_actor <- 15
  neighborhood <- matrix(1, nrow = n_actor, ncol = n_actor)
  diag(neighborhood) <- 0

  obs_x <- rep(c(1, 0), length.out = n_actor)
  obs_y <- rep(c(0, 1), length.out = n_actor)
  obs_z <- matrix(c(1, 2, 2, 3, 4, 5, 6, 7), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    x_attribute = obs_x,
    y_attribute = obs_y,
    z_network = obs_z,
    neighborhood = neighborhood,
    directed = FALSE,
    fix_z = TRUE
  )

  sampler_empty <- sampler.iglm(
    init_empty = TRUE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 10)
  )

  mod <- iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(-1, 0.5, 0.5),
    sampler = sampler_empty
  )

  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()
  obs_z_canon <- obs_z
  for (i in seq_len(nrow(obs_z_canon))) {
    if (obs_z_canon[i, 1] > obs_z_canon[i, 2]) {
      obs_z_canon[i, ] <- obs_z_canon[i, 2:1]
    }
  }
  obs_z_canon <- obs_z_canon[order(obs_z_canon[, 1], obs_z_canon[, 2]), , drop = FALSE]

  for (s in samples) {
    s_z <- s$z_network
    for (i in seq_len(nrow(s_z))) {
      if (s_z[i, 1] > s_z[i, 2]) s_z[i, ] <- s_z[i, 2:1]
    }
    s_z <- s_z[order(s_z[, 1], s_z[, 2]), , drop = FALSE]
    expect_equal(s_z, obs_z_canon)
  }
})

test_that("Simulation with both fix_x and fix_z preserves both", {
  n_actor <- 15
  neighborhood <- matrix(1, nrow = n_actor, ncol = n_actor)
  diag(neighborhood) <- 0

  obs_x <- rep(c(1, 0), length.out = n_actor)
  obs_y <- rep(c(0, 1), length.out = n_actor)
  obs_z <- matrix(c(1, 2, 2, 3, 4, 5), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    x_attribute = obs_x,
    y_attribute = obs_y,
    z_network = obs_z,
    neighborhood = neighborhood,
    directed = FALSE,
    fix_x = TRUE,
    fix_z = TRUE
  )

  sampler_empty <- sampler.iglm(
    init_empty = TRUE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 10)
  )

  mod <- iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(-1, 0.5, 0.5),
    sampler = sampler_empty
  )

  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()

  obs_z_canon <- obs_z
  for (i in seq_len(nrow(obs_z_canon))) {
    if (obs_z_canon[i, 1] > obs_z_canon[i, 2]) {
      obs_z_canon[i, ] <- obs_z_canon[i, 2:1]
    }
  }
  obs_z_canon <- obs_z_canon[order(obs_z_canon[, 1], obs_z_canon[, 2]), , drop = FALSE]

  for (s in samples) {
    expect_equal(s$x_attribute, obs_x)
    s_z <- s$z_network
    for (i in seq_len(nrow(s_z))) {
      if (s_z[i, 1] > s_z[i, 2]) s_z[i, ] <- s_z[i, 2:1]
    }
    s_z <- s_z[order(s_z[, 1], s_z[, 2]), , drop = FALSE]
    expect_equal(s_z, obs_z_canon)
  }
})

test_that("Simulation with fix_z_alocal preserves non-overlapping edges", {
  n_actor <- 10
  neighborhood <- matrix(0, nrow = n_actor, ncol = n_actor)
  neighborhood[1:5, 1:5] <- 1
  neighborhood[6:10, 6:10] <- 1
  diag(neighborhood) <- 0

  obs_x <- rep(c(1, 0), length.out = n_actor)
  obs_y <- rep(c(0, 1), length.out = n_actor)
  obs_z <- matrix(c(1, 2, 6, 7, 1, 8, 2, 9), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    x_attribute = obs_x,
    y_attribute = obs_y,
    z_network = obs_z,
    neighborhood = neighborhood,
    directed = FALSE,
    fix_z_alocal = TRUE
  )

  sampler_empty <- sampler.iglm(
    init_empty = TRUE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 10)
  )

  mod <- iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(-1, 0.5, 0.5),
    sampler = sampler_empty
  )

  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()

  overlap_mat <- data_obj$overlap
  for (s in samples) {
    s_z <- s$z_network
    non_overlap_edges <- list()
    for (row in seq_len(nrow(s_z))) {
      u <- s_z[row, 1]
      v <- s_z[row, 2]
      is_overlap <- FALSE
      for (o_row in seq_len(nrow(overlap_mat))) {
        if ((overlap_mat[o_row, 1] == u && overlap_mat[o_row, 2] == v) ||
            (overlap_mat[o_row, 1] == v && overlap_mat[o_row, 2] == u)) {
          is_overlap <- TRUE
          break
        }
      }
      if (!is_overlap) {
        non_overlap_edges[[length(non_overlap_edges) + 1]] <- sort(c(u, v))
      }
    }
    non_overlap_mat <- do.call(rbind, non_overlap_edges)
    non_overlap_mat <- non_overlap_mat[order(non_overlap_mat[, 1], non_overlap_mat[, 2]), , drop = FALSE]
    expected_non_overlap <- matrix(c(1, 8, 2, 9), ncol = 2, byrow = TRUE)
    expect_equal(non_overlap_mat, expected_non_overlap)
  }
})

test_that("Estimation with fix_x works and produces valid estimates and SEs", {
  n_actor <- 50
  block <- matrix(nrow = 25, ncol = 25, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(2, block, simplify = FALSE)))
  diag(neighborhood) <- 0

  set.seed(42)
  data_x <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    fix_x = TRUE,
    x_attribute = rbinom(50, 1, 0.5)
  )
  sampler_obj <- sampler.iglm(n_burn_in = 10, n_simulation = 5, init_empty = FALSE)
  m_x <- iglm(
    formula = data_x ~ edges(mode = "local") + attribute_y + attribute_xy,
    coef = c(-2, 0.5, 0.5),
    sampler = sampler_obj,
    control = control.iglm(max_it = 10, display_progress = FALSE)
  )
  m_x$simulate()
  m_x$set_target(m_x$get_samples()[[1]])
  m_x$estimate()

  expect_true(m_x$results$estimated)
  expect_equal(nrow(m_x$coef), 3)
  expect_true(all(!is.na(m_x$coef)))
  expect_true(all(!is.na(m_x$results$var)))

  # Post-estimation simulate should also maintain fix_x
  m_x$simulate()
  for (s in m_x$get_samples()) {
    expect_equal(s$x_attribute, data_x$x_attribute)
  }
})

test_that("Estimation with fix_z works (autologistic model with fixed network)", {
  n_actor <- 50
  block <- matrix(nrow = 25, ncol = 25, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(2, block, simplify = FALSE)))
  diag(neighborhood) <- 0

  set.seed(42)
  sim_net <- matrix(c(1, 2, 2, 3, 3, 4, 26, 27, 27, 28), ncol = 2, byrow = TRUE)
  data_z <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    fix_z = TRUE,
    z_network = sim_net
  )
  sampler_obj <- sampler.iglm(n_burn_in = 10, n_simulation = 5, init_empty = FALSE)
  m_z <- iglm(
    formula = data_z ~ attribute_y + attribute_x + attribute_xy,
    coef = c(-0.5, 0.5, 0.5),
    sampler = sampler_obj,
    control = control.iglm(max_it = 10, display_progress = FALSE)
  )
  m_z$simulate()
  m_z$set_target(m_z$get_samples()[[1]])
  m_z$estimate()

  expect_true(m_z$results$estimated)
  expect_equal(nrow(m_z$coef), 3)
  expect_true(all(!is.na(m_z$coef)))
  expect_true(all(!is.na(m_z$results$var)))

  # Post-estimation simulate should also maintain fix_z
  m_z$simulate()
  for (s in m_z$get_samples()) {
    expect_equal(nrow(s$z_network), nrow(sim_net))
  }
})

test_that("Estimation with both fix_x and fix_z works", {
  n_actor <- 50
  block <- matrix(nrow = 25, ncol = 25, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(2, block, simplify = FALSE)))
  diag(neighborhood) <- 0

  set.seed(42)
  sim_net <- matrix(c(1, 2, 2, 3, 3, 4, 26, 27, 27, 28), ncol = 2, byrow = TRUE)
  data_xz <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    fix_x = TRUE,
    fix_z = TRUE,
    x_attribute = rbinom(50, 1, 0.5),
    z_network = sim_net
  )
  sampler_obj <- sampler.iglm(n_burn_in = 10, n_simulation = 5, init_empty = FALSE)
  m_xz <- iglm(
    formula = data_xz ~ attribute_y + attribute_xy + attribute_yz,
    coef = c(-1, 0.5, 0.2),
    sampler = sampler_obj,
    control = control.iglm(max_it = 10, display_progress = FALSE)
  )
  m_xz$simulate()
  m_xz$set_target(m_xz$get_samples()[[1]])
  m_xz$estimate()

  expect_true(m_xz$results$estimated)
  expect_equal(nrow(m_xz$coef), 3)
  expect_true(all(!is.na(m_xz$coef)))
  expect_true(all(!is.na(m_xz$results$var)))

  # Post-estimation simulate should maintain both
  m_xz$simulate()
  for (s in m_xz$get_samples()) {
    expect_equal(s$x_attribute, data_xz$x_attribute)
    expect_equal(nrow(s$z_network), nrow(sim_net))
  }
})

test_that("Estimation with fix_z_alocal works", {
  n_actor <- 50
  block <- matrix(nrow = 25, ncol = 25, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(2, block, simplify = FALSE)))
  diag(neighborhood) <- 0

  set.seed(42)
  data_alocal <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    fix_z_alocal = TRUE
  )
  sampler_obj <- sampler.iglm(n_burn_in = 10, n_simulation = 5, init_empty = FALSE)
  m_alocal <- iglm(
    formula = data_alocal ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(-2, 0.5, 0.5),
    sampler = sampler_obj,
    control = control.iglm(max_it = 10, display_progress = FALSE)
  )
  m_alocal$simulate()
  m_alocal$set_target(m_alocal$get_samples()[[1]])
  m_alocal$estimate()

  expect_true(m_alocal$results$estimated)
  expect_equal(nrow(m_alocal$coef), 3)
  expect_true(all(!is.na(m_alocal$coef)))
  expect_true(all(!is.na(m_alocal$results$var)))
})

test_that("simulate_iglm validates basis consistency", {
  n_actor <- 10
  data_undir <- iglm.data(n_actor = n_actor, directed = FALSE, type_x = "binomial", type_y = "binomial")
  
  # Non-iglm.data basis
  expect_error(
    simulate_iglm(data_undir ~ edges() + attribute_y, coef = c(-1, 0.5), basis = list()),
    "basis must be an object of class iglm.data"
  )

  # Directed basis vs undirected model
  basis_dir <- iglm.data(n_actor = n_actor, directed = TRUE, type_x = "binomial", type_y = "binomial")
  expect_error(
    simulate_iglm(data_undir ~ edges() + attribute_y, coef = c(-1, 0.5), basis = basis_dir),
    "The 'basis' object must have the same directedness as the model data"
  )

  # Different number of actors
  basis_diff_n <- iglm.data(n_actor = 15, directed = FALSE, type_x = "binomial", type_y = "binomial")
  expect_error(
    simulate_iglm(data_undir ~ edges() + attribute_y, coef = c(-1, 0.5), basis = basis_diff_n),
    "The 'basis' object must have the same number of actors as the model data"
  )

  # Different type_x
  basis_diff_x <- iglm.data(n_actor = n_actor, directed = FALSE, type_x = "normal", type_y = "binomial")
  expect_error(
    simulate_iglm(data_undir ~ edges() + attribute_y, coef = c(-1, 0.5), basis = basis_diff_x),
    "The 'basis' object must have the same type_x as the model data"
  )

  # Different type_y
  basis_diff_y <- iglm.data(n_actor = n_actor, directed = FALSE, type_x = "binomial", type_y = "normal")
  expect_error(
    simulate_iglm(data_undir ~ edges() + attribute_y, coef = c(-1, 0.5), basis = basis_diff_y),
    "The 'basis' object must have the same type_y as the model data"
  )
})

test_that("Simulation with init_empty = TRUE and fix_z_alocal = TRUE preserves non-overlapping ties", {
  n_actor <- 20
  neighborhood <- matrix(0, nrow = n_actor, ncol = n_actor)
  neighborhood[1:10, 1:10] <- 1
  neighborhood[11:20, 11:20] <- 1
  diag(neighborhood) <- 0

  obs_z <- matrix(c(
    1, 2,
    2, 3,
    11, 12,
    1, 11,
    5, 15,
    8, 18
  ), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    z_network = obs_z,
    neighborhood = neighborhood,
    n_actor = n_actor,
    directed = FALSE,
    fix_z_alocal = TRUE
  )

  sampler_empty <- sampler.iglm(
    init_empty = TRUE,
    n_burn_in = 5,
    n_simulation = 3,
    sampler_z = sampler.net.attr(n_proposals = 10)
  )

  mod <- iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y,
    coef = c(-1, 0.5),
    sampler = sampler_empty
  )

  mod$simulate(display_progress = FALSE)
  samples <- mod$get_samples()
  expect_equal(length(samples), 3)

  for (s in samples) {
    s_edges <- s$z_network
    for (pair in list(c(1, 11), c(5, 15), c(8, 18))) {
      has_edge <- any((s_edges[, 1] == pair[1] & s_edges[, 2] == pair[2]) |
                      (s_edges[, 1] == pair[2] & s_edges[, 2] == pair[1]))
      expect_true(has_edge)
    }
  }
})

test_that("simulate_iglm preserves custom label_x, label_y, label_z and supports fix_x, fix_z overrides", {
  n_actor <- 10
  neighborhood <- matrix(1, nrow = n_actor, ncol = n_actor)
  diag(neighborhood) <- 0

  obs_x <- rep(c(1, 0), length.out = n_actor)
  obs_y <- rep(c(0, 1), length.out = n_actor)
  obs_z <- matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE)

  data_obj <- iglm.data(
    x_attribute = obs_x,
    y_attribute = obs_y,
    z_network = obs_z,
    neighborhood = neighborhood,
    directed = FALSE,
    label_x = "covariate_age",
    label_y = "outcome_smoking",
    label_z = "friendship_network"
  )

  sim_res <- simulate_iglm(
    formula = data_obj ~ edges(mode = "local") + attribute_y,
    coef = c(-1, 0.5),
    sampler = sampler.iglm(n_burn_in = 2, n_simulation = 2),
    only_stats = FALSE,
    fix_x = TRUE
  )

  expect_equal(length(sim_res$samples), 2)
  for (s in sim_res$samples) {
    expect_equal(s$label_x, "covariate_age")
    expect_equal(s$label_y, "outcome_smoking")
    expect_equal(s$label_z, "friendship_network")
    expect_equal(s$x_attribute, obs_x)
  }
})



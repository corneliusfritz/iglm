test_that("Define a iglm.data object and check all functions", {
  # debugonce(iglm.data)
  tmp <- iglm.data(
    neighborhood = matrix(c(
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 0, 0, 1,
      0, 1, 1, 0
    ), nrow = 4, byrow = TRUE),
    directed = FALSE,
    type_x = "binomial",
    type_y = "binomial"
  )

  expect_equal(inherits(tmp, "iglm.data"), expected = TRUE)
  expect_equal(tmp$degree()$degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(tmp$mean_x(), expected = 0)
  expect_equal(tmp$mean_y(), expected = 0)
  expect_equal(tmp$mean_z(), expected = 0)
  tmp_name <- paste(tempfile(), ".RDS")
  tmp$save(file = tmp_name)
  rm(tmp)

  loaded_tmp <- iglm.data(file = tmp_name)
  expect_equal(inherits(loaded_tmp, "iglm.data"), expected = TRUE)
  expect_equal(loaded_tmp$degree()$degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(loaded_tmp$mean_x(), expected = 0)
  expect_equal(loaded_tmp$mean_y(), expected = 0)
  expect_equal(loaded_tmp$mean_z(), expected = 0)
  rm(loaded_tmp)
  tmp <- iglm.data(
    z_network = matrix(c(
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 0, 0, 1,
      0, 1, 1, 0
    ), nrow = 4, byrow = TRUE),
    directed = FALSE,
    n_actor = 4, x_attribute = c(0, 0, 1, 0),
    y_attribute = c(0, 1, 0, 1),
    type_x = "binomial",
    type_y = "binomial"
  )
  # debugonce(tmp$degree)
  expect_equal(tmp$degree()$degree_seq, expected = c(2, 2, 2, 2))
  expect_equal(tmp$mean_z(), expected = 4 / 6)
  expect_equal(tmp$mean_x(), expected = 1 / 4)
  expect_equal(tmp$mean_y(), expected = 2 / 4)
  expect_equal(nrow(tmp$overlap) == 12,
    expected = nrow(tmp$neighborhood) == 12
  )
  tmp$save(file = tmp_name)
  rm(tmp)
  loaded_tmp <- iglm.data(file = tmp_name)

  expect_equal(loaded_tmp$degree()$degree_seq, expected = c(2, 2, 2, 2))
  expect_equal(loaded_tmp$mean_z(), expected = 4 / 6)
  expect_equal(loaded_tmp$mean_x(), expected = 1 / 4)
  expect_equal(loaded_tmp$mean_y(), expected = 2 / 4)
  expect_equal(nrow(loaded_tmp$overlap) == 12,
    expected = nrow(loaded_tmp$neighborhood) == 12
  )

  file.remove(tmp_name)
})

test_that("Define a directed iglm.data object and check all functions", {
  tmp <- iglm.data(
    neighborhood = matrix(c(
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 0, 0, 1,
      0, 1, 1, 0
    ), nrow = 4, byrow = TRUE),
    directed = TRUE,
    type_x = "binomial",
    type_y = "binomial"
  )

  expect_equal(inherits(tmp, "iglm.data"), expected = TRUE)
  expect_equal(tmp$degree()$in_degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(tmp$degree()$out_degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(tmp$mean_x(), expected = 0)
  expect_equal(tmp$mean_y(), expected = 0)
  expect_equal(tmp$mean_z(), expected = 0)

  tmp_name <- paste(tempfile(), ".RDS")
  tmp$save(file = tmp_name)
  rm(tmp)

  loaded_tmp <- iglm.data(file = tmp_name)
  expect_equal(inherits(loaded_tmp, "iglm.data"), expected = TRUE)
  expect_equal(loaded_tmp$degree()$in_degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(loaded_tmp$degree()$out_degree_seq, expected = c(0, 0, 0, 0))
  expect_equal(loaded_tmp$mean_x(), expected = 0)
  expect_equal(loaded_tmp$mean_y(), expected = 0)
  expect_equal(loaded_tmp$mean_z(), expected = 0)
  rm(loaded_tmp)

  tmp <- iglm.data(
    z_network = matrix(c(
      0, 1, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 1, 0, 0
    ), nrow = 4, byrow = TRUE),
    directed = TRUE,
    n_actor = 4, x_attribute = c(0, 0, 1, 0),
    y_attribute = c(0, 1, 0, 1),
    type_x = "binomial",
    type_y = "binomial"
  )

  expect_equal(tmp$mean_z(), expected = 5 / 12)
  expect_equal(tmp$mean_x(), expected = 1 / 4)
  expect_equal(tmp$mean_y(), expected = 2 / 4)
  network_tmp <- matrix(c(
    0, 1, 1, 0,
    0, 0, 0, 1,
    0, 0, 0, 1,
    0, 1, 0, 0
  ), nrow = 4, byrow = TRUE)
  expect_equal(tmp$degree()$in_degree_seq, expected = colSums(network_tmp))
  expect_equal(tmp$degree()$out_degree_seq, expected = rowSums(network_tmp))
  expect_equal(nrow(tmp$overlap) == 12,
    expected = nrow(tmp$neighborhood) == 12
  )
  tmp$save(file = tmp_name)
  rm(tmp)
  loaded_tmp <- iglm.data(file = tmp_name)
  expect_equal(loaded_tmp$mean_z(), expected = 5 / 12)
  expect_equal(loaded_tmp$mean_x(), expected = 1 / 4)
  expect_equal(loaded_tmp$mean_y(), expected = 2 / 4)
  expect_equal(loaded_tmp$degree()$in_degree_seq, expected = colSums(network_tmp))
  expect_equal(loaded_tmp$degree()$out_degree_seq, expected = rowSums(network_tmp))
  expect_equal(nrow(loaded_tmp$overlap) == 12,
    expected = nrow(loaded_tmp$neighborhood) == 12
  )
  file.remove(tmp_name)
})

test_that("iglm.data validation throws error when attributes or networks contain NA", {
  x_na <- c(0, 1, NA, 0)
  x_clean <- c(0, 1, 0, 0)
  y_na <- c(1, NA, 1, 0)
  y_clean <- c(1, 0, 1, 0)
  z_clean <- matrix(0, 4, 4)
  z_na <- matrix(c(0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4)

  # Check x_attribute NA
  expect_error(
    iglm.data(x_attribute = x_na, y_attribute = y_clean, z_network = z_clean, n_actor = 4, type_x = "normal", type_y = "binomial"),
    pattern = "'x_attribute' contains missing \\(NA/NaN\\) values."
  )

  # Check y_attribute NA
  expect_error(
    iglm.data(x_attribute = x_clean, y_attribute = y_na, z_network = z_clean, n_actor = 4, type_x = "normal", type_y = "binomial"),
    pattern = "'y_attribute' contains missing \\(NA/NaN\\) values."
  )

  # Check z_network matrix NA
  expect_error(
    iglm.data(x_attribute = x_clean, y_attribute = y_clean, z_network = z_na, n_actor = 4, type_x = "normal", type_y = "binomial"),
    pattern = "'z_network' contains missing \\(NA/NaN\\) values."
  )

  # Check neighborhood matrix NA
  expect_error(
    iglm.data(x_attribute = x_clean, y_attribute = y_clean, z_network = z_clean, neighborhood = z_na, n_actor = 4, type_x = "normal", type_y = "binomial"),
    pattern = "'neighborhood' contains missing \\(NA/NaN\\) values."
  )
})

test_that("spillover_degree_distribution works with custom x_i, x_j, y_i, y_j parameters", {
  n_actor <- 6
  z <- matrix(c(
    0, 1, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 1, 0, 1, 0, 0,
    0, 0, 1, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 1, 0
  ), nrow = 6, byrow = TRUE)

  x <- c(1, 1, 0, 0, 1, 0)
  y <- c(0, 1, 1, 0, 0, 1)

  data_obj <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial")

  # Default (backward compatibility)
  res_def <- data_obj$spillover_degree_distribution(plot = FALSE)
  expect_true(is.list(res_def))
  expect_true("out_spillover_degree" %in% names(res_def))
  expect_true("in_spillover_degree" %in% names(res_def))

  # Custom x_i = 1, y_j = 0
  res_x1_y0 <- data_obj$spillover_degree_distribution(x_i = 1, y_j = 0, plot = FALSE)
  expect_true(is.list(res_x1_y0))
  expect_equal(sum(res_x1_y0$out_spillover_degree), 1)

  # Custom x_i = 0, x_j = 1
  res_x0_x1 <- data_obj$spillover_degree_distribution(x_i = 0, x_j = 1, plot = FALSE)
  expect_true(is.list(res_x0_x1))

  # Custom y_i = 1, y_j = 1
  res_yy <- data_obj$spillover_degree_distribution(y_i = 1, y_j = 1, plot = FALSE)
  expect_true(is.list(res_yy))

  # Custom predicate function
  res_pred <- data_obj$spillover_degree_distribution(x_i = function(val) val == 1, y_j = function(val) val == 0, plot = FALSE)
  expect_equal(res_pred, res_x1_y0)

  # When constraints produce 0 matching senders
  res_zero_senders <- data_obj$spillover_degree_distribution(x_i = 999, y_j = 0, prob = FALSE, plot = FALSE)
  expect_equal(sum(res_zero_senders$out_spillover_degree), 0)
  expect_equal(sum(res_zero_senders$in_spillover_degree), 3) # 3 matching receivers, all degree 0

  # When constraints produce 0 matching senders and receivers
  res_zero_both <- data_obj$spillover_degree_distribution(x_i = 999, y_j = 999, prob = FALSE, plot = FALSE)
  expect_equal(sum(res_zero_both$out_spillover_degree), 0)
  expect_equal(sum(res_zero_both$in_spillover_degree), 0)
})

test_that("degree_distribution works with custom x_i, x_j, y_i, y_j parameters", {
  n_actor <- 6
  z <- matrix(c(
    0, 1, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 1, 0, 1, 0, 0,
    0, 0, 1, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 1, 0
  ), nrow = 6, byrow = TRUE)

  x <- c(1, 1, 0, 0, 1, 0)
  y <- c(0, 1, 1, 0, 0, 1)

  data_dir <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = TRUE)
  data_undir <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = FALSE)

  # Directed default vs constrained
  deg_dir_def <- data_dir$degree_distribution(plot = FALSE)
  expect_true(is.list(deg_dir_def))
  expect_true("in_degree" %in% names(deg_dir_def))
  expect_true("out_degree" %in% names(deg_dir_def))

  deg_dir_sub <- data_dir$degree_distribution(x_i = 1, y_j = 1, plot = FALSE)
  expect_true(is.list(deg_dir_sub))

  # Undirected default vs constrained
  deg_undir_def <- data_undir$degree_distribution(plot = FALSE)
  expect_true(inherits(deg_undir_def, "table"))

  deg_undir_sub <- data_undir$degree_distribution(x_i = 1, y_j = 1, plot = FALSE)
  expect_true(inherits(deg_undir_sub, "table"))
})

test_that("spillover_degree_distribution does not overwrite cached descriptives when constraints are passed", {
  n_actor <- 6
  z <- matrix(c(
    0, 1, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 1, 0, 1, 0, 0,
    0, 0, 1, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 1, 0
  ), nrow = 6, byrow = TRUE)

  x <- c(1, 1, 0, 0, 1, 0)
  y <- c(0, 1, 1, 0, 0, 1)

  data_obj <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = FALSE)

  # 1. Unconstrained call populates descriptives
  res_def <- data_obj$spillover_degree_distribution(plot = FALSE)
  expect_equal(data_obj$descriptives$spillover_degree_distribution, res_def)

  # 2. Constrained call does not overwrite cached descriptives
  res_constrained <- data_obj$spillover_degree_distribution(x_i = 1, y_j = 0, plot = FALSE)
  expect_equal(data_obj$descriptives$spillover_degree_distribution, res_def)
})

test_that("Positional arguments work for distribution methods", {
  n_actor <- 5
  z <- matrix(c(
    0, 1, 1, 0, 0,
    1, 0, 1, 0, 0,
    1, 1, 0, 1, 0,
    0, 0, 1, 0, 1,
    0, 0, 0, 1, 0
  ), nrow = 5, byrow = TRUE)
  x <- c(1, 0, 1, 0, 1)
  y <- c(0, 1, 0, 1, 0)
  data_obj <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = TRUE)

  expect_no_error(data_obj$edgewise_shared_partner_distribution("OTP", c(0, 3), FALSE, FALSE))
  expect_no_error(data_obj$dyadwise_shared_partner_distribution("OTP", c(0, 3), FALSE, FALSE))
  expect_no_error(data_obj$degree_distribution(c(0, 4), FALSE, FALSE))
  expect_no_error(data_obj$spillover_degree_distribution(FALSE, c(0, 3), FALSE))
})

test_that("Undirected constrained degree accounts for canonical edge orientation", {
  n_actor <- 4
  # Only edge is between actor 1 and actor 4 (stored as 1 -> 4)
  z <- matrix(0, nrow = 4, ncol = 4)
  z[1, 4] <- z[4, 1] <- 1
  x <- c(1, 0, 0, 0) # actor 1 has x=1
  y <- c(0, 0, 0, 1) # actor 4 has y=1
  data_obj <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = FALSE)

  # Degree of actor 1 connecting to y=1 receivers (actor 4)
  deg_1 <- data_obj$degree(x_i = 1, y_j = 1)
  expect_equal(unname(deg_1$degree_seq), 1)

  # Degree of actor 4 (x=0, y=1) connecting to x=1 receivers (actor 1)
  deg_4 <- data_obj$degree(y_i = 1, x_j = 1)
  expect_equal(unname(deg_4$degree_seq), 1)
})

test_that("Undirected degree and degree_distribution are equivalent for x_i = 1 vs x_j = 1", {
  n_actor <- 5
  z <- matrix(c(
    0, 1, 1, 0, 0,
    1, 0, 1, 0, 1,
    1, 1, 0, 1, 0,
    0, 0, 1, 0, 1,
    0, 1, 0, 1, 0
  ), nrow = 5, byrow = TRUE)
  x <- c(1, 1, 0, 0, 1) # nodes 1, 2, 5 have x=1
  y <- c(0, 1, 1, 0, 0)

  data_undir <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_actor = n_actor, type_x = "binomial", type_y = "binomial", directed = FALSE)

  # Single constraint equivalence
  deg_xi <- data_undir$degree(x_i = 1)
  deg_xj <- data_undir$degree(x_j = 1)
  expect_equal(deg_xi, deg_xj)

  # Single constraint degree distribution equivalence
  dist_xi <- data_undir$degree_distribution(x_i = 1, plot = FALSE)
  dist_xj <- data_undir$degree_distribution(x_j = 1, plot = FALSE)
  expect_equal(dist_xi, dist_xj)

  # Double constraint x_i = 1, x_j = 1 implies both ends have x = 1 (subgraph on nodes 1, 2, 5)
  # Edges among {1, 2, 5}: (1,2), (2,5) -> node 1 deg=1, node 2 deg=2, node 5 deg=1
  deg_both <- data_undir$degree(x_i = 1, x_j = 1)
  expect_equal(unname(deg_both$degree_seq), c(1, 2, 1))

  # Single constraint counts full degrees for actors with x=1:
  # Node 1: edges to 2, 3 -> deg = 2
  # Node 2: edges to 1, 3, 5 -> deg = 3
  # Node 5: edges to 2, 4 -> deg = 2
  expect_equal(unname(deg_xi$degree_seq), c(2, 3, 2))
})





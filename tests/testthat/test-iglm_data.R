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


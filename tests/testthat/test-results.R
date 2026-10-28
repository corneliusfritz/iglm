test_that("Define a iglm object and check all the results information", {
  n_units <- 20
  block <- matrix(nrow = 5, ncol = 5, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_units / 5, block, simplify = FALSE)))

  overlapping_degree <- 0.5
  neighborhood <- matrix(nrow = n_units, ncol = n_units, data = 0)
  block <- matrix(nrow = 5, ncol = 5, data = 0)
  size_neighborhood <- 5
  size_overlap <- ceiling(size_neighborhood * overlapping_degree)

  end <- floor((n_units - size_neighborhood) / size_overlap)
  for (i in 0:end) {
    neighborhood[(1 + size_overlap * i):(size_neighborhood + size_overlap * i), (1 + size_overlap * i):(size_neighborhood + size_overlap * i)] <- 1
  }
  neighborhood[(n_units - size_neighborhood + 1):(n_units), (n_units - size_neighborhood + 1):(n_units)] <- 1
  family_x <- "binomial"
  family_y <- "binomial"

  xyz_obj_new <- iglm.data(neighborhood = neighborhood, directed = FALSE, family_x = family_x, family_y = family_y)
  xyz_obj_new$overlap
  gt_coef <- c(3, -1, -1)
  gt_coef_pop <- c(rnorm(n = n_units, -2, 1))

  sampler_new <- sampler.iglm(
    n_burn_in = 10, n_simulation = 1,
    sampler_x = sampler.net.attr(n_proposals = n_units * 10),
    sampler_y = sampler.net.attr(n_proposals = n_units * 10),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 10),
    init_empty = F
  )
  # xyz_obj_new$neighborhood
  # xyz_obj_new$z_network <- xyz_obj_new$neighborhood

  model_tmp_new <- iglm(
    formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x + degrees,
    coef = gt_coef, coef_degrees = gt_coef_pop, sampler = sampler_new,
    control = control.iglm(
      accelerated = FALSE, max_it = 200,
      display_progress = TRUE
    )
  )

  tmp_name <- paste(tempfile(), ".RDS")
  model_tmp_new$results$save(file = tmp_name)
  loaded_results <- results(file = tmp_name)
  model_tmp_new$simulate()
  expect_equal(length(loaded_results$samples), expected = 0)
  expect_equal(length(model_tmp_new$results$samples), expected = 1)
  model_tmp_new$results$save(file = tmp_name)
  loaded_results <- results(file = tmp_name)
  expect_equal(length(loaded_results$samples), expected = 1)
  file.remove(tmp_name)
})

test_that("generate_colors does not modify the global RNG seed", {
  set.seed(42)
  seed_before <- .Random.seed
  cols <- generate_colors(5)
  seed_after <- .Random.seed
  expect_identical(seed_before, seed_after)
})

test_that("results class exposes prediction active binding and validates inputs", {
  res <- results(size_coef = 2, size_coef_degrees = 0)
  expect_identical(res$prediction, list())
  expect_error(res$prediction <- list(), regexp = "'prediction' is read-only.", fixed = TRUE)

  # Non-iglm.prediction object should error
  expect_error(res$set_prediction("not_a_prediction"), regexp = "'prediction' must be of class 'iglm.prediction'.", fixed = TRUE)

  pred_obj <- structure(list(x = data.frame(a = 1)), class = "iglm.prediction")
  res$set_prediction(pred_obj)
  expect_identical(res$prediction, pred_obj)

  # Check gather contains prediction
  gathered <- res$gather()
  expect_identical(gathered$prediction, pred_obj)

  # Check save/load
  tmp_file <- tempfile(fileext = ".rds")
  res$save(file = tmp_file)
  res_loaded <- results(file = tmp_file)
  expect_identical(res_loaded$prediction, pred_obj)
  unlink(tmp_file)
})



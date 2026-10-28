test_that("Define a iglm object, simulate, estimate, assess", {
  n_actor <- 100
  block <- matrix(nrow = 50, ncol = 50, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 50, block, simplify = FALSE)))

  overlapping_degree <- 0.5
  neighborhood <- matrix(nrow = n_actor, ncol = n_actor, data = 0)
  block <- matrix(nrow = 5, ncol = 5, data = 0)
  size_neighborhood <- 5
  size_overlap <- ceiling(size_neighborhood * overlapping_degree)

  end <- floor((n_actor - size_neighborhood) / size_overlap)
  for (i in 0:end) {
    neighborhood[(1 + size_overlap * i):(size_neighborhood + size_overlap * i), (1 + size_overlap * i):(size_neighborhood + size_overlap * i)] <- 1
  }
  neighborhood[(n_actor - size_neighborhood + 1):(n_actor), (n_actor - size_neighborhood + 1):(n_actor)] <- 1
  family_x <- "binomial"
  family_y <- "binomial"

  xyz_obj_new <- iglm.data(neighborhood = neighborhood, directed = FALSE, family_x = family_x, family_y = family_y)
  gt_coef <- c(3, -1, -1)
  gt_coef_pop <- c(rnorm(n = n_actor, -2, 1))

  sampler_new <- sampler.iglm(
    n_burn_in = 10, n_simulation = 1,
    sampler_x = sampler.net.attr(n_proposals = n_actor * 10),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 10),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 10),
    init_empty = F
  )

  expect_equal(inherits(sampler_new, "sampler.iglm"), expected = TRUE)

  model_tmp_new <- iglm(
    formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x + degrees,
    coef = gt_coef, coef_degrees = gt_coef_pop, sampler = sampler_new,
    control = control.iglm(accelerated = F, max_it = 200, display_progress = F)
  )


  tmp_name <- paste0(tempfile(), ".rds")
  # debugonce(model_tmp_new$save)
  model_tmp_new$save(file = tmp_name)
  model_tmp_loaded <- iglm(file = tmp_name)

  expect_equal(inherits(model_tmp_loaded, "iglm.object"), expected = TRUE)
  expect_equal(length(model_tmp_loaded$results$samples), expected = 0)
  expect_equal(model_tmp_loaded$results$stats, expected = NULL)

  expect_equal(inherits(model_tmp_new, "iglm.object"), expected = TRUE)
  expect_equal(length(model_tmp_new$results$samples), expected = 0)
  expect_equal(model_tmp_new$results$stats, expected = NULL)
  expect_error(model_tmp_new$get_samples(), "No samples available. Please run 'simulate()' first.", fixed = TRUE)

  model_tmp_new$simulate()

  expect_equal(length(model_tmp_new$results$samples), expected = 1)
  expect_equal(nrow(model_tmp_new$results$stats), expected = 1)
  expect_equal(model_tmp_new$iglm.data$mean_z(), expected = 0)


  expect_equal(inherits(model_tmp_loaded, "iglm.object"), expected = TRUE)
  expect_equal(length(model_tmp_loaded$results$samples), expected = 0)
  expect_equal(model_tmp_loaded$results$stats, expected = NULL)

  samples <- model_tmp_new$get_samples()
  model_tmp_new$set_target(samples[[1]])
  expect_equal(model_tmp_new$iglm.data$mean_z(),
    expected = nrow(samples[[1]]$z_network) / (n_actor * (n_actor - 1) / 2)
  )
  # debugonce(model_tmp_new$estimate)
  expect_error(model_tmp_new$estimate())


  sampler_est <- sampler.iglm(
    n_burn_in = 1, n_simulation = 10,
    sampler_x = sampler.net.attr(n_proposals = n_actor * 10),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 10),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 10),
    init_empty = F
  )

  model_tmp_new$set_sampler(sampler_est)
  expect_equal(model_tmp_new$sampler$n_burn_in, 1)
  # debugonce(model_tmp_new$estimate)
  model_tmp_new$estimate()
  expect_no_warning(model_tmp_new$estimate())
  # expect_equal(as.vector(round(model_tmp_new$coef)), round(gt_coef))
  expect_equal(length(model_tmp_new$results$model_assessment$observed), 0)
  model_tmp_new$assess(formula = ~degree_distribution, plot = FALSE)
  expect_equal(length(model_tmp_new$results$model_assessment$observed), 1)

  model_tmp_new$save(file = tmp_name)
  model_tmp_loaded <- iglm(file = tmp_name)
  model_tmp_new$results$model_assessment$observed

  # expect_equal(as.vector(round(model_tmp_loaded$coef)), round(gt_coef))
  expect_equal(length(model_tmp_loaded$results$model_assessment$observed), 1)

  file.remove(tmp_name)
})

test_that("iglm throws error when covariate object does not exist", {
  n_actor <- 4
  z <- matrix(0, n_actor, n_actor)
  x <- c(0, 1, 1, 0)
  y <- c(1, 0, 1, 0)
  data_obj <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_units = 4, family_x = "binomial", family_y = "binomial")

  expect_error(
    iglm(
      formula = data_obj ~ cov_x(data = non_existent_covariate_var),
      coef = c(1),
      sampler = sampler.iglm(n_burn_in = 2, n_simulation = 1, init_empty = FALSE)
    ),
    regexp = "Could not evaluate argument 'data' in term 'cov_x': object 'non_existent_covariate_var' not found",
    fixed = TRUE
  )

  expect_error(
    iglm(
      formula = data_obj ~ cov_x(non_existent_covariate_var),
      coef = c(1),
      sampler = sampler.iglm(n_burn_in = 2, n_simulation = 1, init_empty = FALSE)
    ),
    regexp = "Could not evaluate argument '..1' in term 'cov_x': object 'non_existent_covariate_var' not found",
    fixed = TRUE
  )

  expect_error(
    iglm(
      formula = data_obj ~ cov_x(data = stop("custom error message")),
      coef = c(1),
      sampler = sampler.iglm(n_burn_in = 2, n_simulation = 1, init_empty = FALSE)
    ),
    regexp = "Could not evaluate argument 'data' in term 'cov_x': custom error message",
    fixed = TRUE
  )

  # Finite coefficient checks
  expect_error(
    iglm(
      formula = data_obj ~ edges(mode = "local"),
      coef = c(NA),
      sampler = sampler.iglm(n_burn_in = 2, n_simulation = 1, init_empty = FALSE)
    ),
    regexp = "coef.*must contain finite numeric values"
  )

  # Empty formula RHS check
  expect_error(
    iglm(
      formula = data_obj ~ 1,
      coef = c(1)
    ),
    regexp = "Formula must contain at least one term"
  )

  # Control parameter bounds & type checks
  expect_error(
    control.iglm(max_it = 10.5),
    regexp = "'max_it' must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    control.iglm(max_it = Inf),
    regexp = "'max_it' must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    control.iglm(tol = Inf),
    regexp = "'tol' must be a positive number",
    fixed = TRUE
  )
  expect_error(
    control.iglm(offset_nonoverlap = Inf),
    regexp = "'offset_nonoverlap' must be a single numeric value",
    fixed = TRUE
  )

  # Pure degree model check
  m_degrees_only <- iglm(
    formula = data_obj ~ degrees,
    sampler = sampler.iglm(n_burn_in = 2, n_simulation = 2, init_empty = FALSE),
    control = control.iglm(max_it = 2, display_progress = FALSE, var_method = "Godambe")
  )
  expect_equal(inherits(m_degrees_only, "iglm.object"), TRUE)
  expect_no_error(m_degrees_only$estimate())
  expect_equal(dim(m_degrees_only$results$var), c(0, 0))

  # Pure degree model with Mean-value variance method
  m_degrees_updated <- iglm(
    formula = data_obj ~ degrees,
    sampler = sampler.iglm(n_burn_in = 2, n_simulation = 2, init_empty = FALSE),
    control = control.iglm(max_it = 2, display_progress = FALSE, var_method = "Mean-value")
  )
  expect_no_error(m_degrees_updated$estimate())
  expect_equal(dim(m_degrees_updated$results$var), c(0, 0))

  # Pure degree model on directed network
  data_obj_dir <- iglm.data(x_attribute = x, y_attribute = y, z_network = z, n_units = 4, directed = TRUE, family_x = "binomial", family_y = "binomial")
  m_degrees_dir <- iglm(
    formula = data_obj_dir ~ degrees,
    sampler = sampler.iglm(n_burn_in = 2, n_simulation = 2, init_empty = FALSE),
    control = control.iglm(max_it = 2, display_progress = FALSE, var_method = "Godambe")
  )
  expect_no_error(m_degrees_dir$estimate())
  expect_equal(dim(m_degrees_dir$results$var), c(0, 0))

  # Test assess removes non-distribution terms like geodesic_distances, dsp, geo, geodist, esp
  expect_warning(
    m_degrees_updated$assess(formula = ~ degree_distribution + geodesic_distances, plot = FALSE),
    regexp = "Unrecognized terms deleted: geodesic_distances",
    fixed = TRUE
  )
  expect_warning(
    m_degrees_updated$assess(formula = ~ degree_distribution + dsp + geo + geodist + esp, plot = FALSE),
    regexp = "Unrecognized terms deleted: dsp, geo, geodist, esp",
    fixed = TRUE
  )
})

test_that("update_formula_remove_terms handles dropping all terms", {
  f <- y ~ x1 + x2
  f_updated <- update_formula_remove_terms(f, c("x1", "x2"))
  expect_equal(deparse(f_updated), "y ~ 0")

  f_one <- update_formula_remove_terms(f, "x1")
  expect_equal(deparse(f_one), "y ~ x2")
})

test_that("iglm accepts NULL sampler without atomic vector error", {
  n <- 4
  dat <- iglm.data(n_units = n, directed = FALSE)
  m <- iglm(dat ~ edges(), sampler = NULL, control = control.iglm(max_it = 1))
  expect_true(inherits(m, "iglm.object"))
  m$estimate()
  expect_equal(colnames(m$results$coefficients_path), "edges()")
})

test_that("rhs_terms_as_list supports one-sided formulas", {
  res <- rhs_terms_as_list(~ edges() + attribute_x)
  expect_equal(length(res), 2)
  expect_true("edges" %in% names(res))
  expect_true("attribute_x" %in% names(res))
})

test_that("iglm loads data flags and custom labels correctly from file", {
  n <- 5
  neighborhood <- matrix(1, n, n)
  dat <- iglm.data(
    n_units = n, directed = FALSE,
    neighborhood = neighborhood,
    fix_x = TRUE, fix_z = TRUE, fix_z_alocal = FALSE,
    label_x = "custom_x", label_y = "custom_y", label_z = "custom_z"
  )
  m <- iglm(dat ~ edges(mode = "local"), sampler = sampler.iglm(n_simulation = 2, init_empty = FALSE), control = control.iglm(max_it = 1))
  tmp_file <- tempfile(fileext = ".rds")
  m$save(file = tmp_file)

  m_loaded <- iglm(file = tmp_file)
  expect_true(m_loaded$iglm.data$fix_x)
  expect_true(m_loaded$iglm.data$fix_z)
  expect_false(m_loaded$iglm.data$fix_z_alocal)
  expect_equal(m_loaded$iglm.data$label_x, "custom_x")
  expect_equal(m_loaded$iglm.data$label_y, "custom_y")
  expect_equal(m_loaded$iglm.data$label_z, "custom_z")
  unlink(tmp_file)
})

test_that("simulate_iglm rejects invalid sampler object", {
  n <- 4
  dat <- iglm.data(n_units = n, directed = FALSE)
  expect_error(
    simulate_iglm(dat ~ edges(), coef = 0, sampler = "not_a_sampler"),
    regexp = "'sampler' must be a 'sampler.iglm' object created with 'sampler.iglm()'.",
    fixed = TRUE
  )
})

test_that("statistics rejects one-sided formulas with informative message", {
  expect_error(
    statistics(~ edges()),
    regexp = "Formula must be two-sided with an 'iglm.data' object on the LHS.",
    fixed = TRUE
  )
})

test_that("predict() handles default arguments, restores control state, and sets results$prediction", {
  n_actor <- 6
  block <- matrix(1, n_actor, n_actor)
  dat <- iglm.data(neighborhood = block, directed = FALSE, family_x = "binomial", family_y = "binomial")

  sampler_obj <- sampler.iglm(
    n_burn_in = 2, n_simulation = 2,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 20),
    init_empty = FALSE
  )

  m <- iglm(
    formula = dat ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(1, -0.5, 0.5),
    sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 5, display_progress = FALSE)
  )

  # Check that predict() errors before estimation
  expect_error(m$predict(), regexp = "Model has not been estimated yet", fixed = TRUE)

  # Estimate model
  m$estimate()
  expect_true(m$results$estimated)

  # Check item 2: predict() with default arguments should not error
  pred <- expect_no_error(m$predict())
  expect_s3_class(pred, "iglm.prediction")
  expect_true(all(c("x", "y", "z") %in% names(pred)))

  # Check item 4: results$prediction active binding returns prediction and is read-only
  expect_identical(m$results$prediction, pred)
  expect_error(m$results$prediction <- list(), regexp = "'prediction' is read-only.", fixed = TRUE)

  # Check item 3: predict(variant = 'conditional') restores control and allows future estimation
  expect_true(m$control$estimate_model)
  pred_cond <- m$predict(variant = "conditional")
  expect_true(m$control$estimate_model)
  expect_no_error(m$estimate())
  expect_true(m$results$estimated)

  # Check item 4: saving and restoring preserves results$prediction
  tmp_file <- tempfile(fileext = ".rds")
  m$save(file = tmp_file)
  m_loaded <- iglm(file = tmp_file)
  expect_identical(m_loaded$results$prediction, pred_cond)
  unlink(tmp_file)
})

test_that("simulate() returns samples invisibly and supports parallel cluster from sampler", {
  n_actor <- 6
  block <- matrix(1, n_actor, n_actor)
  dat <- iglm.data(neighborhood = block, directed = FALSE, family_x = "binomial", family_y = "binomial")

  sampler_obj <- sampler.iglm(
    n_burn_in = 2, n_simulation = 2,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 10),
    sampler_z = sampler.net.attr(n_proposals = 20),
    init_empty = FALSE
  )

  m <- iglm(
    formula = dat ~ edges(mode = "local") + attribute_y,
    coef = c(0.5, -0.5),
    sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )

  # Check item 6: simulate() returns samples invisibly
  sim_samples <- m$simulate()
  expect_true(is.list(sim_samples))
  expect_equal(length(sim_samples), 2)
  expect_identical(sim_samples, m$results$samples)

  # Check item 5: simulate() works with parallel cluster attached to sampler
  m_parallel <- iglm(
    formula = dat ~ edges(mode = "local") + attribute_y,
    coef = c(0.5, -0.5),
    sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )
  cl <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  m_parallel$sampler$set_cluster(cl)
  expect_false(is.null(m_parallel$sampler$cluster))

  sim_parallel <- expect_no_error(m_parallel$simulate())
  expect_equal(length(sim_parallel), 2)
  expect_identical(sim_parallel, m_parallel$results$samples)
})

test_that("predict(variant = 'marginal') automatically triggers simulation when samples are empty", {
  n_actor <- 4
  block <- matrix(1, n_actor, n_actor)
  dat <- iglm.data(neighborhood = block, directed = FALSE)
  m <- iglm(
    formula = dat ~ edges(mode = "local"), coef = -1,
    sampler = sampler.iglm(n_burn_in = 2, n_simulation = 3, init_empty = FALSE),
    control = control.iglm(max_it = 1, display_progress = FALSE)
  )
  m$results$update(estimated = TRUE)
  expect_equal(length(m$results$samples), 0)

  pred <- expect_no_error(m$predict(variant = "marginal"))
  expect_equal(length(m$results$samples), 3)
  expect_s3_class(pred, "iglm.prediction")
})

test_that("predict(variant = 'conditional') correctly uses sender degree coefficients and handles single term", {
  n_actor <- 4
  block <- matrix(1, n_actor, n_actor)
  dat <- iglm.data(neighborhood = block, directed = TRUE)
  # Node 1 has positive expansiveness, node 2 has negative expansiveness
  coef_deg <- c(5, -5, 0, 0, 0, 0, 0, 0)
  m <- iglm(
    formula = dat ~ edges(mode = "local") + degrees,
    coef = 0, coef_degrees = coef_deg,
    control = control.iglm(max_it = 1, display_progress = FALSE)
  )
  m$results$update(estimated = TRUE)

  pred <- expect_no_error(m$predict(variant = "conditional", type = "z"))
  pz <- pred$z
  p_1_3 <- pz$prediction[pz$sender == 1 & pz$receiver == 3]
  p_2_3 <- pz$prediction[pz$sender == 2 & pz$receiver == 3]
  # Sender 1 must have higher predicted tie probability to node 3 than sender 2
  expect_gt(p_1_3, p_2_3)
})









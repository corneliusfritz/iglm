test_that("results$plot error handling works as expected", {
  res <- results(size_coef = 2, size_coef_degrees = 0)

  # At least one must be TRUE
  expect_error(res$plot(), "At least one of `stats`, `trace`, or `model_assessment` must be TRUE.")

  # Stats with no samples
  expect_error(res$plot(stats = TRUE), "No samples available to plot.")

  # Trace when not estimated
  expect_error(res$plot(trace = TRUE), "Model has not been estimated yet. Cannot plot results.")

  # Model assessment when not available
  expect_error(res$plot(model_assessment = TRUE), "No model assessment available to plot.")
})

test_that("results$plot works for trace, stats, and model assessment on undirected model", {
  n_actor <- 20
  block <- matrix(nrow = 5, ncol = 5, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 5, block, simplify = FALSE)))

  xyz_obj <- iglm.data(neighborhood = neighborhood, directed = FALSE, type_x = "binomial", type_y = "binomial")
  gt_coef <- c(1, -1, -1)
  gt_coef_pop <- rnorm(n = n_actor, -2, 1)

  sampler_obj <- sampler.iglm(
    n_burn_in = 2, n_simulation = 5,
    sampler_x = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 2),
    init_empty = FALSE
  )

  model_fit <- iglm(
    formula = xyz_obj ~ edges(mode = "local") + attribute_y + attribute_x + degrees,
    coef = gt_coef, coef_degrees = gt_coef_pop, sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )

  model_fit$simulate()
  model_fit$set_target(model_fit$get_samples()[[1]])
  model_fit$estimate()

  # Test trace plot
  pdf(NULL)
  expect_silent(model_fit$results$plot(trace = TRUE))
  dev.off()

  # Test stats plot
  pdf(NULL)
  expect_silent(model_fit$results$plot(stats = TRUE))
  dev.off()

  # Test single-model assessment plot (undirected)
  model_fit$assess(
    formula = ~ degree_distribution + geodesic_distances_distribution +
      edgewise_shared_partner_distribution + dyadwise_shared_partner_distribution +
      y_distribution + x_distribution,
    plot = FALSE
  )

  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE))
  dev.off()

  # Test multi-model comparison plot
  assess_copy <- model_fit$results$model_assessment
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE, assess_copy))
  dev.off()

  # Test MCMC diagnostics
  model_fit$assess(formula = ~ degree_distribution + mcmc_diagnostics, plot = FALSE)
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE))
  dev.off()
})

test_that("results$plot works for directed model with in/out degrees and continuous attributes", {
  n_actor <- 15
  block <- matrix(nrow = 5, ncol = 5, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 5, block, simplify = FALSE)))

  xyz_obj <- iglm.data(neighborhood = neighborhood, directed = TRUE, type_x = "normal", type_y = "normal")
  gt_coef <- c(1, -0.5, -0.5)

  sampler_obj <- sampler.iglm(
    n_burn_in = 2, n_simulation = 2,
    sampler_x = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 2),
    init_empty = FALSE
  )

  model_fit <- iglm(
    formula = xyz_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = gt_coef, sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )

  model_fit$simulate()
  model_fit$set_target(model_fit$get_samples()[[1]])
  model_fit$estimate()

  # Trace without degrees
  pdf(NULL)
  expect_silent(model_fit$results$plot(trace = TRUE))
  dev.off()

  model_fit$assess(
    formula = ~ degree_distribution + y_distribution + x_distribution +
      spillover_degree_distribution,
    plot = FALSE
  )

  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE))
  dev.off()

  # Multi-model comparison with normal distributions
  assess_copy <- model_fit$results$model_assessment
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE, assess_copy))
  dev.off()

  # Custom spillover parameters in assess formula
  model_fit$assess(
    formula = ~ spillover_degree_distribution(x_i = 1, y_j = 0) +
      spillover_degree_distribution(x_i = 0, y_j = 1),
    plot = FALSE
  )
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE))
  dev.off()

  # Multi-model comparison with multiple custom spillover terms
  assess_custom_copy <- model_fit$results$model_assessment
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE, assess_custom_copy))
  dev.off()

  # Multi-model comparison with distinct second model fit
  model_fit2 <- iglm(
    formula = xyz_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = c(0.2, -0.1, -0.1), sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )
  model_fit2$simulate()
  model_fit2$set_target(model_fit$get_samples()[[1]])
  model_fit2$estimate()

  res2 <- model_fit2$assess(
    formula = ~ degree_distribution + geodesic_distances_distribution + y_distribution + spillover_degree_distribution,
    plot = FALSE
  )
  res1 <- model_fit$assess(
    formula = ~ degree_distribution + geodesic_distances_distribution + y_distribution + spillover_degree_distribution,
    plot = FALSE
  )
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE, res2))
  dev.off()
})




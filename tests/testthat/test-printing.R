test_that("iglm.object print and summary works as expected", {
  n_actor <- 20
  neighborhood <- matrix(1, n_actor, n_actor)
  diag(neighborhood) <- 0
  
  xyz_obj <- iglm.data(neighborhood = neighborhood, directed = FALSE)
  
  # Mock estimation results to trigger printing
  model <- iglm(
    formula = xyz_obj ~ edges(mode = "local") + attribute_y,
    coef = c(-2, 0.5),
    control = control.iglm(estimate_model = FALSE)
  )
  
  # Manually inject mock results since we skipped estimation
  # Use the results$update() method
  model$results$update(
    coefficients_path = matrix(c(-2, 0.5), nrow = 1),
    var = diag(c(0.1, 0.05)),
    estimated = TRUE
  )
  
  # We also need to set .coef in the iglm.object itself
  # Since it's private, we use the environment trick
  assign(".coef", matrix(c(-2, 0.5), ncol = 1), envir = model$.__enclos_env__$private)
  rownames(model$.__enclos_env__$private$.coef) <- c("edges(mode = 'local')", "attribute_y")
  
  assign(".time_estimation", structure(1.23456, units = "secs", class = "difftime"), 
         envir = model$.__enclos_env__$private)
  
  # Test default summary
  out_sum <- capture.output(model$summary())
  expect_true(any(grepl("Results:", out_sum)))
  expect_true(any(grepl("edges", out_sum)))
  expect_false(any(grepl("Formula:", out_sum))) # summary() sets print.formula = FALSE by default
  
  # Test print with formula
  out_print <- capture.output(model$print(print.formula = TRUE))
  expect_true(any(grepl("Formula:", out_print)))
  
  # Test digits argument
  out_digits <- capture.output(model$summary(digits = 1))
  # Check if "1.2" (from 1.23456) is in the output for Time
  expect_true(any(grepl("Time for estimation: 1.2 secs", out_digits)))
  
  # Test disabling results
  out_no_coef <- capture.output(model$print(print.coefmat = FALSE))
  expect_false(any(grepl("Estimate", out_no_coef)))
  
  # Test eps.Pvalue
  # Inject a very small p-value
  model$.__enclos_env__$private$.coef <- matrix(c(-200, 0.5), ncol = 1)
  rownames(model$.__enclos_env__$private$.coef) <- c("edges(mode = 'local')", "attribute_y")
  out_eps <- capture.output(model$summary(eps.Pvalue = 0.5))
  # With eps.Pvalue = 0.5, most small p-values should be shown as < 0.5
  expect_true(any(grepl("< 0.5", out_eps)))
})

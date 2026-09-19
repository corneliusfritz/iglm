test_that("glm arguments passed to iglm() trigger informative stop signals", {
  n_units <- 10
  neighborhood <- matrix(1, nrow = n_units, ncol = n_units)
  xyz_obj <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    type_x = "binomial",
    type_y = "binomial"
  )

  # 1. 'data' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), data = data.frame(a = 1)),
    regexp = "iglm does not accept a 'data' argument in the standard glm sense",
    fixed = TRUE
  )

  # 2. 'family' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), family = binomial()),
    regexp = "'family' is not an argument to iglm()",
    fixed = TRUE
  )

  # 3. 'subset' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), subset = 1:5),
    regexp = "Subsetting via 'subset' is not permitted in iglm",
    fixed = TRUE
  )

  # 4. 'weights' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), weights = rep(1, n_units)),
    regexp = "'weights' is not supported in iglm",
    fixed = TRUE
  )

  # 5. 'offset' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), offset = rep(0, n_units)),
    regexp = "'offset' is not supported in iglm",
    fixed = TRUE
  )

  # 6. 'na.action' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), na.action = na.omit),
    regexp = "'na.action' is not supported in iglm",
    fixed = TRUE
  )

  # 7. 'contrasts' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), contrasts = list()),
    regexp = "'contrasts' is not supported in iglm",
    fixed = TRUE
  )

  # 8. 'method' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), method = "glm.fit"),
    regexp = "'method' (such as 'glm.fit') is not supported in iglm",
    fixed = TRUE
  )

  # 9. 'start' argument trap
  expect_error(
    iglm(xyz_obj ~ edges(), start = c(0, 0)),
    regexp = "'start' is not an argument to iglm()",
    fixed = TRUE
  )

  # 10. Arbitrary unknown named arguments
  expect_error(
    iglm(xyz_obj ~ edges(), custom_unknown_param = 123),
    regexp = "Unrecognized argument(s) passed to iglm(): 'custom_unknown_param'",
    fixed = TRUE
  )

  # 11. Unrecognized positional arguments
  expect_error(
    iglm(xyz_obj ~ edges(), NULL, NULL, NULL, NULL, NULL, NULL, 999),
    regexp = "Unrecognized positional argument(s) passed to iglm()",
    fixed = TRUE
  )
})

test_that("glm arguments passed directly to iglm.object.generator$new trigger informative stop signals", {
  n_units <- 10
  neighborhood <- matrix(1, nrow = n_units, ncol = n_units)
  xyz_obj <- iglm.data(
    neighborhood = neighborhood,
    directed = FALSE,
    type_x = "binomial",
    type_y = "binomial"
  )

  expect_error(
    iglm.object.generator$new(formula = xyz_obj ~ edges(), data = data.frame(a = 1)),
    regexp = "iglm does not accept a 'data' argument in the standard glm sense",
    fixed = TRUE
  )

  expect_error(
    iglm.object.generator$new(formula = xyz_obj ~ edges(), family = binomial()),
    regexp = "'family' is not an argument to iglm()",
    fixed = TRUE
  )
})

test_that("LHS in formula not being an iglm.data object gives informative error", {
  # Data frame column
  df <- data.frame(outcome = c(0, 1, 0, 1))
  expect_error(
    iglm(df$outcome ~ edges()),
    regexp = "The LHS of the formula ('df$outcome') is of class 'numeric', but iglm requires an 'iglm.data' object",
    fixed = TRUE
  )

  # Non-existent column / object (common glm habit)
  expect_error(
    iglm(unassigned_column_name ~ edges()),
    regexp = "The LHS of the formula ('unassigned_column_name') could not be found",
    fixed = TRUE
  )

  # One-sided formula
  expect_error(
    iglm(~ edges()),
    regexp = "Formula must be two-sided with an 'iglm.data' object on the LHS",
    fixed = TRUE
  )
})


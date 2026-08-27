test_that("check.IglmTerm generates informative error messages with term name", {
  # Create a dummy data_object
  data_obj_undirected <- list(directed = FALSE)
  class(data_obj_undirected) <- "iglm.data"

  data_obj_directed <- list(directed = TRUE)
  class(data_obj_directed) <- "iglm.data"

  # Test directed error
  arglist_mutual <- list(base_name = "mutual", label = "mutual(mode = 'global')")
  expect_error(
    check.IglmTerm(data_obj_undirected, arglist_mutual, directed = TRUE),
    pattern = "Term 'mutual' is only for directed networks."
  )

  # Test directed error fallback to label if base_name is missing
  arglist_mutual_no_base <- list(label = "mutual(mode = 'global')")
  expect_error(
    check.IglmTerm(data_obj_undirected, arglist_mutual_no_base, directed = TRUE),
    pattern = "Term 'mutual' is only for directed networks."
  )

  # Test directed error fallback to caller function name
  # We define a function starting with InitIglmTerm.
  InitIglmTerm.test_term <- function(data_obj, arglist) {
    check.IglmTerm(data_obj, arglist, directed = TRUE)
  }
  expect_error(
    InitIglmTerm.test_term(data_obj_undirected, list()),
    pattern = "Term 'test_term' is only for directed networks."
  )

  # Test mandatory argument error
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_mutual, mandatory = "mode"),
    pattern = "Argument 'mode' is mandatory for term 'mutual'."
  )

  # Test expected argument value error
  arglist_invalid_mode <- list(base_name = "mutual", mode = "invalid_mode")
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_invalid_mode, expected = list(mode = c("global", "local"))),
    pattern = "Argument 'mode' of term 'mutual' must be one of: global, local"
  )

  # Test expected type error (numeric)
  arglist_non_numeric <- list(base_name = "cov_z", z_matrix = "not_numeric")
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_non_numeric, expected = list(z_matrix = "numeric")),
    pattern = "Argument 'z_matrix' of term 'cov_z' must be numeric."
  )

  # Test expected type error (matrix)
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_non_numeric, expected = list(z_matrix = "matrix")),
    pattern = "Argument 'z_matrix' of term 'cov_z' must be a matrix or numeric vector."
  )

  # Test expected NA/NaN error for numeric vector
  arglist_na_numeric <- list(base_name = "cov_z", z_matrix = c(1, NA, 3))
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_na_numeric, expected = list(z_matrix = "numeric")),
    pattern = "Argument 'z_matrix' of term 'cov_z' contains missing \\(NA/NaN\\) values."
  )

  # Test expected type error (scalar_numeric)
  arglist_vector_numeric <- list(base_name = "gwdegree", decay = c(0.1, 0.2))
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_vector_numeric, expected = list(decay = "scalar_numeric")),
    pattern = "Argument 'decay' of term 'gwdegree' must be a single numeric value."
  )

  arglist_vector_xi <- list(base_name = "gwdegree", x_i = c(1, 2))
  expect_error(
    check.IglmTerm(data_obj_directed, arglist_vector_xi, expected = list(x_i = "scalar_numeric")),
    pattern = "Argument 'x_i' of term 'gwdegree' must be a single numeric value."
  )
})

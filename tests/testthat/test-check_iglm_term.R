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
})

test_that("check.IglmTerm warns about arguments the term does not use", {
  data_obj <- list(directed = FALSE)
  class(data_obj) <- "iglm.data"

  # An argument the term never declares is dropped by the InitIglmTerm.* method.
  # It should not disappear quietly.
  expect_warning(
    iglm:::check.IglmTerm(data_obj,
      list(base_name = "spillover_yy", mode = "local", data = matrix(1)),
      expected = list(mode = "local"), defaults = list(mode = "local")
    ),
    "Term 'spillover_yy' does not use argument"
  )

  # Declared arguments must not warn, whether supplied or defaulted.
  expect_silent(
    iglm:::check.IglmTerm(data_obj,
      list(base_name = "cov_z", mode = "local", data = matrix(1), type = 1),
      expected = list(mode = c("global", "local", "alocal"), data = "matrix", type = "numeric"),
      defaults = list(mode = "global", data = matrix(1), type = 1)
    )
  )
  expect_silent(
    iglm:::check.IglmTerm(data_obj, list(base_name = "spillover_yy"),
      expected = list(mode = "local"), defaults = list(mode = "local")
    )
  )

  # Mandatory arguments count as known.
  expect_silent(
    iglm:::check.IglmTerm(data_obj, list(base_name = "cov_y", data = 1:3),
      mandatory = "data"
    )
  )

  # Metadata attached by InitIglmTerm is never reported.
  expect_silent(
    iglm:::check.IglmTerm(data_obj,
      list(base_name = "spillover_yy", term_name = "spillover_yy",
           label = "spillover_yy(mode = 'local')", mode = "local"),
      expected = list(mode = "local"), defaults = list(mode = "local")
    )
  )

  # Several unknown arguments are reported together.
  expect_warning(
    iglm:::check.IglmTerm(data_obj,
      list(base_name = "spillover_yy", mode = "local", data = matrix(1), decay = 0.5),
      expected = list(mode = "local"), defaults = list(mode = "local")
    ),
    "data, decay"
  )

  # Falls back to a term-less message when no name can be recovered.
  expect_warning(
    iglm:::check.IglmTerm(data_obj, list(bogus = 1)),
    "^Term does not use argument"
  )
})

test_that("Define a sampler object and test all options", {
  tmp_sampler <- sampler.iglm(
    init_empty = TRUE, n_burn_in = 10, n_simulation = 20,
    seed = 42,
    sampler_x = sampler.net.attr(n_proposals = 10),
    sampler_y = sampler.net.attr(n_proposals = 20),
    sampler_z = sampler.net.attr(n_proposals = 30)
  )
  expect_equal(tmp_sampler$n_burn_in, 10)
  expect_equal(tmp_sampler$n_simulation, 20)
  expect_equal(tmp_sampler$seed, 42)
  expect_equal(inherits(tmp_sampler$sampler_x, "sampler.net.attr"), TRUE)
  expect_equal(inherits(tmp_sampler$sampler_y, "sampler.net.attr"), TRUE)
  expect_equal(inherits(tmp_sampler$sampler_z, "sampler.net.attr"), TRUE)
  expect_equal(tmp_sampler$sampler_x$n_proposals, 10)
  expect_equal(tmp_sampler$sampler_y$n_proposals, 20)
  expect_equal(tmp_sampler$sampler_z$n_proposals, 30)
  tmp_name <- paste(tempfile(), ".RDS")
  tmp_sampler$save(file = tmp_name)
  rm(tmp_sampler)

  loaded_sampler <- sampler.iglm(file = tmp_name)
  expect_equal(loaded_sampler$n_burn_in, 10)
  expect_equal(loaded_sampler$n_simulation, 20)
  expect_equal(loaded_sampler$seed, 42)
  expect_equal(inherits(loaded_sampler$sampler_x, "sampler.net.attr"), TRUE)
  expect_equal(inherits(loaded_sampler$sampler_y, "sampler.net.attr"), TRUE)
  expect_equal(inherits(loaded_sampler$sampler_z, "sampler.net.attr"), TRUE)
  expect_equal(loaded_sampler$sampler_x$n_proposals, 10)
  expect_equal(loaded_sampler$sampler_y$n_proposals, 20)
  expect_equal(loaded_sampler$sampler_z$n_proposals, 30)
  file.remove(tmp_name)
})

test_that("sampler.iglm defaults init_empty to FALSE and uses fixed default seed when seed = NA", {
  s1 <- sampler.iglm()
  s2 <- sampler.iglm()
  expect_false(s1$init_empty)
  expect_false(s2$init_empty)
  expect_equal(s1$seed, 123456789)
  expect_equal(s2$seed, 123456789)

  gen1 <- sampler.iglm.generator$new()
  gen2 <- sampler.iglm.generator$new()
  expect_false(gen1$init_empty)
  expect_false(gen2$init_empty)
  expect_equal(gen1$seed, 123456789)
  expect_equal(gen2$seed, 123456789)
})

test_that("set_x_sampler, set_y_sampler, and set_z_sampler accept sampler.net.attr objects and reject invalid types", {
  s <- sampler.iglm()
  sub_x <- sampler.net.attr(n_proposals = 111)
  sub_y <- sampler.net.attr(n_proposals = 222)
  sub_z <- sampler.net.attr(n_proposals = 333)

  expect_no_error(s$set_x_sampler(sub_x))
  expect_no_error(s$set_y_sampler(sub_y))
  expect_no_error(s$set_z_sampler(sub_z))

  expect_equal(s$sampler_x$n_proposals, 111)
  expect_equal(s$sampler_y$n_proposals, 222)
  expect_equal(s$sampler_z$n_proposals, 333)

  # Test active bindings setter
  sub_x2 <- sampler.net.attr(n_proposals = 444)
  expect_no_error(s$sampler_x <- sub_x2)
  expect_equal(s$sampler_x$n_proposals, 444)

  # Test error on invalid types
  expect_error(s$set_x_sampler("not_a_sampler"), "'sampler_x' must be created with 'sampler.net.attr()'.", fixed = TRUE)
  expect_error(s$set_y_sampler(123), "'sampler_y' must be created with 'sampler.net.attr()'.", fixed = TRUE)
  expect_error(s$set_z_sampler(list()), "'sampler_z' must be created with 'sampler.net.attr()'.", fixed = TRUE)
})

test_that("sampler.iglm active binding init_empty validates logical type with correct punctuation", {
  s <- sampler.iglm()
  expect_error(s$init_empty <- "yes", regexp = "'init_empty' must be logical.", fixed = TRUE)
  expect_error(s$init_empty <- 123, regexp = "'init_empty' must be logical.", fixed = TRUE)
})



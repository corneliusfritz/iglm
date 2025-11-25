
test_that('Define a sampler object and test all options', {
  tmp_sampler <- sampler.iglm(init_empty = TRUE, n_burn_in = 10, n_simulation = 20, 
                                sampler_x = sampler.net.attr(n_proposals = 10, seed = 1),
                                sampler_y = sampler.net.attr(n_proposals = 20, seed = 2),
                                sampler_z = sampler.net.attr(n_proposals = 30, seed = 3))
  expect_equal(tmp_sampler$n_burn_in, 10)
  expect_equal(tmp_sampler$n_simulation, 20)
  expect_equal(inherits(tmp_sampler$sampler_x, "sampler.net.attr"), T)
  expect_equal(inherits(tmp_sampler$sampler_y, "sampler.net.attr"), T)
  expect_equal(inherits(tmp_sampler$sampler_z, "sampler.net.attr"), T)
  expect_equal(tmp_sampler$sampler_x$seed, 1)
  expect_equal(tmp_sampler$sampler_y$seed, 2)
  expect_equal(tmp_sampler$sampler_z$seed, 3)
  expect_equal(tmp_sampler$sampler_x$n_proposals, 10)
  expect_equal(tmp_sampler$sampler_y$n_proposals, 20)
  expect_equal(tmp_sampler$sampler_z$n_proposals, 30)
  tmp_name <- paste(tempfile(), ".RDS")
  tmp_sampler$save(file = tmp_name)
  rm(tmp_sampler)
  
  loaded_sampler <- sampler.iglm(file = tmp_name)
  expect_equal(loaded_sampler$n_burn_in, 10)
  expect_equal(loaded_sampler$n_simulation, 20)
  expect_equal(inherits(loaded_sampler$sampler_x, "sampler.net.attr"), T)
  expect_equal(inherits(loaded_sampler$sampler_y, "sampler.net.attr"), T)
  expect_equal(inherits(loaded_sampler$sampler_z, "sampler.net.attr"), T)
  expect_equal(loaded_sampler$sampler_x$seed, 1)
  expect_equal(loaded_sampler$sampler_y$seed, 2)
  expect_equal(loaded_sampler$sampler_z$seed, 3)
  file.remove(tmp_name)
  
})

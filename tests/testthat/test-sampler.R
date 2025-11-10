
test_that('Define a sampler object and test all options', {
  tmp_sampler <- sampler.iglm(init_empty = TRUE, n_burn_in = 10, n_simulation = 20, 
                                sampler.x = sampler.net_attr(n_proposals = 10, seed = 1),
                                sampler.y = sampler.net_attr(n_proposals = 20, seed = 2),
                                sampler.z = sampler.net_attr(n_proposals = 30, seed = 3))
  expect_equal(tmp_sampler$n_burn_in, 10)
  expect_equal(tmp_sampler$n_simulation, 20)
  expect_equal(inherits(tmp_sampler$sampler.x, "sampler_net_attr"), T)
  expect_equal(inherits(tmp_sampler$sampler.y, "sampler_net_attr"), T)
  expect_equal(inherits(tmp_sampler$sampler.z, "sampler_net_attr"), T)
  expect_equal(tmp_sampler$sampler.x$seed, 1)
  expect_equal(tmp_sampler$sampler.y$seed, 2)
  expect_equal(tmp_sampler$sampler.z$seed, 3)
  expect_equal(tmp_sampler$sampler.x$n_proposals, 10)
  expect_equal(tmp_sampler$sampler.y$n_proposals, 20)
  expect_equal(tmp_sampler$sampler.z$n_proposals, 30)
  tmp_name <- paste(tempfile(), ".RDS")
  tmp_sampler$save(file = tmp_name)
  rm(tmp_sampler)
  
  loaded_sampler <- sampler.iglm(file = tmp_name)
  expect_equal(loaded_sampler$n_burn_in, 10)
  expect_equal(loaded_sampler$n_simulation, 20)
  expect_equal(inherits(loaded_sampler$sampler.x, "sampler_net_attr"), T)
  expect_equal(inherits(loaded_sampler$sampler.y, "sampler_net_attr"), T)
  expect_equal(inherits(loaded_sampler$sampler.z, "sampler_net_attr"), T)
  expect_equal(loaded_sampler$sampler.x$seed, 1)
  expect_equal(loaded_sampler$sampler.y$seed, 2)
  expect_equal(loaded_sampler$sampler.z$seed, 3)
  file.remove(tmp_name)
  
})

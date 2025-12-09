
test_that('Define a iglm object, simulate, estimate, assess', {
  n_actor =100
  block <- matrix(nrow = 50, ncol = 50, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor/50, block, simplify=FALSE)))
  
  overlapping_degree = 0.5
  neighborhood = matrix(nrow = n_actor, ncol = n_actor, data = 0)
  block <- matrix(nrow = 5, ncol = 5, data = 0)
  size_neighborhood = 5
  size_overlap = ceiling(size_neighborhood*overlapping_degree)
  
  end = floor((n_actor-size_neighborhood)/size_overlap)
  for(i in 0:end){
    neighborhood[(1+size_overlap*i):(size_neighborhood+size_overlap*i), (1+size_overlap*i):(size_neighborhood+size_overlap*i)] = 1
  }
  neighborhood[(n_actor-size_neighborhood+1):(n_actor), (n_actor-size_neighborhood+1):(n_actor)] = 1
  type_x <- "binomial"
  type_y <- "binomial"
  
  xyz_obj_new = iglm.data(neighborhood = neighborhood, directed = FALSE, type_x = type_x, type_y = type_y)
  gt_coef = c(3, -1,-1)
  gt_coef_pop =  c(rnorm(n = n_actor, -2, 1))
  
  sampler_new = sampler.iglm(n_burn_in = 10, n_simulation = 1,
                               sampler_x = sampler.net.attr(n_proposals =  n_actor*10,seed = 13),
                               sampler_y = sampler.net.attr(n_proposals =  n_actor*10, seed = 32),
                               sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                               init_empty = F)
  
  expect_equal(inherits(sampler_new,"sampler.iglm"),expected = TRUE)
  
  model_tmp_new <- iglm(formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x + degrees,
                          coef = gt_coef,  coef_degrees = gt_coef_pop, sampler = sampler_new, 
                          control = control.iglm(accelerated = F,max_it = 200, display_progress = F, var = T))

  
  tmp_name <- paste0(tempfile(),".rds")
  # debugonce(model_tmp_new$save)
  model_tmp_new$save(file = tmp_name)
  model_tmp_loaded <- iglm(file = tmp_name)
  
  expect_equal(inherits(model_tmp_loaded,"iglm.object"),expected = TRUE)
  expect_equal(length(model_tmp_loaded$results$samples),expected = 0)
  expect_equal(model_tmp_loaded$results$stats,expected = NULL)
  
  expect_equal(inherits(model_tmp_new,"iglm.object"),expected = TRUE)
  expect_equal(length(model_tmp_new$results$samples),expected = 0)
  expect_equal(model_tmp_new$results$stats,expected = NULL)
  
  model_tmp_new$simulate()
  
  expect_equal(length(model_tmp_new$results$samples),expected = 1)
  expect_equal(nrow(model_tmp_new$results$stats),expected = 1)
  expect_equal(model_tmp_new$iglm.data$mean_z(),expected = 0)
  
  
  expect_equal(inherits(model_tmp_loaded,"iglm.object"),expected = TRUE)
  expect_equal(length(model_tmp_loaded$results$samples),expected = 0)
  expect_equal(model_tmp_loaded$results$stats,expected = NULL)
  
  samples <- model_tmp_new$get_samples()
  model_tmp_new$set_target(samples[[1]])
  expect_equal(model_tmp_new$iglm.data$mean_z(),
                         expected = nrow(samples[[1]]$z_network)/(n_actor *(n_actor-1)/2))
  # debugonce(model_tmp_new$estimate)
  expect_error(model_tmp_new$estimate())
  
  
  sampler_est = sampler.iglm(n_burn_in = 1, n_simulation = 10,
                               sampler_x = sampler.net.attr(n_proposals =  n_actor*10,seed = 1),
                               sampler_y = sampler.net.attr(n_proposals =  n_actor*10, seed = 3),
                               sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10, seed = 13),
                               init_empty = F)
  
  model_tmp_new$set_sampler(sampler_est)
  expect_equal(model_tmp_new$sampler$n_burn_in, 1)
  # debugonce(model_tmp_new$estimate)
  model_tmp_new$estimate()
  expect_no_warning(model_tmp_new$estimate())
  # expect_equal(as.vector(round(model_tmp_new$coef)), round(gt_coef))
  expect_equal(length(model_tmp_new$results$model_assessment$observed), 0)
  model_tmp_new$assess(formula = ~  degree_distribution, plot = FALSE)
  expect_equal(length(model_tmp_new$results$model_assessment$observed), 1)
  
  model_tmp_new$save(file = tmp_name)
  model_tmp_loaded <- iglm(file = tmp_name)
  model_tmp_new$results$model_assessment$observed
  
  # expect_equal(as.vector(round(model_tmp_loaded$coef)), round(gt_coef))
  expect_equal(length(model_tmp_loaded$results$model_assessment$observed), 1)
  
  file.remove(tmp_name)
  
  
})

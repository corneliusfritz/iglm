
test_that('Define a iglm object and check all the results information', {
  n_actors =20
  block <- matrix(nrow = 5, ncol = 5, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actors/5, block, simplify=FALSE)))
  
  overlapping_degree = 0.5
  neighborhood = matrix(nrow = n_actors, ncol = n_actors, data = 0)
  block <- matrix(nrow = 5, ncol = 5, data = 0)
  size_neighborhood = 5
  size_overlap = ceiling(size_neighborhood*overlapping_degree)
  
  end = floor((n_actors-size_neighborhood)/size_overlap)
  for(i in 0:end){
    neighborhood[(1+size_overlap*i):(size_neighborhood+size_overlap*i), (1+size_overlap*i):(size_neighborhood+size_overlap*i)] = 1
  }
  neighborhood[(n_actors-size_neighborhood+1):(n_actors), (n_actors-size_neighborhood+1):(n_actors)] = 1
  type_x <- "binomial"
  type_y <- "binomial"
  
  xyz_obj_new = iglm.data(neighborhood = neighborhood, directed = FALSE, type_x = type_x, type_y = type_y)
  
  gt_coef = c(3, -1,-1)
  gt_coef_pop =  c(rnorm(n = n_actors, -2, 1))
  
  sampler_new = sampler.iglm(n_burn_in = 10, n_simulation = 1,
                               sampler.x = sampler.net_attr(n_proposals =  n_actors*10,seed = 13),
                               sampler.y = sampler.net_attr(n_proposals =  n_actors*10, seed = 32),
                               sampler.z = sampler.net_attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                               init_empty = F)
  model_tmp_new <- iglm(formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x + popularity,
                          coef = gt_coef,  coef_popularity = gt_coef_pop, sampler = sampler_new, 
                          control = control.iglm(accelerated = F,max_it = 200, display_progress = T, var = T))
  
  tmp_name <- paste(tempfile(), ".RDS")
  model_tmp_new$results$save(file = tmp_name)
  loaded_results <- results(file = tmp_name) 
  model_tmp_new$simulate()
  expect_equal(length(loaded_results$samples),expected = 0)
  expect_equal(length(model_tmp_new$results$samples),expected = 1)
  model_tmp_new$results$save(file = tmp_name)
  loaded_results <- results(file = tmp_name) 
  expect_equal(length(loaded_results$samples),expected = 1)
  file.remove(tmp_name)
  })

test_that('Test some sufficient statistics for undirected networks', {
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
  type_x <- "poisson"
  type_y <- "poisson"
  
  xyz_obj_new = iglm.data(neighborhood = neighborhood, directed = FALSE, type_x = type_x, type_y = type_y)
  gt_coef = c(3, -1,-1)
  gt_coef_pop =  c(rnorm(n = n_actor, -2, 1))
  
  sampler_new = sampler.iglm(n_burn_in = 1, n_simulation = 5,
                             sampler_x = sampler.net.attr(n_proposals =  n_actor*100,seed = 13),
                             sampler_y = sampler.net.attr(n_proposals =  n_actor*100, seed = 32),
                             sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                             init_empty = F)
  
  
  model_tmp_new <- iglm(formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x + degrees,
                        coef = gt_coef,  coef_degrees = gt_coef_pop, sampler = sampler_new, 
                        control = control.iglm(accelerated = F,max_it = 200, display_progress = F))
  # debugonce(model_tmp_new$simulate)
  model_tmp_new$print()
  model_tmp_new$simulate()
  
  
  count_values_iglm <- count_statistics(model_tmp_new$results$samples[[1]] ~ 
                                          spillover_xx_scaled(mode = "local")+ 
                                        spillover_yy_scaled(mode = "local") +
                                          spillover_xy_scaled(mode = "local") +
                                          spillover_yx_scaled(mode = "local"))
  # Count the statistics by hand
  tmp <- model_tmp_new$get_samples()
  z_network <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  # Undirected network
  z_network[tmp[[1]]$z_network] <- 1
  z_network[cbind(tmp[[1]]$z_network[,2], tmp[[1]]$z_network[,1])] <- 1
  
  overlap <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  overlap[tmp[[1]]$overlap] <- 1
  val_xx <- c("spillover_xx_scaled(mode = 'local')" =0)
  val_yy <- c("spillover_yy_scaled(mode = 'local')" =0)
  val_xy <- c("spillover_xy_scaled(mode = 'local')" =0)
  val_yx <- c("spillover_yx_scaled(mode = 'local')" =0)
  network_nb <- z_network*overlap
  for(i in 1:tmp[[1]]$n_actor){
    if(sum(network_nb[i,]) ==0){
      next
    }
    val_xx <- val_xx + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_yy <- val_yy + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_xy <- val_xy + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_yx <- val_yx + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
  }
  expect_equal(count_values_iglm[1], val_xx)
  expect_equal(count_values_iglm[2], val_yy)
  expect_equal(count_values_iglm[3], val_xy)
  expect_equal(count_values_iglm[4], val_yx)
})


test_that('Test some sufficient statistics for directed networks', {
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
  type_x <- "poisson"
  type_y <- "poisson"
  
  xyz_obj_new = iglm.data(neighborhood = neighborhood, directed = TRUE, type_x = type_x, type_y = type_y)
  gt_coef = c(3, -1,-1)
  gt_coef_pop =  c(rnorm(n = n_actor, -2, 1))
  
  sampler_new = sampler.iglm(n_burn_in = 10, n_simulation = 1,
                             sampler_x = sampler.net.attr(n_proposals =  n_actor*10,seed = 13),
                             sampler_y = sampler.net.attr(n_proposals =  n_actor*10, seed = 32),
                             sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                             init_empty = F)
  
  model_tmp_new <- iglm(formula = xyz_obj_new ~ edges(mode = "local") + attribute_y + attribute_x,
                        coef = gt_coef, sampler = sampler_new, 
                        control = control.iglm(accelerated = F,max_it = 200, display_progress = F))
  # debugonce(model_tmp_new$simulate)
  model_tmp_new$simulate()
  
  count_values_iglm <- count_statistics(model_tmp_new$results$samples[[1]] ~ 
                                          spillover_xx_scaled(mode = "local") +
                                          spillover_yy_scaled(mode = "local") +
                                          spillover_xy_scaled(mode = "local")+
                                          spillover_yx_scaled(mode = "local") )
  # Count the statistics by hand
  tmp <- model_tmp_new$get_samples()
  z_network <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  # Directed network
  z_network[tmp[[1]]$z_network] <- 1
  
  overlap <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  overlap[tmp[[1]]$overlap] <- 1
  val_xx <- c("spillover_xx_scaled(mode = 'local')" =0)
  val_yy <- c("spillover_yy_scaled(mode = 'local')" =0)
  val_xy <- c("spillover_xy_scaled(mode = 'local')" =0)
  val_yx <- c("spillover_yx_scaled(mode = 'local')" =0)
  network_nb <- z_network*overlap
  for(i in 1:tmp[[1]]$n_actor){
    if(sum(network_nb[i,]) ==0){
      next
    }
    val_xx <- val_xx + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_yy <- val_yy + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_xy <- val_xy + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
    val_yx <- val_yx + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])/sum(network_nb[i,])
  }
  expect_equal(count_values_iglm[1], val_xx)
  expect_equal(count_values_iglm[2], val_yy)
  expect_equal(count_values_iglm[3], val_xy)
  expect_equal(count_values_iglm[4], val_yx)
})


test_that('Test the spillover effects', {
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
  
  xyz_obj_new = iglm.data(neighborhood = neighborhood, directed = TRUE, type_x = type_x, type_y = type_y)
  gt_coef = c(5, -1,-1)
  gt_coef_pop =  c(rnorm(n = n_actor, -2, 1))
  
  sampler_new = sampler.iglm(n_burn_in = 1, n_simulation = 3,
                             sampler_x = sampler.net.attr(n_proposals =  n_actor*10,seed = 13),
                             sampler_y = sampler.net.attr(n_proposals =  n_actor*10, seed = 32),
                             sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                             init_empty = T)
  
  model_tmp_new <- iglm(formula = xyz_obj_new ~ edges(mode = "local") + spillover_yx +
                          spillover_yy,
                        coef = gt_coef, sampler = sampler_new, 
                        control = control.iglm(accelerated = FALSE,max_it = 200, display_progress = FALSE))
  # debugonce(model_tmp_new$simulate)
  model_tmp_new$simulate()
  count_values_iglm <- count_statistics(model_tmp_new$results$samples[[1]] ~ 
                                          spillover_xx+ spillover_yy +spillover_xy+spillover_yx)
  
  # Count the statistics by hand
  tmp <- model_tmp_new$get_samples()
  z_network <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  # Directed network
  z_network[tmp[[1]]$z_network] <- 1
  
  overlap <- matrix(0, nrow = tmp[[1]]$n_actor, ncol = tmp[[1]]$n_actor)
  overlap[tmp[[1]]$overlap] <- 1
  val_xx <- c("spillover_xx" =0)
  val_yy <- c("spillover_yy" =0)
  val_xy <- c("spillover_xy" =0)
  val_yx <- c("spillover_yx" =0)
  network_nb <- z_network*overlap
  for(i in 1:tmp[[1]]$n_actor){
    if(sum(network_nb[i,]) ==0){
      next
    }
    val_xx <- val_xx + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])
    val_yy <- val_yy + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])
    val_xy <- val_xy + sum((tmp[[1]]$x_attribute[i])*tmp[[1]]$y_attribute[network_nb[i,] ==1])
    val_yx <- val_yx + sum((tmp[[1]]$y_attribute[i])*tmp[[1]]$x_attribute[network_nb[i,] ==1])
  }
  expect_equal(count_values_iglm[1], val_xx)
  expect_equal(count_values_iglm[2], val_yy)
  expect_equal(count_values_iglm[2], model_tmp_new$results$stats[1,3])
  
  expect_equal(count_values_iglm[3], val_xy)
  expect_equal(count_values_iglm[4], val_yx)
  expect_equal(count_values_iglm[4], model_tmp_new$results$stats[1,2])
})


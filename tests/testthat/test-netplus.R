test_that('Define a netplus object and check all functions', {
  
  tmp <- netplus(neighborhood = matrix(c(0,1,1,0,
                                         1,0,0,1,
                                         1,0,0,1,
                                         0,1,1,0), nrow=4, byrow=TRUE),
                 directed = FALSE,
                 type_x = "binomial",
                 type_y = "binomial")
  
  expect_equal(inherits(tmp,"netplus"),expected = TRUE)
  expect_equal(tmp$degree()$degree_seq,expected = c(0,0,0,0))
  expect_equal(tmp$density_x(),expected = 0)
  expect_equal(tmp$density_y(),expected = 0)
  expect_equal(tmp$density_z(),expected = 0)
  tmp_name <- paste(tempfile(), ".RDS")
  tmp$save(file = tmp_name)
  rm(tmp)
  
  loaded_tmp <- netplus(file = tmp_name)
  expect_equal(inherits(loaded_tmp,"netplus"),expected = TRUE)
  expect_equal(loaded_tmp$degree()$degree_seq,expected = c(0,0,0,0))
  expect_equal(loaded_tmp$density_x(),expected = 0)
  expect_equal(loaded_tmp$density_y(),expected = 0)
  expect_equal(loaded_tmp$density_z(),expected = 0)
  rm(loaded_tmp)
  
  tmp <- netplus(z_network =  matrix(c(0,1,1,0,
                                       1,0,0,1,
                                       1,0,0,1,
                                       0,1,1,0), nrow=4, byrow=TRUE),
                 directed = FALSE,
                 n_actor = 4,x_attribute = c(0,0,1,0),
                 y_attribute = c(0,1,0,1),
                 type_x = "binomial",
                 type_y = "binomial")
  # debugonce(tmp$degree)
  expect_equal(tmp$degree()$degree_seq,expected = c(2,2,2,2))
  expect_equal(tmp$density_z(),expected =4/6)
  expect_equal(tmp$density_x(),expected =1/4)
  expect_equal(tmp$density_y(),expected =2/4)
  expect_equal(nrow(tmp$overlap) == 12,
               expected =nrow(tmp$neighborhood) == 12)
  tmp$save(file = tmp_name)
  rm(tmp)
  loaded_tmp <- netplus(file = tmp_name)
  
  expect_equal(loaded_tmp$degree()$degree_seq,expected = c(2,2,2,2))
  expect_equal(loaded_tmp$density_z(),expected =4/6)
  expect_equal(loaded_tmp$density_x(),expected =1/4)
  expect_equal(loaded_tmp$density_y(),expected =2/4)
  expect_equal(nrow(loaded_tmp$overlap) == 12,
               expected =nrow(loaded_tmp$neighborhood) == 12)
  
  file.remove(tmp_name)
  
  })

test_that('Define a directed netplus object and check all functions', {
  
  tmp <- netplus(neighborhood =matrix(c(0,1,1,0,
                                        1,0,0,1,
                                        1,0,0,1,
                                        0,1,1,0), nrow=4, byrow=TRUE),
          directed = TRUE,
          type_x = "binomial",
          type_y = "binomial")
  
  expect_equal(inherits(tmp,"netplus"),expected = TRUE)
  expect_equal(tmp$degree()$in_degree_seq,expected = c(0,0,0,0))
  expect_equal(tmp$degree()$out_degree_seq,expected = c(0,0,0,0))
  expect_equal(tmp$density_x(),expected = 0)
  expect_equal(tmp$density_y(),expected = 0)
  expect_equal(tmp$density_z(),expected = 0)
  
  tmp_name <- paste(tempfile(), ".RDS")
  tmp$save(file = tmp_name)
  rm(tmp)
  
  loaded_tmp <- netplus(file = tmp_name)
  expect_equal(inherits(loaded_tmp,"netplus"),expected = TRUE)
  expect_equal(loaded_tmp$degree()$in_degree_seq,expected = c(0,0,0,0))
  expect_equal(loaded_tmp$degree()$out_degree_seq,expected = c(0,0,0,0))
  expect_equal(loaded_tmp$density_x(),expected = 0)
  expect_equal(loaded_tmp$density_y(),expected = 0)
  expect_equal(loaded_tmp$density_z(),expected = 0)
  rm(loaded_tmp)
  
  tmp <- netplus(z_network =  matrix(c(0,1,1,0,
                                       0,0,0,1,
                                       0,0,0,1,
                                       0,1,0,0), nrow=4, byrow=TRUE),
                 directed = TRUE,
                 n_actor = 4,x_attribute = c(0,0,1,0),
                 y_attribute = c(0,1,0,1),
                 type_x = "binomial",
                 type_y = "binomial")
  
  expect_equal(tmp$density_z(),expected =5/12)
  expect_equal(tmp$density_x(),expected =1/4)
  expect_equal(tmp$density_y(),expected =2/4)
  network_tmp <- matrix(c(0,1,1,0,
           0,0,0,1,
           0,0,0,1,
           0,1,0,0), nrow=4, byrow=TRUE)
  expect_equal(tmp$degree()$in_degree_seq,expected = colSums(network_tmp))
  expect_equal(tmp$degree()$out_degree_seq,expected = rowSums(network_tmp))
  expect_equal(nrow(tmp$overlap) == 12,
                         expected =nrow(tmp$neighborhood) == 12)
  tmp$save(file = tmp_name)
  rm(tmp)
  loaded_tmp <- netplus(file = tmp_name)
  expect_equal(loaded_tmp$density_z(),expected =5/12)
  expect_equal(loaded_tmp$density_x(),expected =1/4)
  expect_equal(loaded_tmp$density_y(),expected =2/4)
  expect_equal(loaded_tmp$degree()$in_degree_seq,expected = colSums(network_tmp))
  expect_equal(loaded_tmp$degree()$out_degree_seq,expected = rowSums(network_tmp))
  expect_equal(nrow(loaded_tmp$overlap) == 12,
               expected =nrow(loaded_tmp$neighborhood) == 12)
  file.remove(tmp_name)
  
})
  
test_that("dyadwise_shared_partner and edgewise_shared_partner work with mode in undirected networks", {
  # 4-node network: edges (1,2), (1,3), (2,3), (3,4)
  z <- matrix(c(
    0, 1, 1, 0,
    1, 0, 1, 0,
    1, 1, 0, 1,
    0, 0, 1, 0
  ), nrow = 4, byrow = TRUE)

  # Neighborhood: block 1: {1, 2, 3}, block 2: {4}
  nb <- matrix(c(
    1, 1, 1, 0,
    1, 1, 1, 0,
    1, 1, 1, 0,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)

  data_obj <- iglm.data(
    z_network = z,
    neighborhood = nb,
    directed = FALSE,
    type_x = "binomial",
    type_y = "binomial"
  )

  # Hand-calculated dyadwise shared partners (upper triangular):
  # (1,2) -> common neighbor 3 -> 1
  # (1,3) -> common neighbor 2 -> 1
  # (1,4) -> common neighbor 3 -> 1
  # (2,3) -> common neighbor 1 -> 1
  # (2,4) -> common neighbor 3 -> 1
  # (3,4) -> common neighbor none -> 0

  # Global dyadwise shared partners (default)
  dsp_global <- as.matrix(data_obj$dyadwise_shared_partner(mode = "global"))
  expect_equal(dsp_global[1, 2], 1)
  expect_equal(dsp_global[1, 3], 1)
  expect_equal(dsp_global[1, 4], 1)
  expect_equal(dsp_global[2, 3], 1)
  expect_equal(dsp_global[2, 4], 1)
  expect_equal(dsp_global[3, 4], 0)
  expect_true(is.na(dsp_global[1, 1]))
  expect_true(is.na(dsp_global[2, 1]))

  # Local dyadwise shared partners (only pairs within {1, 2, 3} are kept)
  dsp_local <- as.matrix(data_obj$dyadwise_shared_partner(mode = "local"))
  expect_equal(dsp_local[1, 2], 1)
  expect_equal(dsp_local[1, 3], 1)
  expect_equal(dsp_local[2, 3], 1)
  expect_true(is.na(dsp_local[1, 4]))
  expect_true(is.na(dsp_local[2, 4]))
  expect_true(is.na(dsp_local[3, 4]))

  # Dyadwise distributions
  dist_global <- data_obj$dyadwise_shared_partner_distribution(mode = "global", value_range = c(0, 1), prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(dist_global["0"]), 1)
  expect_equal(as.numeric(dist_global["1"]), 5)

  dist_local <- data_obj$dyadwise_shared_partner_distribution(mode = "local", value_range = c(0, 1), prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(dist_local["0"]), 0)
  expect_equal(as.numeric(dist_local["1"]), 3)

  # Global edgewise shared partners: edges are (1,2), (1,3), (2,3), (3,4) -> counts 1, 1, 1, 0
  esp_global <- data_obj$edgewise_shared_partner(mode = "global")
  expect_equal(sort(esp_global), c(0, 1, 1, 1))

  # Local edgewise shared partners: only edges (1,2), (1,3), (2,3) -> counts 1, 1, 1
  esp_local <- data_obj$edgewise_shared_partner(mode = "local")
  expect_equal(sort(esp_local), c(1, 1, 1))

  # Edgewise distributions
  edist_global <- data_obj$edgewise_shared_partner_distribution(mode = "global", value_range = c(0, 1), prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(edist_global["0"]), 1)
  expect_equal(as.numeric(edist_global["1"]), 3)

  edist_local <- data_obj$edgewise_shared_partner_distribution(mode = "local", value_range = c(0, 1), prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(edist_local["0"]), 0)
  expect_equal(as.numeric(edist_local["1"]), 3)
})

test_that("dyadwise_shared_partner and edgewise_shared_partner work with mode in directed networks", {
  # 4-node directed network: edges 1->2, 1->3, 2->3, 3->4
  z <- matrix(c(
    0, 1, 1, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE)

  # Neighborhood: block 1: {1, 2, 3}, block 2: {4}
  nb <- matrix(c(
    1, 1, 1, 0,
    1, 1, 1, 0,
    1, 1, 1, 0,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)

  data_obj <- iglm.data(
    z_network = z,
    neighborhood = nb,
    directed = TRUE,
    type_x = "binomial",
    type_y = "binomial"
  )

  # 1. OTP (Outgoing Two-Path: z_{i,h} * z_{j,h})
  # (1,2) and (2,1) both send to 3 -> count 1
  # All other pairs -> 0
  otp_global <- as.matrix(data_obj$dyadwise_shared_partner(type = "OTP", mode = "global"))
  expect_equal(otp_global[1, 2], 1)
  expect_equal(otp_global[2, 1], 1)
  expect_equal(otp_global[1, 3], 0)
  expect_true(is.na(otp_global[1, 1]))

  otp_local <- as.matrix(data_obj$dyadwise_shared_partner(type = "OTP", mode = "local"))
  expect_equal(otp_local[1, 2], 1)
  expect_equal(otp_local[2, 1], 1)
  expect_equal(otp_local[1, 3], 0)
  expect_equal(otp_local[3, 1], 0)
  expect_equal(otp_local[2, 3], 0)
  expect_equal(otp_local[3, 2], 0)
  expect_true(is.na(otp_local[1, 4]))
  expect_true(is.na(otp_local[4, 1]))

  otp_dist_local <- data_obj$dyadwise_shared_partner_distribution(type = "OTP", mode = "local", prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(otp_dist_local["0"]), 4)
  expect_equal(as.numeric(otp_dist_local["1"]), 2)

  # 2. ISP (Ingoing Shared Partner / Two-Path: z_{i,h} * z_{h,j} -> path i -> h -> j)
  # (1,3) via 2: 1->2->3 (count 1)
  # (1,4) via 3: 1->3->4 (count 1)
  # (2,4) via 3: 2->3->4 (count 1)
  isp_global <- as.matrix(data_obj$dyadwise_shared_partner(type = "ISP", mode = "global"))
  expect_equal(isp_global[1, 3], 1)
  expect_equal(isp_global[1, 4], 1)
  expect_equal(isp_global[2, 4], 1)
  expect_equal(isp_global[1, 2], 0)

  isp_local <- as.matrix(data_obj$dyadwise_shared_partner(type = "ISP", mode = "local"))
  expect_equal(isp_local[1, 3], 1)
  expect_equal(isp_local[1, 2], 0)
  expect_equal(isp_local[2, 1], 0)
  expect_equal(isp_local[3, 1], 0)
  expect_equal(isp_local[2, 3], 0)
  expect_equal(isp_local[3, 2], 0)
  expect_true(is.na(isp_local[1, 4]))
  expect_true(is.na(isp_local[2, 4]))

  isp_dist_local <- data_obj$dyadwise_shared_partner_distribution(type = "ISP", mode = "local", prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(isp_dist_local["0"]), 5)
  expect_equal(as.numeric(isp_dist_local["1"]), 1)

  # 3. ITP (Incoming Two-Path: z_{h,i} * z_{h,j} -> h sends to both i and j)
  # 1 sends to 2 and 3 -> (2,3) = 1, (3,2) = 1
  itp_global <- as.matrix(data_obj$dyadwise_shared_partner(type = "ITP", mode = "global"))
  expect_equal(itp_global[2, 3], 1)
  expect_equal(itp_global[3, 2], 1)
  expect_equal(itp_global[1, 2], 0)

  itp_local <- as.matrix(data_obj$dyadwise_shared_partner(type = "ITP", mode = "local"))
  expect_equal(itp_local[2, 3], 1)
  expect_equal(itp_local[3, 2], 1)
  expect_equal(itp_local[1, 2], 0)
  expect_equal(itp_local[2, 1], 0)
  expect_equal(itp_local[1, 3], 0)
  expect_equal(itp_local[3, 1], 0)
  expect_true(is.na(itp_local[1, 4]))

  itp_dist_local <- data_obj$dyadwise_shared_partner_distribution(type = "ITP", mode = "local", prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(itp_dist_local["0"]), 4)
  expect_equal(as.numeric(itp_dist_local["1"]), 2)

  # 4. OSP (Outgoing Shared Partner: z_{h,i} * z_{j,h} -> j -> h -> i)
  # (3,1) via 2: 1->2->3 (count 1)
  # (4,1) via 3: 1->3->4 (count 1)
  # (4,2) via 3: 2->3->4 (count 1)
  osp_global <- as.matrix(data_obj$dyadwise_shared_partner(type = "OSP", mode = "global"))
  expect_equal(osp_global[3, 1], 1)
  expect_equal(osp_global[4, 1], 1)
  expect_equal(osp_global[4, 2], 1)

  osp_local <- as.matrix(data_obj$dyadwise_shared_partner(type = "OSP", mode = "local"))
  expect_equal(osp_local[3, 1], 1)
  expect_equal(osp_local[1, 3], 0)
  expect_true(is.na(osp_local[4, 1]))

  osp_dist_local <- data_obj$dyadwise_shared_partner_distribution(type = "OSP", mode = "local", prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(osp_dist_local["0"]), 5)
  expect_equal(as.numeric(osp_dist_local["1"]), 1)

  # Directed edgewise shared partners:
  # Edges: 1->2, 1->3, 2->3 (in overlap), 3->4 (outside overlap)
  # OTP on local edges:
  # (1,2) count 1; (1,3) count 0; (2,3) count 0
  esp_otp_local <- data_obj$edgewise_shared_partner(type = "OTP", mode = "local")
  expect_equal(sort(esp_otp_local), c(0, 0, 1))

  edist_otp_local <- data_obj$edgewise_shared_partner_distribution(type = "OTP", mode = "local", prob = FALSE, plot = FALSE)
  expect_equal(as.numeric(edist_otp_local["0"]), 2)
  expect_equal(as.numeric(edist_otp_local["1"]), 1)
})

test_that("dyadwise and edgewise shared partner functions validate inputs properly", {
  n_actor <- 4
  nb <- matrix(1, 4, 4)
  diag(nb) <- 0
  data_obj <- iglm.data(neighborhood = nb, directed = FALSE, type_x = "binomial", type_y = "binomial")

  # Invalid mode argument
  expect_error(data_obj$dyadwise_shared_partner(mode = "invalid"), "'mode' must be either 'global' or 'local'.")
  expect_error(data_obj$dyadwise_shared_partner(mode = TRUE), "'mode' must be either 'global' or 'local'.")
  expect_error(data_obj$dyadwise_shared_partner(mode = 123), "'mode' must be either 'global' or 'local'.")

  expect_error(data_obj$edgewise_shared_partner(mode = "invalid"), "'mode' must be either 'global' or 'local'.")
  expect_error(data_obj$edgewise_shared_partner(mode = 123), "'mode' must be either 'global' or 'local'.")

  expect_error(data_obj$dyadwise_shared_partner_distribution(mode = "foo"), "'mode' must be either 'global' or 'local'.")
  expect_error(data_obj$edgewise_shared_partner_distribution(mode = "foo"), "'mode' must be either 'global' or 'local'.")

  # Invalid type
  expect_error(data_obj$dyadwise_shared_partner(type = "INVALID"), "type must be one of 'OTP', 'ISP', 'OSP', 'ITP', 'ALL', or 'symm'.")
  expect_error(data_obj$edgewise_shared_partner(type = "INVALID"), "type must be one of 'OTP', 'ISP', 'OSP', 'ITP', 'ALL', or 'symm'.")
  expect_error(data_obj$dyadwise_shared_partner_distribution(type = "INVALID"), "type must be one of 'OTP', 'ISP', 'OSP', 'ITP', 'ALL', or 'symm'.")
  expect_error(data_obj$edgewise_shared_partner_distribution(type = "INVALID"), "type must be one of 'OTP', 'ISP', 'OSP', 'ITP', 'ALL', or 'symm'.")

  # Directed types on undirected network
  expect_error(data_obj$dyadwise_shared_partner(type = "OTP"), "Type 'OTP' is only for directed networks.")
  expect_error(data_obj$edgewise_shared_partner(type = "ISP"), "Type 'ISP' is only for directed networks.")
  expect_error(data_obj$dyadwise_shared_partner_distribution(type = "OSP"), "Type 'OSP' is only for directed networks.")
  expect_error(data_obj$edgewise_shared_partner_distribution(type = "ITP"), "Type 'ITP' is only for directed networks.")

  # Undirected type (symm) on directed network
  data_dir <- iglm.data(neighborhood = nb, directed = TRUE, type_x = "binomial", type_y = "binomial")
  expect_error(data_dir$dyadwise_shared_partner(type = "symm"), "Type 'symm' is only for undirected networks.")
  expect_error(data_dir$edgewise_shared_partner(type = "symm"), "Type 'symm' is only for undirected networks.")
  expect_error(data_dir$dyadwise_shared_partner_distribution(type = "symm"), "Type 'symm' is only for undirected networks.")
  expect_error(data_dir$edgewise_shared_partner_distribution(type = "symm"), "Type 'symm' is only for undirected networks.")

  # Invalid value_range
  expect_error(data_obj$dyadwise_shared_partner_distribution(value_range = c(-1, 2)), "'value_range' values must be non-negative.")
  expect_error(data_obj$edgewise_shared_partner_distribution(value_range = c(-1, 2)), "'value_range' values must be non-negative.")
  expect_error(data_obj$dyadwise_shared_partner_distribution(value_range = 1:3), "'value_range' must be a numeric vector of length 2.")
  expect_error(data_obj$edgewise_shared_partner_distribution(value_range = 1:3), "'value_range' must be a numeric vector of length 2.")
})

test_that("assess and results$plot work with mode in shared partner distributions", {
  n_actor <- 15
  block <- matrix(nrow = 5, ncol = 5, data = 1)
  neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 5, block, simplify = FALSE)))

  xyz_obj <- iglm.data(neighborhood = neighborhood, directed = FALSE, type_x = "binomial", type_y = "binomial")
  gt_coef <- c(1, -0.5, -0.5)

  sampler_obj <- sampler.iglm(
    n_burn_in = 2, n_simulation = 2,
    sampler_x = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_y = sampler.net.attr(n_proposals = n_actor * 2),
    sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 2)
  )

  model_fit <- iglm(
    formula = xyz_obj ~ edges(mode = "local") + attribute_y + attribute_x,
    coef = gt_coef, sampler = sampler_obj,
    control = control.iglm(accelerated = FALSE, max_it = 2, display_progress = FALSE)
  )

  model_fit$simulate()
  model_fit$set_target(model_fit$get_samples()[[1]])
  model_fit$estimate()

  # Run assess with mode in dyadwise and edgewise shared partner distributions
  assessment <- model_fit$assess(
    formula = ~ dyadwise_shared_partner_distribution(mode = "local") +
      edgewise_shared_partner_distribution(mode = "local") +
      dyadwise_shared_partner_distribution(mode = "global") +
      edgewise_shared_partner_distribution(mode = "global"),
    plot = FALSE
  )

  expect_true(inherits(assessment, "iglm_model_assessment"))
  expect_true("dyadwise_shared_partner_distribution_mode_local" %in% names(assessment$observed))
  expect_true("edgewise_shared_partner_distribution_mode_local" %in% names(assessment$observed))
  expect_true("dyadwise_shared_partner_distribution_mode_global" %in% names(assessment$observed))
  expect_true("edgewise_shared_partner_distribution_mode_global" %in% names(assessment$observed))

  # Test plotting without errors
  pdf(NULL)
  expect_silent(model_fit$results$plot(model_assessment = TRUE))
  dev.off()
})

library(iglm)

test_that("spillover_yx_scaled_local gives correct statistics for an empty to one-edge network", {
    n_actor <- 3
    adj <- matrix(0, n_actor, n_actor)
    x <- c(1, 2, 3)
    y <- c(1, 4, 1)

    data_obj <- iglm:::iglm.data(
        x_attribute = x,
        y_attribute = y,
        type_x = "normal",
        type_y = "normal",
        z_network = adj,
        directed = TRUE,
        n_actor = n_actor
    )

    S_with <- x[2]
    d_with <- 1
    A_with <- S_with / d_with
    hand_global_empty <- 0
    hand_global_edge12 <- y[1] * A_with

    res_empty <- iglm:::xyz_count_global(
        z_network = adj,
        x_attribute = x,
        y_attribute = y,
        neighborhood = matrix(1, n_actor, n_actor) - diag(n_actor),
        overlap = matrix(1, n_actor, n_actor) - diag(n_actor),
        directed = TRUE,
        terms = "spillover_yx_scaled_local",
        n_actor = n_actor,
        data_list = list(matrix(0)),
        type_list = c(0),
        type_x = "normal",
        type_y = "normal",
        attr_x_scale = 1.0,
        attr_y_scale = 1.0
    )

    adj2 <- adj
    adj2[1, 2] <- 1
    res_edge12 <- iglm:::xyz_count_global(
        z_network = adj2,
        x_attribute = x,
        y_attribute = y,
        neighborhood = matrix(1, n_actor, n_actor) - diag(n_actor),
        overlap = matrix(1, n_actor, n_actor) - diag(n_actor),
        directed = TRUE,
        terms = "spillover_yx_scaled_local",
        n_actor = n_actor,
        data_list = list(matrix(0)),
        type_list = c(0),
        type_x = "normal",
        type_y = "normal",
        attr_x_scale = 1.0,
        attr_y_scale = 1.0
    )

    expect_equal(as.numeric(res_empty), hand_global_empty)
    expect_equal(as.numeric(res_edge12), hand_global_edge12)
})

test_that("spillover_yx_scaled_local gives correct statistics for a full network", {
    # 1. Setup Data for a FULL network
    n_actor <- 3
    adj_full <- matrix(1, n_actor, n_actor) - diag(n_actor)
    x <- c(1, 2, 3)
    y <- c(1, 4, 1)

    data_obj_full <- iglm:::iglm.data(
        x_attribute = x,
        y_attribute = y,
        type_x = "normal",
        type_y = "normal",
        z_network = adj_full,
        directed = TRUE,
        n_actor = n_actor
    )

    # 2. Hand calculation for global statistic on FULL network
    # Formula: sum_i Y_i * (sum_{j in N(i)} X_j) / d_i
    # Since it's a full network, N(i) is everyone except i, d_i = 2
    sum_x <- sum(x)
    hand_global_full <- sum(y * (sum_x - x) / (n_actor - 1))

    # 3. Calculate via internal count_global
    res_full <- iglm:::xyz_count_global(
        z_network = adj_full,
        x_attribute = x,
        y_attribute = y,
        neighborhood = matrix(1, n_actor, n_actor) - diag(n_actor),
        overlap = matrix(1, n_actor, n_actor) - diag(n_actor),
        directed = TRUE,
        terms = "spillover_yx_scaled_local",
        n_actor = n_actor,
        data_list = list(matrix(0)),
        type_list = c(0),
        type_x = "normal",
        type_y = "normal",
        attr_x_scale = 1.0,
        attr_y_scale = 1.0
    )
    expect_equal(as.numeric(res_full), hand_global_full)

    # 4. Check via simulate_iglm
    formula_full <- data_obj_full ~ spillover_yx_scaled_local
    sampler <- iglm:::sampler.iglm(
        n_simulation = 1,
        n_burn_in = 0,
        init_empty = FALSE,
        sampler_x = iglm:::sampler.net.attr(n_proposals = 0),
        sampler_y = iglm:::sampler.net.attr(n_proposals = 0),
        sampler_z = iglm:::sampler.net.attr(n_proposals = 0)
    )
    res_sim_full <- simulate_iglm(formula = formula_full, coef = c(0), sampler = sampler, only_stats = TRUE)
    expect_equal(as.numeric(res_sim_full$stats[1]), hand_global_full)
})


test_that("spillover_xy_scaled_local gives correct statistics for a full network", {
    # 1. Setup Data for a FULL network
    n_actor <- 3
    adj_full <- matrix(1, n_actor, n_actor) - diag(n_actor)
    x <- c(1, 2, 3)
    y <- c(1, 4, 1)

    data_obj_full <- iglm:::iglm.data(
        x_attribute = x,
        y_attribute = y,
        type_x = "normal",
        type_y = "normal",
        z_network = adj_full,
        directed = TRUE,
        n_actor = n_actor
    )

    # 2. Hand calculation for global statistic on FULL network
    # Formula: sum_i X_i * (sum_{j in N(i)} Y_j) / d_i
    # Since it's a full network, N(i) is everyone except i, d_i = 2
    sum_y <- sum(y)
    hand_global_full <- sum(x * (sum_y - y) / (n_actor - 1))

    # 3. Calculate via internal count_global
    res_full <- iglm:::xyz_count_global(
        z_network = adj_full,
        x_attribute = x,
        y_attribute = y,
        neighborhood = matrix(1, n_actor, n_actor) - diag(n_actor),
        overlap = matrix(1, n_actor, n_actor) - diag(n_actor),
        directed = TRUE,
        terms = "spillover_xy_scaled_local",
        n_actor = n_actor,
        data_list = list(matrix(0)),
        type_list = c(0),
        type_x = "normal",
        type_y = "normal",
        attr_x_scale = 1.0,
        attr_y_scale = 1.0
    )
    expect_equal(as.numeric(res_full), hand_global_full)

    # 4. Check via simulate_iglm
    formula_full <- data_obj_full ~ spillover_xy_scaled_local
    sampler <- iglm:::sampler.iglm(
        n_simulation = 1,
        n_burn_in = 0,
        init_empty = FALSE,
        sampler_x = iglm:::sampler.net.attr(n_proposals = 0),
        sampler_y = iglm:::sampler.net.attr(n_proposals = 0),
        sampler_z = iglm:::sampler.net.attr(n_proposals = 0)
    )
    res_sim_full <- simulate_iglm(formula = formula_full, coef = c(0), sampler = sampler, only_stats = TRUE)
    expect_equal(as.numeric(res_sim_full$stats[1]), hand_global_full)
})

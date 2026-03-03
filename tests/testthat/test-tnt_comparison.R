library(testthat)
library(iglm)

test_that("TNT and normal sampling lead to similar results", {
    skip_on_cran()

    # Set up a small network
    n_actor <- 30
    neighborhood <- matrix(1, n_actor, n_actor)
    diag(neighborhood) <- 0

    # Create iglm.data object
    xyz_obj <- iglm.data(neighborhood = neighborhood, directed = FALSE)

    # Define a simple model
    formula <- xyz_obj ~ edges(mode = "local")

    # Create samplers
    # TNT Sampler
    sampler_tnt <- sampler.iglm(
        n_simulation = 50,
        n_burn_in = 20,
        sampler_z = sampler.net.attr(n_proposals = 2000, tnt = TRUE, seed = 1)
    )

    # Normal Sampler
    sampler_normal <- sampler.iglm(
        n_simulation = 50,
        n_burn_in = 20,
        sampler_z = sampler.net.attr(n_proposals = 2000, tnt = FALSE, seed = 1)
    )

    coef <- c(edges = -2)

    # Run simulations
    sim_tnt <- simulate_iglm(formula = formula, coef = coef, sampler = sampler_tnt)
    sim_normal <- simulate_iglm(formula = formula, coef = coef, sampler = sampler_normal)

    # Compare statistics
    stats_tnt <- colMeans(sim_tnt$stats)
    stats_normal <- colMeans(sim_normal$stats)

    cat("\nMean stats (TNT): ", stats_tnt, "\n")
    cat("Mean stats (Normal): ", stats_normal, "\n")

    # Check if they are reasonably close
    expect_true(abs(stats_tnt - stats_normal) < 1)
})

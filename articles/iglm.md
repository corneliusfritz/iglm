# An Introduction to Estimating Joint Probability Models with \`iglm\`

## Overview

This vignette provides an introduction to the `iglm` package, which is
designed for estimating joint probability models that incorporate
network structures. The package allows users to analyze how individual
attributes and network connections jointly influence outcomes of
interest.

## Basic Usage

To use the `iglm` package, you first need to load it into your R session

``` r
library(iglm)
```

Next, you can create a `iglm` object by specifying the network structure
and the attributes of interest. Here is a simple example:

``` r
n_actor <- 100

attribute_info <- rnorm(n_actor)
attribute_cov <- diag(attribute_info)
edge_cov <- outer(attribute_info, attribute_info, FUN = function(x, y) {
  abs(x - y)
})
set.seed(123)

alpha <- 0.3
block <- matrix(nrow = 50, ncol = 50, data = 1)
neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actor / 50, block, simplify = FALSE)))

overlapping_degree <- 0.5
neighborhood <- matrix(nrow = n_actor, ncol = n_actor, data = 0)
block <- matrix(nrow = 5, ncol = 5, data = 0)
size_neighborhood <- 5
size_overlap <- ceiling(size_neighborhood * overlapping_degree)

end <- floor((n_actor - size_neighborhood) / size_overlap)
for (i in 0:end) {
  neighborhood[(1 + size_overlap * i):(size_neighborhood + size_overlap * i), (1 + size_overlap * i):(size_neighborhood + size_overlap * i)] <- 1
}
neighborhood[(n_actor - size_neighborhood + 1):(n_actor), (n_actor - size_neighborhood + 1):(n_actor)] <- 1

type_x <- "binomial"
type_y <- "binomial"
formula_beg <- as.formula("xyz_obj ~ 1 ")
formula_model <- as.formula("xyz_object ~ 1 ")

object <- iglm.data(neighborhood = neighborhood, directed = F, type_x = type_x, type_y = type_y)
```

## Model Specification

You can specify a model formula that includes various network statistics
and attribute effects. For example:

``` r
formula <- object ~ edges + attribute_y + attribute_x + degrees
```

To fully define the model, you need to set up a sampler for the MCMC
estimation and set all necessary parameters:

``` r
# Parameters of edges(mode = "local"), attribute_y, and attribute_x
gt_coef <- c(3, -1, -1)
# Parameters for popularity effect
gt_coef_pop <- c(rnorm(n = n_actor, -2, 1))
# Define the sampler
sampler_tmp <- sampler.iglm(
  n_burn_in = 100, n_simulation = 10,
  sampler_x = sampler.net.attr(n_proposals = n_actor * 10, seed = 13),
  sampler_y = sampler.net.attr(n_proposals = n_actor * 10, seed = 32),
  sampler_z = sampler.net.attr(n_proposals = sum(neighborhood > 0) * 10, seed = 134),
  init_empty = F
)

model_tmp_new <- iglm(
  formula = formula,
  coef = gt_coef, coef_degrees = gt_coef_pop, sampler = sampler_tmp,
  control = control.iglm(accelerated = F, max_it = 200, display_progress = F)
)
```

## Model Simulation

Once you have specified a model, you can simulate new data based on the
fitted parameters:

``` r
# Simulate new networks
model_tmp_new$simulate()
# Get the samples
tmp <- model_tmp_new$get_samples()
```

## Model Estimation

You can estimate the model parameters using the `estimate` method:

``` r
# First set the first simulated network as the target for estimation
model_tmp_new$set_target(tmp[[1]])
model_tmp_new$estimate()
model_tmp_new$iglm.data$degree_distribution(plot = TRUE)
```

![](iglm_files/figure-html/unnamed-chunk-7-1.png)

## Model Assessment

After estimation, you can assess the model fit using various
diagnostics:

``` r
model_tmp_new$assess(formula = ~ degree_distribution +
  geodesic_distances_distribution + edgewise_shared_partner_distribution + mcmc_diagnostics)
```

![](iglm_files/figure-html/unnamed-chunk-8-1.png)![](iglm_files/figure-html/unnamed-chunk-8-2.png)![](iglm_files/figure-html/unnamed-chunk-8-3.png)![](iglm_files/figure-html/unnamed-chunk-8-4.png)![](iglm_files/figure-html/unnamed-chunk-8-5.png)![](iglm_files/figure-html/unnamed-chunk-8-6.png)![](iglm_files/figure-html/unnamed-chunk-8-7.png)![](iglm_files/figure-html/unnamed-chunk-8-8.png)![](iglm_files/figure-html/unnamed-chunk-8-9.png)![](iglm_files/figure-html/unnamed-chunk-8-10.png)![](iglm_files/figure-html/unnamed-chunk-8-11.png)![](iglm_files/figure-html/unnamed-chunk-8-12.png)![](iglm_files/figure-html/unnamed-chunk-8-13.png)![](iglm_files/figure-html/unnamed-chunk-8-14.png)![](iglm_files/figure-html/unnamed-chunk-8-15.png)![](iglm_files/figure-html/unnamed-chunk-8-16.png)![](iglm_files/figure-html/unnamed-chunk-8-17.png)![](iglm_files/figure-html/unnamed-chunk-8-18.png)![](iglm_files/figure-html/unnamed-chunk-8-19.png)![](iglm_files/figure-html/unnamed-chunk-8-20.png)![](iglm_files/figure-html/unnamed-chunk-8-21.png)![](iglm_files/figure-html/unnamed-chunk-8-22.png)![](iglm_files/figure-html/unnamed-chunk-8-23.png)![](iglm_files/figure-html/unnamed-chunk-8-24.png)![](iglm_files/figure-html/unnamed-chunk-8-25.png)![](iglm_files/figure-html/unnamed-chunk-8-26.png)![](iglm_files/figure-html/unnamed-chunk-8-27.png)![](iglm_files/figure-html/unnamed-chunk-8-28.png)![](iglm_files/figure-html/unnamed-chunk-8-29.png)![](iglm_files/figure-html/unnamed-chunk-8-30.png)![](iglm_files/figure-html/unnamed-chunk-8-31.png)![](iglm_files/figure-html/unnamed-chunk-8-32.png)![](iglm_files/figure-html/unnamed-chunk-8-33.png)![](iglm_files/figure-html/unnamed-chunk-8-34.png)![](iglm_files/figure-html/unnamed-chunk-8-35.png)![](iglm_files/figure-html/unnamed-chunk-8-36.png)![](iglm_files/figure-html/unnamed-chunk-8-37.png)![](iglm_files/figure-html/unnamed-chunk-8-38.png)![](iglm_files/figure-html/unnamed-chunk-8-39.png)

``` r
model_tmp_new$results$plot(model_assessment = T)
```

![](iglm_files/figure-html/unnamed-chunk-8-40.png)![](iglm_files/figure-html/unnamed-chunk-8-41.png)![](iglm_files/figure-html/unnamed-chunk-8-42.png)![](iglm_files/figure-html/unnamed-chunk-8-43.png)![](iglm_files/figure-html/unnamed-chunk-8-44.png)![](iglm_files/figure-html/unnamed-chunk-8-45.png)

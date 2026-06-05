
# iglm <img src="man/figures/logo.png" align="right" height="139" alt="iglm logo" />

# An Introduction to Estimating Joint Probability Models with `iglm`

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
#> Loading required package: RcppProgress
```

Next, you can create a `iglm` object by specifying the network structure
and the attributes of interest. Here is a simple example:

``` r
n_actors = 100

attribute_info = rnorm(n_actors)
attribute_cov = diag(attribute_info)
edge_cov = outer(attribute_info, attribute_info, FUN = function(x,y){abs(x-y)})
set.seed(123)

alpha = 0.3
block <- matrix(nrow = 50, ncol = 50, data = 1)
neighborhood <- as.matrix(Matrix::bdiag(replicate(n_actors/50, block, simplify=FALSE)))

overlapping_degree = 0.5
neighborhood = matrix(nrow = n_actors, ncol = n_actors, data = 0)
block <- matrix(nrow = 5, ncol = 5, data = 0)
size_neighborhood <- 5
size_overlap <-  ceiling(size_neighborhood*overlapping_degree)

end <- floor((n_actors-size_neighborhood)/size_overlap)
for(i in 0:end){
  neighborhood[(1+size_overlap*i):(size_neighborhood+size_overlap*i), (1+size_overlap*i):(size_neighborhood+size_overlap*i)] = 1
}
neighborhood[(n_actors-size_neighborhood+1):(n_actors), (n_actors-size_neighborhood+1):(n_actors)] = 1

type_x <- "binomial"
type_y <- "binomial"
object = iglm.data(neighborhood = neighborhood, directed = F, type_x = type_x, type_y = type_y, n_actor = n_actors)
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
gt_coef = c(3,-1,-1)
# Parameters for degree effect
gt_coef_degrees =  c(rnorm(n = n_actors, -2, 1))
# Define the sampler
sampler_tmp = sampler.iglm(n_burn_in = 100, n_simulation = 10,
                               sampler_x = sampler.net.attr(n_proposals =  n_actors*10),
                               sampler_y = sampler.net.attr(n_proposals =  n_actors*10),
                               sampler_z = sampler.net.attr(n_proposals = sum(neighborhood>0)*10),
                               init_empty = F)

model_tmp_new <- iglm(formula = formula,
                           coef = gt_coef,  coef_degrees = gt_coef_degrees, sampler = sampler_tmp, 
                          control = control.iglm(accelerated = F,max_it = 200, display_progress = F))
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

<img src="man/figures/model_est-1.png" alt="" width="100%" />

## Model Assessment

After estimation, you can assess the model fit using various
diagnostics:

``` r
model_tmp_new$assess(formula = ~  degree_distribution + 
                                 geodesic_distances_distribution + edgewise_shared_partner_distribution + mcmc_diagnostics)
```

<img src="man/figures/model_assess-1.png" alt="" width="100%" /><img src="man/figures/model_assess-2.png" alt="" width="100%" /><img src="man/figures/model_assess-3.png" alt="" width="100%" /><img src="man/figures/model_assess-4.png" alt="" width="100%" /><img src="man/figures/model_assess-5.png" alt="" width="100%" /><img src="man/figures/model_assess-6.png" alt="" width="100%" />

``` r
model_tmp_new$results$plot(model_assessment = T)
```

<img src="man/figures/model_assess-7.png" alt="" width="100%" /><img src="man/figures/model_assess-8.png" alt="" width="100%" /><img src="man/figures/model_assess-9.png" alt="" width="100%" /><img src="man/figures/model_assess-10.png" alt="" width="100%" /><img src="man/figures/model_assess-11.png" alt="" width="100%" /><img src="man/figures/model_assess-12.png" alt="" width="100%" />

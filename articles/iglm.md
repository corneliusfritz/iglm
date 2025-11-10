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
n_actors =100

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
size_neighborhood = 5
size_overlap = ceiling(size_neighborhood*overlapping_degree)

end = floor((n_actors-size_neighborhood)/size_overlap)
for(i in 0:end){
  neighborhood[(1+size_overlap*i):(size_neighborhood+size_overlap*i), (1+size_overlap*i):(size_neighborhood+size_overlap*i)] = 1
}
neighborhood[(n_actors-size_neighborhood+1):(n_actors), (n_actors-size_neighborhood+1):(n_actors)] = 1

type_x <- "binomial"
type_y <- "binomial"
formula_beg = as.formula("xyz_obj ~ 1 ")
formula_model = as.formula("xyz_object ~ 1 ")

object = netplus(neighborhood = neighborhood, directed = F, type_x = type_x, type_y = type_y)
```

## Model Specification

You can specify a model formula that includes various network statistics
and attribute effects. For example:

``` r
formula <- object ~ edges(mode = "local") + attribute_y + attribute_x + popularity
```

To fully define the model, you need to set up a sampler for the MCMC
estimation and set all necessary parameters:

``` r
# Parameters of edges(mode = "local"), attribute_y, and attribute_x
gt_coef = c(3,-1,-1)
# Parameters for popularity effect
gt_coef_pop =  c(rnorm(n = n_actors, -2, 1))
# Define the sampler
sampler_tmp = sampler.iglm(n_burn_in = 100, n_simulation = 10,
                               sampler.x = sampler.net_attr(n_proposals =  n_actors*10,seed = 13),
                               sampler.y = sampler.net_attr(n_proposals =  n_actors*10, seed = 32),
                               sampler.z = sampler.net_attr(n_proposals = sum(neighborhood>0)*10, seed = 134),
                               init_empty = F)

model_tmp_new <- iglm(formula = formula,
                           coef = gt_coef,  coef_popularity = gt_coef_pop, sampler = sampler_tmp, 
                          control = control.iglm(accelerated = F,max_it = 200, display_progress = T, var = T))
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
#> Target netplus object has been set successfully.
model_tmp_new$estimate()
#> Starting with the preprocessing
#> Starting with the estimation
#> Iteration = 1Iteration = 2Iteration = 3Iteration = 4Iteration = 5Iteration = 6Iteration = 7Iteration = 8Iteration = 9Iteration = 10Iteration = 11Iteration = 12Iteration = 13Iteration = 14Iteration = 15Iteration = 16Iteration = 17Iteration = 18Iteration = 19Iteration = 20Iteration = 21Iteration = 22Iteration = 23Iteration = 24Iteration = 25Iteration = 26Iteration = 27Iteration = 28Iteration = 29Iteration = 30Iteration = 31Iteration = 32Iteration = 33Iteration = 34Iteration = 35Iteration = 36Iteration = 37Iteration = 38Iteration = 39Iteration = 40Iteration = 41Iteration = 42Iteration = 43Iteration = 44Iteration = 45Iteration = 46Iteration = 47Iteration = 48Iteration = 49Iteration = 50Iteration = 51Iteration = 52Iteration = 53Iteration = 54Iteration = 55Iteration = 56Iteration = 57Iteration = 58Iteration = 59Iteration = 60Iteration = 61Iteration = 62Iteration = 63Iteration = 64Iteration = 65Iteration = 66Iteration = 67Iteration = 68Iteration = 69Iteration = 70Iteration = 71Iteration = 72Iteration = 73Iteration = 74Iteration = 75Iteration = 76Iteration = 77Iteration = 78Iteration = 79Iteration = 80Iteration = 81Iteration = 82Iteration = 83Iteration = 84Iteration = 85Iteration = 86Iteration = 87Iteration = 88Iteration = 89Iteration = 90Iteration = 91Iteration = 92Iteration = 93Iteration = 94Iteration = 95Iteration = 96Iteration = 97Iteration = 98Iteration = 99Iteration = 100Iteration = 101Iteration = 102Iteration = 103Iteration = 104Iteration = 105Iteration = 106Iteration = 107Iteration = 108Iteration = 109Iteration = 110Iteration = 111Iteration = 112Iteration = 113Iteration = 114Iteration = 115Iteration = 116Iteration = 117Iteration = 118Iteration = 119Iteration = 120Iteration = 121Iteration = 122Iteration = 123Iteration = 124Iteration = 125Iteration = 126Iteration = 127Iteration = 128Iteration = 129Iteration = 130Iteration = 131Iteration = 132Iteration = 133Iteration = 134Iteration = 135Iteration = 136Iteration = 137Iteration = 138Iteration = 139Iteration = 140Iteration = 141Iteration = 142Iteration = 143Iteration = 144Iteration = 145Iteration = 146Iteration = 147Iteration = 148Iteration = 149Iteration = 150Iteration = 151Iteration = 152Iteration = 153Iteration = 154Iteration = 155Iteration = 156Iteration = 157Iteration = 158Iteration = 159Iteration = 160Iteration = 161Iteration = 162Iteration = 163Iteration = 164Iteration = 165Iteration = 166Iteration = 167Iteration = 168Iteration = 169Iteration = 170Iteration = 171Iteration = 172Iteration = 173Iteration = 174Iteration = 175Iteration = 176Iteration = 177Iteration = 178Iteration = 179Iteration = 180Iteration = 181Iteration = 182Iteration = 183Iteration = 184Iteration = 185Iteration = 186Iteration = 187Iteration = 188Iteration = 189Iteration = 190Iteration = 191Iteration = 192Iteration = 193Iteration = 194Iteration = 195Iteration = 196Iteration = 197Iteration = 198Iteration = 199Iteration = 200Done with the estimation
#> Starting with samples to estimate uncertainty 
#> Sample: 1Sample: 2Sample: 3Sample: 4Sample: 5Sample: 6Sample: 7Sample: 8Sample: 9Sample: 10Sample: 11Sample: 12Sample: 13Sample: 14Sample: 15Sample: 16Sample: 17Sample: 18Sample: 19Sample: 20Sample: 21Sample: 22Sample: 23Sample: 24Sample: 25Sample: 26Sample: 27Sample: 28Sample: 29Sample: 30Sample: 31Sample: 32Sample: 33Sample: 34Sample: 35Sample: 36Sample: 37Sample: 38Sample: 39Sample: 40Sample: 41Sample: 42Sample: 43Sample: 44Sample: 45Sample: 46Sample: 47Sample: 48Sample: 49Sample: 50Sample: 51Sample: 52Sample: 53Sample: 54Sample: 55Sample: 56Sample: 57Sample: 58Sample: 59Sample: 60Sample: 61Sample: 62Sample: 63Sample: 64Sample: 65Sample: 66Sample: 67Sample: 68Sample: 69Sample: 70Sample: 71Sample: 72Sample: 73Sample: 74Sample: 75Sample: 76Sample: 77Sample: 78Sample: 79Sample: 80Sample: 81Sample: 82Sample: 83Sample: 84Sample: 85Sample: 86Sample: 87Sample: 88Sample: 89Sample: 90Sample: 91Sample: 92Sample: 93Sample: 94Sample: 95Sample: 96Sample: 97Sample: 98Sample: 99Sample: 100Sample: 101Sample: 102Sample: 103Sample: 104Sample: 105Sample: 106Sample: 107Sample: 108Sample: 109Sample: 110
#> Results: 
#> 
#>                       Estimate Std. Error
#> edges(mode = 'local')    2.806      0.510
#> attribute_y             -1.099      0.171
#> attribute_x             -1.046      0.208
```

## Model Assessment

After estimation, you can assess the model fit using various
diagnostics:

``` r
model_tmp_new$model_assessment(formula = ~  degree_distribution + 
                                 geodesic_distances_distribution + edgewise_shared_partner_distribution)
model_tmp_new$results$plot(model_assessment = T)
```

![](iglm_files/figure-html/unnamed-chunk-8-1.png)![](iglm_files/figure-html/unnamed-chunk-8-2.png)![](iglm_files/figure-html/unnamed-chunk-8-3.png)

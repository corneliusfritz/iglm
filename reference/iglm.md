# Construct a iglm Model Specification Object

The `iglm` package implements a comprehensive regression framework
introduced in Fritz et al. (2025) for studying relationships among
attributes \\(X, Y)\\ under network interference \\(Z)\\. It is based on
a joint probability model for dependent outcomes (\\Y\\) and network
connections \\(Z)\\, conditional on a fixed set of predictors (X). This
approach generalizes standard Generalized Linear Models (GLMs) to
settings where the responses and connections of units are
interdependent. The framework is designed to be interpretable by
representing conditional distributions as GLMs, scalable to large
networks via pseudo-likelihood and convex optimization, and provides
insight into outcome-connection dependencies (i.e., spillover effects)
that are missed by conditional models.

The joint probability density is specified as an exponential-family
model of the form: \$\$f\_{\theta}(y,z,x) \propto \Big\[\prod\_{i=1}^{N}
a_y(y_i) \exp(\theta_g^T g_i(x_i, y_i^\*)) \Big\] \times \Big\[\prod\_{i
\ne j} a_z(z\_{i,j}) \exp(\theta_h^T h\_{i,j}(x, y_i^\*, y_j^\*, z))
\Big\],\$\$ which is defined by two distinct sets of user-specified
features:

- **\\g_i(x,y,z)\\**: A vector of actor-level functions (or "g-terms")
  that describe the relationship between an individual actor \\i\\'s
  predictors (\\x_i\\) and their own response (\\y_i\\).

- **\\h\_{i,j}(x,y,z)\\**: A vector of pair-level functions (or
  "h-terms") that specify how the connections (\\z\\) and responses
  (\\y_i, y_j\\) of a pair of units \\\\i,j\\\\ depend on each other and
  the wider network structure.

This separation allows the model to simultaneously capture
individual-level behavior (via \\g_i\\) and dyadic, network-based
dependencies (via \\h\_{i,j}\\), including local dependence limited to
overlapping neighborhoods (see, Fritz et al., 2025). This help page
documents the various statistics available in 'iglm', corresponding to
the \\g_i\\ (attribute-level) and \\h\_{i,j}\\ (pair-level) components
of the joint model. This is a user-facing constructor for creating a
`iglm_object`. This `R6` object encompasses the complete model
specification, linking the formula, data ([`netplus`](netplus.md)
object), initial coefficients, MCMC sampler settings, and estimation
controls. It serves as the primary input for subsequent methods like
`$estimate()` and `$simulate()`.

## Usage

``` r
iglm(
  formula = NULL,
  coef = NULL,
  coef_popularity = NULL,
  sampler = NULL,
  control = NULL,
  file = NULL
)
```

## Arguments

- formula:

  A model \`formula\` object. The left-hand side should be the name of a
  \`netplus\` object available in the calling environment. See
  [`model_terms`](model_terms.md) for details on specifying the
  right-hand side terms.

- coef:

  Optional numeric vector of initial coefficients for the structural
  (non-popularity) terms in \`formula\`. If \`NULL\`, coefficients are
  initialized to zero. Length must match the number of terms.

- coef_popularity:

  Optional numeric vector specifying the initial popularity
  coefficients. Required if \`formula\` includes popularity terms,
  otherwise should be \`NULL\`. Length must match \`n_actor\` (for
  undirected) or \`2 \* n_actor\` (for directed).

- sampler:

  An object of class [`sampler_iglm`](sampler_iglm_generator.md),
  controlling the MCMC sampling scheme. If \`NULL\`, default sampler
  settings will be used.

- control:

  An object of class [`control.iglm`](control.iglm.md), specifying
  parameters for the estimation algorithm. If \`NULL\`, default control
  settings will be used.

- file:

  Optional character string specifying a file path to load a previously
  saved `iglm_object` from disk (in RDS format). If provided, other
  arguments are ignored and the object is loaded from the file.

## Value

An object of class `iglm_object`.

## References

Fritz, C., Schweinberger, M., Bhadra, S., and D.R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

Schweinberger, M. and M.S. Handcock (2015). Local Dependence in Random
Graph Models: Characterization, Properties, and Statistical Inference.
Journal of the Royal Statistical Society, Series B (Statistical
Methodology), 7, 647-676.

Schweinberger, M. and J.R. Stewart (2020). Concentration and Consistency
Results for Canonical and Curved Exponential-Family Models of Random
Graphs. The Annals of Statistics, 48, 374-396.

Stewart, J.R. and M. Schweinberger (2025). Pseudo-Likelihood-Based
M-Estimation of Random Graphs with Dependent Edges and Parameter Vectors
of Increasing Dimension. The Annals of Statistics, to appear.

## Examples

``` r
# Example usage:
library(iglm)
# Create a netplus data object (example)
n_actors <- 50
neighborhood <- matrix(1, nrow = n_actors, ncol = n_actors)
xyz_obj <- netplus(neighborhood = neighborhood, directed = FALSE,
                   type_x = "binomial", type_y = "binomial")
# Define ground truth coefficients
gt_coef <- c("edges_local" = 3, "attribute_y" = -1, "attribute_x" = -1)
gt_coef_pop <- rnorm(n = n_actors, -2, 1)
# Define MCMC sampler
sampler_new <- sampler.iglm(n_burn_in = 100, n_simulation = 10,
                               sampler.x = sampler.net_attr(n_proposals = n_actors * 10, seed = 13),
                               sampler.y = sampler.net_attr(n_proposals = n_actors * 10, seed = 32),
                               sampler.z = sampler.net_attr(n_proposals = sum(neighborhood > 0
                               ) * 10, seed = 134),
                               init_empty = FALSE)
# Create iglm model specification
model_tmp_new <- iglm(formula = xyz_obj ~ edges(mode = "local") +
                          attribute_y + attribute_x + popularity,
                          coef = gt_coef,
                          coef_popularity = gt_coef_pop,
                          sampler = sampler_new,
                          control = control.iglm(accelerated = FALSE,
                          max_it = 200, display_progress = TRUE, var = TRUE))
# Simulate from the model
model_tmp_new$simulate()
model_tmp_new$set_target(model_tmp_new$get_samples()[[1]])
#> Target netplus object has been set successfully.

# Estimate model parameters
model_tmp_new$estimate()
#> Starting with the preprocessing
#> Starting with the estimation
#> Iteration = 1Iteration = 2Iteration = 3Iteration = 4Iteration = 5Iteration = 6Iteration = 7Iteration = 8Iteration = 9Iteration = 10Iteration = 11Iteration = 12Iteration = 13Iteration = 14Iteration = 15Iteration = 16Iteration = 17Iteration = 18Iteration = 19Iteration = 20Iteration = 21Iteration = 22Done with the estimation
#> Starting with samples to estimate uncertainty 
#> Sample: 1Sample: 2Sample: 3Sample: 4Sample: 5Sample: 6Sample: 7Sample: 8Sample: 9Sample: 10Sample: 11Sample: 12Sample: 13Sample: 14Sample: 15Sample: 16Sample: 17Sample: 18Sample: 19Sample: 20Sample: 21Sample: 22Sample: 23Sample: 24Sample: 25Sample: 26Sample: 27Sample: 28Sample: 29Sample: 30Sample: 31Sample: 32Sample: 33Sample: 34Sample: 35Sample: 36Sample: 37Sample: 38Sample: 39Sample: 40Sample: 41Sample: 42Sample: 43Sample: 44Sample: 45Sample: 46Sample: 47Sample: 48Sample: 49Sample: 50Sample: 51Sample: 52Sample: 53Sample: 54Sample: 55Sample: 56Sample: 57Sample: 58Sample: 59Sample: 60Sample: 61Sample: 62Sample: 63Sample: 64Sample: 65Sample: 66Sample: 67Sample: 68Sample: 69Sample: 70Sample: 71Sample: 72Sample: 73Sample: 74Sample: 75Sample: 76Sample: 77Sample: 78Sample: 79Sample: 80Sample: 81Sample: 82Sample: 83Sample: 84Sample: 85Sample: 86Sample: 87Sample: 88Sample: 89Sample: 90Sample: 91Sample: 92Sample: 93Sample: 94Sample: 95Sample: 96Sample: 97Sample: 98Sample: 99Sample: 100Sample: 101Sample: 102Sample: 103Sample: 104Sample: 105Sample: 106Sample: 107Sample: 108Sample: 109Sample: 110
#> Results: 
#> 
#>                       Estimate Std. Error
#> edges(mode = 'local')    2.989      0.073
#> attribute_y             -0.847      0.340
#> attribute_x             -1.386      0.324

# Model Assessment
model_tmp_new$model_assessment(formula = ~  degree_distribution )
# model_tmp_new$results$plot(model_assessment = TRUE)
                                                   
```

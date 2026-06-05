# Constructor for a iglm Sampler

Creates an object of class \`sampler.iglm\` (and \`R6\`) which holds all
parameters controlling the MCMC sampling process for \`iglm\` models.
This includes global settings like the number of simulations and
burn-in, as well as references to specific samplers for the network
(\`z\`) and attribute (\`x\`, \`y\`) components.

This function provides a convenient way to specify these settings before
passing them to the \`iglm\` constructor or simulation functions.

## Usage

``` r
sampler.iglm(
  sampler_x = NULL,
  sampler_y = NULL,
  sampler_z = NULL,
  n_simulation = 100,
  n_burn_in = 10,
  init_empty = TRUE,
  seed = NA,
  cluster = NULL,
  file = NULL
)
```

## Arguments

- sampler_x:

  An object of class \`sampler.net.attr\` (created by
  \`sampler.net.attr()\`) specifying how to sample the \`x_attribute\`.
  If \`NULL\` (default), default \`sampler.net.attr()\` settings are
  used.

- sampler_y:

  An object of class \`sampler.net.attr\` specifying how to sample the
  \`y_attribute\`. If \`NULL\` (default), default settings are used.

- sampler_z:

  An object of class \`sampler.net.attr\` specifying how to sample the
  \`z_network\` ties \*within\* the defined neighborhood/overlap region.
  If \`NULL\` (default), default settings are used.

- n_simulation:

  (integer) The number of independent samples to generate after the
  burn-in period. Default: 100. Must be non-negative.

- n_burn_in:

  (integer) The number of MCMC iterations to discard at the start for
  burn-in. Default: 10. Must be non-negative.

- init_empty:

  (logical) If \`TRUE\` (default), initialize the MCMC chain from an
  empty state.

- seed:

  (integer or \`NA\`) A single integer seed set once before sampling
  begins to ensure reproducibility. If \`NA\` (default), a random seed
  is generated automatically.

- cluster:

  A parallel cluster object (e.g., from \`parallel::makeCluster()\`) for
  parallel simulations. If \`NULL\` (default), simulations run
  sequentially.

- file:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

## Value

An object of class \`sampler.iglm\` (and \`R6\`).

## See also

\`sampler.net.attr\`, \`iglm\`, \`control.iglm\`

## Examples

``` r
n_actor <- 50
sampler_new <- sampler.iglm(
  n_burn_in = 100, n_simulation = 10,
  seed = 42,
  sampler_x = sampler.net.attr(n_proposals = n_actor * 10),
  sampler_y = sampler.net.attr(n_proposals = n_actor * 10),
  sampler_z = sampler.net.attr(n_proposals = n_actor^2, tnt = TRUE),
  init_empty = FALSE
)
sampler_new
#> Sampler settings
#> ------------------------------------------------------------
#> Core parameters
#>   n_simulation :10
#>   n_burn_in    :100
#>   init_empty   :FALSE
#>   seed         :42
#> 
#> Sub-samplers
#>   sampler_x:
#>     Number of proposals : 500
#>     TNT sampling        : TRUE
#>   sampler_y:
#>     Number of proposals : 500
#>     TNT sampling        : TRUE
#>   sampler_z:
#>     Number of proposals : 2500
#>     TNT sampling        : TRUE
sampler_new$seed
#> [1] 42
sampler_new$set_n_simulation(100)
sampler_new$n_simulation
#> [1] 100
```

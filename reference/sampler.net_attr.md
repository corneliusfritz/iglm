# Constructor for Single Component Sampler Settings

Creates an object of class \`sampler_net_attr\` (and \`R6\`). This
object specifies the MCMC sampling parameters for a single component
(like an attribute vector or a network structure) within the larger
\`iglm\` simulation framework. It is typically used as input when
creating a \`sampler_iglm\` object.

## Usage

``` r
sampler.net_attr(n_proposals = 10000, seed = NA, file = NULL)
```

## Arguments

- n_proposals:

  (integer) The number of MCMC proposals (iterations) to perform for
  this specific component during each sampling update. Default: 10000.

- seed:

  (integer or \`NA\`) An integer seed for the random number generator to
  ensure reproducibility for this component's sampling process. If
  \`NA\` (default), a random seed will be generated automatically.

- file:

  (character or \`NULL\`) If provided, loads the sampler state from '
  the specified .rds file instead of initializing from parameters.

## Value

An object of class \`sampler_net_attr\` (and \`R6\`).

## See also

\`sampler.iglm\`

## Examples

``` r
# Default settings
sampler_comp_default <- sampler.net_attr()
sampler_comp_default
#>   Number of proposals : 10000
#>   Random seed         : 167759

# Custom settings
sampler_comp_custom <- sampler.net_attr(n_proposals = 50000, seed = 123)
sampler_comp_custom
#>   Number of proposals : 50000
#>   Random seed         : 123
```

# Constructor for Single Component Sampler Settings

Creates an object of class \`sampler_net_attr\` (and \`R6\`). Specifies
MCMC sampling parameters for one component (attribute or network) within
the \`iglm\` simulation framework. Used as input to \`sampler.iglm()\`.

## Usage

``` r
sampler.net.attr(n_proposals = 10000, file = NULL, tnt = TRUE)
```

## Arguments

- n_proposals:

  (integer) Number of MCMC proposals per sampling update. Default:
  10000.

- file:

  (character or \`NULL\`) If provided, loads state from an .rds file.

- tnt:

  (logical) If \`TRUE\` (default), use Tie-No-Tie sampling.

## Value

An object of class \`sampler_net_attr\` (and \`R6\`).

## See also

\`sampler.iglm\`

## Examples

``` r
sampler_comp_default <- sampler.net.attr()
sampler_comp_default
#>   Number of proposals : 10000
#>   TNT sampling        : TRUE

sampler_comp_custom <- sampler.net.attr(n_proposals = 50000, tnt = FALSE)
sampler_comp_custom
#>   Number of proposals : 50000
#>   TNT sampling        : FALSE
```

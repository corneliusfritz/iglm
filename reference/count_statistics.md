# Compute Statistics

Computes statistics.

## Usage

``` r
count_statistics(formula)
```

## Arguments

- formula:

  A model \`formula\` object. The left-hand side should be the name of a
  [`iglm.data`](iglm.data.md) object available in the calling
  environment. See [`model_terms`](model_terms.md) for details on
  specifying the right-hand side terms.

## Value

A named numeric vector. Each element corresponds to a term in the
\`formula\`, and its value is the calculated observed feature for that
term based on the data in the [`iglm.data`](iglm.data.md) object. The
names of the vector match the coefficient names derived from the formula
terms.

## Examples

``` r
# Create a iglm.data object
n_actors = 10
neighborhood = matrix(1, nrow = n_actors, ncol = n_actors)
type_x <- "binomial"
type_y <- "binomial"
x_attr_data <- rbinom(n_actors, 1, 0.5)
y_attr_data <- rbinom(n_actors, 1, 0.5)
z_net_data  <- matrix(0, nrow = n_actors, ncol = n_actors) 
object = iglm.data(z_network = z_net_data, x_attribute = x_attr_data,
 y_attribute = y_attr_data, neighborhood = neighborhood,
 directed = FALSE,type_x = type_x,type_y = type_y)
count_statistics(object ~ edges(mode = "local") + attribute_y + attribute_x)
#> edges(mode = 'local')           attribute_y           attribute_x 
#>                     0                     2                     4 
```

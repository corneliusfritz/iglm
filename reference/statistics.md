# Compute Statistics

Computes observed sufficient statistics for the terms specified in a
formula.

## Usage

``` r
statistics(formula, canonical_names = FALSE)
```

## Arguments

- formula:

  A model \`formula\` object. The left-hand side should be the name of a
  [`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
  object available in the calling environment. Alternatively, the
  left-hand side can be a `iglm.data.list` object to compute statistics
  for multiple
  [`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
  objects at once (e.g., the output of
  [`simulate_iglm`](https://corneliusfritz.github.io/iglm/reference/simulate_iglm.md)).
  See
  [`iglm.terms`](https://corneliusfritz.github.io/iglm/reference/iglm.terms.md)
  for details on specifying the right-hand side terms.

- canonical_names:

  (logical) If \`TRUE\`, returns canonical term names without label
  substitution. Default is \`FALSE\`.

## Value

A named numeric vector (or matrix if a list of data objects is
provided). Each element corresponds to a term in the \`formula\`, and
its value is the calculated observed feature for that term based on the
data in the
[`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
object. The names of the vector match the (optionally labeled) term
names derived from the formula terms.

## See also

[`iglm`](https://corneliusfritz.github.io/iglm/reference/iglm.md),
[`iglm.terms`](https://corneliusfritz.github.io/iglm/reference/iglm.terms.md),
[`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)

## Examples

``` r
# Create an iglm.data object
n_units <- 10
neighborhood <- matrix(1, nrow = n_units, ncol = n_units)
family_x <- "binomial"
family_y <- "binomial"
x_attr_data <- rbinom(n_units, 1, 0.5)
y_attr_data <- rbinom(n_units, 1, 0.5)
z_net_data <- matrix(0, nrow = n_units, ncol = n_units)
object <- iglm.data(
  z_network = z_net_data, x_attribute = x_attr_data,
  y_attribute = y_attr_data, neighborhood = neighborhood,
  directed = FALSE, family_x = family_x, family_y = family_y,
  n_units = n_units,
  label_x = "age", label_y = "income"
)
statistics(object ~ edges(mode = "local") + attribute_y + attribute_x)
#> edges(mode = 'local')      attribute_income         attribute_age 
#>                     0                     4                     4 
statistics(object ~ edges(mode = "local") + attribute_y + attribute_x, canonical_names = TRUE)
#> edges(mode = 'local')           attribute_y           attribute_x 
#>                     0                     4                     4 
```

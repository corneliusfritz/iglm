# Constructor for the iglm.data R6 object

Creates an \`iglm.data\` object, which stores network and attribute
data. This function acts as a user-friendly interface to the
\`iglm.data\` R6 class generator. It handles data input, infers
parameters like the number of units (\`n_units\`) and network
directedness (\`directed\`) if not explicitly provided, processes
network data into a consistent edgelist format, calculates the overlap
relation based on an optional neighborhood definition, and performs
extensive validation of all inputs.

## Usage

``` r
iglm.data(
  x_attribute = NULL,
  y_attribute = NULL,
  z_network = NULL,
  neighborhood = NULL,
  directed = TRUE,
  n_units = NA,
  family_x = "binomial",
  family_y = "binomial",
  scale_x = 1,
  scale_y = 1,
  fix_x = FALSE,
  fix_z = FALSE,
  fix_z_alocal = FALSE,
  return_neighborhood = TRUE,
  file = NULL,
  label_x = "x",
  label_y = "y",
  label_z = "z",
  n_actor = NULL
)
```

## Arguments

- x_attribute:

  A numeric vector for the first unit-level attribute.

- y_attribute:

  A numeric vector for the second unit-level attribute.

- z_network:

  A matrix representing the network. Can be a 2-column edgelist or a
  square adjacency matrix.

- neighborhood:

  An optional matrix for the neighborhood representing local dependence.
  Can be a 2-column edgelist or a square adjacency matrix. A tie in
  \`neighborhood\` between unit i and j indicates that j is in the
  neighborhood of i, implying dependence between the respective units.

- directed:

  A logical value indicating if \`z_network\` is directed. If \`NA\`
  (default), directedness is inferred from the symmetry of
  \`z_network\`.

- n_units:

  An integer for the number of units in the system. If \`NA\` (default),
  \`n_units\` is inferred from the attributes or network matrices.

- family_x:

  Character string for the family of \`x_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- family_y:

  Character string for the family of \`y_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- scale_x:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- scale_y:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- fix_x:

  (logical) If \`TRUE\`, the 'x' predictor is held fixed during
  estimation/simulation (fixed design in regression). Default is
  \`FALSE\`.

- fix_z:

  (logical) If \`TRUE\`, the 'z' network is held fixed during
  estimation/simulation (fixed network design). Default is \`FALSE\`.
  Setting this to TRUE, allows practitioners to estimate autologistic
  unit attribute models, which were introduced in binary settings in
  Daraganova, G., & Robins, G. (2013).

- fix_z_alocal:

  (logical) If \`TRUE\`, edges outside the overlap region are fixed,
  else they are random (default).

- return_neighborhood:

  Logical. If \`TRUE\` (default) and \`neighborhood\` is \`NULL\`, a
  full neighborhood (all dyads) is generated implying global dependence.
  If \`FALSE\`, no neighborhood is set.

- file:

  (character) Optional file path to load a saved \`iglm.data\` object
  state.

- label_x:

  Character string for the label/name of \`x_attribute\`. Default is
  \`"x"\`.

- label_y:

  Character string for the label/name of \`y_attribute\`. Default is
  \`"y"\`.

- label_z:

  Character string for the label/name of \`z_network\`. Default is
  \`"z"\`.

- n_actor:

  (Deprecated) An integer for the number of units. Please use `n_units`
  instead.

## Value

An object of class \`iglm.data\` (and \`R6\`).

## References

Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

Daraganova, G., and Robins, G. (2013). Exponential random graph models
for social networks: Theory, methods and applications, 102-114.
Cambridge University Press.

## Examples

``` r
# \donttest{
data("state_twitter")
state_twitter$iglm.data$degree_distribution(prob = FALSE, plot = TRUE)


state_twitter$iglm.data$geodesic_distances_distribution(prob = FALSE, plot = TRUE)

state_twitter$iglm.data$mean_x()
#> [1] 0.3939394
state_twitter$iglm.data$mean_y()
#> [1] 0.4121212
# }

# Generate a small iglm data object either via adjacency matrix or edgelist
tmp_adjacency <- iglm.data(
  z_network = matrix(c(
    0, 1, 1, 0,
    1, 0, 0, 1,
    1, 0, 0, 1,
    0, 1, 1, 0
  ), nrow = 4, byrow = TRUE),
  directed = FALSE,
  n_units = 4,
  family_x = "binomial",
  family_y = "binomial"
)


tmp_edgelist <- iglm.data(
  z_network = tmp_adjacency$z_network,
  directed = FALSE,
  n_units = 4,
  family_x = "binomial",
  family_y = "binomial"
)

tmp_edgelist$mean_z()
#> [1] 0.6666667
tmp_adjacency$mean_z()
#> [1] 0.6666667
```

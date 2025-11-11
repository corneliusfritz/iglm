# Constructor for the iglm.data R6 object

Creates a \`iglm.data\` object, which stores network and attribute data.
This function acts as a user-friendly interface to the \`iglm.data\` R6
class generator. It handles data input, infers parameters like the
number of actors (\`n_actor\`) and network directedness (\`directed\`)
if not explicitly provided, processes network data into a consistent
edgelist format, calculates the overlap relation based on an optional
neighborhood definition, and performs extensive validation of all
inputs.

## Usage

``` r
iglm.data(
  x_attribute = NULL,
  y_attribute = NULL,
  z_network = NULL,
  neighborhood = NULL,
  directed = NA,
  n_actor = NA,
  type_x = "binomial",
  type_y = "binomial",
  scale_x = 1,
  scale_y = 1,
  return_neighborhood = TRUE,
  file = NULL
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
  \`neighborhood\` between actor i and j indicates that j is in the
  neighborhood of i, implying dependence between the respective actors.

- directed:

  A logical value indicating if \`z_network\` is directed. If \`NA\`
  (default), directedness is inferred from the symmetry of
  \`z_network\`.

- n_actor:

  An integer for the number of actors in the system. If \`NA\`
  (default), \`n_actor\` is inferred from the attributes or network
  matrices.

- type_x:

  Character string for the type of \`x_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- type_y:

  Character string for the type of \`y_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- scale_x:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- scale_y:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- return_neighborhood:

  Logical. If \`TRUE\` (default) and \`neighborhood\` is \`NULL\`, a
  full neighborhood (all dyads) is generated implying global dependence.
  If \`FALSE\`, no neighborhood is set.

- file:

  (character) Optional file path to load a saved \`iglm.data\` object
  state.

## Value

An object of class \`iglm.data\` (and \`R6\`).

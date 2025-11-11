# A R6 class to represent networks with unit-level attributes

The \`iglm_data\` class is a container for storing, validating, and
analyzing unit-level attributes (x_attribute, y_attribute) and
connections (z_network).

## Active bindings

- `x_attribute`:

  (\`numeric\`) Read-only. The vector for the first unit-level
  attribute.

- `y_attribute`:

  (\`numeric\`) Read-only. The vector for the second unit-level
  attribute.

- `z_network`:

  (\`matrix\`) Read-only. The primary network structure as a 2-column
  integer edgelist.

- `neighborhood`:

  (\`matrix\` or \`NULL\`) Read-only. The secondary/neighborhood
  structure as a 2-column integer edgelist. \`NULL\` if not provided.

- `overlap`:

  (\`matrix\`) Read-only. The calculated overlap relation (dyads with
  shared neighbors in \`neighborhood\`) as a 2-column integer edgelist.

- `directed`:

  (\`logical\`) Read-only. Indicates if the \`z_network\` is treated as
  directed.

- `n_actor`:

  (\`integer\`) Read-only. The total number of actors (nodes) in the
  network.

- `type_x`:

  (\`character\`) Read-only. The specified distribution type for the
  \`x_attribute\`.

- `type_y`:

  (\`character\`) Read-only. The specified distribution type for the
  \`y_attribute\`.

- `scale_x`:

  (\`numeric\`) Read-only. The scale parameter associated with the
  \`x_attribute\`.

- `scale_y`:

  (\`numeric\`) Read-only. The scale parameter associated with the
  \`y_attribute\`.

## Methods

### Public methods

- [`iglm_data_generator$new()`](#method-iglm_data-new)

- [`iglm_data_generator$gather()`](#method-iglm_data-gather)

- [`iglm_data_generator$save()`](#method-iglm_data-save)

- [`iglm_data_generator$density_z()`](#method-iglm_data-density_z)

- [`iglm_data_generator$density_x()`](#method-iglm_data-density_x)

- [`iglm_data_generator$density_y()`](#method-iglm_data-density_y)

- [`iglm_data_generator$edgewise_shared_partner()`](#method-iglm_data-edgewise_shared_partner)

- [`iglm_data_generator$set_neighborhood_overlap()`](#method-iglm_data-set_neighborhood_overlap)

- [`iglm_data_generator$dyadwise_shared_partner()`](#method-iglm_data-dyadwise_shared_partner)

- [`iglm_data_generator$geodesic_distances_distribution()`](#method-iglm_data-geodesic_distances_distribution)

- [`iglm_data_generator$geodesic_distances()`](#method-iglm_data-geodesic_distances)

- [`iglm_data_generator$edgewise_shared_partner_distribution()`](#method-iglm_data-edgewise_shared_partner_distribution)

- [`iglm_data_generator$dyadwise_shared_partner_distribution()`](#method-iglm_data-dyadwise_shared_partner_distribution)

- [`iglm_data_generator$degree_distribution()`](#method-iglm_data-degree_distribution)

- [`iglm_data_generator$degree()`](#method-iglm_data-degree)

- [`iglm_data_generator$spillover_degree_distribution()`](#method-iglm_data-spillover_degree_distribution)

- [`iglm_data_generator$plot()`](#method-iglm_data-plot)

- [`iglm_data_generator$print()`](#method-iglm_data-print)

- [`iglm_data_generator$clone()`](#method-iglm_data-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \`iglm_data\` object, that includes data on two attributes
and one network.

#### Usage

    iglm_data_generator$new(
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

#### Arguments

- `x_attribute`:

  A numeric vector for the first unit-level attribute.

- `y_attribute`:

  A numeric vector for the second unit-level attribute.

- `z_network`:

  A matrix representing the network. Can be a 2-column edgelist or a
  square adjacency matrix.

- `neighborhood`:

  An optional matrix for the neighborhood representing local dependence.
  Can be a 2-column edgelist or a square adjacency matrix. A tie in
  \`neighborhood\` between actor i and j indicates that j is in the
  neighborhood of i, implying dependence between the respective actors.

- `directed`:

  A logical value indicating if \`z_network\` is directed. If \`NA\`
  (default), directedness is inferred from the symmetry of
  \`z_network\`.

- `n_actor`:

  An integer for the number of actors in the system. If \`NA\`
  (default), \`n_actor\` is inferred from the attributes or network
  matrices.

- `type_x`:

  Character string for the type of \`x_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- `type_y`:

  Character string for the type of \`y_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. Default is
  \`"binomial"\`.

- `scale_x`:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- `scale_y`:

  A positive numeric value for scaling (e.g., variance for "normal"
  type). Default is 1.

- `return_neighborhood`:

  Logical. If \`TRUE\` (default) and \`neighborhood\` is \`NULL\`, a
  full neighborhood (all dyads) is generated implying global dependence.
  If \`FALSE\`, no neighborhood is set.

- `file`:

  (character) Optional file path to load a saved \`iglm_data\` object
  state.

#### Returns

A new \`iglm_data\` object.

------------------------------------------------------------------------

### Method `gather()`

Gathers the current state of the \`iglm_data\` object into a list. This
includes all attributes, network, and configuration details necessary to
reconstruct the object later.

#### Usage

    iglm_data_generator$gather()

#### Returns

A list containing the current state of the \`iglm_data\` object.

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Saves the current state of the \`iglm_data\` object to a specified file
path in RDS format. This includes all attributes, network, and
configuration details necessary to reconstruct the object later.

#### Usage

    iglm_data_generator$save(file)

#### Arguments

- `file`:

  (character) The file where the object state should be saved.

#### Returns

The \`iglm_data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### Method `density_z()`

Calculates the density of the \`z_network\`.

#### Usage

    iglm_data_generator$density_z()

#### Returns

A numeric value for the network density.

------------------------------------------------------------------------

### Method `density_x()`

Calculates the mean of the \`x_attribute\`.

#### Usage

    iglm_data_generator$density_x()

#### Returns

A numeric value for the mean of \`x_attribute\`.

------------------------------------------------------------------------

### Method `density_y()`

Calculates the mean of the \`y_attribute\`.

#### Usage

    iglm_data_generator$density_y()

#### Returns

A numeric value for the mean of \`y_attribute\`.

------------------------------------------------------------------------

### Method `edgewise_shared_partner()`

Calculates the matrix of edgewise shared partners. This is a two-path
matrix (e.g., \$A A^T\$ or \$A^T A\$).

#### Usage

    iglm_data_generator$edgewise_shared_partner(type = "ALL")

#### Arguments

- `type`:

  (character) The type of two-path to calculate for directed networks.
  Ignored if network is undirected. Must be one of: \`"OTP"\` (Outgoing
  Two-Path), \`"ISP"\` (In-Star), \`"OSP"\` (Out-Star), \`"ITP"\`
  (Incoming Two-Path), \`"ALL"\` (Symmetric all-partner). Default is
  \`"ALL"\`.

#### Returns

A sparse matrix (\`dgCMatrix\`) of shared partner counts.

------------------------------------------------------------------------

### Method `set_neighborhood_overlap()`

Sets the neighborhood and overlap matrices.

#### Usage

    iglm_data_generator$set_neighborhood_overlap(neighborhood, overlap)

#### Arguments

- `neighborhood`:

  A matrix for a secondary neighborhood. Can be a 2-column edgelist or a
  square adjacency matrix.

- `overlap`:

  A matrix for the overlap network. Can be a 2-column edgelist or a
  square adjacency matrix.

#### Returns

None. Updates the internal neighborhood and overlap matrices.

------------------------------------------------------------------------

### Method `dyadwise_shared_partner()`

Calculates the matrix of edgewise shared partners. This is a two-path
matrix (e.g., \$A A^T\$ or \$A^T A\$).

#### Usage

    iglm_data_generator$dyadwise_shared_partner(type = "ALL")

#### Arguments

- `type`:

  (character) The type of two-path to calculate for directed networks.
  Ignored if network is undirected. Must be one of: \`"OTP"\` (Outgoing
  Two-Path, z_i,j\*z_j,h ), \`"ISP"\` (In-Star), \`"OSP"\` (Out-Star),
  \`"ITP"\` (Incoming Two-Path), \`"ALL"\` (Symmetric all-partner).
  Default is \`"ALL"\`.

#### Returns

A sparse matrix (\`dgCMatrix\`) of shared partner counts.

------------------------------------------------------------------------

### Method `geodesic_distances_distribution()`

Calculates the geodesic distance distribution of the symmetrized
\`z_network\`.

#### Usage

    iglm_data_generator$geodesic_distances_distribution(
      value_range = NULL,
      prob = TRUE
    )

#### Arguments

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  distances to tabulate. If \`NULL\` (default), the range is inferred
  from the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

#### Returns

A named vector (a \`table\` object) with the distribution of geodesic
distances. Includes \`Inf\` for unreachable pairs.

------------------------------------------------------------------------

### Method `geodesic_distances()`

Calculates the all-pairs geodesic distance matrix for the symmetrized
\`z_network\` using a matrix-based BFS algorithm.

#### Usage

    iglm_data_generator$geodesic_distances()

#### Returns

A sparse matrix (\`dgCMatrix\`) where \`D\[i, j\]\` is the shortest path
distance from i to j. \`Inf\` indicates no path.

------------------------------------------------------------------------

### Method `edgewise_shared_partner_distribution()`

Calculates the distribution of edgewise shared partners.

#### Usage

    iglm_data_generator$edgewise_shared_partner_distribution(
      type = "ALL",
      value_range = NULL,
      prob = TRUE
    )

#### Arguments

- `type`:

  (character) The type of shared partner matrix to use. See
  \`edgewise_shared_partner\` for details. Default is \`"ALL"\`.

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  counts to tabulate. If \`NULL\` (default), the range is inferred from
  the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

#### Returns

A named vector (a \`table\` object) with the distribution of shared
partner counts.

------------------------------------------------------------------------

### Method `dyadwise_shared_partner_distribution()`

Calculates the distribution of edgewise shared partners.

#### Usage

    iglm_data_generator$dyadwise_shared_partner_distribution(
      type = "ALL",
      value_range = NULL,
      prob = TRUE
    )

#### Arguments

- `type`:

  (character) The type of shared partner matrix to use. See
  \`edgewise_shared_partner\` for details. Default is \`"ALL"\`.

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  counts to tabulate. If \`NULL\` (default), the range is inferred from
  the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

#### Returns

A named vector (a \`table\` object) with the distribution of shared
partner counts.

------------------------------------------------------------------------

### Method `degree_distribution()`

Calculates the degree distribution of the \`z_network\`.

#### Usage

    iglm_data_generator$degree_distribution(value_range = NULL, prob = TRUE)

#### Arguments

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  degrees to tabulate. If \`NULL\` (default), the range is inferred from
  the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

#### Returns

If the network is directed, a list containing two \`table\` objects:
\`in_degree\` and \`out_degree\`. If undirected, a single \`table\`
object with the degree distribution.

------------------------------------------------------------------------

### Method `degree()`

Calculates the degree sequence(s) of the \`z_network\`.

#### Usage

    iglm_data_generator$degree()

#### Returns

If the network is directed, a list containing two vectors:
\`in_degree_seq\` and \`out_degree_seq\`. If undirected, a single list
containing the vector \`degree_seq\`.

------------------------------------------------------------------------

### Method `spillover_degree_distribution()`

Calculates the spillover degree distribution between actors with
\`x_attribute == 1\` and actors with \`y_attribute == 1\`.

#### Usage

    iglm_data_generator$spillover_degree_distribution(
      prob = TRUE,
      value_range = NULL
    )

#### Arguments

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  degrees to tabulate. If \`NULL\` (default), the range is inferred from
  the data.

#### Returns

A list containing two \`table\` objects: \`out_spillover_degree\` (from
x_i=1 to y_j=1) and \`in_spillover_degree\` (from y_i=1 to x_j=1).

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot the network using \`igraph\`.

Visualizes the \`z_network\` using the \`igraph\` package. Nodes can be
colored by \`x_attribute\` and sized by \`y_attribute\`.
\`neighborhood\` edges can be plotted as a background layer.

#### Usage

    iglm_data_generator$plot(
      node_color = "x",
      node_size = "y",
      show_overlap = TRUE,
      layout = igraph::layout_with_fr,
      network_edges_col = "grey60",
      neighborhood_edges_col = "orange",
      main = "",
      legend_col_n_levels = NULL,
      legend_size_n_levels = NULL,
      legend_pos = "right",
      alpha_neighborhood = 0.2,
      edge.width = 1,
      edge.arrow.size = 1,
      vertex.frame.width = 0.5,
      coords = NULL,
      ...
    )

#### Arguments

- `node_color`:

  (character) Attribute to map to node color. One of \`"x"\` (default),
  \`"y"\`, or \`"none"\`.

- `node_size`:

  (character) Attribute to map to node size. One of \`"y"\` (default),
  \`"x"\`, or \`"constant"\`.

- `show_overlap`:

  (logical) If \`TRUE\` (default), plot the \`neighborhood\` edges as a
  background layer.

- `layout`:

  An \`igraph\` layout function (e.g., \`igraph::layout_with_fr\`).

- `network_edges_col`:

  (character) Color for the \`z_network\` edges.

- `neighborhood_edges_col`:

  (character) Color for the \`neighborhood\` edges.

- `main`:

  (character) The main title for the plot.

- `legend_col_n_levels`:

  (integer) Number of levels for the color legend.

- `legend_size_n_levels`:

  (integer) Number of levels for the size legend.

- `legend_pos`:

  (character) Position of the legend (e.g., \`"right"\`).

- `alpha_neighborhood`:

  (numeric) Alpha transparency for neighborhood edges.

- `edge.width`:

  (numeric) Width of the network edges.

- `edge.arrow.size`:

  (numeric) Size of the arrowheads for directed edges.

- `vertex.frame.width`:

  (numeric) Width of the vertex frame.

- `coords`:

  (matrix) Optional matrix of x-y coordinates for node layout.

- `...`:

  Additional arguments passed to \`plot.igraph\`.

#### Returns

A list containing the \`igraph\` object (\`graph\`) and the layout
coordinates (\`coords\`), invisibly.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a summary of the \`iglm_data\` object to the console.

#### Usage

    iglm_data_generator$print(digits = 3, ...)

#### Arguments

- `digits`:

  (integer) Number of digits to round numeric output to.

- `...`:

  Additional arguments (not used).

#### Returns

The object's private environment, invisibly.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    iglm_data_generator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

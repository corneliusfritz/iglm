# Networks with Unit-Level Attributes (R6 Class)

The \`iglm.data\` class is a container for storing, validating, and
analyzing unit-level attributes (x_attribute, y_attribute) and
connections (z_network).

## Active bindings

- `label_x`:

  (\`character\`) Label/name for \`x_attribute\`.

- `label_y`:

  (\`character\`) Label/name for \`y_attribute\`.

- `label_z`:

  (\`character\`) Label/name for \`z_network\`.

- `x_attribute`:

  (\`numeric\`) The vector for the first unit-level attribute.

- `y_attribute`:

  (\`numeric\`) The vector for the second unit-level attribute.

- `z_network`:

  (\`matrix\`) The primary network structure as a 2-column integer
  edgelist.

- `neighborhood`:

  (\`matrix\`) Read-only. The secondary/neighborhood structure as a
  2-column integer edgelist. An empty matrix if not provided.

- `overlap`:

  (\`matrix\`) Read-only. The calculated overlap relation (dyads with
  shared neighbors in \`neighborhood\`) as a 2-column integer edgelist.
  An empty matrix if overlap hasn't been computed or is not available.

- `directed`:

  (\`logical\`) Indicates if the \`z_network\` is treated as directed.

- `n_actor`:

  (\`integer\`) The total number of actors (nodes) in the network.

- `type_x`:

  (\`character\`) The specified distribution type for the
  \`x_attribute\`.

- `type_y`:

  (\`character\`) The specified distribution type for the
  \`y_attribute\`.

- `scale_x`:

  (\`numeric\`) The scale parameter associated with the \`x_attribute\`.

- `scale_y`:

  (\`numeric\`) The scale parameter associated with the \`y_attribute\`.

- `fix_x`:

  (\`logical\`) Indicates if the \`x_attribute\` is fixed during
  estimation/simulation.

- `fix_z`:

  (\`logical\`) RIndicates if the \`z_network\` is fixed during
  estimation/simulation.

- `descriptives`:

  (\`list\`)A list storing computed descriptive statistics for the
  network and attributes.

- `fix_z_alocal`:

  (\`logical\`) Flag indicating whether nonoverlap edges are treated as
  random.

## Methods

### Public methods

- [`iglm.data$new()`](#method-iglm.data-initialize)

- [`iglm.data$set_z_network()`](#method-iglm.data-set_z_network)

- [`iglm.data$set_type_x()`](#method-iglm.data-set_type_x)

- [`iglm.data$set_type_y()`](#method-iglm.data-set_type_y)

- [`iglm.data$set_scale_x()`](#method-iglm.data-set_scale_x)

- [`iglm.data$set_scale_y()`](#method-iglm.data-set_scale_y)

- [`iglm.data$set_x_attribute()`](#method-iglm.data-set_x_attribute)

- [`iglm.data$set_y_attribute()`](#method-iglm.data-set_y_attribute)

- [`iglm.data$set_label_x()`](#method-iglm.data-set_label_x)

- [`iglm.data$set_label_y()`](#method-iglm.data-set_label_y)

- [`iglm.data$set_label_z()`](#method-iglm.data-set_label_z)

- [`iglm.data$gather()`](#method-iglm.data-gather)

- [`iglm.data$set_fix_z_alocal()`](#method-iglm.data-set_fix_z_alocal)

- [`iglm.data$delete_isolates()`](#method-iglm.data-delete_isolates)

- [`iglm.data$save()`](#method-iglm.data-save)

- [`iglm.data$set_fix_x()`](#method-iglm.data-set_fix_x)

- [`iglm.data$set_fix_z()`](#method-iglm.data-set_fix_z)

- [`iglm.data$mean_z()`](#method-iglm.data-mean_z)

- [`iglm.data$mean_x()`](#method-iglm.data-mean_x)

- [`iglm.data$mean_y()`](#method-iglm.data-mean_y)

- [`iglm.data$x_distribution()`](#method-iglm.data-x_distribution)

- [`iglm.data$x_dist()`](#method-iglm.data-x_dist)

- [`iglm.data$y_distribution()`](#method-iglm.data-y_distribution)

- [`iglm.data$y_dist()`](#method-iglm.data-y_dist)

- [`iglm.data$edgewise_shared_partner()`](#method-iglm.data-edgewise_shared_partner)

- [`iglm.data$set_neighborhood_overlap()`](#method-iglm.data-set_neighborhood_overlap)

- [`iglm.data$dyadwise_shared_partner()`](#method-iglm.data-dyadwise_shared_partner)

- [`iglm.data$geodesic_distances_distribution()`](#method-iglm.data-geodesic_distances_distribution)

- [`iglm.data$geodesic_distances()`](#method-iglm.data-geodesic_distances)

- [`iglm.data$esp()`](#method-iglm.data-esp)

- [`iglm.data$esp_dist()`](#method-iglm.data-esp_dist)

- [`iglm.data$dsp()`](#method-iglm.data-dsp)

- [`iglm.data$dsp_dist()`](#method-iglm.data-dsp_dist)

- [`iglm.data$geo()`](#method-iglm.data-geo)

- [`iglm.data$geo_dist()`](#method-iglm.data-geo_dist)

- [`iglm.data$edgewise_shared_partner_distribution()`](#method-iglm.data-edgewise_shared_partner_distribution)

- [`iglm.data$dyadwise_shared_partner_distribution()`](#method-iglm.data-dyadwise_shared_partner_distribution)

- [`iglm.data$degree_distribution()`](#method-iglm.data-degree_distribution)

- [`iglm.data$degree()`](#method-iglm.data-degree)

- [`iglm.data$deg()`](#method-iglm.data-deg)

- [`iglm.data$deg_dist()`](#method-iglm.data-deg_dist)

- [`iglm.data$plot()`](#method-iglm.data-plot)

- [`iglm.data$print()`](#method-iglm.data-print)

- [`iglm.data$clone()`](#method-iglm.data-clone)

------------------------------------------------------------------------

### `iglm.data$new()`

Create a new \`iglm.data\` object, that includes data on two attributes
and one network.

#### Usage

    iglm.data$new(
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
      fix_x = FALSE,
      fix_z = FALSE,
      fix_z_alocal = TRUE,
      return_neighborhood = TRUE,
      file = NULL,
      label_x = "x",
      label_y = "y",
      label_z = "z"
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

- `fix_x`:

  Logical. If \`TRUE\`, the \`x_attribute\` is treated as fixed during
  model estimation and simulation. Default is \`FALSE\`.

- `fix_z`:

  Logical. If \`TRUE\`, the \`z_network\` is treated as fixed during
  model estimation and simulation. Default is \`FALSE\`.

- `fix_z_alocal`:

  Logical. If \`TRUE\` (default), alocal dyads in the neighborhood are
  fixed.

- `return_neighborhood`:

  Logical. If \`TRUE\` (default) and \`neighborhood\` is \`NULL\`, a
  full neighborhood (all dyads) is generated implying global dependence.
  If \`FALSE\`, no neighborhood is set.

- `file`:

  (character) Optional file path to load a saved \`iglm.data\` object
  state.

- `label_x`:

  Character string for the label/name of \`x_attribute\`. Default is
  \`"x"\`.

- `label_y`:

  Character string for the label/name of \`y_attribute\`. Default is
  \`"y"\`.

- `label_z`:

  Character string for the label/name of \`z_network\`. Default is
  \`"z"\`.

#### Returns

A new \`iglm.data\` object.

------------------------------------------------------------------------

### `iglm.data$set_z_network()`

Sets the \`z_network\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_z_network(z_network)

#### Arguments

- `z_network`:

  A matrix representing the network. Can be a 2-column edgelist or a
  square adjacency matrix. @return The \`iglm.data\` object itself
  (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_type_x()`

Sets the \`type_x\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_type_x(type_x)

#### Arguments

- `type_x`:

  A character string for the type of \`x_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`. @return The
  \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_type_y()`

Sets the \`type_y\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_type_y(type_y)

#### Arguments

- `type_y`:

  A character string for the type of \`y_attribute\`. Must be one of
  \`"binomial"\`, \`"poisson"\`, or \`"normal"\`.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_scale_x()`

Sets the \`scale_x\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_scale_x(scale_x)

#### Arguments

- `scale_x`:

  A positive numeric value for scaling (e.g., variance for "normal"
  type).

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_scale_y()`

Sets the \`scale_y\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_scale_y(scale_y)

#### Arguments

- `scale_y`:

  A positive numeric value for scaling (e.g., variance for "normal"
  type).

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_x_attribute()`

Sets the \`x_attribute\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_x_attribute(x_attribute)

#### Arguments

- `x_attribute`:

  A numeric vector for the first unit-level attribute.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_y_attribute()`

Sets the \`y_attribute\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_y_attribute(y_attribute)

#### Arguments

- `y_attribute`:

  A numeric vector for the first unit-level attribute.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_label_x()`

Sets the label for the \`x_attribute\`.

#### Usage

    iglm.data$set_label_x(label_x)

#### Arguments

- `label_x`:

  A character string for the label of \`x_attribute\`.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_label_y()`

Sets the label for the \`y_attribute\`.

#### Usage

    iglm.data$set_label_y(label_y)

#### Arguments

- `label_y`:

  A character string for the label of \`y_attribute\`.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_label_z()`

Sets the label for the \`z_network\`.

#### Usage

    iglm.data$set_label_z(label_z)

#### Arguments

- `label_z`:

  A character string for the label of \`z_network\`.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$gather()`

Gathers the current state of the \`iglm.data\` object into a list. This
includes all attributes, network, and configuration details necessary to
reconstruct the object later.

#### Usage

    iglm.data$gather()

#### Returns

A list containing the current state of the \`iglm.data\` object.

------------------------------------------------------------------------

### `iglm.data$set_fix_z_alocal()`

Sets the option whether alocal edges are fixed or not.

#### Usage

    iglm.data$set_fix_z_alocal(fix_z_alocal)

#### Arguments

- `fix_z_alocal`:

  A logical value indicating whether alocal edges should be treated as
  fixed or not.

------------------------------------------------------------------------

### `iglm.data$delete_isolates()`

Deletes isolates from the \`z_network\` and updates the attributes and
neighborhood accordingly. Isolates are actors that do not have any
connections in the \`z_network\`. This method identifies such actors,
removes them from the attributes and neighborhood, and updates the
\`z_network\` to reflect the new actor indices.

#### Usage

    iglm.data$delete_isolates()

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$save()`

Saves the current state of the \`iglm.data\` object to a specified file
path in RDS format. This includes all attributes, network, and
configuration details necessary to reconstruct the object later.

#### Usage

    iglm.data$save(file)

#### Arguments

- `file`:

  (character) The file where the object state should be saved. Must have
  a .rds extension.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_fix_x()`

Sets the \`fix_x\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_fix_x(fix_x)

#### Arguments

- `fix_x`:

  A logical value indicating if \`x_attribute\` is fixed or random.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$set_fix_z()`

Sets the \`fix_z\` of the \`iglm.data\` object.

#### Usage

    iglm.data$set_fix_z(fix_z)

#### Arguments

- `fix_z`:

  A logical value indicating if \`z_network\` is fixed or random.

#### Returns

The \`iglm.data\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$mean_z()`

Calculates the density of the \`z_network\`.

#### Usage

    iglm.data$mean_z()

#### Returns

A numeric value for the network density.

------------------------------------------------------------------------

### `iglm.data$mean_x()`

Calculates the mean of the \`x_attribute\`.

#### Usage

    iglm.data$mean_x()

#### Returns

A numeric value for the mean of \`x_attribute\`.

------------------------------------------------------------------------

### `iglm.data$mean_y()`

Calculates the mean of the \`y_attribute\`.

#### Usage

    iglm.data$mean_y()

#### Returns

A numeric value for the mean of \`y_attribute\`.

------------------------------------------------------------------------

### `iglm.data$x_distribution()`

Calculates the distribution of the \`x_attribute\`.

#### Usage

    iglm.data$x_distribution(value_range = NULL, prob = TRUE, plot = TRUE)

#### Arguments

- `value_range`:

  (numeric vector) Optional range of values to consider for the
  distribution. If \`NULL\` (default), the range is inferred from the
  data.

- `prob`:

  (logical) If \`TRUE\` (default), returns probabilities; if \`FALSE\`,
  returns frequencies.

- `plot`:

  (logical) If \`TRUE\` (default), plots the distribution using a
  density plot for continuous data or a bar plot for discrete data.

#### Returns

A numeric vector representing the distribution of \`x_attribute\`
(invisible).

------------------------------------------------------------------------

### `iglm.data$x_dist()`

Short alias for \`x_distribution\`.

#### Usage

    iglm.data$x_dist(value_range = NULL, prob = TRUE, plot = TRUE)

#### Arguments

- `value_range`:

  (numeric vector) Optional range of values to consider for the
  distribution. If \`NULL\` (default), the range is inferred from the
  data.

- `prob`:

  (logical) If \`TRUE\` (default), returns probabilities; if \`FALSE\`,
  returns frequencies.

- `plot`:

  (logical) If \`TRUE\` (default), plots the distribution.

#### Returns

A numeric vector representing the distribution of \`x_attribute\`
(invisible).

------------------------------------------------------------------------

### `iglm.data$y_distribution()`

Calculates the distribution of the \`y_attribute\`.

#### Usage

    iglm.data$y_distribution(value_range = NULL, prob = TRUE, plot = TRUE)

#### Arguments

- `value_range`:

  (numeric vector) Optional range of values to consider for the
  distribution. If \`NULL\` (default), the range is inferred from the
  data.

- `prob`:

  (logical) If \`TRUE\` (default), returns probabilities; if \`FALSE\`,
  returns frequencies.

- `plot`:

  (logical) If \`TRUE\` (default), plots the distribution using a
  density plot for continuous data or a bar plot for discrete data.

#### Returns

A numeric vector representing the distribution of \`y_attribute\`
(invisible).

------------------------------------------------------------------------

### `iglm.data$y_dist()`

Short alias for \`y_distribution\`.

#### Usage

    iglm.data$y_dist(value_range = NULL, prob = TRUE, plot = TRUE)

#### Arguments

- `value_range`:

  (numeric vector) Optional range of values to consider for the
  distribution. If \`NULL\` (default), the range is inferred from the
  data.

- `prob`:

  (logical) If \`TRUE\` (default), returns probabilities; if \`FALSE\`,
  returns frequencies.

- `plot`:

  (logical) If \`TRUE\` (default), plots the distribution.

#### Returns

A numeric vector representing the distribution of \`y_attribute\`
(invisible).

------------------------------------------------------------------------

### `iglm.data$edgewise_shared_partner()`

Calculates the matrix of edgewise shared partners. This is a two-path
matrix (e.g., \$A A^T\$ or \$A^T A\$).

#### Usage

    iglm.data$edgewise_shared_partner(type = "ALL", mode = "global")

#### Arguments

- `type`:

  (character) The type of two-path to calculate for directed networks.
  Ignored if network is undirected. Must be one of: \`"OTP"\` (Outgoing
  Two-Path, \\z\_{i,j}\\ z\_{i,h} \\ z\_{j,h}\\ ), \`"ISP"\` (Ingoing
  Shared Partner, \\z\_{i,j}\\ z\_{h,i} \\ z\_{j,h}\\), \`"OSP"\`
  (Outgoing Shared Partner, \\z\_{i,j}\\ z\_{i,h} \\ z\_{j,h}\\),
  \`"ITP"\` (Incoming Two-Path, \\z\_{i,j}\\ z\_{h,i} \\ z\_{j,h}\\),
  \`"ALL"\` (Any one of the above). Default is \`"ALL"\`.

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all
  edges, or \`"local"\` to evaluate only edges with overlapping
  neighborhoods (from \`overlap\`).

#### Returns

A numeric vector of shared partner counts for edges.

------------------------------------------------------------------------

### `iglm.data$set_neighborhood_overlap()`

Sets the neighborhood and overlap matrices.

#### Usage

    iglm.data$set_neighborhood_overlap(neighborhood, overlap)

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

### `iglm.data$dyadwise_shared_partner()`

Calculates the matrix of dyadwise shared partners.

#### Usage

    iglm.data$dyadwise_shared_partner(type = "ALL", mode = "global")

#### Arguments

- `type`:

  (character) The type of two-path to calculate for directed networks.
  Ignored if network is undirected. Must be one of: \`"OTP"\` (Outgoing
  Two-Path, \\z\_{i,h} \\ z\_{j,h}\\ ), \`"ISP"\` (Ingoing Shared
  Partner, \\z\_{h,i} \\ z\_{j,h}\\), \`"OSP"\` (Outgoing Shared
  Partner, \\z\_{i,h} \\ z\_{j,h}\\), \`"ITP"\` (Incoming Two-Path,
  \\z\_{h,i} \\ z\_{j,h}\\), \`"ALL"\` (Any one of the above). Default
  is \`"ALL"\`.

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all
  dyads, or \`"local"\` to evaluate only dyads with overlapping
  neighborhoods.

#### Returns

A sparse matrix (\`dgCMatrix\`) of shared partner counts.

------------------------------------------------------------------------

### `iglm.data$geodesic_distances_distribution()`

Calculates the geodesic distance distribution of the symmetrized
\`z_network\`.

#### Usage

    iglm.data$geodesic_distances_distribution(
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
    )

#### Arguments

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  distances to tabulate. If \`NULL\` (default), the range is inferred
  from the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all node
  pairs, or \`"local"\` to evaluate only pairs with overlapping
  neighborhoods (from \`overlap\`).

#### Returns

A named vector (a \`table\` object) with the distribution of geodesic
distances. Includes \`Inf\` for unreachable pairs.

------------------------------------------------------------------------

### `iglm.data$geodesic_distances()`

Calculates the all-pairs geodesic distance matrix for the symmetrized
\`z_network\` using a matrix-based BFS algorithm.

#### Usage

    iglm.data$geodesic_distances(mode = "global")

#### Arguments

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all
  pairs, or \`"local"\` to evaluate only pairs with overlapping
  neighborhoods.

#### Returns

A sparse matrix (\`dgCMatrix\`) where \`D\[i, j\]\` is the shortest path
distance from i to j. \`Inf\` indicates no path.

------------------------------------------------------------------------

### `iglm.data$esp()`

Short alias for \`edgewise_shared_partner\`.

#### Usage

    iglm.data$esp(type = "ALL", mode = "global")

#### Arguments

- `type`:

  (character) The type of two-path to calculate. Default is \`"ALL"\`.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A numeric vector of shared partner counts for edges.

------------------------------------------------------------------------

### `iglm.data$esp_dist()`

Short alias for \`edgewise_shared_partner_distribution\`.

#### Usage

    iglm.data$esp_dist(
      type = "ALL",
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
    )

#### Arguments

- `type`:

  (character) The type of two-path to calculate. Default is \`"ALL"\`.

- `value_range`:

  (numeric vector) Range of counts to tabulate.

- `prob`:

  (logical) If \`TRUE\` (default), returns proportions.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A named vector with the distribution.

------------------------------------------------------------------------

### `iglm.data$dsp()`

Short alias for \`dyadwise_shared_partner\`.

#### Usage

    iglm.data$dsp(type = "ALL", mode = "global")

#### Arguments

- `type`:

  (character) The type of two-path to calculate. Default is \`"ALL"\`.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A sparse matrix (\`dgCMatrix\`) of shared partner counts.

------------------------------------------------------------------------

### `iglm.data$dsp_dist()`

Short alias for \`dyadwise_shared_partner_distribution\`.

#### Usage

    iglm.data$dsp_dist(
      type = "ALL",
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
    )

#### Arguments

- `type`:

  (character) The type of two-path to calculate. Default is \`"ALL"\`.

- `value_range`:

  (numeric vector) Range of counts to tabulate.

- `prob`:

  (logical) If \`TRUE\` (default), returns proportions.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A named vector with the distribution.

------------------------------------------------------------------------

### `iglm.data$geo()`

Short alias for \`geodesic_distances\`.

#### Usage

    iglm.data$geo(mode = "global")

#### Arguments

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A sparse matrix (\`dgCMatrix\`) of geodesic distances.

------------------------------------------------------------------------

### `iglm.data$geo_dist()`

Short alias for \`geodesic_distances_distribution\`.

#### Usage

    iglm.data$geo_dist(
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
    )

#### Arguments

- `value_range`:

  (numeric vector) Range of distances to tabulate.

- `prob`:

  (logical) If \`TRUE\` (default), returns proportions.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

A named vector with the distribution.

------------------------------------------------------------------------

### `iglm.data$edgewise_shared_partner_distribution()`

Calculates the distribution of edgewise shared partners.

#### Usage

    iglm.data$edgewise_shared_partner_distribution(
      type = "ALL",
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
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

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all
  edges, or \`"local"\` to evaluate only edges with overlapping
  neighborhoods.

#### Returns

A named vector (a \`table\` object) with the distribution of shared
partner counts.

------------------------------------------------------------------------

### `iglm.data$dyadwise_shared_partner_distribution()`

Calculates the distribution of dyadwise shared partners.

#### Usage

    iglm.data$dyadwise_shared_partner_distribution(
      type = "ALL",
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      mode = "global"
    )

#### Arguments

- `type`:

  (character) The type of shared partner matrix to use. See
  \`dyadwise_shared_partner\` for details. Default is \`"ALL"\`.

- `value_range`:

  (numeric vector) A vector \`c(min, max)\` specifying the range of
  counts to tabulate. If \`NULL\` (default), the range is inferred from
  the data.

- `prob`:

  (logical) If \`TRUE\` (default), returns a probability distribution
  (proportions). If \`FALSE\`, returns raw counts.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `mode`:

  (character) Either \`"global"\` (default) to evaluate across all
  dyads, or \`"local"\` to evaluate only dyads with overlapping
  neighborhoods.

#### Returns

A named vector (a \`table\` object) with the distribution of shared
partner counts.

------------------------------------------------------------------------

### `iglm.data$degree_distribution()`

Calculates the degree distribution of the \`z_network\`.

A flexible, general function for evaluating network connectivity across
global networks, local neighborhoods, attribute-defined subgroups, and
cross-group spillover pathways.

#### Topological Scope (`mode`)

- `"global"` (default): Evaluates degree distributions across all dyads
  in the network.

- `"local"`: Evaluates local degree distributions restricted strictly to
  actor pairs that share an overlapping neighborhood (`overlap`).

#### Directionality and Bipartite Subgroups

- **Directed networks**: Always returns both `out_degree` (ties sent)
  and `in_degree` (ties received).

- **Undirected networks**: Returns a single overall `degree`
  distribution when unconstrained or single-side constrained. When
  bilateral constraints are supplied (e.g., sender \\i\\ and receiver
  \\j\\), ties are evaluated directionally (\\i \to j\\), returning both
  `out_degree` and `in_degree`.

#### Spillover and Subgroup Conditioning

Any combination of sender attributes (`x_i`, `y_i`) and receiver
attributes (`x_j`, `y_j`) can be specified to measure spillover
dynamics:

- `out_degree`: Distribution of ties sent from matching senders \\i\\ to
  matching receivers \\j\\ (spillover sending capacity).

- `in_degree`: Distribution of ties received by matching receivers \\j\\
  from matching senders \\i\\ (spillover exposure).

#### Supported Constraint Formats & Internal Handling

Attribute constraints (`x_i`, `x_j`, `y_i`, `y_j`) accept:

- **Exact scalar values**: For binary attributes, matches actors with
  that exact value (e.g., `x_i = 1`).

- **Continuous / Count shortcuts**: When an attribute is continuous
  (`"normal"`) or count (`"poisson"`), setting `1` internally selects
  above-mean actors (\\x_i \> \bar{x}\\), and `0` selects
  below-or-equal-to-mean actors (\\x_i \le \bar{x}\\). Any other numeric
  value \\v\\ matches actors with exact value \\v\\.

- **Discrete value sets**: Vectors such as `x_i = c(1, 2)` match actors
  with any value in that set.

- **Filtering functions**: Custom functions (vectorized or scalar),
  e.g., `x_i = function(x) x > 0.5` or
  `y_j = \(y) if (y > 2) TRUE else FALSE`.

#### Plotting and Axis Labels

When `plot = TRUE`, mathematical expressions are formatted automatically
for the x-axis:

- Exact values and sets display as \\x_i == 1\\ or \\x_i == \text{c(1,
  2)}\\.

- Continuous shortcuts display with sample mean bars as \\x_i \>
  \bar{x}\\ or \\x_i \le \bar{x}\\.

- Filtering functions display as \\x_i == \text{"fn"}\\.

#### Usage

    iglm.data$degree_distribution(
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      x_i = NULL,
      x_j = NULL,
      y_i = NULL,
      y_j = NULL,
      mode = "global"
    )

#### Arguments

- `value_range`:

  (numeric vector or list) A vector `c(min, max)` specifying the range
  of degrees to tabulate, or a list with `in_degree` and `out_degree`.
  If `NULL` (default), ranges are inferred from the data.

- `prob`:

  (logical) If `TRUE` (default), returns a probability distribution
  (proportions). If `FALSE`, returns raw counts.

- `plot`:

  (logical) If `TRUE`, plots the degree distribution barplot(s).

- `x_i`:

  (optional) Exact value, vector, or filtering function for attribute
  `x` of sender actor \\i\\.

- `x_j`:

  (optional) Exact value, vector, or filtering function for attribute
  `x` of receiver actor \\j\\.

- `y_i`:

  (optional) Exact value, vector, or filtering function for attribute
  `y` of sender actor \\i\\.

- `y_j`:

  (optional) Exact value, vector, or filtering function for attribute
  `y` of receiver actor \\j\\.

- `mode`:

  (character) Either `"global"` (default) to evaluate across all dyads,
  or `"local"` to evaluate only ties within overlapping neighborhoods
  (`overlap`).

#### Returns

If the network is directed or if bilateral constraints are provided, a
list containing two `table` objects: `out_degree` and `in_degree`. If
undirected without bilateral constraints, a single `table` object with
the degree distribution.

#### Examples

    data(copenhagen)

    # 1. Standard global degree distribution
    copenhagen$degree_distribution(plot = FALSE)

    # 2. Local degree distribution restricted to overlapping neighborhoods
    copenhagen$degree_distribution(mode = "local", plot = FALSE)

    # 3. Spillover degree using exact attribute values
    copenhagen$deg_dist(x_i = 1, y_j = 1, mode = "local", plot = FALSE)

    # 4. Spillover degree using filtering functions
    copenhagen$deg_dist(
      x_i = function(x) x > mean(x),
      y_j = function(y) y > mean(y),
      mode = "local",
      plot = FALSE
    )

------------------------------------------------------------------------

### `iglm.data$degree()`

Calculates the degree sequence(s) of the \`z_network\`.

General function for calculating actor-level degree sequences across
global topologies, local neighborhoods, attribute-defined subsets, or
directional spillover pathways.

#### Usage

    iglm.data$degree(
      x_i = NULL,
      x_j = NULL,
      y_i = NULL,
      y_j = NULL,
      mode = "global"
    )

#### Arguments

- `x_i`:

  (optional) Exact value, vector, or filtering function for attribute
  `x` of sender actor \\i\\.

- `x_j`:

  (optional) Exact value, vector, or filtering function for attribute
  `x` of receiver actor \\j\\.

- `y_i`:

  (optional) Exact value, vector, or filtering function for attribute
  `y` of sender actor \\i\\.

- `y_j`:

  (optional) Exact value, vector, or filtering function for attribute
  `y` of receiver actor \\j\\.

- `mode`:

  (character) `"global"` (default) or `"local"`.

#### Returns

If the network is directed or if bilateral constraints are given, a list
containing two numeric vectors: `out_degree_seq` and `in_degree_seq`. If
undirected without bilateral constraints, a list containing the vector
`degree_seq`.

#### Examples

    data(copenhagen)

    # Global degree sequence
    copenhagen$degree()

    # Local spillover degree sequence with filtering functions
    copenhagen$deg(
      x_i = function(x) x > 2,
      y_j = function(y) y > 2,
      mode = "local"
    )

------------------------------------------------------------------------

### `iglm.data$deg()`

Short alias for \`degree\`.

#### Usage

    iglm.data$deg(x_i = NULL, x_j = NULL, y_i = NULL, y_j = NULL, mode = "global")

#### Arguments

- `x_i`:

  Optional sender attribute constraint.

- `x_j`:

  Optional receiver attribute constraint.

- `y_i`:

  Optional sender attribute constraint.

- `y_j`:

  Optional receiver attribute constraint.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

Degree sequence(s).

------------------------------------------------------------------------

### `iglm.data$deg_dist()`

Short alias for \`degree_distribution\`. Supports standard, local, and
attribute-constrained (spillover) degree distributions.

#### Usage

    iglm.data$deg_dist(
      value_range = NULL,
      prob = TRUE,
      plot = TRUE,
      x_i = NULL,
      x_j = NULL,
      y_i = NULL,
      y_j = NULL,
      mode = "global"
    )

#### Arguments

- `value_range`:

  Optional range of degrees to tabulate.

- `prob`:

  (logical) If \`TRUE\`, returns proportions.

- `plot`:

  (logical) If \`TRUE\`, plots the distribution.

- `x_i`:

  Optional sender attribute constraint.

- `x_j`:

  Optional receiver attribute constraint.

- `y_i`:

  Optional sender attribute constraint.

- `y_j`:

  Optional receiver attribute constraint.

- `mode`:

  (character) \`"global"\` (default) or \`"local"\`.

#### Returns

Degree distribution table(s).

------------------------------------------------------------------------

### `iglm.data$plot()`

Plot the network using \`igraph\`.

Visualizes the \`z_network\` using the \`igraph\` package. Nodes can be
colored by \`x_attribute\` and sized by \`y_attribute\`.
\`neighborhood\` edges can be plotted as a background layer.

#### Usage

    iglm.data$plot(
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
      legend_size = 0.5,
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

- `legend_size`:

  (numeric) Scaling factor for the size legend.

- `...`:

  Additional arguments passed to \`plot.igraph\`.

#### Returns

A list containing the \`igraph\` object (\`graph\`) and the layout
coordinates (\`coords\`), invisibly.

------------------------------------------------------------------------

### `iglm.data$print()`

Print a summary of the \`iglm.data\` object to the console.

#### Usage

    iglm.data$print(digits = 3, ...)

#### Arguments

- `digits`:

  (integer) Number of digits to round numeric output to.

- `...`:

  Additional arguments (not used).

#### Returns

The object's private environment, invisibly.

------------------------------------------------------------------------

### `iglm.data$clone()`

The objects of this class are cloneable with this method.

#### Usage

    iglm.data$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `iglm.data$degree_distribution()`
## ------------------------------------------------

data(copenhagen)

# 1. Standard global degree distribution
copenhagen$degree_distribution(plot = FALSE)

# 2. Local degree distribution restricted to overlapping neighborhoods
copenhagen$degree_distribution(mode = "local", plot = FALSE)

# 3. Spillover degree using exact attribute values
copenhagen$deg_dist(x_i = 1, y_j = 1, mode = "local", plot = FALSE)

# 4. Spillover degree using filtering functions
copenhagen$deg_dist(
  x_i = function(x) x > mean(x),
  y_j = function(y) y > mean(y),
  mode = "local",
  plot = FALSE
)

## ------------------------------------------------
## Method `iglm.data$degree()`
## ------------------------------------------------

data(copenhagen)

# Global degree sequence
copenhagen$degree()
#> $degree_seq
#>   [1] 13 28  4 27  5 11 26  4 14  1 15  5  4 11  2  9 20  4  7  5 56 14 14 13 16
#>  [26] 24 11 11 13  1 18  8 38 24 10  3  6  9 22  5  6  9  6  6  4  9 21  3 26  8
#>  [51] 47 15 15  7 13 14  8 40  5  5 34 14 23 12 12  4 12  4  6  4 13 24  7  5  4
#>  [76]  2 34 12  8 11 44  7  4  6 10  7  7 11 20  8 13  8  8  3  6  2  7  6  4  6
#> [101]  3 10  7  5  9 11  8 22 54  8 15 17 31 13  9 21  7 13  5  7 32 12 13 12 23
#> [126]  1  7 11  5 24 10 16 15 13 45 26 12 10  9  8  3  8 11 50  5  6  3 13 20  9
#> [151]  5  1 11  2  6 36  4 14 10 10  1 13 19 15  7  6 12 14 30 14 19  9  6  4  5
#> [176] 29 29 18 40 20  5 15  5 21 14  4  5  6 10 20 11  3 11  8 10 14  5 21 38 13
#> [201] 15 21 19  5 11 16 15  8  8 11  9 20  3  2  9  2 12  9 12 11  5 25  5 10 17
#> [226]  9  4  7 17  1 14  5 15  5 21  5  3  7 16 10 53  5 18  9  8 15  7 12  5 65
#> [251]  1 30  9 13  8  1 17  5  1  9 33  7 15  8 11  5  8  9 29  6  3  4  9 39 39
#> [276]  6 32  6 10  7 27 11  3 11  7  9 11 21 14  2  3 16  6  9  1 12  6  9  9  7
#> [301] 14  9  5  8 10 11 10 34 13 11  9 21  5  5  7 14 16  3  3  6  5 24 17 12 10
#> [326]  8  3 27 12 11 12  7  7  2  6  7  5 10 13  8  4 13  6  8 11 12 27  4 18 17
#> [351]  4 11  5  7 28  4  3 24  5 10  8 15  5  6  1 19  7  7 28 15 12 11  9  8  9
#> [376] 13 25 25  6  3 12 11  5 10  8 22 47 10  8  5 47  1 12 18 14 18  8  4  7 10
#> [401] 41 11  1  8  2 11 14  8 13
#> 

# Local spillover degree sequence with filtering functions
copenhagen$deg(
  x_i = function(x) x > 2,
  y_j = function(y) y > 2,
  mode = "local"
)
#> $out_degree_seq
#> numeric(0)
#> 
#> $in_degree_seq
#>   [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#>  [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#>  [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
```

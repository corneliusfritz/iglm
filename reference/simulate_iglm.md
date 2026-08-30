# Simulate Responses and Connections

Simulate responses and connections.

## Usage

``` r
simulate_iglm(
  formula,
  basis = NULL,
  coef,
  coef_degrees = NULL,
  sampler = NULL,
  only_stats = TRUE,
  display_progress = FALSE,
  offset_nonoverlap = 0,
  cluster = NULL,
  fix_x = NULL,
  fix_z = NULL
)
```

## Arguments

- formula:

  A model \`formula\` object. The left-hand side should be the name of a
  \`iglm.data\` object available in the calling environment. See
  [`iglm-terms`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  for details on specifying the right-hand side terms.

- basis:

  An optional
  [`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
  object to serve as the basis (initial state) for the simulation. If
  provided, the simulation starts from the state (`x_attribute`,
  `y_attribute`, `z_network`) defined in this object. The \`basis\`
  object must be consistent with the model data referenced in
  \`formula\` (same number of actors, directedness, and attribute
  types). All structural specifications (terms, neighborhood, overlap,
  and fixed variable flags) are taken from the \`formula\` data object.
  If \`NULL\` (default), the initial state is taken from the
  \`iglm.data\` object referenced in the \`formula\`.

- coef:

  Numeric vector containing the coefficient values for the structural
  (non-degrees) terms defined in the \`formula\`.

- coef_degrees:

  Numeric vector specifying the degrees coefficient values
  (expansiveness/attractiveness). This is required **only if** the
  \`formula\` includes degrees terms. Its length must be \`n_actor\`
  (for undirected networks) or \`2 \* n_actor\` (for directed networks),
  where \`n_actor\` is determined from the \`iglm.data\` object in the
  formula.

- sampler:

  An object of class \`sampler.iglm\` (created by \`sampler.iglm()\`)
  specifying the MCMC sampling parameters. This includes the number of
  simulations (\`n_simulation\`), burn-in iterations (\`n_burn_in\`),
  initialization settings (\`init_empty\`), and component sampler
  settings (\`sampler_x\`, \`sampler_y\`, etc.). If \`NULL\` (default),
  default settings from \`sampler.iglm()\` are used.

- only_stats:

  (logical). If `TRUE` (default, consistent with the usage signature),
  the function returns only the matrix of features calculated for each
  simulation. The full simulated `iglm.data` objects are discarded to
  minimize memory usage. If `FALSE`, the complete simulated `iglm.data`
  objects are created and returned within the `samples` component of the
  output list.

- display_progress:

  Logical. If \`TRUE\`, progress messages or a progress bar (depending
  on the backend implementation) are displayed during simulation.
  Default is \`FALSE\`.

- offset_nonoverlap:

  Numeric scalar value passed to the C++ simulator. This value is
  typically added to the linear predictor for dyads that are **not**
  part of the 'overlap' set defined in the \`iglm.data\` object,
  potentially modifying tie probabilities outside the primary
  neighborhood. Default is \`0\`.

- cluster:

  Optional parallel cluster object created, for example, by
  [`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html).
  If provided and valid, the function performs a single burn-in
  simulation on the main R process, then distributes the remaining
  \`n_simulation\` tasks across the cluster workers using
  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html).
  The master seed is offset for each worker to ensure different random
  streams. If \`NULL\` (default), all simulations are run sequentially
  in the main R process.

- fix_x:

  Optional logical override indicating if \`x_attribute\` should be held
  fixed during simulation. If \`NULL\` (default), the setting from the
  model \`iglm.data\` object is used.

- fix_z:

  Optional logical override indicating if \`z_network\` should be held
  fixed during simulation. If \`NULL\` (default), the setting from the
  model \`iglm.data\` object is used.

## Value

A list containing one or two components (depending on \`only_stats\`):

- \`samples\`:

  If \`only_stats = FALSE\`, this is a list of length
  \`sampler\$n_simulation\` where each element is an
  [`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
  object representing one simulated draw from the model. The list has
  the S3 class \`"iglm.data.list"\`. If \`only_stats = TRUE\`, this
  component is omitted.

- \`stats\`:

  A numeric matrix with \`sampler\$n_simulation\` rows and
  \`length(coef)\` columns. Each row contains the features
  (corresponding to the model terms in \`formula\`) calculated for one
  simulation draw. Column names are set to match the term names.

## Details

**Parallel Execution:** When a \`cluster\` object is provided, the
simulation process is adapted:

1.  A single simulation run (including burn-in specified by
    \`sampler\$n_burn_in\`) is performed on the master node to obtain a
    starting state for the parallel chains.

2.  The total number of requested simulations
    (\`sampler\$n_simulation\`) is divided among the cluster workers.

3.  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html)
    is used to run simulations on each worker. Each worker starts from
    the state obtained after the initial burn-in, performs **zero**
    additional burn-in (\`n_burn_in = 0\` passed to workers), and
    generates its assigned share of the simulations. Component sampler
    seeds are offset based on the worker ID to ensure pseudo-independent
    random number streams.

4.  Results (simulated objects or statistics) from all workers are
    collected and combined.

This approach ensures that the initial burn-in phase happens only once,
saving time.

## Errors

The function stops with an error if:

- The length of \`coef\` does not match the number of terms derived from
  \`formula\`.

- Formula preprocessing fails.

- The \`sampler\` object is not of class
  [`sampler.iglm`](https://corneliusfritz.github.io/iglm/reference/sampler.iglm.md).

- The C++ backend \`xyz_simulate_cpp\` encounters an error.

## See also

[`iglm`](https://corneliusfritz.github.io/iglm/reference/iglm.md) for
creating the model object,
[`sampler.iglm`](https://corneliusfritz.github.io/iglm/reference/sampler.iglm.md)
for creating the sampler object,
[`iglm.data`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
for the data object structure.

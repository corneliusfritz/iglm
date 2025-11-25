# Simulate responses and connections

Simulate responses and connections.

## Usage

``` r
# S3 method for class 'iglm'
simulate(
  formula,
  coef,
  coef_popularity = NULL,
  sampler = NULL,
  only_stats = TRUE,
  display_progress = FALSE,
  offset_nonoverlap = 0,
  cluster = NULL,
  fix_x = FALSE
)
```

## Arguments

- formula:

  A model \`formula\` object. The left-hand side should be the name of a
  \`iglm.data\` object available in the calling environment. See
  `model_terms` for details on specifying the right-hand side terms.

- coef:

  Numeric vector containing the coefficient values for the structural
  (non-popularity) terms defined in the \`formula\`.

- coef_popularity:

  Numeric vector specifying the popularity coefficient values
  (expansiveness/attractiveness). This is required **only if** the
  \`formula\` includes popularity terms. Its length must be \`n_actor\`
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
  “parallel::makeCluster“. If provided and valid, the function performs
  a single burn-in simulation on the main R process, then distributes
  the remaining \`n_simulation\` tasks across the cluster workers using
  “parallel::parLapply“. Seeds for component samplers are offset for
  each worker to ensure different random streams. If \`NULL\` (default),
  all simulations are run sequentially in the main R process.

- fix_x:

  Logical. If \`TRUE\`, the simulation holds the \`x_attribute\` fixed
  at its initial state (from the [`iglm.data`](iglm.data.md) object) and
  only simulates the \`y_attribute\` and \`z_network\`. If \`FALSE\`
  (default), all components (x, y, z) are simulated according to the
  model and sampler settings.

## Value

A list containing two components:

- \`samples\`:

  If \`only_stats = FALSE\`, this is a list of length
  \`sampler\$n_simulation\` where each element is a \`iglm.data\` object
  representing one simulated draw from the model. The list has the S3
  class \`"iglm.data.list"\`. If \`only_stats = TRUE\`, this is
  typically an empty list.

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

3.  “parallel::parLapply“ is used to run simulations on each worker.
    Each worker starts from the state obtained after the initial
    burn-in, performs **zero** additional burn-in (\`n_burn_in = 0\`
    passed to workers), and generates its assigned share of the
    simulations. Component sampler seeds are offset based on the worker
    ID to ensure pseudo-independent random number streams.

4.  Results (simulated objects or statistics) from all workers are
    collected and combined.

This approach ensures that the initial burn-in phase happens only once,
saving time.

## Errors

The function stops with an error if:

- The length of \`coef\` does not match the number of terms derived from
  \`formula\`.

- \`formula_preprocess\` fails.

- The \`sampler\` object is not of class \`sampler.iglm\`.

- The C++ backend \`xyz_simulate_cpp\` encounters an error.

- Helper functions like \`XYZ_to_R\` or \`is_cluster_active\` are not
  found.

Warnings may be issued if default sampler settings are used.

## See also

`iglm` for creating the model object, `sampler.iglm` for creating the
sampler object, `iglm.data` for the data object structure.

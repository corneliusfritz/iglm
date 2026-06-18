# iglm Sampler Settings (R6 Class)

The \`sampler.iglm\` class is an R6 container for specifying and storing
the parameters that control the MCMC (Markov Chain Monte Carlo) sampling
process used in
[`iglm`](https://corneliusfritz.github.io/iglm/reference/iglm.md)
simulations and potentially during estimation. It includes settings for
the number of simulations, burn-in period, initialization, and
parallelization options. It also holds references to component samplers
([`sampler.net.attr`](https://corneliusfritz.github.io/iglm/reference/sampler.net.attr.md)
objects) responsible for sampling individual parts (attributes x, y,
network z).

## Active bindings

- `sampler_x`:

  (\`sampler_net_attr\`) The sampler configuration object for the x
  attribute.

- `sampler_y`:

  (\`sampler_net_attr\`) The sampler configuration object for the y
  attribute.

- `sampler_z`:

  (\`sampler_net_attr\`) The sampler configuration object for the z
  network (overlap region).

- `n_simulation`:

  (\`integer\`) The number of configurations to simulate.

- `n_burn_in`:

  (\`integer\`) The number of initial MCMC iterations to discard.

- `init_empty`:

  (\`logical\`) Whether to initialize simulations from an empty state.

- `seed`:

  (\`integer\`) Read-only. The random seed used for sampling.

- `cluster`:

  (\`cluster\` object or \`NULL\`) The parallel cluster object being
  used, or \`NULL\`.

## Methods

### Public methods

- [`sampler.iglm$new()`](#method-sampler.iglm-initialize)

- [`sampler.iglm$set_cluster()`](#method-sampler.iglm-set_cluster)

- [`sampler.iglm$deactive_cluster()`](#method-sampler.iglm-deactive_cluster)

- [`sampler.iglm$set_n_simulation()`](#method-sampler.iglm-set_n_simulation)

- [`sampler.iglm$set_n_burn_in()`](#method-sampler.iglm-set_n_burn_in)

- [`sampler.iglm$set_init_empty()`](#method-sampler.iglm-set_init_empty)

- [`sampler.iglm$set_x_sampler()`](#method-sampler.iglm-set_x_sampler)

- [`sampler.iglm$set_y_sampler()`](#method-sampler.iglm-set_y_sampler)

- [`sampler.iglm$set_z_sampler()`](#method-sampler.iglm-set_z_sampler)

- [`sampler.iglm$set_seed()`](#method-sampler.iglm-set_seed)

- [`sampler.iglm$print()`](#method-sampler.iglm-print)

- [`sampler.iglm$gather()`](#method-sampler.iglm-gather)

- [`sampler.iglm$save()`](#method-sampler.iglm-save)

- [`sampler.iglm$clone()`](#method-sampler.iglm-clone)

------------------------------------------------------------------------

### `sampler.iglm$new()`

Create a new \`sampler.iglm\` object. Initializes all sampler settings,
using defaults for component samplers (\`sampler.net.attr\`) if not
provided, and validates inputs.

#### Usage

    sampler.iglm$new(
      sampler_x = NULL,
      sampler_y = NULL,
      sampler_z = NULL,
      n_simulation = 100,
      n_burn_in = 10,
      init_empty = TRUE,
      seed = NA,
      cluster = NULL,
      file = NULL
    )

#### Arguments

- `sampler_x`:

  An object of class \`sampler.net.attr\` controlling sampling for the x
  attribute. If \`NULL\`, defaults from \`sampler.net.attr()\` are used.

- `sampler_y`:

  An object of class \`sampler.net.attr\` controlling sampling for the y
  attribute. If \`NULL\`, defaults from \`sampler.net.attr()\` are used.

- `sampler_z`:

  An object of class \`sampler.net.attr\` controlling sampling for the z
  network (within the defined neighborhood/overlap). If \`NULL\`,
  defaults from \`sampler.net.attr()\` are used.

- `n_simulation`:

  (integer) The number of network/attribute configurations to simulate
  and store after the burn-in period. Default is 100. Must be
  non-negative.

- `n_burn_in`:

  (integer) The number of initial MCMC iterations to discard (burn-in)
  before starting to collect simulations. Default is 10. Must be
  non-negative.

- `init_empty`:

  (logical) If \`TRUE\` (default), the MCMC chain is initialized from an
  empty state (e.g., empty network, attributes at mean). If \`FALSE\`,
  initialization might depend on the specific sampler implementation
  (e.g., starting from observed data).

- `seed`:

  (integer or \`NA\`) A single integer seed for the random number
  generator, set once before sampling begins. If \`NA\` (default), a
  random seed is generated automatically.

- `cluster`:

  A parallel cluster object (e.g., from the \`parallel\` package) to use
  for running simulations in parallel. If \`NULL\` (default),
  simulations are run sequentially.

- `file`:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

#### Returns

A new \`sampler.iglm\` object.

------------------------------------------------------------------------

### `sampler.iglm$set_cluster()`

Sets the parallel cluster object to be used for simulations.

#### Usage

    sampler.iglm$set_cluster(cluster)

#### Arguments

- `cluster`:

  A parallel cluster object from the \`parallel\` package.

------------------------------------------------------------------------

### `sampler.iglm$deactive_cluster()`

Deactivates parallel processing for this sampler instance by setting the
internal cluster object reference to \`NULL\`.

#### Usage

    sampler.iglm$deactive_cluster()

#### Returns

The \`sampler.iglm\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `sampler.iglm$set_n_simulation()`

Sets the number of simulations to generate after burn-in.

#### Usage

    sampler.iglm$set_n_simulation(n_simulation)

#### Arguments

- `n_simulation`:

  (integer) The number of simulations to set.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_n_burn_in()`

Sets the number of burn-in iterations.

#### Usage

    sampler.iglm$set_n_burn_in(n_burn_in)

#### Arguments

- `n_burn_in`:

  (integer) The number of burn-in iterations to set.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_init_empty()`

Sets whether to initialize simulations from an empty state.

#### Usage

    sampler.iglm$set_init_empty(init_empty)

#### Arguments

- `init_empty`:

  (logical) \`TRUE\` to initialize from empty, \`FALSE\` otherwise.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_x_sampler()`

Sets the sampler configuration for the x attribute.

#### Usage

    sampler.iglm$set_x_sampler(sampler_x)

#### Arguments

- `sampler_x`:

  An object of class \`sampler_net_attr\`.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_y_sampler()`

Sets the sampler configuration for the y attribute.

#### Usage

    sampler.iglm$set_y_sampler(sampler_y)

#### Arguments

- `sampler_y`:

  An object of class \`sampler_net_attr\`.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_z_sampler()`

Sets the sampler configuration for the z attribute.

#### Usage

    sampler.iglm$set_z_sampler(sampler_z)

#### Arguments

- `sampler_z`:

  An object of class \`sampler_net_attr\`.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$set_seed()`

Sets the random seed for this sampler.

#### Usage

    sampler.iglm$set_seed(seed)

#### Arguments

- `seed`:

  (integer) The random seed to set.

#### Returns

None.

------------------------------------------------------------------------

### `sampler.iglm$print()`

Prints a formatted summary of the sampler configuration to the console.

#### Usage

    sampler.iglm$print(digits = 3, ...)

#### Arguments

- `digits`:

  (integer) Number of digits for formatting numeric values. Default: 3.

- `...`:

  Additional arguments (currently ignored).

#### Returns

The \`sampler.iglm\` object itself (\`self\`), invisibly.

------------------------------------------------------------------------

### `sampler.iglm$gather()`

Gathers all data from private fields into a list.

#### Usage

    sampler.iglm$gather()

#### Returns

A list containing all information of the sampler.

------------------------------------------------------------------------

### `sampler.iglm$save()`

Save the object's complete state to a directory. This will save the main
sampler's settings to a file named 'sampler.iglm_state.rds' within the
specified directory, and will also call the \`save()\` method for each
nested sampler (.x, .y, .z), saving them into the same directory.

#### Usage

    sampler.iglm$save(file)

#### Arguments

- `file`:

  (character) The file to a directory where the state files will be
  saved. The directory will be created if it does not exist.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### `sampler.iglm$clone()`

The objects of this class are cloneable with this method.

#### Usage

    sampler.iglm$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

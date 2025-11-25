# R6 Class for Single Component Sampler Settings

The \`sampler_net_attr\` class is a simple R6 container used within the
\`sampler.iglm\` class. It holds the MCMC sampling parameters for a
single component of the \`iglm\` model, such as one attribute (e.g.,
\`x_attribute\`) or a part of the network (e.g., \`z_network\` within
the overlap). It primarily stores the number of proposals and a random
seed.

## Active bindings

- `n_proposals`:

  (\`integer\`) Read-only. The number of MCMC proposals per sampling
  step.

- `seed`:

  (\`integer\`) Read-only. The random seed used for this component's
  sampler.

## Methods

### Public methods

- [`sampler.net.attr.generator$new()`](#method-sampler.net.attr-new)

- [`sampler.net.attr.generator$print()`](#method-sampler.net.attr-print)

- [`sampler.net.attr.generator$gather()`](#method-sampler.net.attr-gather)

- [`sampler.net.attr.generator$set_n_proposals()`](#method-sampler.net.attr-set_n_proposals)

- [`sampler.net.attr.generator$set_seed()`](#method-sampler.net.attr-set_seed)

- [`sampler.net.attr.generator$save()`](#method-sampler.net.attr-save)

- [`sampler.net.attr.generator$clone()`](#method-sampler.net.attr-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \`sampler_net_attr\` object. Validates inputs and sets a
random seed if none is provided.

#### Usage

    sampler.net.attr.generator$new(n_proposals = 10000, seed = NA, file = NULL)

#### Arguments

- `n_proposals`:

  (integer) The number of MCMC proposals (iterations) to perform for
  this specific component during each sampling step. Default is 10000.
  Must be a non-negative integer.

- `seed`:

  (integer or \`NA\`) An integer seed for the random number generator to
  ensure reproducibility for this component's sampling. If \`NA\`
  (default), a random seed is generated automatically.

- `file`:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

#### Returns

A new \`sampler_net_attr\` object.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a summary of the sampler settings for this component.

#### Usage

    sampler.net.attr.generator$print(indent = "  ")

#### Arguments

- `indent`:

  (character) A string used for indentation (e.g., spaces) when
  printing, useful for nested structures. Default is " ".

#### Returns

The object itself, invisibly. Called for side effect.

------------------------------------------------------------------------

### Method `gather()`

Gathers all data from private fields into a list.

#### Usage

    sampler.net.attr.generator$gather()

#### Returns

A list containing all information of the sampler.

------------------------------------------------------------------------

### Method `set_n_proposals()`

Sets the number of MCMC proposals for this component.

#### Usage

    sampler.net.attr.generator$set_n_proposals(n_proposals)

#### Arguments

- `n_proposals`:

  (integer) The number of proposals to set.

#### Returns

None.

------------------------------------------------------------------------

### Method `set_seed()`

Sets the random seed for this component's sampler.

#### Usage

    sampler.net.attr.generator$set_seed(seed)

#### Arguments

- `seed`:

  (integer) The random seed to set.

#### Returns

None.

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Save the object's state to an .rds file.

#### Usage

    sampler.net.attr.generator$save(file)

#### Arguments

- `file`:

  (character) The file file where the state will be saved.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    sampler.net.attr.generator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

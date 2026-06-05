# Single Component Sampler Settings (R6 Class)

The \`sampler_net_attr\` class is a simple R6 container used within the
\`sampler.iglm\` class. It holds the MCMC sampling parameters for a
single component of the \`iglm\` model, such as one attribute (e.g.,
\`x_attribute\`) or a part of the network (e.g., \`z_network\` within
the overlap). It stores the number of proposals and the TNT flag. The
random seed is managed centrally by the parent \`sampler.iglm\` object.

## Active bindings

- `n_proposals`:

  (\`integer\`) Read-only. Number of MCMC proposals per step.

- `tnt`:

  (\`logical\`) Read-only. Whether TNT sampling is used.

## Methods

### Public methods

- [`sampler.net.attr.generator$new()`](#method-sampler.net.attr-new)

- [`sampler.net.attr.generator$print()`](#method-sampler.net.attr-print)

- [`sampler.net.attr.generator$gather()`](#method-sampler.net.attr-gather)

- [`sampler.net.attr.generator$set_n_proposals()`](#method-sampler.net.attr-set_n_proposals)

- [`sampler.net.attr.generator$set_tnt()`](#method-sampler.net.attr-set_tnt)

- [`sampler.net.attr.generator$save()`](#method-sampler.net.attr-save)

- [`sampler.net.attr.generator$clone()`](#method-sampler.net.attr-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \`sampler_net_attr\` object.

#### Usage

    sampler.net.attr.generator$new(n_proposals = 10000, file = NULL, tnt = TRUE)

#### Arguments

- `n_proposals`:

  (integer) The number of MCMC proposals (iterations) to perform for
  this specific component during each sampling step. Default is 10000.
  Must be a non-negative integer.

- `file`:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

- `tnt`:

  (logical) If \`TRUE\` (default), use Tie-No-Tie sampling (only if used
  for networks).

#### Returns

A new \`sampler_net_attr\` object.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a summary of the sampler settings for this component.

#### Usage

    sampler.net.attr.generator$print(indent = "  ")

#### Arguments

- `indent`:

  (character) Indentation string. Default is " ".

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `gather()`

Gathers all data into a list.

#### Usage

    sampler.net.attr.generator$gather()

#### Returns

A list with \`n_proposals\` and \`tnt\`.

------------------------------------------------------------------------

### Method `set_n_proposals()`

Sets the number of MCMC proposals.

#### Usage

    sampler.net.attr.generator$set_n_proposals(n_proposals)

#### Arguments

- `n_proposals`:

  (integer) Number of proposals.

------------------------------------------------------------------------

### Method `set_tnt()`

Sets whether to use TNT sampling.

#### Usage

    sampler.net.attr.generator$set_tnt(tnt)

#### Arguments

- `tnt`:

  (logical) \`TRUE\` to use TNT sampling.

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Save state to an .rds file.

#### Usage

    sampler.net.attr.generator$save(file)

#### Arguments

- `file`:

  (character) File path.

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

# An R6 class for Network GLM (Generalized Linear Model) Objects

The \`iglm_object\` class encapsulates all components required to
define, estimate, and simulate from a network generalized linear model.
This includes the model formula, coefficients, the underlying network
and attribute data (via a \`iglm.data\` object), sampler controls,
estimation controls, and storage for results.

## References

Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

Stewart, J. R. and M. Schweinberger (2025). Pseudo-Likelihood-Based
M-Estimation of Random Graphs with Dependent Edges and Parameter Vectors
of Increasing Dimension. Annals of Statistics, to appear.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random
graph models: characterization, properties, and statistical inference.
Journal of the Royal Statistical Society, Series B (Statistical
Methodology), 7, 647-676.

## Active bindings

- `formula`:

  (\`formula\`) Read-only. The model formula specifying terms and data
  object.

- `coef`:

  (\`numeric\`) Read-only. The current vector of non-popularity
  coefficient estimates or initial values.

- `coef_popularity`:

  (\`numeric\` or \`NULL\`) Read-only. The current vector of popularity
  coefficient estimates or initial values, or \`NULL\` if not
  applicable.

- `results`:

  (\`results\`) Read-only. The [`results`](results.md) R6 object
  containing all estimation and simulation outputs.

- `iglm.data`:

  (\`iglm.data\`) Read-only. The associated [`iglm.data`](iglm.data.md)
  R6 object containing the network and attribute data.

- `control`:

  (\`control.iglm\`) Read-only. The [`control.iglm`](control.iglm.md)
  object specifying estimation parameters.

- `sampler`:

  (\`sampler_iglm\`) Read-only. The
  [`sampler_iglm`](sampler_iglm_generator.md) object specifying MCMC
  sampling parameters.

- `count_statistics`:

  (\`numeric\`) Read-only. A named vector of the observed network
  statistics corresponding to the model terms, calculated on the current
  \`iglm.data\` data.

## Methods

### Public methods

- [`iglm_object_generator$new()`](#method-iglm_object-new)

- [`iglm_object_generator$model_assessment()`](#method-iglm_object-model_assessment)

- [`iglm_object_generator$print()`](#method-iglm_object-print)

- [`iglm_object_generator$plot()`](#method-iglm_object-plot)

- [`iglm_object_generator$gather()`](#method-iglm_object-gather)

- [`iglm_object_generator$save()`](#method-iglm_object-save)

- [`iglm_object_generator$estimate()`](#method-iglm_object-estimate)

- [`iglm_object_generator$summary()`](#method-iglm_object-summary)

- [`iglm_object_generator$simulate()`](#method-iglm_object-simulate)

- [`iglm_object_generator$get_samples()`](#method-iglm_object-get_samples)

- [`iglm_object_generator$set_sampler()`](#method-iglm_object-set_sampler)

- [`iglm_object_generator$set_target()`](#method-iglm_object-set_target)

- [`iglm_object_generator$clone()`](#method-iglm_object-clone)

------------------------------------------------------------------------

### Method `new()`

Internal method to calculate the observed count statistics based on the
model formula and the data in the \`iglm.data\` object. Populates the
\`private\$.count_statistics\` field.

Internal validation method. Checks the consistency and validity of all
components of the \`iglm_object\`. Stops with an error if any check
fails.

Creates a new \`iglm_object\`. This involves parsing the formula,
linking the data object, initializing coefficients, setting up sampler
and control objects, calculating initial statistics, and validating.

#### Usage

    iglm_object_generator$new(
      formula = NULL,
      coef = NULL,
      coef_popularity = NULL,
      sampler = NULL,
      control = NULL,
      file = NULL
    )

#### Arguments

- `formula`:

  A model \`formula\` object. The left-hand side should be the name of a
  [`iglm.data`](iglm.data.md) object available in the calling
  environment. See [`model_terms`](model_terms.md) for details on
  specifying the right-hand side terms.

- `coef`:

  A numeric vector of initial coefficients for the terms in the formula
  (excluding popularity). If \`NULL\`, coefficients are initialized to
  zero.

- `coef_popularity`:

  An optional numeric vector of initial popularity coefficients. Should
  be \`NULL\` if the formula does not include popularity terms.

- `sampler`:

  A [`sampler_iglm`](sampler_iglm_generator.md) object specifying the
  MCMC sampler settings. If \`NULL\`, default settings are used.

- `control`:

  A [`control.iglm`](control.iglm.md) object specifying estimation
  control parameters. If \`NULL\`, default settings are used.

- `file`:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

#### Returns

A new \`iglm_object\`.

------------------------------------------------------------------------

### Method `model_assessment()`

Performs model assessment by calculating specified network statistics on
the observed network and comparing their distribution to the
distribution obtained from simulated networks based on the current model
parameters. Requires simulations to have been run first (via
`iglm_object$simulate` or `iglm_object_generator$estimate`).

#### Usage

    iglm_object_generator$model_assessment(formula)

#### Arguments

- `formula`:

  A formula specifying the network statistics to assess (e.g., \`~
  degree_distribution() + geodesic_distances_distribution()\`). The
  terms should correspond to methods available in the
  [`iglm.data`](iglm.data.md) object that end with \`distributions\`.

#### Returns

An object of class \`iglm_model_assessment\` containing the observed
statistics and the distribution of simulated statistics. The result is
also stored internally.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a summary of the \`iglm_object\`. If estimation results are
available, they are printed in a standard coefficient table format.

#### Usage

    iglm_object_generator$print(digits = 4, ...)

#### Arguments

- `digits`:

  (integer) Number of digits for rounding numeric output.

- `...`:

  Additional arguments (not used).

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot the estimation results, including coefficient convergence paths and
model assessment diagnostics if available.

#### Usage

    iglm_object_generator$plot(
      stats = FALSE,
      trace = FALSE,
      model_assessment = FALSE
    )

#### Arguments

- `stats`:

  (logical) If \`TRUE\`, plot the observed vs. simulated statistics from
  model assessment. Default is \`FALSE\`.

- `trace`:

  (logical) If \`TRUE\`, plot the coefficient convergence paths. Default
  is \`FALSE\`.

- `model_assessment`:

  (logical) If \`TRUE\`, plot diagnostics from the model assessment (if
  already carried out). Default is \`FALSE\`.

------------------------------------------------------------------------

### Method `gather()`

Gathers all components of the [`iglm_object`](iglm.md) into a single
list for easy saving or inspection.

#### Usage

    iglm_object_generator$gather()

#### Returns

A list containing all key components of the [`iglm_object`](iglm.md).
This includes the formula, coefficients, sampler, control settings,
preprocessing info, time taken for estimation, count statistics,
results, and the underlying [`iglm.data`](iglm.data.md) data object.

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Save the [`iglm_object`](iglm.md) to a file in RDS format.

#### Usage

    iglm_object_generator$save(file = NULL)

#### Arguments

- `file`:

  (character) File path to save the object to.

#### Returns

Invisibly returns \`NULL\`.

------------------------------------------------------------------------

### Method `estimate()`

Estimate the model parameters using the specified control settings.
Stores the results internally and updates the coefficient fields.

#### Usage

    iglm_object_generator$estimate()

#### Returns

If the no preprocessing should be returned (as per control settings),
this function returns a list containing detailed estimation results,
invisibly. Includes final coefficients, variance-covariance matrix,
convergence path, Fisher information, score vector, log-likelihood, and
any simulations performed during estimation. Else, the function returns
a list of the desired preprocessed data (as a data.frame) and needed
time.

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Provides a summary of the estimation results. Requires the model to have
been estimated first.

#### Usage

    iglm_object_generator$summary(digits = 3)

#### Arguments

- `digits`:

  (integer) Number of digits for rounding numeric output.

#### Returns

Prints the summary to the console and returns \`NULL\` invisibly.

------------------------------------------------------------------------

### Method [`simulate()`](https://rdrr.io/r/stats/simulate.html)

Simulate networks from the fitted model or a specified model. Stores the
simulations and/or summary statistics internally. The simulation is
carried out using the internal MCMC sampler described in
[`simulate_iglm`](simulate_iglm.md).

#### Usage

    iglm_object_generator$simulate(
      nsim = 1,
      only_stats = FALSE,
      display_progress = TRUE,
      offset_nonoverlap = 0
    )

#### Arguments

- `nsim`:

  (integer) Number of networks to simulate. Default is 1.

- `only_stats`:

  (logical) If \`TRUE\`, only calculate and store summary statistics for
  each simulation, discarding the network object itself. Default is
  \`FALSE\`.

- `display_progress`:

  (logical) If \`TRUE\` (default), display a progress bar during
  simulation.

- `offset_nonoverlap`:

  (numeric) Offset to apply for non-overlapping dyads during simulation
  (if applicable to the sampler). This option is useful if the sparsity
  of edges of units with non-overlapping neighborhoods is known. Default
  is 0.

#### Returns

A list containing the simulated networks (\`samples\`, as a
\`iglm.data.list\` if \`only_stats = FALSE\`) and/or their summary
statistics (\`stats\`), invisibly.

------------------------------------------------------------------------

### Method `get_samples()`

Retrieve the simulated networks stored in the object. Requires
`simulate` or `estimate` to have been run first.

#### Usage

    iglm_object_generator$get_samples()

#### Returns

A list of [`iglm.data`](iglm.data.md) objects representing the simulated
networks, invisibly. Returns an error if no samples are available.

------------------------------------------------------------------------

### Method `set_sampler()`

Replace the internal MCMC sampler with a new one. This is useful for
changing the sampling scheme without redefining the entire model.

#### Usage

    iglm_object_generator$set_sampler(sampler)

#### Arguments

- `sampler`:

  A [`sampler_iglm`](sampler_iglm_generator.md) object. @return The
  [`iglm_object`](iglm.md) itself, invisibly.

------------------------------------------------------------------------

### Method `set_target()`

Replace the internal \`iglm.data\` data object with a new one. This is
useful for applying a fitted model to new observed data. Recalculates
count statistics and re-validates the object.

#### Usage

    iglm_object_generator$set_target(x)

#### Arguments

- `x`:

  A [`iglm.data`](iglm.data.md) “ object containing the new observed
  data.

#### Returns

The [`iglm_object`](iglm.md) itself, invisibly.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    iglm_object_generator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

# An R6 class for Network GLM (Generalized Linear Model) Objects

The \`iglm.object\` class encapsulates all components required to
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

  (\`sampler.iglm\`) Read-only. The [`sampler.iglm`](sampler.iglm.md)
  object specifying MCMC sampling parameters.

- `sufficient_statistics`:

  (\`numeric\`) Read-only. A named vector of the observed network
  statistics corresponding to the model terms, calculated on the current
  \`iglm.data\` data.

## Methods

### Public methods

- [`iglm.object.generator$new()`](#method-iglm.object-new)

- [`iglm.object.generator$model_assessment()`](#method-iglm.object-model_assessment)

- [`iglm.object.generator$print()`](#method-iglm.object-print)

- [`iglm.object.generator$plot()`](#method-iglm.object-plot)

- [`iglm.object.generator$gather()`](#method-iglm.object-gather)

- [`iglm.object.generator$save()`](#method-iglm.object-save)

- [`iglm.object.generator$estimate()`](#method-iglm.object-estimate)

- [`iglm.object.generator$summary()`](#method-iglm.object-summary)

- [`iglm.object.generator$simulate()`](#method-iglm.object-simulate)

- [`iglm.object.generator$predict()`](#method-iglm.object-predict)

- [`iglm.object.generator$set_coefficients()`](#method-iglm.object-set_coefficients)

- [`iglm.object.generator$get_samples()`](#method-iglm.object-get_samples)

- [`iglm.object.generator$set_sampler()`](#method-iglm.object-set_sampler)

- [`iglm.object.generator$set_target()`](#method-iglm.object-set_target)

- [`iglm.object.generator$clone()`](#method-iglm.object-clone)

------------------------------------------------------------------------

### Method `new()`

Internal method to calculate the observed count statistics based on the
model formula and the data in the \`iglm.data\` object. Populates the
\`private\$.sufficient_statistics\` field.

Internal validation method. Checks the consistency and validity of all
components of the \`iglm.object\`. Stops with an error if any check
fails.

Creates a new \`iglm.object\`. This involves parsing the formula,
linking the data object, initializing coefficients, setting up sampler
and control objects, calculating initial statistics, and validating.

#### Usage

    iglm.object.generator$new(
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
  environment. See [`model.terms`](model.terms.md) for details on
  specifying the right-hand side terms.

- `coef`:

  A numeric vector of initial coefficients for the terms in the formula
  (excluding popularity). If \`NULL\`, coefficients are initialized to
  zero.

- `coef_popularity`:

  An optional numeric vector of initial popularity coefficients. Should
  be \`NULL\` if the formula does not include popularity terms.

- `sampler`:

  A [`sampler.iglm`](sampler.iglm.md) object specifying the MCMC sampler
  settings. If \`NULL\`, default settings are used.

- `control`:

  A [`control.iglm`](control.iglm.md) object specifying estimation
  control parameters. If \`NULL\`, default settings are used.

- `file`:

  (character or \`NULL\`) If provided, loads the sampler state from the
  specified .rds file instead of initializing from parameters.

#### Returns

A new \`iglm.object\`.

------------------------------------------------------------------------

### Method `model_assessment()`

Performs model assessment by calculating specified network statistics on
the observed network and comparing their distribution to the
distribution obtained from simulated networks based on the current model
parameters. Requires simulations to have been run first (via
`iglm.object$simulate` or `iglm.object_generator$estimate`).

#### Usage

    iglm.object.generator$model_assessment(formula, plot = TRUE)

#### Arguments

- `formula`:

  A formula specifying the network statistics to assess (e.g., \`~
  degree_distribution() + geodesic_distances_distribution()\`). The
  terms should correspond to methods available in the
  [`iglm.data`](iglm.data.md) object that end with \`distributions\`. If
  the term mcmc_diagnostics is included, MCMC diagnostics will also be
  computed.

- `plot`:

  (logical) If \`TRUE\`, generates plots comparing observed and
  simulated statistics. Default is \`TRUE\`.

#### Returns

An object of class \`iglm_model_assessment\` containing the observed
statistics and the distribution of simulated statistics. The result is
also stored internally.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a summary of the \`iglm.object\`. If estimation results are
available, they are printed in a standard coefficient table format.

#### Usage

    iglm.object.generator$print(digits = 4, ...)

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

    iglm.object.generator$plot(
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

Gathers all components of the [`iglm.object`](iglm.md) into a single
list for easy saving or inspection.

#### Usage

    iglm.object.generator$gather()

#### Returns

A list containing all key components of the [`iglm.object`](iglm.md).
This includes the formula, coefficients, sampler, control settings,
preprocessing info, time taken for estimation, count statistics,
results, and the underlying [`iglm.data`](iglm.data.md) data object.

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Save the [`iglm.object`](iglm.md) to a file in RDS format.

#### Usage

    iglm.object.generator$save(file = NULL)

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

    iglm.object.generator$estimate()

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

    iglm.object.generator$summary(digits = 3)

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

    iglm.object.generator$simulate(
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

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Calculates predicted values for the nodal covariates (`x`), the outcome
variable (`y`), and the network structure (`z`). The function supports
two prediction modes: *marginal* (based on Monte Carlo integration over
simulated samples) and *conditional* (based on the analytical linear
predictor and point estimates).

#### Usage

    iglm.object.generator$predict(
      variant = c("conditional", "marginal"),
      type = c("x", "y", "z")
    )

#### Arguments

- `variant`:

  A character string specifying the type of prediction to generate. Must
  be one of:

  - `"marginal"`: Computes predictions by aggregating over the MCMC
    samples stored in the internal results. If samples do not exist,
    `self$simulate()` is triggered automatically. This represents the
    expectation integrated over the uncertainty of the latent process.

  - `"conditional"`: Computes predictions using the systematic component
    of the Generalized Linear Model (GLM). It calculates the linear
    predictor \\\eta = X\beta\\ (plus offset and popularity terms for
    the network) and applies the inverse link function \\\mu =
    g^{-1}(\eta)\\.

  Defaults to `c("conditional", "marginal")`.

- `type`:

  A character vector indicating which components to predict. Options
  are:

  - `"x"`: Nodal covariates.

  - `"y"`: Nodal outcome variable.

  - `"z"`: Dyadic network structure (interaction probabilities).

  Defaults to `c("x", "y", "z")`.

#### Details

**Marginal Predictions:** When `variant = "marginal"`, the function
approximates the expected value via Monte Carlo integration:
\$\$\hat{\mu} = \frac{1}{S} \sum\_{s=1}^{S} y^{(s)}\$\$ where
\\y^{(s)}\\ are the realized values from the \\s\\-th simulation sample.
For the network `z`, this results in an edge probability matrix averaged
over all sampled networks.

**Conditional Predictions:** When `variant = "conditional"`, the
function calculates the theoretical mean \\\mu\\ based on the estimated
coefficients \\\hat{\theta}\\:

- For **Binomial** families: \\\mu = (1 + \exp(-\eta))^{-1}\\
  (Logistic).

- For **Poisson** families: \\\mu = \exp(\eta)\\ (Exponential).

- For **Gaussian** families: \\\mu = \eta\\ (Identity).

For the network component `z`, the linear predictor includes dyadic
covariates, popularity effects (sender/receiver variances), and
structural offsets.

#### Returns

A list containing the requested predictions:

- `x`, `y`:

  A matrix or data frame where the first column is the actor ID and
  subsequent columns represent the predicted mean values.

- `z`:

  A data frame containing the edgelist with columns: `sender`,
  `receiver`, and `prediction` (probability or intensity).

The results are also invisibly stored in the internal state
`private$.results`.

------------------------------------------------------------------------

### Method `set_coefficients()`

Manually set the model coefficients to new values. This is useful for
sensitivity analyses or applying the model to different scenarios.

#### Usage

    iglm.object.generator$set_coefficients(coef, coef_popularity = NULL)

#### Arguments

- `coef`:

  A numeric vector of new coefficient values for the non-popularity
  terms.

- `coef_popularity`:

  A numeric vector of new coefficient values for the popularity terms,
  if applicable. Must be provided if the model includes popularity
  effects.

#### Returns

The [`iglm.object`](iglm.md) itself, invisibly.

------------------------------------------------------------------------

### Method `get_samples()`

Retrieve the simulated networks stored in the object. Requires
`simulate` or `estimate` to have been run first.

#### Usage

    iglm.object.generator$get_samples()

#### Returns

A list of [`iglm.data`](iglm.data.md) objects representing the simulated
networks, invisibly. Returns an error if no samples are available.

------------------------------------------------------------------------

### Method `set_sampler()`

Replace the internal MCMC sampler with a new one. This is useful for
changing the sampling scheme without redefining the entire model.

#### Usage

    iglm.object.generator$set_sampler(sampler)

#### Arguments

- `sampler`:

  A [`sampler.iglm`](sampler.iglm.md) object. @return The
  [`iglm.object`](iglm.md) itself, invisibly.

------------------------------------------------------------------------

### Method `set_target()`

Replace the internal \`iglm.data\` data object with a new one. This is
useful for applying a fitted model to new observed data. Recalculates
count statistics and re-validates the object.

#### Usage

    iglm.object.generator$set_target(x)

#### Arguments

- `x`:

  A [`iglm.data`](iglm.data.md) “ object containing the new observed
  data.

#### Returns

The [`iglm.object`](iglm.md) itself, invisibly.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    iglm.object.generator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

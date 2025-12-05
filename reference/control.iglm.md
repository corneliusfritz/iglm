# Set Control Parameters for iglm Estimation

Create a list of control parameters for the \`iglm\` estimation
algorithm.

## Usage

``` r
control.iglm(
  estimate_model = TRUE,
  fix_x = FALSE,
  fix_z = FALSE,
  display_progress = FALSE,
  return_samples = TRUE,
  offset_nonoverlap = 0,
  var = FALSE,
  non_stop = FALSE,
  tol = 0.001,
  max_it = 100,
  return_x = FALSE,
  return_y = FALSE,
  return_z = FALSE,
  accelerated = TRUE,
  cluster = NULL,
  exact = FALSE,
  updated_uncertainty = TRUE
)
```

## Arguments

- estimate_model:

  (logical) If \`TRUE\` (default), the main model parameters are
  estimated. If \`FALSE\`, estimation is skipped and only the
  preprocessing is done.

- fix_x:

  (logical) If \`TRUE\`, the 'x' predictor is held fixed during
  estimation/simulation (fixed design in regression). Default is
  \`FALSE\`.

- fix_z:

  (logical) If \`TRUE\`, the 'z' network is held fixed during
  estimation/simulation (fixed network design). Default is \`FALSE\`.
  Setting this to TRUE, allows practicioners to estimate autologistic
  actor attribute models, which were introduced in binary settings in
  Daraganova, G., & Robins, G. (2013).

- display_progress:

  (logical) If \`TRUE\`, display progress messages or a progress bar
  during estimation. Default is \`FALSE\`.

- return_samples:

  (logical). If `TRUE` (default), return simulated network/attribute
  samples (i.e., `iglm.data` objects) generated during estimation (if
  applicable).

- offset_nonoverlap:

  (numeric) A value added to the linear predictor for dyads not in the
  'overlap' set. Default is \`0\`.

- var:

  (logical) If \`TRUE\`, attempt to calculate and return the
  variance-covariance matrix of the estimated parameters. Default is
  \`FALSE\`.

- non_stop:

  (logical) If \`TRUE\`, the estimation algorithm continues until
  \`max_it\` iterations, ignoring the \`tol\` convergence criterion.
  Default is \`FALSE\`.

- tol:

  (numeric) The tolerance level for convergence. The estimation stops
  when the change in coefficients between iterations is less than
  \`tol\`. Default is \`0.001\`.

- max_it:

  (integer) The maximum number of iterations for the estimation
  algorithm. Default is \`100\`.

- return_x:

  (logical). If `TRUE`, return the change statistics for the `x`
  attribute Default is `FALSE`. from samples. Default is \`FALSE\`.
  (Note: \`return_samples=TRUE\` likely implies this).

- return_y:

  (logical). If `TRUE`, return the change statistics for the `y`
  attribute Default is `FALSE`.

- return_z:

  (logical). If `TRUE`, return the change statistics for the `z`
  network. Default is `FALSE`.

- accelerated:

  (logical) If \`TRUE\` (default), an accelerated MM algorithm is used
  based on a Quasi Newton scheme described in the Supplemental Material
  of Fritz et al (2025).

- cluster:

  A parallel cluster object (e.g., from the \`parallel\` package) to use
  for potentially parallelizing parts of the estimation or simulation.
  Default is \`NULL\` (no parallelization).

- exact:

  (logical) If \`TRUE\`, potentially use an exact calculation method of
  the pseudo Fisher information for assessing the uncertainty of the
  estimates. Default is \`FALSE\`.

- updated_uncertainty:

  (logical) If \`TRUE\` (default), potentially use an updated method for
  calculating uncertainty estimates (based on the mean-value theorem as
  opposed to the Godambe Information).

## Value

A list object of class \`"control.iglm"\` containing the specified
control parameters.

## References

Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear. Daraganova, G., and Robins, G. (2013). Exponential random
graph models for social networks: Theory, methods and applications,
102-114. Cambridge University Press.

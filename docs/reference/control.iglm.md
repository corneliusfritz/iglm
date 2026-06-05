# Set Control Parameters for iglm Estimation

Create a list of control parameters for the \`iglm\` estimation
algorithm.

## Usage

``` r
control.iglm(
  estimate_model = TRUE,
  display_progress = FALSE,
  return_samples = TRUE,
  offset_nonoverlap = 0,
  var_method = "Mean-value",
  non_stop = FALSE,
  tol = 0.001,
  max_it = 100,
  return_x = FALSE,
  return_y = FALSE,
  return_z = FALSE,
  accelerated = TRUE,
  exact = TRUE
)
```

## Arguments

- estimate_model:

  (logical) If \`TRUE\` (default), the main model parameters are
  estimated. If \`FALSE\`, estimation is skipped and only the
  preprocessing is done.

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

- var_method:

  (string) Method for variance estimation. Options are "Mean-value"
  (default), "Godambe", and "Hessian". The mean-value version is
  described in Section 3.3 of Fritz et al. (2025), the Godambe method is
  described in Schmid and Hunter (2023), and the "Hessian" option just
  assumes that the pseudo likelihood is the correct likelihood.

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

- exact:

  (logical) If \`TRUE\`, the pseudo Fisher information is calculated
  exact for assessing the uncertainty of the estimates. Default is
  \`FALSE\`.

## Value

A list object of class \`"control.iglm"\` containing the specified
control parameters.

## References

Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

Schmid, C.S. and D. R. Hunter (2023). Computing Pseudolikelihood
Estimators for Exponential-Family Random Graph Models. Journal of Data
Science

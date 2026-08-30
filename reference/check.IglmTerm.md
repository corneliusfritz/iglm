# Check Arguments for iglm Model Terms

This is an internal helper function used to validate and set defaults
for arguments passed to iglm model terms.

## Usage

``` r
check.IglmTerm(
  data_object,
  arglist,
  mandatory = character(0),
  expected = list(),
  defaults = list(),
  directed = NULL
)
```

## Arguments

- data_object:

  The iglm.data object.

- arglist:

  The list of arguments passed to the term.

- mandatory:

  Character vector of mandatory argument names.

- expected:

  A list where keys are expected argument names and values are either a
  character vector of allowed values, or a type string ("numeric",
  "matrix").

- defaults:

  a list of default values for arguments.

- directed:

  Logical indicating if the term is only for directed (TRUE) or
  undirected (FALSE) networks.

## Value

A modified `arglist` with defaults applied and validated values.

## Details

`check.IglmTerm` normalizes and validates arguments passed to model
terms in formulas:

1.  **Positional Argument Normalization**: Arguments passed without
    parameter names are sorted by position and mapped sequentially to
    available candidate parameters. Candidate parameters are identified
    from `mandatory`, `defaults`, and `expected` (excluding metadata and
    arguments already supplied by name).

2.  **Excess Positional Arguments**: Unmatched positional arguments
    (where position index exceeds the number of candidate parameters)
    remain in `arglist` and trigger an unexpected-argument error during
    validation.

3.  **Type and Value Validation**: Validates mandatory arguments,
    allowed categorical values, numeric/matrix types, scalar
    constraints, and ensures no `NA`/`NaN` values are present in numeric
    inputs.

4.  **Default Injection**: Injects default values for any omitted
    optional parameters.

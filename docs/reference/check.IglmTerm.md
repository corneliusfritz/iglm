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

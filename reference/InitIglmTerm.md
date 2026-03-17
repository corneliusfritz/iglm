# Initialize iglm Model Terms

This is an internal generic function used to initialize mapping and data
for terms in a [`iglm.object`](iglm.md) formula.

## Usage

``` r
InitIglmTerm(data_object, arglist, ...)
```

## Arguments

- data_object:

  An [`iglm.data`](iglm.data.md) object.

- arglist:

  A list of arguments passed to the term in the formula.

- ...:

  Additional arguments.

## Value

A list containing information for the C++ backend.

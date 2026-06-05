# Constructor for the results R6 Object

Creates a new instance of the \`results\` R6 class. This class is
designed to store various outputs from \`iglm\` model estimation and
simulation. Users typically do not need to call this constructor
directly; it is used internally by the \`iglm_object\`.

## Usage

``` r
results(size_coef, size_coef_degrees, file = NULL)
```

## Arguments

- size_coef:

  (integer) The number of non-degrees coefficients the object should be
  initialized to accommodate.

- size_coef_degrees:

  (integer) The number of degrees coefficients the object should be
  initialized to accommodate.

- file:

  (character or NULL) Optional file path to load a previously saved
  \`results\` object. If provided, the object will be initialized by
  loading from this file.

## Value

An object of class \`results\` (and \`R6\`), initialized with empty or
NA structures appropriately sized based on the input dimensions.

# Generate the Skeleton for an R package to implement additional iglm terms

This function generates the directory structure and source files for a
new R package named `iglm.userterms` (or whatever name is provided in
the parameter `pkg_name`). This auxiliary package serves as a template
for extending the `iglm` framework to user-defined sufficient
statistics. By compiling this package, users can link custom C++
implementations of change statistics directly with the `iglm` package,
enabling seamless integration of new model terms.

## Usage

``` r
create_userterms_skeleton(path = ".", pkg_name = "iglm.userterms")
```

## Arguments

- path:

  A character string specifying the path where the package directory
  should be created. Defaults to the current working directory (`"."`).

- pkg_name:

  A character string specifying the name of the package to be created.

## Details

The function creates a directory with the name specified in `pkg_name`
at the specified location. As an example for a possible statistic, the
statistic counting mutual connections in the network is implemented.
After defining all possible change-statistics in the c++ function (this
has to include a change for `z_ij` (network), `x_i` (attribute x), and
`y_i` (attribute y) all toggling from 0 to 1), the function has to be
registered using the `EFFECT_REGISTER` macro. After compiling the
package, users have to load the package using
[`library(pkg_name)`](https://rdrr.io/r/base/library.html) before using
it in `iglm`.

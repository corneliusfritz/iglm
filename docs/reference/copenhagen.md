# Copenhagen Network Study

A preprocessed dataset containing social ties, physical proximity, and
nodal attributes for a subset of participants in the Copenhagen Networks
Study. The object is provided as an `iglm.data` class.

## Usage

``` r
data(copenhagen)
```

## Format

The `iglm.data` provides the following information:

- z_network:

  A \\E \times 2\\ integer matrix representing the undirected friendship
  network (\$Z\$).

- x_attribute:

  A logical/binomial vector of length \\N\\ indicating gender (1 for
  female, 0 for male).

- y_attribute:

  A numeric vector of length \\N\\ representing the log-transformed
  total call duration in minutes: \\y_i =
  \log(\frac{\text{seconds}}{60})\\.

- neighborhood:

  A matrix defining the proximity-based constraint space. Pairs are
  included if their cumulative physical proximity exceeded 24 hours
  during the observation period.

- fix_x:

  Boolean `TRUE`, indicating that the \\x\\ attribute is treated as
  exogenous.

## Details

The following preprocessing steps were carried out:

- **Temporal Aggregation:** Proximity data (Bluetooth scans) were
  aggregated into sessions. A session break was defined by any temporal
  gap exceeding 10 minutes.

- **Recursive Pruning:** A recursive filter removed actors with missing
  gender information or isolated actors in either the friendship
  (`z_network`) or proximity (`neighborhood`) networks,

## References

Sapiezynski, P., Stopczynski, A., Lassen, D. D. and Lehmann, S. (2019),
Interaction data from the Copenhagen Networks Study. Scientific Data
6(1), 315.

# Twitter (X) data list for U.S. state legislators (10-state subset)

This data object is data derived from the Twitter (X) interactions
between U.S. state legislators, which is a subset of the data analyzed
in Fritz et al. (2025).' The data is filtered to include only
legislators from 10 states (NY, CA, TX, FL, IL, PA, OH, GA, NC, MI) and
is further subset to the largest connected component based on mention or
retweet activity.

This object contains the main [`iglm_data`](iglm_data.md) object and 5
pre-computed dyadic covariates.

## Usage

``` r
data(state_twitter)
```

## Format

A `list` object containing 6 components. Let N be the number of
legislators in the filtered 10-state subset.

- iglm_data:

  A [`iglm_data`](iglm_data.md) object (which is also a `list`)
  parameterized as follows:

  - `x_attribute`: A binary numeric vector of length N. Value is `1` if
    the legislator's party is 'Republican', `0` otherwise.

  - `y_attribute`: A Poisson numeric vector of length N. Represents the
    count of hatespeech incidents (`actors_data$number_hatespeech`) for
    each legislator.

  - `z_network`: A directed edgelist (2-column matrix) of size
    `n_edges x 2`. A tie `(i, j)` exists if legislator `i` either
    mentioned or retweeted legislator `j`.

  - `neighborhood`: A directed edgelist (2-column matrix). Represents
    the follower network, where a tie `(i, j)` exists if legislator `i`
    follows legislator `j`. Self-loops (diagonal) are included.

- match_gender:

  An N x N `matrix`. `matrix[i, j] = 1` if legislator `i` and legislator
  `j` have the same gender, `0` otherwise.

- match_race:

  An N x N `matrix`. `matrix[i, j] = 1` if legislator `i` and legislator
  `j` have the same race, `0` otherwise.

- match_state:

  An N x N `matrix`. `matrix[i, j] = 1` if legislator `i` and legislator
  `j` are from the same state, `0` otherwise.

- white_attribute:

  A 1 x N `matrix` (a row vector). `matrix[1, i] = 1` if legislator `i`
  is 'White', `0` otherwise.

- gender_attribute:

  A 1 x N `matrix` (a row vector). `matrix[1, i] = 1` if legislator `i`
  is 'female', `0` otherwise.

## References

Gopal, Kim, Nakka, Boehmke, Harden, Desmarais. The National Network of
U.S. State Legislators on Twitter. Political Science Research & Methods,
Forthcoming.

Kim, Nakka, Gopal, Desmarais,Mancinelli, Harden, Ko, and Boehmke (2022).
Attention to the COVID-19 pandemic on Twitter: Partisan differences
among U.S. state legislators. Legislative Studies Quarterly 47,
1023–1041.

Fritz, C., Schweinberger, M. , Bhadra S., and D. R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

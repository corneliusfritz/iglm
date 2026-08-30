# Model Terms in iglm

## Overview

This vignette describes all model terms available in `iglm` (version
1.2.6) for specifying the sufficient statistics of joint
network-attribute models. Terms are passed on the right-hand side of the
`formula` argument in
[`iglm()`](https://corneliusfritz.github.io/iglm/reference/iglm.md) and
govern how individual attributes and network connections jointly
determine the log-linear probabilities of the model.

A model in `iglm` decomposes its sufficient statistics into two
families:

- **Unit-level terms** \\g_i(x_i, y_i)\\: depend only on unit \\i\\’s
  own attributes.
- **Pair-level terms** \\h\_{i,j}(x, y, z)\\: depend on the connection
  \\z\_{i,j}\\ and the attributes of units \\i\\ and \\j\\ as well as
  the wider network.

The total sufficient statistic of the model is then \\ S(x, y, z) =
\sum_i g_i(x_i, y_i) + \sum\_{i \ne j} h\_{i,j}(x, y, z). \\

------------------------------------------------------------------------

## Key Definitions

Before stating all statistics, we introduce the formal notation and
definitions used throughout this vignette:

- **Population and Dyads:**
  - \\𝒫 = \\1, \ldots, N\\\\ denotes the population of \\N\\ units.
  - \\𝒟\\ denotes the set of dyads (pairs of distinct units): \\𝒟 =
    \\(i,j) : 1 \le i \neq j \le N\\\\ for directed connections and \\𝒟
    = \\(i,j) : 1 \le i \< j \le N\\\\ for undirected connections.
- **Variables and Attributes:**
  - \\x_i\\: Exogenous (or secondary) predictor attribute of unit \\i
    \in 𝒫\\.
  - \\y_i\\: Endogenous outcome attribute of unit \\i \in 𝒫\\.
  - \\z\_{i,j} \in \\0, 1\\\\: Binary connection indicator from unit
    \\i\\ to unit \\j\\ for \\(i,j) \in 𝒟\\, collected in the connection
    matrix \\\mathbf{z}\\.
  - \\v_i\\: Optional unit-level exogenous covariate.
  - \\w\_{i,j}\\: Optional dyadic exogenous covariate.
- **Neighbourhoods and Local Structure:**
  - \\𝒩_i \subset 𝒫\\ denotes the local neighbourhood of unit \\i\\
    (with \\i \in 𝒩_i\\).
  - \\c\_{i,j} \in \\0, 1\\\\ is the neighbourhood overlap indicator,
    taking the value 1 if \\𝒩_i \cap 𝒩_j \neq \emptyset\\, and 0
    otherwise.
- **Connections:** Different types of indicators for connections:
  - Overlapping: \\u\_{i,j} = c\_{i,j} z\_{i,j}\\, a connection between
    units \\i\\ and \\j\\ where \\𝒩_i \cap 𝒩_j \neq \emptyset\\.
  - Non-overlapping: \\k\_{i,j} = (1-c\_{i,j}) z\_{i,j}\\, a connection
    between units \\i\\ and \\j\\ where \\𝒩_i \cap 𝒩_j = \emptyset\\.
  - \\e\_{i,j}^{(\mathtt{s})}\\ for \\\mathtt{s} \in \\\mathtt{global},
    \mathtt{local}, \mathtt{alocal}\\\\ is defined by: \\
    e\_{i,j}^{(\mathtt{s})} = \begin{cases} z\_{i,j} & \text{if }
    \mathtt{s} = \mathtt{global}\\ u\_{i,j} & \text{if } \mathtt{s} =
    \mathtt{local} \\ k\_{i,j} & \text{if } \mathtt{s} = \mathtt{alocal}
    \end{cases} \\ The mode parameter \\\mathtt{s}\\ is generally
    defined as \\\mathtt{s} \in \\\mathtt{global}, \mathtt{local},
    \mathtt{alocal}\\\\, but note that for the terms `gwesp`, `gwdsp`,
    `gwodegree`, `gwidegree`, `edges_x_match`, and `edges_y_match`
    (defined in the summary table), only the options \\\mathtt{s} \in
    \\\mathtt{global}, \mathtt{local}\\\\ are implemented as their
    \\\mathtt{alocal}\\ version is not very useful.
- **Degree Statistics:** For unit \\i \in 𝒫\\ and mode \\\mathtt{s} \in
  \\\mathtt{global}, \mathtt{local}\\\\:
  - Out-degree: \\\operatorname{deg}(i, \mathtt{s}) = \sum\_{j \in 𝒫
    \setminus \\i\\} e\_{i,j}^{(\mathtt{s})}\\ with
    \\\operatorname{deg}(i) = \operatorname{deg}(i, \mathtt{global})\\.
  - In-degree: \\\operatorname{ideg}(i, \mathtt{s}) = \sum\_{j \in 𝒫
    \setminus \\i\\} e\_{j,i}^{(\mathtt{s})}\\ with
    \\\operatorname{ideg}(i) = \operatorname{ideg}(i,
    \mathtt{global})\\.
- **Common Partners (CP):** For a dyad \\(i,j) \in 𝒟\\ and mode
  \\\mathtt{s} \in \\\mathtt{global}, \mathtt{local}\\\\, the number of
  shared partners via distinct path structures is defined as:
  - Outgoing Two-Paths (OTP): \\\operatorname{CP}(i, j, \mathtt{s},
    \mathtt{OTP}) = \sum\_{h \in 𝒫 \setminus \\i,j\\}
    e\_{i,h}^{(\mathtt{s})}\\ e\_{h,j}^{(\mathtt{s})}\\.
  - Incoming Shared Partners (ISP): \\\operatorname{CP}(i, j,
    \mathtt{s}, \mathtt{ISP}) = \sum\_{h \in 𝒫 \setminus \\i,j\\}
    e\_{h,i}^{(\mathtt{s})}\\ e\_{h,j}^{(\mathtt{s})}\\.
  - Outgoing Shared Partners (OSP): \\\operatorname{CP}(i, j,
    \mathtt{s}, \mathtt{OSP}) = \sum\_{h \in 𝒫 \setminus \\i,j\\}
    e\_{i,h}^{(\mathtt{s})}\\ e\_{j,h}^{(\mathtt{s})}\\.
  - Incoming Two-Paths (ITP): \\\operatorname{CP}(i, j, \mathtt{s},
    \mathtt{ITP}) = \sum\_{h \in 𝒫 \setminus \\i,j\\}
    e\_{h,i}^{(\mathtt{s})}\\ e\_{j,h}^{(\mathtt{s})}\\.
  - Undirected Version: \\\operatorname{CP}(i, j, \mathtt{s}) = \sum\_{h
    \in 𝒫 \setminus \\i,j\\} e\_{i,h}^{(\mathtt{s})}\\
    e\_{h,j}^{(\mathtt{s})}\\.
- **Miscellaneous:**
  - Geometrically-weighted weight: \\w_k(\alpha) = \exp(\alpha) \left\[
    1 - (1 - \exp(-\alpha))^k \right\]\\.
  - Indicator for directionality: \\\mathbb{I}\_U(\mathbf{z})\\, taking
    the value 1 if connections in \\\mathbf{z}\\ are undirected, and 0
    otherwise.
  - Indicator for transitive connection: \\d\_{i,j}(\mathbf{z}) =
    \mathbb{I}(\exists\\ k \in 𝒩_i \cap 𝒩_j: z\_{i,k} = z\_{k,j} = 1)\\.

The sections below and the summary table list all implemented terms as
of `iglm` version 1.2.6 and will be extended in future releases.

------------------------------------------------------------------------

## Category 1: Attribute Dependence Terms (\\g_i\\ Terms)

These terms capture how individual predictors \\x_i\\ (exogenous) and
\\y_i\\ (endogenous) relate to each other, without reference to the
network.

### `attribute_x`

**Description:** Intercept for the endogenous \\x\\-attribute.

\\ g_i(x_i, y_i) = x_i \\

``` r

formula <- object ~ attribute_x
```

------------------------------------------------------------------------

### `attribute_y`

**Description:** Intercept for the endogenous \\y\\-attribute.

\\ g_i(x_i, y_i) = y_i \\

``` r

formula <- object ~ attribute_y
```

------------------------------------------------------------------------

### `cov_x(data = v)`

**Description:** Effect of a unit-level exogenous covariate \\v_i\\ on
attribute \\x_i\\.

\\ g_i(x_i, y_i) = v_i\\ x_i \\

``` r

formula <- object ~ cov_x(data = v)
```

------------------------------------------------------------------------

### `cov_y(data = v)`

**Description:** Effect of a unit-level exogenous covariate \\v_i\\ on
attribute \\y_i\\.

\\ g_i(x_i, y_i) = v_i\\ y_i \\

``` r

formula <- object ~ cov_y(data = v)
```

------------------------------------------------------------------------

### `attribute_xy(mode = "global" | "local" | "alocal")`

**Description:** Interaction between the two attributes \\x_i\\ and
\\y_i\\, optionally mediated by the neighbourhood structure.

| Mode     | Formula                                                         |
|----------|-----------------------------------------------------------------|
| `global` | \\x_i\\ y_i\\                                                   |
| `local`  | \\x_i \sum\_{j \in 𝒩_i} y_j + y_i \sum\_{j \in 𝒩_i} x_j\\       |
| `alocal` | \\x_i \sum\_{j \notin 𝒩_i} y_j + y_i \sum\_{j \notin 𝒩_i} x_j\\ |

``` r

formula <- object ~ attribute_xy(mode = "local")
```

------------------------------------------------------------------------

## Category 2: Network Dependence Terms (\\h\_{i,j}\\ Terms)

These terms capture how the network topology \\z\\ drives edge
formation. All are pair-level statistics.

### `degrees`

**Description:** Node-level degree fixed effects. One parameter per
unit, capturing heterogeneity in activity not explained by other terms.
Estimation relies on an MM algorithm constraint.

``` r

formula <- object ~ degrees
```

------------------------------------------------------------------------

### `edges(mode = "global" | "local" | "alocal")`

**Description:** Baseline propensity for a tie \\z\_{i,j}\\ to form; the
network analogue of an intercept.

\\ h\_{i,j}(x, y, z) = e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ edges(mode = "global")
formula <- object ~ edges(mode = "local")
formula <- object ~ edges(mode = "alocal")
```

------------------------------------------------------------------------

### `mutual(mode = "global" | "local" | "alocal")`

**Description:** Reciprocity in directed networks. Counts pairs where
\\i \to j\\ and \\j \to i\\ both exist (counted once per unordered pair,
hence the factor \\1/2\\).

\\ h\_{i,j}(x, y, z) = \frac{e\_{i,j}^{(\mathtt{s})}\\
e\_{j,i}^{(\mathtt{s})}}{2} \\

Only valid for **directed** networks.

``` r

formula <- object ~ mutual(mode = "global")
```

------------------------------------------------------------------------

### `cov_z(data = w, mode = "global" | "local" | "alocal")`

**Description:** Dyadic covariate — exogenous edge-level covariate
\\w\_{i,j}\\ influences tie formation.

\\ h\_{i,j}(x, y, z) = w\_{i,j}\\ e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ cov_z(data = W, mode = "global")
```

------------------------------------------------------------------------

### `cov_z_out(data = v, mode = "global" | "local" | "alocal")`

**Description:** Sender covariate — exogenous nodal attribute \\v_i\\
influences the propensity to *send* a tie.

\\ h\_{i,j}(x, y, z) = v_i\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ cov_z_out(data = v, mode = "global")
```

------------------------------------------------------------------------

### `cov_z_in(data = v, mode = "global" | "local" | "alocal")`

**Description:** Receiver covariate — exogenous nodal attribute \\v_j\\
influences the propensity to *receive* a tie.

\\ h\_{i,j}(x, y, z) = v_j\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ cov_z_in(data = v, mode = "global")
```

------------------------------------------------------------------------

### `isolates`

**Description:** Captures the proportion of units with no connections at
all (total degree zero).

\\ h\_{i,j}(x, y, z) = \mathbb{I}\\\left(\sum\_{j \in 𝒫 \setminus \\i\\}
z\_{i,j} + z\_{j,i} = 0\right) \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ isolates
```

------------------------------------------------------------------------

### `nonisolates`

**Description:** Captures the proportion of units that have at least one
connection.

\\ h\_{i,j}(x, y, z) = \mathbb{I}\\\left(\sum\_{j \in 𝒫 \setminus \\i\\}
z\_{i,j} + z\_{j,i} \ne 0\right) \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ nonisolates
```

------------------------------------------------------------------------

### `gwdegree(mode = "global" | "local", decay = α)`

**Description:** Geometrically Weighted Degree — captures the overall
degree distribution with exponential decay parameter \\\alpha\\.

\\ h\_{i,j}(x, y, z) = w\_{\operatorname{deg}(i)}(\alpha) +
w\_{\operatorname{deg}(j)}(\alpha) \\

Suitable for both directed and undirected networks. Only
`mode %in% c("global", "local")` is available.

``` r

formula <- object ~ gwdegree(mode = "global", decay = 0.5)
```

------------------------------------------------------------------------

### `gwodegree(mode = "global" | "local", decay = α)`

**Description:** Geometrically Weighted Out-Degree — captures the
out-degree distribution in directed networks.

\\ h\_{i,j}(x, y, z) = w\_{\operatorname{deg}(i,\\\mathtt{s})}(\alpha)
\\

Only valid for **directed** networks. Only
`mode %in% c("global", "local")` is available.

``` r

formula <- object ~ gwodegree(mode = "global", decay = 0.5)
```

------------------------------------------------------------------------

### `gwidegree(mode = "global" | "local", decay = α)`

**Description:** Geometrically Weighted In-Degree — captures the
in-degree distribution in directed networks.

\\ h\_{i,j}(x, y, z) = w\_{\operatorname{ideg}(i,\\\mathtt{s})}(\alpha)
\\

Only valid for **directed** networks. Only
`mode %in% c("global", "local")` is available.

``` r

formula <- object ~ gwidegree(mode = "global", decay = 0.5)
```

------------------------------------------------------------------------

### `transitive`

**Description:** Transitivity indicator — rewards edges that close a
locally transitive triple.

\\ h\_{i,j}(x, y, z) = d\_{i,j}(\mathbf{z})\\ z\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ transitive
```

------------------------------------------------------------------------

### `gwesp_symm(mode = "global" | "local", decay = α)`

**Description:** Geometrically Weighted Edgewise Shared Partners
(undirected) — the classic GWESP statistic for undirected networks.

\\ h\_{i,j}(x, y, z) = e\_{i,j}^{(\mathtt{s})}\\
w\_{\operatorname{CP}(i,j,\mathtt{s})}(\alpha) \\

Suitable for undirected networks only.

``` r

formula <- object ~ gwesp_symm(mode = "global", decay = 0.5)
```

------------------------------------------------------------------------

### `gwesp(mode = "global" | "local", type = "OTP" | "ISP" | "OSP" | "ITP", decay = α)`

**Description:** Geometrically Weighted Edgewise Shared Partners
(directed) — conditions shared partners on a specific path type.

\\ h\_{i,j}(x, y, z) = e\_{i,j}^{(\mathtt{s})}\\
w\_{\operatorname{CP}(i,j,\mathtt{s},\mathtt{type})}(\alpha) \\

Only valid for **directed** networks. Only
`mode %in% c("global", "local")` is available.

``` r

formula <- object ~ gwesp(mode = "global", type = "OTP", decay = 0.5)
```

------------------------------------------------------------------------

### `gwdsp_symm(mode = "local", decay = α)`

**Description:** Geometrically Weighted Dyadwise Shared Partners
(undirected) — models triadic potential irrespective of the closing
edge.

\\ h\_{i,j}(x, y, z) =
w\_{\operatorname{CP}(i,j,\mathtt{local})}(\alpha) \\

Suitable for undirected networks only.

``` r

formula <- object ~ gwdsp_symm(mode = "local", decay = 0.5)
```

------------------------------------------------------------------------

### `gwdsp(mode = "global" | "local", type = "OTP" | "ISP" | "OSP" | "ITP", decay = α)`

**Description:** Geometrically Weighted Dyadwise Shared Partners
(directed) — models directed triadic potential irrespective of the
closing edge.

\\ h\_{i,j}(x, y, z) =
w\_{\operatorname{CP}(i,j,\mathtt{s},\mathtt{type})}(\alpha) \\

Only valid for **directed** networks. Only
`mode %in% c("global", "local")` is available.

``` r

formula <- object ~ gwdsp(mode = "global", type = "OTP", decay = 0.5)
```

------------------------------------------------------------------------

## Category 3: Joint Attribute/Network Dependence Terms (\\h\_{i,j}\\ Terms)

These terms capture the interplay between nodal attributes and network
position. They are the key building blocks for studying spillover
effects.

### `attribute_xz(mode = "local")`

**Description:** Additive effect of \\x_i\\ and \\x_j\\ on local edge
formation.

\\ h\_{i,j}(x, y, z) = (x_i + x_j)\\ u\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ attribute_xz(mode = "local")
```

------------------------------------------------------------------------

### `attribute_yz(mode = "local")`

**Description:** Additive effect of \\y_i\\ and \\y_j\\ on local edge
formation.

\\ h\_{i,j}(x, y, z) = (y_i + y_j)\\ u\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ attribute_yz(mode = "local")
```

------------------------------------------------------------------------

### `edges_x_match(mode = "global" | "local")`

**Description:** Homophily on \\x\\ — rewards edges between units with
equal \\x\\-values.

\\ h\_{i,j}(x, y, z) = \mathbb{I}(x_i = x_j)\\ e\_{i,j}^{(\mathtt{s})}
\\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ edges_x_match(mode = "global")
```

------------------------------------------------------------------------

### `edges_y_match(mode = "global" | "local")`

**Description:** Homophily on \\y\\ — rewards edges between units with
equal \\y\\-values.

\\ h\_{i,j}(x, y, z) = \mathbb{I}(y_i = y_j)\\ e\_{i,j}^{(\mathtt{s})}
\\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ edges_y_match(mode = "global")
```

------------------------------------------------------------------------

### `outedges_x(mode = "global" | "local" | "alocal")`

**Description:** Effect of sender attribute \\x_i\\ on out-degree
formation.

\\ h\_{i,j}(x, y, z) = x_i\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ outedges_x(mode = "global")
```

------------------------------------------------------------------------

### `inedges_x(mode = "global" | "local" | "alocal")`

**Description:** Effect of receiver attribute \\x_j\\ on in-degree
reception.

\\ h\_{i,j}(x, y, z) = x_j\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ inedges_x(mode = "global")
```

------------------------------------------------------------------------

### `outedges_y(mode = "global" | "local" | "alocal")`

**Description:** Effect of sender attribute \\y_i\\ on out-degree
formation.

\\ h\_{i,j}(x, y, z) = y_i\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ outedges_y(mode = "global")
```

------------------------------------------------------------------------

### `inedges_y(mode = "global" | "local" | "alocal")`

**Description:** Effect of receiver attribute \\y_j\\ on in-degree
reception.

\\ h\_{i,j}(x, y, z) = y_j\\ e\_{i,j}^{(\mathtt{s})} \\

Only valid for **directed** networks.

``` r

formula <- object ~ inedges_y(mode = "global")
```

------------------------------------------------------------------------

### `spillover_xx(mode = "local")`

**Description:** Symmetric \\x\\-to-\\x\\ spillover — the product \\x_i
x_j\\ along local connections, capturing peer effects in the \\x\\
attribute.

\\ h\_{i,j}(x, y, z) = x_i\\ x_j\\ u\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_xx(mode = "local")
```

------------------------------------------------------------------------

### `spillover_xx_scaled(mode = "global" | "local")`

**Description:** Degree-normalised \\x\\-to-\\x\\ spillover, accounting
for the number of neighbours.

\\ h\_{i,j}(x, y, z) = \left(\frac{x_i\\
x_j}{\operatorname{deg}(i,\mathtt{s})} + \frac{x_j\\
x_i}{\operatorname{deg}(j,\mathtt{s})}\\\mathbb{I}\_U(\mathbf{z})\right)
e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_xx_scaled(mode = "global")
```

------------------------------------------------------------------------

### `spillover_yy(mode = "local")`

**Description:** Symmetric \\y\\-to-\\y\\ spillover — the product \\y_i
y_j\\ along local connections.

\\ h\_{i,j}(x, y, z) = y_i\\ y_j\\ u\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_yy(mode = "local")
```

------------------------------------------------------------------------

### `spillover_yy_scaled(mode = "global" | "local")`

**Description:** Degree-normalised \\y\\-to-\\y\\ spillover.

\\ h\_{i,j}(x, y, z) = \left(\frac{y_i\\
y_j}{\operatorname{deg}(i,\mathtt{s})} + \frac{y_j\\
y_i}{\operatorname{deg}(j,\mathtt{s})}\\\mathbb{I}\_U(\mathbf{z})\right)
e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_yy_scaled(mode = "global")
```

------------------------------------------------------------------------

### `spillover_xy(mode = "local")`

**Description:** Symmetric cross-attribute spillover — \\x_i \to y_j\\
and \\x_j \to y_i\\ along local connections. For undirected networks
both directions are summed.

\\ h\_{i,j}(x, y, z) = x_i\\ y_j\\ u\_{i,j} + x_j\\ y_i\\ u\_{i,j}\\
\mathbb{I}\_U(\mathbf{z}) \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_xy(mode = "local")
```

------------------------------------------------------------------------

### `spillover_xy_scaled(mode = "global" | "local")`

**Description:** Degree-normalised symmetric cross-attribute spillover
(\\x \to y\\).

\\ h\_{i,j}(x, y, z) = \left(\frac{x_i\\
y_j}{\operatorname{deg}(i,\mathtt{s})} + \frac{x_j\\
y_i}{\operatorname{deg}(j,\mathtt{s})}\\\mathbb{I}\_U(\mathbf{z})\right)
e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_xy_scaled(mode = "global")
```

------------------------------------------------------------------------

### `spillover_yx(mode = "local")`

**Description:** Directed cross-attribute spillover — \\y_i \to x_j\\
only (no symmetrisation). Only for directed networks.

\\ h\_{i,j}(x, y, z) = y_i\\ x_j\\ u\_{i,j} \\

Only valid for **directed** networks.

``` r

formula <- object ~ spillover_yx(mode = "local")
```

------------------------------------------------------------------------

### `spillover_yx_scaled(mode = "global" | "local")`

**Description:** Degree-normalised cross-attribute spillover (\\y \to
x\\), with symmetrisation for undirected networks.

\\ h\_{i,j}(x, y, z) = \left(\frac{y_i\\
x_j}{\operatorname{deg}(i,\mathtt{s})} + \frac{y_j\\
x_i}{\operatorname{deg}(j,\mathtt{s})}\\\mathbb{I}\_U(\mathbf{z})\right)
e\_{i,j}^{(\mathtt{s})} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_yx_scaled(mode = "global")
```

------------------------------------------------------------------------

### `spillover_yc(mode = "local", data = v)`

**Description:** Interaction of endogenous attribute \\y\\ with
exogenous covariate \\v\\ along overlapping connections, with
symmetrisation for undirected networks.

\\ h\_{i,j}(x, y, z) = c\_{i,j}\bigl(v_j\\ y_i +
\mathbb{I}\_U(\mathbf{z})\\ v_i\\ y_j\bigr)\\ z\_{i,j} \\

Suitable for both directed and undirected networks.

``` r

formula <- object ~ spillover_yc(data = v, mode = "local")
```

------------------------------------------------------------------------

## Quick-Reference Table

The table below summarises all implemented terms, indicating which
variables (\\x\\, \\y\\, \\z\\) they involve, and whether they support
undirected networks.

| Term | \\x\\ | \\y\\ | \\z\\ | Undirected |
|:---|:--:|:--:|:--:|:--:|
| [`attribute_x`](#attribute_x) | ✓ |  |  | ✓ |
| [`attribute_y`](#attribute_y) |  | ✓ |  | ✓ |
| [`cov_x`](#cov_x) | ✓ |  |  | ✓ |
| [`cov_y`](#cov_y) |  | ✓ |  | ✓ |
| [`attribute_xy(mode = "s")`](#attribute_xy) | ✓ | ✓ |  | ✓ |
| [`degrees`](#degrees) |  |  | ✓ | ✓ |
| [`edges(mode = "s")`](#edges) |  |  | ✓ | ✓ |
| [`mutual(mode = "s")`](#mutual) |  |  | ✓ | ✗ |
| [`cov_z(mode = "s")`](#cov_z) |  |  | ✓ | ✓ |
| [`cov_z_out(mode = "s")`](#cov_z_out) |  |  | ✓ | ✗ |
| [`cov_z_in(mode = "s")`](#cov_z_in) |  |  | ✓ | ✗ |
| [`isolates`](#isolates) |  |  | ✓ | ✓ |
| [`nonisolates`](#nonisolates) |  |  | ✓ | ✓ |
| [`gwdegree(mode = "global")`](#gwdegree) |  |  | ✓ | ✓ |
| [`gwodegree(mode = "s")`](#gwodegree) |  |  | ✓ | ✗ |
| [`gwidegree(mode = "s")`](#gwidegree) |  |  | ✓ | ✗ |
| [`transitive`](#transitive) |  |  | ✓ | ✓ |
| [`gwesp_symm(mode = "s")`](#gwesp_symm) |  |  | ✓ | ✓ |
| [`gwesp(mode = "s", type = "…")`](#gwesp) |  |  | ✓ | ✗ |
| [`gwdsp_symm(mode = "local")`](#gwdsp_symm) |  |  | ✓ | ✓ |
| [`gwdsp(mode = "s", type = "…")`](#gwdsp) |  |  | ✓ | ✗ |
| [`attribute_xz(mode = "local")`](#attribute_xz) | ✓ |  | ✓ | ✓ |
| [`attribute_yz(mode = "local")`](#attribute_yz) |  | ✓ | ✓ | ✓ |
| [`edges_x_match(mode = "s")`](#edges_x_match) | ✓ |  | ✓ | ✓ |
| [`edges_y_match(mode = "s")`](#edges_y_match) |  | ✓ | ✓ | ✓ |
| [`outedges_x(mode = "s")`](#outedges_x) | ✓ |  | ✓ | ✗ |
| [`inedges_x(mode = "s")`](#inedges_x) | ✓ |  | ✓ | ✗ |
| [`outedges_y(mode = "s")`](#outedges_y) |  | ✓ | ✓ | ✗ |
| [`inedges_y(mode = "s")`](#inedges_y) |  | ✓ | ✓ | ✗ |
| [`spillover_xx(mode = "local")`](#spillover_xx) | ✓ |  | ✓ | ✓ |
| [`spillover_xx_scaled(mode = "s")`](#spillover_xx_scaled) | ✓ |  | ✓ | ✓ |
| [`spillover_yy(mode = "local")`](#spillover_yy) |  | ✓ | ✓ | ✓ |
| [`spillover_yy_scaled(mode = "s")`](#spillover_yy_scaled) |  | ✓ | ✓ | ✓ |
| [`spillover_xy(mode = "local")`](#spillover_xy) | ✓ | ✓ | ✓ | ✓ |
| [`spillover_xy_scaled(mode = "s")`](#spillover_xy_scaled) | ✓ | ✓ | ✓ | ✓ |
| [`spillover_yx(mode = "local")`](#spillover_yx) | ✓ | ✓ | ✓ | ✗ |
| [`spillover_yx_scaled(mode = "s")`](#spillover_yx_scaled) | ✓ | ✓ | ✓ | ✓ |
| [`spillover_yc(mode = "local")`](#spillover_yc) |  | ✓ | ✓ | ✓ |

------------------------------------------------------------------------

## References

Fritz, C., Schweinberger, M., Bhadra, S., and D.R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. *Journal of the American Statistical Association*,
to appear. <doi:10.1080/01621459.2025.2565851>

Schweinberger, M. and M.S. Handcock (2015). Local Dependence in Random
Graph Models: Characterization, Properties, and Statistical Inference.
*Journal of the Royal Statistical Society, Series B*, 7, 647–676.

Schweinberger, M. and J.R. Stewart (2020). Concentration and Consistency
Results for Canonical and Curved Exponential-Family Models of Random
Graphs. *The Annals of Statistics*, 48, 374–396.

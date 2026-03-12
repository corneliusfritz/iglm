# Model specification for a \`iglm' object

The help pages of [`iglm`](iglm.md) describe the model with details on
model fitting and estimation. Generally, a model is specified via it's
sufficient statistics, that can be further decomposed into two parts:

- **\\\mathbf{g}\_i(x_i^\*,y_i^\*) = \mathbf{g}\_i(x_i,y_i)=
  (g_i(x_i,y_i))\\**: A vector of unit-level functions (or "g-terms")
  that describe the relationship between an individual actor \\i\\'s
  predictors (\\x_i\\) and their own response (\\y_i\\).

- **\\\mathbf{h}\_{i,j}(x_i^\*,x_j^\*, y_i^\*, y_j^\*, z) =
  \mathbf{h}\_{i,j}(x,y,z)= (h\_{i,j}(x,y,z))\\**: A vector of
  pair-level functions (or "h-terms") that specify how the connections
  (\\z\\) and responses (\\y_i, y_j\\) of a pair of units \\\\i,j\\\\
  depend on each other and the wider network structure.

Each term defines a component for the model's features, which are a sum
of unit-level components, \\\sum_i g_i(x_i,y_i)\\, and/or pair-level
components, \\\sum\_{i \ne j} h\_{i,j}(x,y,z)\\. The implemented terms
are grouped into three categories:

1.  \\{g}\_i\\ terms for attribute dependence,

2.  \\{h}\_{i,j}\\ terms for network dependence,

3.  \\{h}\_{i,j}\\ Terms for joint attribute/network dependence.

Below is a detailed description of each term that can be specified in
the model formula (see, [`iglm`](iglm.md)):
`iglm.data ~ <term_1> + <term_2> + ... `, where the left-hand side of
the formula has to be a [`iglm.data`](iglm.data.md) object, and the
right-hand side lists some of the terms from the list below, which
should be included in the model. Setting, e.g., `<term_1>` to
`attribute_x` includes the term \\g_i(x,y,z) = x_i\\ in \\g_i\\. Note
that the function
[`create_userterms_skeleton`](create_userterms_skeleton.md) can be used
to create custom terms.

## Notation

Here, \\x_i\\ and \\y_i\\ are the attributes for actor \\i\\, and
\\z\_{i,j}\\ indicates the presence (1) or absence (0) of a tie from
actor \\i\\ to actor \\j\\. The local neighborhood of actor \\i\\ is
denoted \\\mathcal{N}\_i\\, and the indicator for whether actors \\i\\
and \\j\\ share a local neighborhood is given by \\c\_{i,j} =
\mathbb{I}(\mathcal{N}\_i \cap \mathcal{N}\_j \neq \emptyset)\\. The
functions below specify the forms of \\g_i(x_i,y_i)\\ and
\\h\_{i,j}(x,y,z)\\ for each term. Some terms also depend on other
covariates, which are denoted by \\v = (v_1, ..., v_N)\\ (unit-level)
and \\w = (w\_{i,j}) \in \mathbb{R}^{N \times N}\\(dyadic). These
covariates must be provided by the user via the `data` argument of the
terms. Assuming that the matriv `x` exists in the environment associated
with the used formula, `cov_z(data = v)` includes the dyadic covariable
\\v = (v\_{i,j})\\ in the model. Some terms also have a `mode` argument,
which can take values `"global"`, `"local"`, or `"alocal"`. The
`"global"` mode indicates that the statistic is computed over the entire
network, while `"local"` restricts the statistic to local neighborhoods
only (i.e., edges where \\c\_{i,j} = 1\\). The `"alocal"` mode restricts
the statistic to non-local edges only (i.e., edges where \\c\_{i,j} =
0\\). For instance, `edges(mode = "local")` counts the number of edges
that connect actors with overlapping neighborhoods. See underneath for
which options are implemented for each term. See the documentation for
[`iglm`](iglm.md) for details on model fitting and estimation.

## List of Terms

- **1. \\g_i\\ Terms for Attribute Dependence**:

- `attribute_x`:

  **Attribute (X) \[g-term\]:** Intercept for attribute \\x\\.
  \\g_i(x_i,y_i) = x_i\\

- `attribute_y`:

  **Attribute (Y) \[g-term\]:** Intercept for attribute \\y\\.
  \\g_i(x_i,y_i) = y_i\\

- `cov_x`:

  **Nodal Covariate (X) \[g-term\]:** Effect of a unit-level exogenous
  covariate \\v_i\\ on endogenous attribute \\x_i\\. \\g_i(x_i,y_i) =
  v_i x_i\\

- `cov_y(data = v)`:

  **Nodal Covariate (Y) \[g-term\]:** Effect of a unit-level exogenous
  covariate \\v_i\\ on endogenous attribute \\y_i\\. \\g_i(x_i,y_i) =
  v_i y_i\\

- `attribute_xy(mode = "global")`:

  **Nodal Attribute Interaction (X-Y) \[g-term\]:** Interaction of
  attributes \\x_i\\ and \\y_i\\. Spatial or structural boundaries
  determine the interaction scale.

  - `global`: \\g_i(x_i,y_i) = x_i y_i\\

  - `local`: \\g_i(x_i,y_i) = x_i \sum\_{j \in \mathcal{N}\_i} y_j + y_i
    \sum\_{j \in \mathcal{N}\_i} x_j\\

  - `alocal`: \\g_i(x_i,y_i) = x_i \sum\_{j \notin \mathcal{N}\_i} y_j +
    y_i \sum\_{j \notin \mathcal{N}\_i} x_j\\

- **2. \\h\_{i,j}\\ Terms for Network Dependence**:

- `degrees`:

  **Degrees \[h-term\]:** Specifies node-level fixed effects. Estimation
  requires an MM algorithm constraint. Directed networks assign distinct
  sender and receiver effects (constraining node N's out-effect to 0 for
  identifiability). Undirected networks assign a singular degree effect
  per node.

- `edges(mode = "global")`:

  **Edges \[h-term\]:** Captures the baseline propensity of tie
  formation \\z\_{i,j}\\, partitioned by structural boundary
  \\c\_{i,j}\\.

  - `global`: \\h\_{i,j}(x,y,z) = z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j}\\

- `mutual(mode = "global")`:

  **Mutual Reciprocity \[h-term\]:** Evaluates reciprocal tie formation
  in directed networks.

  - `global`: \\h\_{i,j}(x,y,z) = z\_{i,j} z\_{j,i}\\ (for \\i \< j\\)

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j} z\_{j,i}\\ (for \\i
    \< j\\)

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j} z\_{j,i}\\
    (for \\i \< j\\)

- `cov_z(data, mode = "global")`:

  **Dyadic Covariate \[h-term\]:** Exogenous dyadic covariate
  \\w\_{i,j}\\ influence on edge formation.

  - `global`: \\h\_{i,j}(x,y,z) = w\_{i,j} z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} w\_{i,j} z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) w\_{i,j} z\_{i,j}\\

- `isolates`:

  **Isolates \[z-term\]:** Captures the frequency of nodes with degree
  zero \\\mathbb{I} (\sum\_{j\neq i} z\_{i,j} = 0)\\.

- `nonisolates`:

  **Non-Isolates \[z-term\]:** Captures the frequency of nodes with
  degree strictly greater than zero \\\mathbb{I} (\sum\_{j\neq i}
  z\_{i,j} \neq 0)\\.

- `gwodegree(mode = "global", decay = 2.0)`:

  **Geometrically Weighted Out-Degree \[z-term\]:** Captures the
  out-degree distribution utilizing an exponential decay parameter
  \\\alpha\\.

  - `global`: Evaluated across the entire adjacency matrix.

  - `local`: Evaluated strictly within local neighborhood partitions
    \\c\_{i,j}\\.

- `gwidegree(mode = "global", decay = 2.0)`:

  **Geometrically Weighted In-Degree \[z-term\]:** Captures the
  in-degree distribution utilizing an exponential decay parameter
  \\\alpha\\.

  - `global`: Evaluated across the entire adjacency matrix.

  - `local`: Evaluated strictly within local neighborhood partitions
    \\c\_{i,j}\\.

- `gwesp(data, mode = "global", type = "OTP", decay = 2.0)`:

  **Geometrically Weighted Edgewise-Shared Partners \[h-term\]:** Models
  triadic closure propensity conditioning on existing edges, scaled by
  decay parameter \\\alpha\\. Types dictate the structural path
  constraint: OTP (outgoing two-paths, \\z\_{i,h}\\z\_{h,j}\\
  z\_{i,j}\\), ITP (incoming two-paths, \\z\_{h,i}\\z\_{j,h}\\
  z\_{i,j}\\), OSP (outgoing shared partners, \\z\_{i,h}\\z\_{j,h}\\
  z\_{i,j}\\), ISP (incoming shared partners, \\z\_{h,i}\\z\_{h,j}\\
  z\_{i,j}\\).

  - `global`: Path evaluations span the entire network topology.

  - `local`: Path evaluations are constrained exclusively to
    \\c\_{i,j}\\ overlapping neighborhoods.

- `gwdsp(data, mode = "global", type = "OTP", decay = 2.0)`:

  **Geometrically Weighted Dyadwise-Shared Partners \[h-term\]:** Models
  the baseline triadic potential irrespective of the closing
  \\z\_{i,j}\\ edge, scaled by decay parameter \\\alpha\\. Types dictate
  the structural path constraint: OTP (outgoing two-paths,
  \\z\_{i,h}\\z\_{h,j}\\), ITP (incoming two-paths,
  \\z\_{h,i}\\z\_{j,h}\\), OSP (outgoing shared partners,
  \\z\_{i,h}\\z\_{j,h}\\), ISP (incoming shared partners,
  \\z\_{h,i}\\z\_{h,j}\\).

  - `global`: Path evaluations span the entire network topology.

  - `local`: Path evaluations are constrained exclusively to
    \\c\_{i,j}\\ overlapping neighborhoods.

- `cov_z_out(data, mode = "global")`:

  **Covariate Sender \[h-term\]:** Exogenous monadic covariate
  \\v\_{i}\\ influence on generating an outgoing tie.

  - `global`: \\h\_{i,j}(x,y,z) = v_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_i z\_{i,j}\\

- `cov_z_in(data, mode = "global")`:

  **Covariate Receiver \[h-term\]:** Exogenous monadic covariate
  \\v\_{j}\\ influence on receiving an incoming tie.

  - `global`: \\h\_{i,j}(x,y,z) = v_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_j z\_{i,j}\\

- `transitive`:

  **Transitivity (Local) \[Joint\]:** Indicator evaluating the presence
  of a local transitive triad configuration. An actor \\h \neq i,j\\
  must exist fulfilling \\h\in \mathcal{N}\_i, h\in \mathcal{N}\_j\\
  mapping to \\z\_{i,j} = z\_{i,h} = z\_{h,j}\\: \\h\_{i,j} = c\_{i,j}
  z\_{i,j} \mathbb{I}(\sum\_{k} c\_{i,k} c\_{j,k} z\_{i,k}
  z\_{k,j}\>1)\\.

- **3. \\h\_{i,j}\\ Terms for Joint Attribute/Network Dependence**:

- `attribute_xz(mode = "local")`:

  **Attribute Sum (X-Z) \[Joint h-term\]:** Models the additive effect
  of \\x_i\\ and \\x_j\\ on edge formation within local neighborhoods.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} (x_i + x_j) z\_{i,j}\\

- `attribute_yz(mode = "local")`:

  **Attribute Sum (Y-Z) \[Joint h-term\]:** Models the additive effect
  of \\y_i\\ and \\y_j\\ on edge formation within local neighborhoods.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} (y_i + y_j) z\_{i,j}\\

- `edges_x_match(mode = "global")`:

  **Attribute Match (X-Z) \[Joint h-term\]:** Models homophily/matching
  on the binary attribute \\x\\.

  - `global`: \\h\_{i,j}(x,y,z) = \mathbb{I}(x_i = x_j) z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} \mathbb{I}(x_i = x_j)
    z\_{i,j}\\

- `edges_y_match(mode = "global")`:

  **Attribute Match (Y-Z) \[Joint h-term\]:** Models homophily/matching
  on the binary attribute \\y\\.

  - `global`: \\h\_{i,j}(x,y,z) = \mathbb{I}(y_i = y_j) z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} \mathbb{I}(y_i = y_j)
    z\_{i,j}\\

- `outedges_x(mode = "global")`:

  **Attribute Out-Degree (X-Z) \[Joint h-term\]:** Influence of
  endogenous \\x_i\\ on out-degree formation.

  - `global`: \\h\_{i,j}(x,y,z) = x_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_i z\_{i,j}\\

- `inedges_x(mode = "global")`:

  **Attribute In-Degree (X-Z) \[Joint h-term\]:** Influence of
  endogenous \\x_j\\ on in-degree reception.

  - `global`: \\h\_{i,j}(x,y,z) = x_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_j z\_{i,j}\\

- `outedges_y(mode = "global")`:

  **Attribute Out-Degree (Y-Z) \[Joint h-term\]:** Influence of
  endogenous \\y_i\\ on out-degree formation.

  - `global`: \\h\_{i,j}(x,y,z) = y_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_i z\_{i,j}\\

- `inedges_y(mode = "global")`:

  **Attribute In-Degree (Y-Z) \[Joint h-term\]:** Influence of
  endogenous \\y_j\\ on in-degree reception.

  - `global`: \\h\_{i,j}(x,y,z) = y_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_j z\_{i,j}\\

- `spillover_xx(mode = "local")`:

  **Symmetric X-X-Z Outcome Spillover \[h-term\]:** Propagates
  \\x\\-outcome spillover effects through explicit network ties.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i x_j z\_{i,j}\\

- `spillover_xx_scaled(mode = "local")`:

  **Scaled X-X-Z Outcome Spillover \[h-term\]:** Normalizes the
  \\x\\-outcome spillover influence by the relevant out-degree topology
  of actor \\i\\.

  - `global`: \\h\_{i,j}(x,y,z) = x_i x_j z\_{i,j} / \text{deg}(i)\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i x_j z\_{i,j} /
    \text{deg}(i, \text{local})\\

- `spillover_yy(mode = "local")`:

  **Symmetric Y-Y-Z Outcome Spillover \[h-term\]:** Propagates
  \\y\\-outcome spillover effects through explicit network ties.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i y_j z\_{i,j}\\

- `spillover_yy_scaled(mode = "local")`:

  **Scaled Y-Y-Z Outcome Spillover \[h-term\]:** Normalizes the
  \\y\\-outcome spillover influence by the relevant out-degree topology
  of actor \\i\\.

  - `global`: \\h\_{i,j}(x,y,z) = y_i y_j z\_{i,j} / \text{deg}(i)\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i y_j z\_{i,j} /
    \text{deg}(i, \text{local})\\

- `spillover_xy(mode = "local")`:

  **Directed X-Y-Z Treatment Spillover \[h-term\]:** Maps the
  cross-attribute \\x_i \to y_j\\ treatment assignment across the
  network graph.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j}\\

- `spillover_xy_scaled(mode = "local")`:

  **Scaled X-Y-Z Treatment Spillover \[h-term\]:** Normalizes the
  cross-attribute \\x_i \to y_j\\ spillover influence by the relevant
  out-degree topology of actor \\i\\.

  - `global`: \\h\_{i,j}(x,y,z) = x_i y_j z\_{i,j} / \text{deg}(i)\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j} /
    \text{deg}(i, \text{local})\\

- `spillover_yx(mode = "local")`:

  **Directed Y-X-Z Treatment Spillover \[h-term\]:** Maps the
  cross-attribute \\y_i \to x_j\\ treatment assignment across the
  network graph.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j}\\

- `spillover_yx_scaled(mode = "local")`:

  **Scaled Y-X-Z Treatment Spillover \[h-term\]:** Normalizes the
  cross-attribute \\y_i \to x_j\\ spillover influence by the relevant
  out-degree topology of actor \\i\\.

  - `global`: \\h\_{i,j}(x,y,z) = y_i x_j z\_{i,j} / \text{deg}(i)\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j} /
    \text{deg}(i, \text{local})\\

- `spillover_yc(mode = "local")`:

  **Directed Y-C-Z Treatment Spillover \[h-term\]:** Specifies exogenous
  covariate \\v\\ interacting with endogenous trait \\y\\ over the
  network.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i v_j z\_{i,j}\\

- `spillover_yc_symm(data = v, mode = "local")`:

  **Symmetric Y-C-Z Treatment Spillover \[h-term\]:** Bidirectional
  mapping of the exogenous covariate \\v\\ and endogenous trait \\y\\
  interaction.

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} (v_i y_j + v_j y_i) z\_{i,j}\\

## References

Fritz, C., Schweinberger, M., Bhadra, S., and D.R. Hunter (2025). A
Regression Framework for Studying Relationships among Attributes under
Network Interference. Journal of the American Statistical Association,
to appear.

Schweinberger, M. and M.S. Handcock (2015). Local Dependence in Random
Graph Models: Characterization, Properties, and Statistical Inference.
Journal of the Royal Statistical Society, Series B (Statistical
Methodology), 7, 647-676.

Schweinberger, M. and J.R. Stewart (2020). Concentration and Consistency
Results for Canonical and Curved Exponential-Family Models of Random
Graphs. The Annals of Statistics, 48, 374-396.

Stewart, J.R. and M. Schweinberger (2025). Pseudo-Likelihood-Based
M-Estimation of Random Graphs with Dependent Edges and Parameter Vectors
of Increasing Dimension. The Annals of Statistics, to appear.

# Model specification for a \`iglm' object

`R` package `iglm` implements generalized linear models (GLMs) for
studying relationships among attributes in connected populations, where
responses of connected units can be dependent. It extends GLMs for
independent responses to dependent responses and can be used for
studying spillover in connected populations and other network-mediated
phenomena. It is based on a joint probability model for dependent
responses (\\Y\\) and connections \\(Z)\\ conditional on predictors (X).

The joint probability density is specified as \$\$f\_{\theta}(y,z,x)
\propto \Big\[\prod\_{i=1}^{N} a_y(y_i) \exp(\theta_g^T g_i(x_i,
y_i^\*)) \Big\] \times \Big\[\prod\_{i \ne j} a_z(z\_{i,j})
\exp(\theta_h^T h\_{i,j}(x, y_i^\*, y_j^\*, z)) \Big\],\$\$ which is
defined by two distinct sets of user-specified features:

- **\\g_i(x,y,z)\\**: A vector of unit-level functions (or "g-terms")
  that describe the relationship between an individual actor \\i\\'s
  predictors (\\x_i\\) and their own response (\\y_i\\).

- **\\h\_{i,j}(x,y,z)\\**: A vector of pair-level functions (or
  "h-terms") that specify how the connections (\\z\\) and responses
  (\\y_i, y_j\\) of a pair of units \\\\i,j\\\\ depend on each other and
  the wider network structure.

This separation allows the model to simultaneously capture
individual-level behavior (via \\g_i\\) and dyadic, network-based
dependencies (via \\h\_{i,j}\\), including local dependence limited to
overlapping neighborhoods. This help page documents the various
statistics available in 'iglm', corresponding to the \\g_i\\
(attribute-level) and \\h\_{i,j}\\ (pair-level) components of the joint
model. In the formula interface, these terms can be specified by adding
them in the right-hand side of a model formula, e.g.,
`netplus_object ~ attribute_x + edges(mode = "local") + popularity, ...`.
See the documentation for [`iglm`](iglm.md) for details on model fitting
and estimation.

## Details

Each term defines a component for the model's features, which are a sum
of unit-level components, \\\sum_i g_i(x,y,z)\\, and/or pair-level
components, \\\sum\_{i \ne j} h\_{i,j}(x,y,z)\\. Here, \\x_i\\ and
\\y_i\\ are the attributes for actor \\i\\, and \\z\_{i,j}\\ indicates
the presence (1) or absence (0) of a tie from actor \\i\\ to actor
\\j\\. The local neighborhood of actor \\i\\ is denoted
\\\mathcal{N}\_i\\, and the indicator for whether actors \\i\\ and \\j\\
share a local neighborhood is given by \\c\_{i,j} =
\mathbb{I}(\mathcal{N}\_i \cap \mathcal{N}\_j \neq \emptyset)\\. The
functions below specify the forms of \\g_i(x,y,z)\\ and
\\h\_{i,j}(x,y,z)\\ for each term. Some terms also depend on other
covariates, which are denoted by \\v = (v_1, ..., v_N)\\ (unit-level)
and \\w = (w\_{i,j}) \in \mathbb{R}^{N \times N}\\(dyadic). These
covariates must be provided by the user via the `data` argument. The
implemented terms are grouped into three categories:

1.  \\g_i\\ terms for attribute dependence,

2.  \\h\_{i,j}\\ terms for network dependence,

3.  \\h\_{i,j}\\ Terms for joint attribute/network dependence.

- **1. \\g_i\\ Terms for Attribute Dependence**:

- `attribute_x`:

  **Attribute (X) \[g-term\]:** Intercept for attribute 'x'.
  \\g_i(x,y,z) = x_i\\

- `attribute_y`:

  **Attribute (Y) \[g-term\]:** Intercept for attribute 'y'.
  \\g_i(x,y,z) = y_i\\

- `cov_x`:

  **Nodal Covariate (X) \[g-term\]:** Effect of a unit-level covariate
  \\v_i\\ on attribute \\x_i\\. \\g_i(x,y,z) = v_i x_i\\

- `cov_y(data = v)`:

  **Nodal Covariate (Y) \[g-term\]:** Effect of a unit-level covariate
  \\v_i\\ on attribute \\y_i\\. \\g_i(x,y,z) = v_i y_i\\

- `attribute_xy(mode = "global")`:

  **Nodal Attribute Interaction (X-Y) \[g-term\]:** Interaction of
  attributes \\x_i\\ and \\y_i\\ on the same node. For mode different
  from "global", we count interactions of an actor's attributes with
  their local neighbors' attributes.

  - `global`: \\g_i(x,y,z) = x_i y_i\\

  - `local`: \\g_i(x,y,z) = x_i \sum\_{j \in \mathcal{N}\_i} y_j + y_i
    \sum\_{j \in \mathcal{N}\_i} x_j\\

  - `alocal`: \\g_i(x,y,z) = x_i \sum\_{j \notin \mathcal{N}\_i} y_j +
    y_i \sum\_{j \notin \mathcal{N}\_i} x_j\\

- **2. \\h\_{i,j}\\ Terms for Network Dependence**:

- `popularity`:

  **Popularity \[h-term\]:** Adds fixed effects for all actors in the
  network. Estimation of popularity effects is carried out using a MM
  algorithm. For directed networks, each actors has a sender and
  receiver effect (we assume that the out effect of actor N is 0 for
  identifiability). For undirected networks, each actor has a single
  popularity effect.

- `edges(mode = "global")`:

  **Edges \[h-term\]:** Counts different types of edges.

  - `global`: \\h\_{i,j}(x,y,z) = z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j}\\

- `mutual(mode = "global")`:

  **Mutual Reciprocity \[h-term\]:** Counts whether the reciprocal tie
  between actors \\i\\ and \\j\\ is present. This term should only be
  used for directed networks.

  - `global`: \\h\_{i,j}(x,y,z) = z\_{i,j} z\_{j,i}\\ (for \\i \< j\\)

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j} z\_{j,i}\\ (for \\i
    \< j\\)

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j} z\_{j,i}\\
    (for \\i \< j\\)

- `cov_z(data, mode = "global")`:

  **Dyadic Covariate \[h-term\]:** The effect of a dyadic covariate
  \\w\_{i,j}\\ for directed or undirected networks.

  - `global`: \\h\_{i,j}(x,y,z) = w\_{i,j} z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} w\_{i,j} z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) w\_{i,j} z\_{i,j}\\

- `cov_z_out(data, mode = "global")`:

  **Covariate Sender \[h-term\]:** The effect of a monadic covariate
  \\v\_{i}\\ on being the sender in a directed network.

  - `global`: \\h\_{i,j}(x,y,z) = v_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_i z\_{i,j}\\

- `cov_z_in(data, mode = "global")`:

  **Covariate Receiver \[h-term\]:** The effect of a monadic covariate
  \\v\_{i}\\ on being the receiver in a directed network.

  - `global`: \\h\_{i,j}(x,y,z) = v_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_j z\_{i,j}\\

- `transitive`:

  **Transitivity (Local) \[Joint\]:** A statistic checking whether the
  dyad is a local transitive edge, meaning that there exists an actor
  \\h \neq i,j\\ such that \\h\in \mathcal{N}\_i, h\in \mathcal{N}\_j\\
  with \\z\_{i,j} = z\_{i,h} = z\_{h,j}\\: \\h\_{i,j} = c\_{i,j}
  z\_{i,j} \mathbb{I}(\sum\_{k} c\_{i,k} c\_{j,k} z\_{i,k}
  z\_{k,j}\>1)\\

- **3. \\h\_{i,j}\\ Terms for Joint Attribute/Network Dependence**:

- `outedges_x_global()`:

  **Attribute Out-Degree (X-Z Global) \[h-term\]:** Models \\x_i\\'s
  effect on its out-degree. Corresponds to \\h\_{i,j}(x,y,z) = x_i
  z\_{i,j}\\.

- `outedges_x(mode = "global")`:

  **Attribute Out-Degree (X-Z) \[Joint h-term\]:** Models \\x_i\\'s
  effect on its out-degree.

  - `global`: \\h\_{i,j}(x,y,z) = x_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_i z\_{i,j}\\

- `inedges_x(mode = "global")`:

  **Attribute In-Degree (X-Z) \[Joint h-term\]:** Models \\x_j\\'s
  effect on its in-degree.

  - `global`: \\h\_{i,j}(x,y,z) = x_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_j z\_{i,j}\\

- `outedges_y(mode = "global")`:

  **Attribute Out-Degree (Y-Z) \[Joint h-term\]:** Models \\y_i\\'s
  effect on its out-degree.

  - `global`: \\h\_{i,j}(x,y,z) = y_i z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_i z\_{i,j}\\

- `inedges_y(mode = "global")`:

  **Attribute In-Degree (Y-Z) \[Joint h-term\]:** Models \\y_j\\'s
  effect on its in-degree.

  - `global`: \\h\_{i,j}(x,y,z) = y_j z\_{i,j}\\

  - `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_j z\_{i,j}\\

  - `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_j z\_{i,j}\\

- `spillover_xx`:

  **Symmetric X-X-Z Outcome Spillover \[h-term\]:** Models \\x\\-outcome
  spillover **within** the local neighborhood. Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} x_i x_j z\_{i,j}\\.

- `spillover_yy`:

  **Symmetric Y-Y-Z Outcome Spillover \[h-term\]:** Models \\y\\-outcome
  spillover **within** the local neighborhood. Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} y_i y_j z\_{i,j}\\.

- `spillover_xy`:

  **Directed X-Y-Z Treatment Spillover \[h-term\]:** Models the \\x_i
  \to y_j\\ treatment spillover **within** the local neighborhood.
  Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j}\\.

- `spillover_xy_symm`:

  **Symmetric X-Y-Z Treatment Spillover \[h-term\]:** Models the \\x_i
  \leftrightarrow y_j\\ treatment spillover **within** the local
  neighborhood. Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} (x_i y_j +
  x_j y_i) z\_{i,j}\\.

- `spillover_yx`:

  **Directed Y-X-Z Treatment Spillover \[h-term\]:** Models the \\y_i
  \to x_j\\ treatment spillover **within** the local neighborhood.
  Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j}\\.

- `spillover_yc`:

  **Directed Y-C-Z Treatment Spillover \[h-term\]:** Models \\y\\-treat
  spillover to a covariate \\v\\ **within** the local neighborhood.
  Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} y_i v_j z\_{i,j}\\.

- `spillover_yc_symm(data = v)`:

  **Symmetric Treatment Spillover \[h-term\]:** Models the \\v_i
  \leftrightarrow y_j \\ treatment spillover . Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} (v_i y_j + v_j y_i) z\_{i,j}\\.

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

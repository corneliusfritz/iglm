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

  **Attribute (X) \[g-term\]:** Intercept for attribute 'x'.
  \\g_i(x_i,y_i) = x_i\\

- `attribute_y`:

  **Attribute (Y) \[g-term\]:** Intercept for attribute 'y'.
  \\g_i(x_i,y_i) = y_i\\

- `cov_x`:

  **Nodal Covariate (X) \[g-term\]:** Effect of a unit-level covariate
  \\v_i\\ on attribute \\x_i\\. \\g_i(x_i,y_i) = v_i x_i\\

- `cov_y(data = v)`:

  **Nodal Covariate (Y) \[g-term\]:** Effect of a unit-level covariate
  \\v_i\\ on attribute \\y_i\\. \\g_i(x_i,y_i) = v_i y_i\\

- `attribute_xy(mode = "global")`:

  **Nodal Attribute Interaction (X-Y) \[g-term\]:** Interaction of
  attributes \\x_i\\ and \\y_i\\ on the same node. For mode different
  from "global", we count interactions of an actor's attributes with
  their local neighbors' attributes.

  - `global`: \\g_i(x_i,y_i) = x_i y_i\\

  - `local`: \\g_i(x_i,y_i) = x_i \sum\_{j \in \mathcal{N}\_i} y_j + y_i
    \sum\_{j \in \mathcal{N}\_i} x_j\\

  - `alocal`: \\g_i(x_i,y_i) = x_i \sum\_{j \notin \mathcal{N}\_i} y_j +
    y_i \sum\_{j \notin \mathcal{N}\_i} x_j\\

- **2. \\h\_{i,j}\\ Terms for Network Dependence**:

- `degrees`:

  **degrees \[h-term\]:** Adds fixed effects for all actors in the
  network. Estimation of degrees effects is carried out using a MM
  algorithm. For directed networks, each actors has a sender and
  receiver effect (we assume that the out effect of actor N is 0 for
  identifiability). For undirected networks, each actor has a single
  degrees effect.

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

- `isolates`:

  **Isolates \[z-term\]:** Counts and accounts for the number of
  non-isolated nodes.

- `nonisolates`:

  **Non-Isolates \[z-term\]:** Counts and accounts for the number of
  non-isolated nodes. It is the exact negative of the `isolates`
  statistic.

- `gwodegree(decay)`:

  **Geometrically Weighted Out-Degree \[z-term\]:** The Geometrically
  Weighted Out-Degree statistic is implemented as in the \`ergm\`
  package.

- `gwidegree(decay)`:

  **Geometrically Weighted In-Degree \[z-term\]:** The Geometrically
  Weighted In-Degree (GWIDegree) statistic is implemented as in the
  \`ergm\` package.

- `gwesp(data, mode = "global", variant = "OTP")`:

  **Geometrically Weighted Edegewise-Shared Partners \[h-term\]:**
  Geometrically weighted edgewise shared partners (GWESP) statistic for
  directed networks as implemented in the \`ergm\` package. Variants
  include: OTP (outgoing two-paths, \\z\_{i,h}\\z\_{h,j}\\ z\_{i,j}\\),
  ITP (incoming two-paths, \\z\_{h,i}\\z\_{j,h}\\ z\_{i,j}\\), OSP
  (outgoing shared partners, \\z\_{i,h}\\z\_{j,h}\\ z\_{i,j}\\), ISP
  (incoming shared partners, \\z\_{h,i}\\z\_{h,j}\\ z\_{i,j}\\).

  - `global`: ESP counts are calculated over all edges in the network.

  - `local`: ESP counts are restricted to local edges only (edges with
    non-overlapping neighborhoods).

- `gwdsp(data, mode = "global", variant = "OTP")`:

  **Geometrically Weighted Dyadwise-Shared Partners \[h-term\]:**
  Geometrically weighted dyadwise shared partners (GWDSP) statistic for
  directed networks as implemented in the \`ergm\` package. Variants
  include: OTP (outgoing two-paths, \\z\_{i,h}\\z\_{h,j}\\), ITP
  (incoming two-paths, \\z\_{h,i}\\z\_{j,h}\\), OSP (outgoing shared
  partners, \\z\_{i,h}\\z\_{j,h}\\), ISP (incoming shared partners,
  \\z\_{h,i}\\z\_{h,j}\\).

  - `global`: ESP counts are calculated over all edges in the network.

  - `local`: ESP counts are restricted to local edges only (edges with
    non-overlapping neighborhoods).

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

- `spillover_xx_scaled`:

  **X-X-Z Outcome Spillover \[h-term\]:** Models \\x\\-outcome spillover
  **within** the local neighborhood but weights the influence of \\x_j\\
  on \\x_i\\ by the out-degree of actor \\i\\ with other actors in its
  neighborhood, denoted by \\\text{local\\degree(i)} (for undirected
  networks, the degree is used)\\. Corresponds to \\h\_{i,j}(x,y,z) =
  c\_{i,j} x_i x_j z\_{i,j} / \text{local\\degree(i)}\\.

- `spillover_yy`:

  **Symmetric Y-Y-Z Outcome Spillover \[h-term\]:** Models \\y\\-outcome
  spillover **within** the local neighborhood. Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} y_i y_j z\_{i,j}\\.

- `spillover_yy_scaled`:

  **Y-Y-Z Outcome Spillover \[h-term\]:** Models \\y\\-outcome spillover
  **within** the local neighborhood but weights the influence of \\y_j\\
  on \\y_i\\ by the degree of actor \\i\\ with other actors in its
  neighborhood, defined above. Corresponds to \\h\_{i,j}(x,y,z) =
  c\_{i,j} y_i y_j z\_{i,j} / \text{local\\degree(i)}\\.

- `spillover_xy`:

  **Directed X-Y-Z Treatment Spillover \[h-term\]:** Models the \\x_i
  \to y_j\\ treatment spillover **within** the local neighborhood.
  Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j}\\.

- `spillover_xy_scaled`:

  **X-Y-Z Outcome Spillover \[h-term\]:** Models the \\x_i \to y_j\\
  treatment spillover **within** the local neighborhood but weights the
  influence of \\y_j\\ on \\x_i\\ by the degree of actor \\i\\ with
  other actors in its neighborhood, defined above. Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j} /
  \text{local\\degree(i)}\\.

- `spillover_xy_symm`:

  **Symmetric X-Y-Z Treatment Spillover \[h-term\]:** Models the \\x_i
  \leftrightarrow y_j\\ treatment spillover **within** the local
  neighborhood. Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} (x_i y_j +
  x_j y_i) z\_{i,j}\\.

- `spillover_yx`:

  **Directed Y-X-Z Treatment Spillover \[h-term\]:** Models the \\y_i
  \to x_j\\ treatment spillover **within** the local neighborhood.
  Corresponds to \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j}\\.

- `spillover_yx_scaled`:

  **Y-X-Z Outcome Spillover \[h-term\]:** Models the \\y_i \to x_j\\
  treatment spillover **within** the local neighborhood but weights the
  influence of \\x_j\\ on \\y_i\\ by the degree of actor \\i\\ with
  other actors in its neighborhood, defined above. Corresponds to
  \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j} /
  \text{local\\degree(i)}\\.

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

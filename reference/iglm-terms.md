# Model Specification for iglm Terms

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

1.  **Attribute Terms**: Depend only on individual attributes \\x_i\\ or
    \\y_i\\.

2.  **Network Terms**: Depend only on the connections \\z\_{i,j}\\.

3.  **Joint Attribute/Network Terms**: Depend on both individual
    attributes and connections.

`degrees`: Degrees: Specifies node-level fixed effects. Estimation
requires an MM algorithm constraint.

`edges(mode = "global")`: Edges: Captures the baseline propensity of tie
formation \\z\_{i,j}\\, partitioned by structural boundary \\c\_{i,j}\\.

- `global`: \\h\_{i,j}(x,y,z) = z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j}\\

`mutual(mode = "global")`: Mutual Reciprocity: Evaluates reciprocal tie
formation in directed networks.

- `global`: \\h\_{i,j}(x,y,z) = z\_{i,j} z\_{j,i}\\ (for \\i \< j\\)

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} z\_{i,j} z\_{j,i}\\ (for \\i \<
  j\\)

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) z\_{i,j} z\_{j,i}\\ (for
  \\i \< j\\)

`cov_z(data, mode = "global")`: Dyadic Covariate: Exogenous dyadic
covariate \\w\_{i,j}\\ influence on edge formation.

- `global`: \\h\_{i,j}(x,y,z) = w\_{i,j} z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} w\_{i,j} z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) w\_{i,j} z\_{i,j}\\

`cov_z_out(data, mode = "global")`: Covariate Sender: Exogenous monadic
covariate \\v\_{i}\\ influence on generating an outgoing tie.

- `global`: \\h\_{i,j}(x,y,z) = v_i z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_i z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_i z\_{i,j}\\

`cov_z_in(data, mode = "global")`: Covariate Receiver: Exogenous monadic
covariate \\v\_{j}\\ influence on receiving an incoming tie.

- `global`: \\h\_{i,j}(x,y,z) = v_j z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} v_j z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) v_j z\_{i,j}\\

`cov_x(data = v)`: Nodal Covariate (X): Effect of a unit-level exogenous
covariate \\v_i\\ on endogenous attribute \\x_i\\. \\g_i(x_i,y_i) = v_i
x_i\\

`cov_y(data = v)`: Nodal Covariate (Y): Effect of a unit-level exogenous
covariate \\v_i\\ on endogenous attribute \\y_i\\. \\g_i(x_i,y_i) = v_i
y_i\\

`attribute_xy(mode = "global")`: Nodal Attribute Interaction (X-Y):
Interaction of attributes \\x_i\\ and \\y_i\\.

- `global`: \\g_i(x_i,y_i) = x_i y_i\\

- `local`: \\g_i(x_i,y_i) = x_i \sum\_{j \in \mathcal{N}\_i} y_j + y_i
  \sum\_{j \in \mathcal{N}\_i} x_j\\

- `alocal`: \\g_i(x_i,y_i) = x_i \sum\_{j \notin \mathcal{N}\_i} y_j +
  y_i \sum\_{j \notin \mathcal{N}\_i} x_j\\

`attribute_yz(mode = "local")`: Attribute Sum (Y-Z): Models the additive
effect of \\y_i\\ and \\y_j\\ on edge formation within local
neighborhoods.

`attribute_xz(mode = "local")`: Attribute Sum (X-Z): Models the additive
effect of \\x_i\\ and \\x_j\\ on edge formation within local
neighborhoods.

`inedges_y(mode = "global")`: Attribute In-Degree (Y-Z): Influence of
endogenous \\y_j\\ on in-degree reception.

- `global`: \\h\_{i,j}(x,y,z) = y_j z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_j z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_j z\_{i,j}\\

`outedges_y(mode = "global")`: Attribute Out-Degree (Y-Z): Influence of
endogenous \\y_i\\ on out-degree formation.

- `global`: \\h\_{i,j}(x,y,z) = y_i z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) y_i z\_{i,j}\\

`inedges_x(mode = "global")`: Attribute In-Degree (X-Z): Influence of
endogenous \\x_j\\ on in-degree reception.

- `global`: \\h\_{i,j}(x,y,z) = x_j z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_j z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_j z\_{i,j}\\

`outedges_x(mode = "global")`: Attribute Out-Degree (X-Z): Influence of
endogenous \\x_i\\ on out-degree formation.

- `global`: \\h\_{i,j}(x,y,z) = x_i z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i z\_{i,j}\\

- `alocal`: \\h\_{i,j}(x,y,z) = (1 - c\_{i,j}) x_i z\_{i,j}\\

`attribute_x`: Attribute (X): Intercept for attribute \\x\\.
\\g_i(x_i,y_i) = x_i\\

`attribute_y`: Attribute (Y): Intercept for attribute \\y\\.
\\g_i(x_i,y_i) = y_i\\

`edges_x_match(mode = "global")`: Attribute Match (X-Z): Models
homophily/matching on the binary attribute \\x\\.

- `global`: \\h\_{i,j}(x,y,z) = \mathbb{I}(x_i = x_j) z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} \mathbb{I}(x_i = x_j) z\_{i,j}\\

`edges_y_match(mode = "global")`: Attribute Match (Y-Z): Models
homophily/matching on the binary attribute \\y\\.

- `global`: \\h\_{i,j}(x,y,z) = \mathbb{I}(y_i = y_j) z\_{i,j}\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} \mathbb{I}(y_i = y_j) z\_{i,j}\\

`spillover_yy_scaled(mode = "global")`: Scaled Y-Y-Z Outcome Spillover:
Normalizes the \\y\\-outcome spillover influence by the relevant
out-degree topology.

- `global`: \\h\_{i,j}(x,y,z) = y_i y_j z\_{i,j} / \text{deg}(i)\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i y_j z\_{i,j} / \text{deg}(i,
  \text{local})\\

`spillover_xx_scaled(mode = "global")`: Scaled X-X-Z Outcome Spillover:
Normalizes the \\x\\-outcome spillover influence by the relevant
out-degree topology.

- `global`: \\h\_{i,j}(x,y,z) = x_i x_j z\_{i,j} / \text{deg}(i)\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i x_j z\_{i,j} / \text{deg}(i,
  \text{local})\\

`spillover_yx_scaled(mode = "global")`: Scaled Y-X-Z Treatment
Spillover: Normalizes cross-attribute \\y_i \to x_j\\ spillover
influence.

- `global`: \\h\_{i,j}(x,y,z) = y_i x_j z\_{i,j} / \text{deg}(i)\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} y_i x_j z\_{i,j} / \text{deg}(i,
  \text{local})\\

`spillover_xy_scaled(mode = "global")`: Scaled X-Y-Z Treatment
Spillover: Normalizes cross-attribute \\x_i \to y_j\\ spillover
influence.

- `global`: \\h\_{i,j}(x,y,z) = x_i y_j z\_{i,j} / \text{deg}(i)\\

- `local`: \\h\_{i,j}(x,y,z) = c\_{i,j} x_i y_j z\_{i,j} / \text{deg}(i,
  \text{local})\\

`gwesp(data, mode = "global", variant = "OSP", decay = 0)`:
Geometrically Weighted Edgewise-Shared Partners: Models triadic closure
propensity conditioning on existing edges. Types dictate path
constraint: OTP, ITP, OSP, ISP for directed; symm for undirected.

`gwdsp(data, mode = "global", variant = "OSP", decay = 0)`:
Geometrically Weighted Dyadwise-Shared Partners: Models triadic
potential irrespective of the closing edge. Types dictate path
constraint: OTP, ITP, OSP, ISP for directed; symm for undirected.

`gwdegree(mode = "global", decay = 0)`: Geometrically Weighted Degree:
Captures the degree distribution utilizing an exponential decay
parameter.

`gwidegree(mode = "global", decay = 0)`: Geometrically Weighted
In-Degree: Captures the in-degree distribution utilizing an exponential
decay parameter.

`gwodegree(mode = "global", decay = 0)`: Geometrically Weighted
Out-Degree: Captures the out-degree distribution utilizing an
exponential decay parameter.

`spillover_yc_symm(data = v, mode = "local")`: Symmetric Y-C-Z Treatment
Spillover: Bidirectional mapping of exogenous covariate \\v\\ and
endogenous trait \\y\\ interaction.

`spillover_xy(mode = "local")`: Directed X-Y-Z Treatment Spillover: Maps
cross-attribute \\x_i \to y_j\\ treatment assignment.

`spillover_yc(mode = "local")`: Directed Y-C-Z Treatment Spillover:
Exogenous covariate \\v\\ interacting with endogenous trait \\y\\.

`spillover_yx(mode = "local")`: Directed Y-X-Z Treatment Spillover: Maps
cross-attribute \\y_i \to x_j\\ treatment assignment.

`spillover_yy(mode = "local")`: Symmetric Y-Y-Z Outcome Spillover:
Propagates \\y\\-outcome spillover effects.

`spillover_xx(mode = "local")`: Symmetric X-X-Z Outcome Spillover:
Propagates \\x\\-outcome spillover effects.

`transitive`: Transitivity (Local): Indicator evaluating the presence of
a local transitive triad configuration.

`nonisolates`: Non-Isolates: Captures frequency of nodes with degree
strictly greater than zero.

`isolates`: Isolates: Captures frequency of nodes with degree zero.

## Category 1

Attribute Terms:

Below is a detailed description of terms that depend only on nodal
attributes:

- `attribute_x-term`, `attribute_y-term`

- `cov_x-term`, `cov_y-term`

- `attribute_xy-term`

## Category 2

Network Terms:

Below is a detailed description of terms that depend only on the network
structure:

- `edges-term`, `mutual-term`

- `cov_z-term`, `cov_z_in-term`, `cov_z_out-term`

- `degrees-term`

- `gwdegree-term`, `gwidegree-term`, `gwodegree-term`

- `gwesp-term`, `gwdsp-term`

- `transitive-term`, `nonisolates-term`, `isolates-term`

## Category 3

Joint Attribute/Network Terms:

Below is a detailed description of terms that depend on both attributes
and the network:

- `attribute_xz-term`, `attribute_yz-term`

- `inedges_x-term`, `inedges_y-term`, `outedges_x-term`,
  `outedges_y-term`

- `edges_x_match-term`, `edges_y_match-term`

- `spillover_xx-term`, `spillover_yy-term`

- `spillover_yx-term`, `spillover_xy-term`, `spillover_yc-term`,
  `spillover_yc_symm-term`

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

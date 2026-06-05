# Package index

## Model Estimation

Main estimation functions and controls for fitting joint models of
networks and attributes.

- [`iglm()`](https://corneliusfritz.github.io/iglm/reference/iglm.md) :
  Construct an iglm Model Specification Object
- [`control.iglm()`](https://corneliusfritz.github.io/iglm/reference/control.iglm.md)
  : Set Control Parameters for iglm Estimation
- [`iglm.object.generator`](https://corneliusfritz.github.io/iglm/reference/iglm.object.generator.md)
  : iglm Objects (R6 Class)

## Data Preparation & Setup

Functions and R6 classes to manage network and attribute data.

- [`iglm.data()`](https://corneliusfritz.github.io/iglm/reference/iglm.data.md)
  : Constructor for the iglm.data R6 object
- [`iglm.data_generator`](https://corneliusfritz.github.io/iglm/reference/iglm.data_generator.md)
  : Networks with Unit-Level Attributes (R6 Class)

## Model Simulation

Functions to simulate network and attribute outcomes from fitted models.

- [`simulate_iglm()`](https://corneliusfritz.github.io/iglm/reference/simulate_iglm.md)
  : Simulate Responses and Connections

## MCMC Samplers

Sampler configurations and classes for network and attribute variables.

- [`sampler.iglm()`](https://corneliusfritz.github.io/iglm/reference/sampler.iglm.md)
  : Constructor for a iglm Sampler
- [`sampler.iglm.generator`](https://corneliusfritz.github.io/iglm/reference/sampler.iglm.generator.md)
  : iglm Sampler Settings (R6 Class)
- [`sampler.net.attr()`](https://corneliusfritz.github.io/iglm/reference/sampler.net.attr.md)
  : Constructor for Single Component Sampler Settings
- [`sampler.net.attr.generator`](https://corneliusfritz.github.io/iglm/reference/sampler.net.attr.generator.md)
  : Single Component Sampler Settings (R6 Class)

## Diagnostics & Results

Classes and functions for retrieving results, statistics, and
assessment.

- [`results()`](https://corneliusfritz.github.io/iglm/reference/results.md)
  : Constructor for the results R6 Object
- [`results.generator`](https://corneliusfritz.github.io/iglm/reference/results.generator.md)
  : iglm Estimation and Simulation Results (R6 Class)
- [`statistics()`](https://corneliusfritz.github.io/iglm/reference/statistics.md)
  : Compute Statistics

## Model Specification & Custom Terms

Utilities for defining custom terms and analyzing sufficient statistics.

- [`iglm-terms`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`iglm.terms`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`degrees-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`edges-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`mutual-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`cov_z-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`cov_z_out-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`cov_z_in-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`cov_x-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`cov_y-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`attribute_xy-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`attribute_yz-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`attribute_xz-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`inedges_y-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`outedges_y-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`inedges_x-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`outedges_x-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`attribute_x-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`attribute_y-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`edges_x_match-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`edges_y_match-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yy_scaled-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_xx_scaled-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yx_scaled-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_xy_scaled-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`gwesp-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`gwdsp-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`gwdegree-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`gwidegree-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`gwodegree-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yc_symm-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_xy-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yc-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yx-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_yy-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`spillover_xx-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`transitive-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`nonisolates-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  [`isolates-term`](https://corneliusfritz.github.io/iglm/reference/iglm-terms.md)
  : Model Specification for iglm Terms
- [`create_userterms_skeleton()`](https://corneliusfritz.github.io/iglm/reference/create_userterms_skeleton.md)
  : Generate the Skeleton for an R Package to Implement Additional iglm
  Terms

## Datasets

Example network and attribute datasets included in the package.

- [`copenhagen`](https://corneliusfritz.github.io/iglm/reference/copenhagen.md)
  : Copenhagen Network Study
- [`state_twitter`](https://corneliusfritz.github.io/iglm/reference/state_twitter.md)
  : Twitter (X) data list for U.S. state legislators (10-state subset)

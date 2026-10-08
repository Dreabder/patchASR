
<!-- README.md is generated from README.Rmd. Please edit README.Rmd, not README.md. -->

# patchASR

[![R-CMD-check](https://github.com/Dreabder/patchASR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Dreabder/patchASR/actions/workflows/R-CMD-check.yaml)

`patchASR` implements patch-aware ancestral-state reconstruction for
continuous traits on rooted phylogenetic trees containing one or more
specified branch-localized shifts.

The package is designed for situations in which a localized large trait
change occurs on a particular phylogenetic branch and the downstream
observations may influence ancestral-state reconstruction elsewhere in
the tree.

The central idea is to separate the shift-associated descendant region
from the remainder of the phylogeny, perform ancestral-state
reconstruction on the resulting components, reconstruct the boundary
node connecting those components, and finally merge all estimates back
onto the internal nodes of the original tree.

## Development installation

`patchASR` is currently distributed as a development version from
GitHub.

It can be installed with:

``` r
remotes::install_github("Dreabder/patchASR")
```

or, if `devtools` is already installed:

``` r
devtools::install_github("Dreabder/patchASR")
```

For package development from a local source checkout, the working copy
can be loaded with:

``` r
devtools::load_all(".")
```

## Minimal example

The following example uses a simple toy ancestral-state reconstruction
backend to demonstrate the `patchASR` interface.

The toy backend is used only for illustrating package mechanics; it is
not intended as a biological ancestral-state reconstruction method.

``` r
library(ape)

tree <- read.tree(
  text = "((((A:1,B:1)J:1,C:1)I:1,D:1)H:1,(E:1,F:1)K:1)ROOT;")

states <- c(A = 1, B = 2, C = 3, D = 4, E = 5, F = 6)

I_node <- Ntip(tree) + match("I", tree$node.label)
J_node <- Ntip(tree) + match("J", tree$node.label)

shifts <- data.frame(parent = I_node, child = J_node)

mean_backend <- function(tree, states) {
  stopifnot(identical(names(states), tree$tip.label))

  stats::setNames(
    rep(mean(states), tree$Nnode),
    tree$node.label)}

fit <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = mean_backend)

fit
#> <patch_asr>
#>   Internal nodes: 5
#>   Patches: 1
#>   Boundary nodes: 1
#>   Boundary residual: 0e+00
```

The final ancestral-state estimates can be inspected with:

``` r
fit$ancestral_states
#>      node original_node original_label estimate                  source
#> 1  node_7             7           ROOT     4.50           component_asr
#> 2  node_8             8              H     4.50           component_asr
#> 3  node_9             9              I     3.75 boundary_reconstruction
#> 4 node_10            10              J     1.50           component_asr
#> 5 node_11            11              K     4.50           component_asr
#>   component
#> 1    mother
#> 2    mother
#> 3  boundary
#> 4   patch_1
#> 5    mother
```

The result records both the numerical estimate and its provenance.
Component-level estimates use `source = "component_asr"`, whereas
reconstructed boundary nodes use `source = "boundary_reconstruction"`.

Component provenance is recorded separately, for example as `mother`,
`patch_1`, or `boundary`.

## Core workflow

For a specified shift edge from parent node `I` to child node `J`,
`patchASR` performs the following operations:

1.  validates the phylogeny, observed terminal states, shift edges, and
    ancestral-state reconstruction backend;
2.  assigns stable node identifiers before pruning or subtree
    extraction;
3.  identifies the descendant patch associated with each shift edge;
4.  separates the phylogeny into a mother component and one or more
    patch components;
5.  performs ancestral-state reconstruction independently on eligible
    components;
6.  reconstructs the designated boundary nodes using branch lengths from
    the original phylogeny; and
7.  merges component and boundary estimates back onto the complete set
    of internal nodes of the original tree.

The primary user interface is:

``` r
patch_asr(
  tree,
  states,
  shifts,
  asr_fun
)
```

## Current scope

The current implementation supports:

- rooted phylogenetic trees;
- fully bifurcating topologies;
- finite and strictly positive branch lengths;
- continuous traits;
- one finite observed trait value for every terminal taxon;
- internal or terminal shift children; and
- one or more compatible shift edges.

The shift parent must be an internal node and cannot be the root. Each
shift parent may define at most one patch region.

For multiple shifts, descendant patch-tip sets must be pairwise
disjoint. Nested or overlapping patch regions are not currently
supported and are rejected explicitly rather than resolved
heuristically.

## Ancestral-state reconstruction backend

`patchASR` is not tied to one particular ancestral-state reconstruction
method.

Users provide an ancestral-state reconstruction function with the
interface:

``` r
asr_fun(tree, states)
```

For each phylogenetic component:

- `tree` is the component phylogeny;
- `states` is a named numeric vector ordered exactly as
  `tree$tip.label`; and
- the function must return one finite ancestral-state estimate for every
  internal node.

The returned vector must be named using the stable internal-node
identifiers supplied in `tree$node.label`.

This interface allows different ancestral-state reconstruction methods
to be connected to the same patching framework through method-specific
adapters.

### Built-in ASR backends

`patchASR` currently provides five built-in ancestral-state
reconstruction interfaces:

- `asr_ape_pic()` performs PIC-style reconstruction using
  `ape::ace(method = "pic")`;
- `asr_ape_ml_bm()` performs maximum-likelihood reconstruction under a
  Brownian-motion model using `ape::ace(method = "ML", model = "BM")`;
- `asr_ape_gls_bm()` performs generalized least-squares reconstruction
  under a Brownian covariance model;
- `make_asr_rphylopars_bm()` creates a Brownian-motion reconstruction
  backend based on `Rphylopars::phylopars()`; and
- `make_asr_phytools_bayes()` creates a Bayesian reconstruction backend
  based on `phytools::anc.Bayes()`.

The first three backends use the required dependency `ape`. `Rphylopars`
and `phytools` are optional dependencies and are required only when
their corresponding backends are used.

All built-in backends enforce explicit internal-node identity mapping.
They do not silently substitute a different reconstruction method if the
requested method fails or if reconstructed node identities cannot be
mapped unambiguously.

The Bayesian backend uses deterministic component-specific random seeds.
With the same input and base seed, repeated analyses are reproducible,
and compatible multiple-shift analyses are invariant to the row order of
the supplied shift table. The current `phytools::anc.Bayes()` adapter
does not support components containing only one internal node.

### Using the built-in backends

For example, Brownian-motion maximum-likelihood reconstruction can be
used directly as the component-level ASR backend:

``` r
fit_ml <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_ml_bm
)
```

PIC-style and GLS-BM reconstruction can be supplied in the same way:

``` r
fit_pic <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_pic
)

fit_gls <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_gls_bm
)
```

Backends requiring additional configuration are created first and then
passed to `patch_asr()`. For example:

``` r
rph_backend <- make_asr_rphylopars_bm(
  REML = TRUE
)

fit_rph <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = rph_backend
)
```

A Bayesian backend can be constructed similarly:

``` r
bayes_backend <- make_asr_phytools_bayes(
  ngen = 10000L,
  sample_freq = 100L,
  burnin_frac = 0.20,
  seed = 123L
)

fit_bayes <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = bayes_backend
)
```

The MCMC settings above illustrate the software interface only. Users
should choose MCMC settings and assess convergence as appropriate for
their own scientific analysis.

## Shift-edge specification

Shift edges are supplied as a data frame containing two columns:

``` r
data.frame(
  parent = ...,
  child = ...
)
```

Both values refer to the original `ape` node numbers in the exact input
tree supplied to `patch_asr()`.

For example:

``` r
shifts <- data.frame(
  parent = 9,
  child = 10
)
```

means that the directed edge:

``` text
9 -> 10
```

is treated as the designated branch-localized shift.

Edge-row numbers from `tree$edge` should not be used as node
identifiers.

## Output

`patch_asr()` returns an object of class:

``` text
patch_asr
```

The main result table is:

``` r
fit$ancestral_states
```

and contains:

- stable node identifier;
- original `ape` node number;
- original node label, when available;
- reconstructed ancestral-state estimate;
- computational source of the estimate; and
- component provenance.

Additional information is available in:

``` r
fit$patch_plan
fit$boundary_states
fit$component_states
fit$node_map
fit$diagnostics
```

## Methodological boundaries

`patchASR` deliberately separates the patching algorithm from the
ancestral-state reconstruction method itself.

The core package therefore handles:

- shift specification;
- stable node identity;
- phylogenetic partitioning;
- component management;
- boundary reconstruction; and
- result merging.

Specific ancestral-state reconstruction methods are connected through
the standardized `asr_fun(tree, states)` interface.

This separation makes it possible to evaluate the patching framework
with different reconstruction methods without changing the core
partitioning algorithm.

The package treats reconstruction backends as explicit computational
components rather than interchangeable fallbacks. If a requested backend
cannot return a complete, finite, and unambiguously identified set of
internal-node estimates, the analysis stops with an error instead of
silently switching to another reconstruction method.

StableTraits or StableTrait-based reconstruction is not part of the
patching algorithm implemented by `patchASR`; such methods may instead
be used externally as comparison or benchmark methods.

## Development status

`patchASR` is currently distributed as a development version
(`0.0.0.9000`).

The package includes automated tests covering input validation, stable
node identity, tree partitioning, component-wise reconstruction,
boundary reconstruction, result merging, built-in ASR adapters, Bayesian
reproducibility, and multiple-shift row-order invariance.

Numerical regression tests against the original research implementation
and direct backend calculations have also been performed. Detailed
algorithmic specifications and validation information are provided in
the package vignette.

See the [patching algorithm
specification](vignettes/patching-algorithm-specification.Rmd) for
further details.

#' phyloPatch: Patch-Aware Ancestral State Reconstruction
#'
#' @description
#' `phyloPatch` provides infrastructure for patch-aware ancestral-state
#' reconstruction of continuous traits on phylogenetic trees containing
#' one or more specified branch-localized shifts.
#'
#' The package separates compatible shift-associated descendant regions
#' from the remaining phylogeny, performs ancestral-state reconstruction
#' on the resulting mother and patch components, reconstructs designated
#' boundary nodes using branch-length geometry from the original tree, and
#' merges all estimates back onto the internal nodes of the original
#' phylogeny.
#'
#' @section Core workflow:
#'
#' The main user-facing function is [patch_asr()].
#'
#' The internal workflow consists of:
#'
#' 1. validating the phylogeny, observed terminal states, shift edges, and
#'    ancestral-state reconstruction backend;
#' 2. assigning stable node identifiers before any pruning or subtree
#'    extraction;
#' 3. constructing a patch plan for each specified shift edge;
#' 4. partitioning the original phylogeny into a mother component,
#'    internal patch components, terminal patches, and boundary nodes;
#' 5. applying the user-supplied ancestral-state reconstruction function
#'    independently to eligible phylogenetic components;
#' 6. reconstructing boundary ancestral states from the original
#'    branch-length geometry; and
#' 7. merging component and boundary estimates back onto the complete set
#'    of internal nodes of the original tree.
#'
#' @section Scope and assumptions:
#'
#' The current implementation is designed for rooted, fully bifurcating
#' phylogenies with finite, strictly positive branch lengths and complete finite continuous
#' trait values for all terminal taxa.
#' Shift edges are specified by the original `ape` node numbers of their
#' parent and child nodes. The shift parent must be an internal non-root
#' node. Shift children may be internal or terminal.
#'
#' Multiple shifts are supported when their patch descendant sets are
#' pairwise disjoint. Nested or overlapping patch regions are not supported
#' in the current implementation.
#'
#' @section Ancestral-state reconstruction backend:
#'
#' `phyloPatch` is independent of any particular ancestral-state
#' reconstruction method. Users provide an `asr_fun` function with the
#' standardized interface `asr_fun(tree, states)`.
#'
#' The `states` argument supplied to the backend is a named numeric vector
#' ordered exactly as `tree$tip.label`. The backend must return one finite
#' estimate for every internal node of the supplied component as a named
#' numeric vector whose names match the stable internal-node identifiers
#' carried by that component tree.
#'
#' Terminal patch components contain no internal nodes and therefore bypass
#' the ancestral-state reconstruction backend.
#'
#' @section Boundary reconstruction:
#'
#' Shift-parent nodes are treated as designated boundary nodes rather than
#' ordinary component-level ancestral-state reconstruction nodes.
#'
#' Boundary states are reconstructed from neighboring known states using
#' branch lengths from the original unmodified phylogeny. When multiple
#' unresolved boundary nodes are adjacent, their values are solved jointly
#' as a linear system.
#'
#' @section Node identity and provenance:
#'
#' Original `ape` node numbers are used only as references to the exact
#' input tree. Before tree modification, `phyloPatch` assigns package-level
#' stable node identifiers so that biological node identity can be retained
#' across pruning, subtree extraction, and local node renumbering.
#'
#' Final ancestral-state results retain both the computational source of an
#' estimate and the component from which it was obtained.
#'
#' @section Main function:
#'
#' See [patch_asr()] for the primary user interface.
#'
#' @keywords internal
"_PACKAGE"

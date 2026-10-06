# Patch-plan construction infrastructure for phyloPatch.
#
# A patch plan freezes all topology and branch-length information required
# for later tree partitioning and boundary reconstruction.


#' Return the parent of one node
#'
#' @param tree A rooted `phylo` object.
#' @param node An original-tree ape node number.
#'
#' @return The parent node number, or `NA_integer_` for the root.
#'
#' @noRd
tree_parent <- function(
    tree,
    node
) {

  rows <- which(
    tree$edge[, 2L] ==
      as.integer(node)
  )

  if (length(rows) == 0L) {
    return(NA_integer_)
  }

  if (length(rows) != 1L) {
    stop(
      "A non-root node must have exactly one parent.",
      call. = FALSE
    )
  }

  as.integer(
    tree$edge[rows, 1L]
  )
}


#' Return the original branch length between two adjacent nodes
#'
#' @param tree A `phylo` object with branch lengths.
#' @param parent Parent node number.
#' @param child Child node number.
#'
#' @return A positive finite branch length.
#'
#' @noRd
tree_edge_length <- function(
    tree,
    parent,
    child
) {

  rows <- which(
    tree$edge[, 1L] ==
      as.integer(parent) &
      tree$edge[, 2L] ==
      as.integer(child)
  )

  if (length(rows) != 1L) {
    stop(
      "Requested parent-child edge was not found uniquely in the original tree.",
      call. = FALSE
    )
  }

  value <- tree$edge.length[[rows]]

  if (!is.finite(value) ||
      value <= 0) {
    stop(
      "Patch-plan branch lengths must be finite and strictly positive.",
      call. = FALSE
    )
  }

  as.numeric(value)
}


#' Build one patch-plan entry
#'
#' @param tree Validated original input tree.
#' @param node_map Stable node map for the original tree.
#' @param parent Original-tree node number of the shift parent.
#' @param child Original-tree node number of the shift child.
#' @param patch_index Internal bookkeeping index.
#'
#' @return A list describing one patch.
#'
#' @noRd
build_patch_entry <- function(
    tree,
    node_map,
    parent,
    child,
    patch_index
) {

  parent <- as.integer(parent)
  child <- as.integer(child)

  n_tip <- ape::Ntip(tree)

  ## ----------------------------------------------------------
  ## Identify H = parent(I)
  ## ----------------------------------------------------------

  upstream <- tree_parent(
    tree,
    parent
  )

  if (is.na(upstream)) {
    stop(
      "A shift parent cannot be the root in the current implementation.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Identify J' = the other child of I
  ## ----------------------------------------------------------

  children <- tree_children(
    tree,
    parent
  )

  if (length(children) != 2L) {
    stop(
      "A shift parent must have exactly two children.",
      call. = FALSE
    )
  }

  if (!(child %in% children)) {
    stop(
      "The specified shift child is not a direct child of its parent.",
      call. = FALSE
    )
  }

  sibling <- setdiff(
    children,
    child
  )

  if (length(sibling) != 1L) {
    stop(
      "Could not identify exactly one sibling-side node for the shift.",
      call. = FALSE
    )
  }

  sibling <- as.integer(sibling)

  ## ----------------------------------------------------------
  ## Determine the complete descendant patch rooted at J
  ## ----------------------------------------------------------

  descendant_nodes <-
    descendant_nodes_including_self(
      tree,
      child
    )

  descendant_tip_nodes <-
    descendant_nodes[
      descendant_nodes <= n_tip
    ]

  descendant_internal_nodes <-
    descendant_nodes[
      descendant_nodes > n_tip
    ]

  if (length(descendant_tip_nodes) < 1L) {
    stop(
      "Every patch must contain at least one descendant terminal taxon.",
      call. = FALSE
    )
  }

  child_type <- if (
    child <= n_tip
  ) {
    "terminal"
  } else {
    "internal"
  }

  requires_asr <- identical(
    child_type,
    "internal"
  )

  ## ----------------------------------------------------------
  ## Enforce internal-versus-terminal patch invariants
  ## ----------------------------------------------------------

  if (identical(
    child_type,
    "terminal"
  )) {

    if (length(descendant_tip_nodes) != 1L ||
        descendant_tip_nodes[[1L]] != child) {
      stop(
        "A terminal-child patch must contain exactly its focal terminal taxon.",
        call. = FALSE
      )
    }

    if (length(descendant_internal_nodes) != 0L) {
      stop(
        "A terminal-child patch must contain no internal nodes.",
        call. = FALSE
      )
    }

  } else {

    if (!(child %in%
          descendant_internal_nodes)) {
      stop(
        "An internal-child patch must contain its child root as an internal node.",
        call. = FALSE
      )
    }

    if (length(descendant_internal_nodes) < 1L) {
      stop(
        "An internal-child patch must contain at least one internal node.",
        call. = FALSE
      )
    }
  }

  ## ----------------------------------------------------------
  ## Convert original node numbers to stable IDs
  ## ----------------------------------------------------------

  parent_id <-
    stable_id_from_original_node(
      node_map,
      parent
    )

  child_id <-
    stable_id_from_original_node(
      node_map,
      child
    )

  upstream_id <-
    stable_id_from_original_node(
      node_map,
      upstream
    )

  sibling_id <-
    stable_id_from_original_node(
      node_map,
      sibling
    )

  descendant_tip_ids <-
    stable_id_from_original_node(
      node_map,
      descendant_tip_nodes
    )

  descendant_internal_ids <-
    if (length(
      descendant_internal_nodes
    ) == 0L) {
      character(0L)
    } else {
      stable_id_from_original_node(
        node_map,
        descendant_internal_nodes
      )
    }

  ## ----------------------------------------------------------
  ## Record original branch lengths required later
  ## ----------------------------------------------------------

  d_parent_boundary <-
    tree_edge_length(
      tree,
      upstream,
      parent
    )

  d_boundary_sibling <-
    tree_edge_length(
      tree,
      parent,
      sibling
    )

  d_shift_edge <-
    tree_edge_length(
      tree,
      parent,
      child
    )

  ## ----------------------------------------------------------
  ## Build patch entry
  ## ----------------------------------------------------------

  list(
    patch_index = as.integer(
      patch_index
    ),
    patch_id = paste0(
      "patch_",
      patch_index
    ),

    parent_node = parent,
    parent_id = parent_id,

    child_node = child,
    child_id = child_id,

    upstream_node = upstream,
    upstream_id = upstream_id,

    sibling_node = sibling,
    sibling_id = sibling_id,

    child_type = child_type,

    descendant_tip_nodes =
      as.integer(
        descendant_tip_nodes
      ),

    descendant_tip_ids =
      descendant_tip_ids,

    descendant_tip_labels =
      tree$tip.label[
        descendant_tip_nodes
      ],

    descendant_internal_nodes =
      as.integer(
        descendant_internal_nodes
      ),

    descendant_internal_ids =
      descendant_internal_ids,

    requires_asr =
      requires_asr,

    branch_length_parent_boundary =
      d_parent_boundary,

    branch_length_boundary_sibling =
      d_boundary_sibling,

    branch_length_shift_edge =
      d_shift_edge
  )
}


#' Build a complete patch plan
#'
#' @param tree Original input phylogeny.
#' @param node_map Stable node map produced by `freeze_node_identity()`.
#' @param shifts Validated shift-edge specification.
#'
#' @return An internal patch-plan object.
#'
#' @noRd
build_patch_plan <- function(
    tree,
    node_map,
    shifts
) {

  ## Revalidate shifts defensively.
  shifts <- validate_shifts(
    tree,
    shifts
  )

  required_map_columns <- c(
    "original_node",
    "node_type",
    "original_label",
    "stable_id"
  )

  if (!is.data.frame(node_map) ||
      !all(
        required_map_columns %in%
        names(node_map)
      )) {
    stop(
      "`node_map` has an invalid structure.",
      call. = FALSE
    )
  }

  expected_original_nodes <-
    seq_len(
      ape::Ntip(tree) +
        tree$Nnode
    )

  if (!identical(
    as.integer(node_map$original_node),
    as.integer(expected_original_nodes)
  )) {
    stop(
      "`node_map` does not correspond to the supplied original tree.",
      call. = FALSE
    )
  }

  if (anyNA(node_map$stable_id) ||
      any(node_map$stable_id == "") ||
      anyDuplicated(node_map$stable_id)) {
    stop(
      "`node_map` contains invalid stable node identifiers.",
      call. = FALSE
    )
  }

  patches <- lapply(
    seq_len(nrow(shifts)),
    function(i) {

      build_patch_entry(
        tree = tree,
        node_map = node_map,
        parent =
          shifts$parent[[i]],
        child =
          shifts$child[[i]],
        patch_index = i
      )
    }
  )

  boundary_ids <- vapply(
    patches,
    function(x) {
      x$parent_id
    },
    character(1)
  )

  patch_tip_ids <- unique(
    unlist(
      lapply(
        patches,
        function(x) {
          x$descendant_tip_ids
        }
      ),
      use.names = FALSE
    )
  )

  patch_tip_labels <- unique(
    unlist(
      lapply(
        patches,
        function(x) {
          x$descendant_tip_labels
        }
      ),
      use.names = FALSE
    )
  )

  patch_internal_ids <- unique(
    unlist(
      lapply(
        patches,
        function(x) {
          x$descendant_internal_ids
        }
      ),
      use.names = FALSE
    )
  )

  ## Pairwise-disjoint patches imply that no node identity should
  ## occur in more than one patch. validate_shifts() already checks
  ## this topologically; this is an additional post-construction check.

  all_patch_node_ids <- unlist(
    lapply(
      patches,
      function(x) {
        c(
          x$descendant_tip_ids,
          x$descendant_internal_ids
        )
      }
    ),
    use.names = FALSE
  )

  if (anyDuplicated(
    all_patch_node_ids
  )) {
    stop(
      "Patch-plan construction produced overlapping patch node sets.",
      call. = FALSE
    )
  }

  if (anyDuplicated(
    boundary_ids
  )) {
    stop(
      "Patch-plan construction produced duplicated boundary nodes.",
      call. = FALSE
    )
  }

  structure(
    list(
      patches = patches,
      n_patches = length(patches),
      boundary_ids = boundary_ids,
      patch_tip_ids = patch_tip_ids,
      patch_tip_labels = patch_tip_labels,
      patch_internal_ids = patch_internal_ids
    ),
    class = "phyloPatch_patch_plan"
  )
}

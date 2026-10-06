# Tree-partition infrastructure for phyloPatch.
#
# The original phylogeny is partitioned into one mother component,
# zero or more internal patch trees, terminal patch representations,
# and the designated boundary-node set.


#' Create an identity-carrying working copy of the original tree
#'
#' @param tree Original input `phylo` object.
#' @param node_map Stable original-tree node map.
#'
#' @return A working copy whose internal node labels contain stable IDs.
#'
#' @noRd
make_identity_working_tree <- function(
    tree,
    node_map
) {

  if (!inherits(tree, "phylo")) {
    stop(
      "`tree` must inherit from class `phylo`.",
      call. = FALSE
    )
  }

  required_columns <- c(
    "original_node",
    "node_type",
    "original_label",
    "stable_id"
  )

  if (!is.data.frame(node_map) ||
      !all(required_columns %in%
           names(node_map))) {
    stop(
      "`node_map` has an invalid structure.",
      call. = FALSE
    )
  }

  n_tip <- ape::Ntip(tree)
  n_total <- n_tip + tree$Nnode

  expected_nodes <- seq_len(
    n_total
  )

  if (!identical(
    as.integer(node_map$original_node),
    as.integer(expected_nodes)
  )) {
    stop(
      "`node_map` does not correspond to the supplied original tree.",
      call. = FALSE
    )
  }

  expected_types <- c(
    rep("tip", n_tip),
    rep("internal", tree$Nnode)
  )

  if (!identical(
    as.character(node_map$node_type),
    expected_types
  )) {
    stop(
      "`node_map` contains invalid node-type assignments.",
      call. = FALSE
    )
  }

  tip_rows <- which(
    node_map$node_type == "tip"
  )

  if (!identical(
    as.character(
      node_map$original_label[
        tip_rows
      ]
    ),
    as.character(tree$tip.label)
  )) {
    stop(
      "`node_map` terminal labels do not correspond to the supplied tree.",
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

  internal_rows <- which(
    node_map$node_type == "internal"
  )

  if (length(internal_rows) !=
      tree$Nnode) {
    stop(
      "`node_map` has an incorrect number of internal nodes.",
      call. = FALSE
    )
  }

  working_tree <- tree

  working_tree$node.label <-
    node_map$stable_id[
      internal_rows
    ]

  working_tree
}


#' Construct one partitioned patch component
#'
#' @param identity_tree Original-tree working copy carrying stable IDs.
#' @param node_map Original stable node map.
#' @param patch One patch-plan entry.
#'
#' @return A structured patch component.
#'
#' @noRd
build_partition_patch_component <- function(
    identity_tree,
    node_map,
    patch
) {

  ## ----------------------------------------------------------
  ## Terminal patch
  ## ----------------------------------------------------------

  if (!isTRUE(
    patch$requires_asr
  )) {

    if (!identical(
      patch$child_type,
      "terminal"
    )) {
      stop(
        "A non-ASR patch must be a terminal-child patch.",
        call. = FALSE
      )
    }

    if (length(
      patch$descendant_tip_nodes
    ) != 1L ||
    length(
      patch$descendant_internal_nodes
    ) != 0L) {
      stop(
        "Terminal patch content is inconsistent with the patch plan.",
        call. = FALSE
      )
    }

    component_map <- data.frame(
      local_node = 1L,
      node_type = "tip",
      stable_id =
        patch$child_id,
      original_node =
        patch$child_node,
      original_label =
        identity_tree$tip.label[
          patch$child_node
        ],
      stringsAsFactors = FALSE
    )

    return(
      list(
        patch_id =
          patch$patch_id,
        requires_asr =
          FALSE,
        child_type =
          "terminal",
        tree =
          NULL,
        component_map =
          component_map,
        tip_ids =
          patch$descendant_tip_ids,
        internal_ids =
          character(0L)
      )
    )
  }

  ## ----------------------------------------------------------
  ## Internal patch
  ## ----------------------------------------------------------

  if (!identical(
    patch$child_type,
    "internal"
  )) {
    stop(
      "An ASR-requiring patch must be an internal-child patch.",
      call. = FALSE
    )
  }

  patch_tree <- ape::extract.clade(
    identity_tree,
    node =
      patch$child_node
  )

  component_map <- map_component_nodes(
    patch_tree,
    node_map
  )

  expected_ids <- c(
    patch$descendant_tip_ids,
    patch$descendant_internal_ids
  )

  if (!setequal(
    component_map$stable_id,
    expected_ids
  )) {
    stop(
      "Extracted patch component does not match its patch-plan node set.",
      call. = FALSE
    )
  }

  actual_tip_ids <-
    component_map$stable_id[
      component_map$node_type ==
        "tip"
    ]

  actual_internal_ids <-
    component_map$stable_id[
      component_map$node_type ==
        "internal"
    ]

  if (!setequal(
    actual_tip_ids,
    patch$descendant_tip_ids
  )) {
    stop(
      "Extracted patch component has an incorrect terminal-node set.",
      call. = FALSE
    )
  }

  if (!setequal(
    actual_internal_ids,
    patch$descendant_internal_ids
  )) {
    stop(
      "Extracted patch component has an incorrect internal-node set.",
      call. = FALSE
    )
  }

  list(
    patch_id =
      patch$patch_id,
    requires_asr =
      TRUE,
    child_type =
      "internal",
    tree =
      patch_tree,
    component_map =
      component_map,
    tip_ids =
      actual_tip_ids,
    internal_ids =
      actual_internal_ids
  )
}


#' Partition a phylogeny according to a validated patch plan
#'
#' @param tree Original validated input phylogeny.
#' @param node_map Stable original-tree node map.
#' @param patch_plan Patch plan produced by `build_patch_plan()`.
#'
#' @return An internal `phyloPatch_partition` object.
#'
#' @noRd
partition_phylogeny <- function(
    tree,
    node_map,
    patch_plan
) {

  ## ----------------------------------------------------------
  ## Defensive input checks
  ## ----------------------------------------------------------

  validate_phylo_tree(
    tree
  )

  if (!inherits(
    patch_plan,
    "phyloPatch_patch_plan"
  )) {
    stop(
      "`patch_plan` must be a valid phyloPatch patch-plan object.",
      call. = FALSE
    )
  }

  if (!is.list(
    patch_plan$patches
  ) ||
  length(
    patch_plan$patches
  ) < 1L) {
    stop(
      "`patch_plan` contains no patch entries.",
      call. = FALSE
    )
  }

  if (!identical(
    length(
      patch_plan$patches
    ),
    patch_plan$n_patches
  )) {
    stop(
      "`patch_plan$n_patches` is inconsistent with its patch entries.",
      call. = FALSE
    )
  }

  identity_tree <-
    make_identity_working_tree(
      tree,
      node_map
    )

  ## ----------------------------------------------------------
  ## Build all patch components from the original working tree
  ## ----------------------------------------------------------

  patch_components <- lapply(
    patch_plan$patches,
    function(patch) {

      build_partition_patch_component(
        identity_tree =
          identity_tree,
        node_map =
          node_map,
        patch =
          patch
      )
    }
  )

  ## ----------------------------------------------------------
  ## Construct the mother component by removing all patch tips
  ## simultaneously.
  ## ----------------------------------------------------------

  patch_tip_labels <-
    patch_plan$patch_tip_labels

  if (length(
    patch_tip_labels
  ) < 1L) {
    stop(
      "The patch plan contains no terminal taxa to remove.",
      call. = FALSE
    )
  }

  if (!all(
    patch_tip_labels %in%
    identity_tree$tip.label
  )) {
    stop(
      "The patch plan contains terminal labels absent from the original tree.",
      call. = FALSE
    )
  }

  mother_tree <- ape::drop.tip(
    identity_tree,
    tip =
      patch_tip_labels,
    collapse.singles =
      TRUE
  )

  if (!inherits(
    mother_tree,
    "phylo"
  )) {
    stop(
      "Mother-tree construction did not return a phylogenetic tree.",
      call. = FALSE
    )
  }

  mother_map <- map_component_nodes(
    mother_tree,
    node_map
  )

  ## ----------------------------------------------------------
  ## Collect actual partition node sets
  ## ----------------------------------------------------------

  original_internal_ids <-
    node_map$stable_id[
      node_map$node_type ==
        "internal"
    ]

  original_tip_ids <-
    node_map$stable_id[
      node_map$node_type ==
        "tip"
    ]

  mother_internal_ids <-
    mother_map$stable_id[
      mother_map$node_type ==
        "internal"
    ]

  mother_tip_ids <-
    mother_map$stable_id[
      mother_map$node_type ==
        "tip"
    ]

  patch_internal_ids <-
    as.character(
      unlist(
        lapply(
          patch_components,
          function(x) {
            x$internal_ids
          }
        ),
        use.names = FALSE
      )
    )

  patch_tip_ids <-
    as.character(
      unlist(
        lapply(
          patch_components,
          function(x) {
            x$tip_ids
          }
        ),
        use.names = FALSE
      )
    )

  boundary_ids <-
    as.character(
      patch_plan$boundary_ids
    )

  ## ----------------------------------------------------------
  ## Check component content against the patch plan
  ## ----------------------------------------------------------

  if (!setequal(
    patch_internal_ids,
    patch_plan$patch_internal_ids
  )) {
    stop(
      "Actual patch internal-node set does not match the patch plan.",
      call. = FALSE
    )
  }

  if (!setequal(
    patch_tip_ids,
    patch_plan$patch_tip_ids
  )) {
    stop(
      "Actual patch terminal-node set does not match the patch plan.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Mother-side pruning-collapse invariant
  ##
  ## Remove:
  ##   - internal nodes retained in mother
  ##   - internal nodes assigned to patches
  ##
  ## Whatever remains must be exactly the designated boundaries.
  ## ----------------------------------------------------------

  collapsed_mother_side_ids <-
    setdiff(
      setdiff(
        original_internal_ids,
        mother_internal_ids
      ),
      patch_internal_ids
    )

  if (!setequal(
    collapsed_mother_side_ids,
    boundary_ids
  )) {
    stop(
      "Mother-side internal nodes removed by pruning-induced collapse ",
      "do not match the designated boundary-node set.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Internal-node partition must be disjoint and complete
  ## ----------------------------------------------------------

  assigned_internal_ids <- c(
    mother_internal_ids,
    patch_internal_ids,
    boundary_ids
  )

  if (anyDuplicated(
    assigned_internal_ids
  )) {
    stop(
      "An original internal node was assigned to more than one ",
      "partition category.",
      call. = FALSE
    )
  }

  if (!setequal(
    assigned_internal_ids,
    original_internal_ids
  )) {
    stop(
      "The internal-node partition is incomplete or contains ",
      "unknown node identities.",
      call. = FALSE
    )
  }

  if (length(
    assigned_internal_ids
  ) != length(
    original_internal_ids
  )) {
    stop(
      "The internal-node partition does not contain exactly one ",
      "assignment per original internal node.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Tip partition must also be disjoint and complete
  ## ----------------------------------------------------------

  expected_mother_tip_ids <-
    setdiff(
      original_tip_ids,
      patch_plan$patch_tip_ids
    )

  if (!setequal(
    mother_tip_ids,
    expected_mother_tip_ids
  )) {
    stop(
      "Mother component contains an incorrect terminal-node set.",
      call. = FALSE
    )
  }

  assigned_tip_ids <- c(
    mother_tip_ids,
    patch_tip_ids
  )

  if (anyDuplicated(
    assigned_tip_ids
  )) {
    stop(
      "A terminal node was assigned to both mother and patch components.",
      call. = FALSE
    )
  }

  if (!setequal(
    assigned_tip_ids,
    original_tip_ids
  )) {
    stop(
      "The terminal-node partition is incomplete.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Return structured partition object
  ## ----------------------------------------------------------

  structure(
    list(
      mother_tree =
        mother_tree,

      mother_node_map =
        mother_map,

      patch_components =
        patch_components,

      boundary_ids =
        boundary_ids,

      mother_internal_ids =
        mother_internal_ids,

      patch_internal_ids =
        patch_internal_ids,

      collapsed_mother_side_ids =
        collapsed_mother_side_ids,

      mother_tip_ids =
        mother_tip_ids,

      patch_tip_ids =
        patch_tip_ids
    ),
    class =
      "phyloPatch_partition"
  )
}

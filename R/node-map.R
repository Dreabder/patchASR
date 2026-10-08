# Stable node-identity infrastructure for patchASR.
#
# All biological node identities are defined relative to the original
# input tree. Derived trees may use different local ape node numbers.


#' Freeze stable node identities on a working copy
#'
#' Construct a stable mapping between nodes in the original input tree and
#' package-internal node identifiers. Stable internal-node identifiers are
#' attached only to a working copy of the tree; the user-supplied tree is
#' not modified.
#'
#' @param tree A phylogenetic tree inheriting from class `phylo`.
#'
#' @return A list with elements `tree` and `node_map`.
#'
#' @noRd
freeze_node_identity <- function(tree) {

  if (!inherits(tree, "phylo")) {
    stop(
      "`tree` must inherit from class `phylo`.",
      call. = FALSE
    )
  }

  if (is.null(tree$tip.label) ||
      length(tree$tip.label) < 1L) {
    stop(
      "`tree` must contain terminal labels.",
      call. = FALSE
    )
  }

  if (anyNA(tree$tip.label) ||
      any(tree$tip.label == "") ||
      anyDuplicated(tree$tip.label)) {
    stop(
      "`tree$tip.label` must contain unique, non-empty labels.",
      call. = FALSE
    )
  }

  n_tip <- length(tree$tip.label)
  n_internal <- tree$Nnode

  if (length(n_internal) != 1L ||
      is.na(n_internal) ||
      n_internal < 1L) {
    stop(
      "`tree` must contain at least one internal node.",
      call. = FALSE
    )
  }

  n_total <- n_tip + n_internal

  original_node <- seq_len(n_total)

  node_type <- ifelse(
    original_node <= n_tip,
    "tip",
    "internal"
  )

  original_label <- rep(
    NA_character_,
    n_total
  )

  original_label[
    seq_len(n_tip)
  ] <- tree$tip.label

  if (!is.null(tree$node.label)) {

    if (length(tree$node.label) != n_internal) {
      stop(
        "`tree$node.label` has an invalid length.",
        call. = FALSE
      )
    }

    original_label[
      (n_tip + 1L):n_total
    ] <- tree$node.label
  }

  stable_id <- ifelse(
    node_type == "tip",
    paste0("tip_", original_node),
    paste0("node_", original_node)
  )

  if (anyDuplicated(stable_id)) {
    stop(
      "Internal error: stable node identifiers are not unique.",
      call. = FALSE
    )
  }

  node_map <- data.frame(
    original_node = original_node,
    node_type = node_type,
    original_label = original_label,
    stable_id = stable_id,
    stringsAsFactors = FALSE
  )

  working_tree <- tree

  working_tree$node.label <-
    node_map$stable_id[
      node_map$node_type == "internal"
    ]

  list(
    tree = working_tree,
    node_map = node_map
  )
}


#' Map derived component nodes back to the original tree
#'
#' @param component_tree A derived phylogenetic component whose internal
#'   node labels contain stable `patchASR` identifiers.
#' @param node_map The original-tree node map produced by
#'   `freeze_node_identity()`.
#'
#' @return A data frame mapping local component nodes to stable and original
#'   node identities.
#'
#' @noRd
map_component_nodes <- function(
    component_tree,
    node_map
) {

  if (!inherits(component_tree, "phylo")) {
    stop(
      "`component_tree` must inherit from class `phylo`.",
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
      !all(required_columns %in% names(node_map))) {
    stop(
      "`node_map` has an invalid structure.",
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

  n_tip <- length(component_tree$tip.label)
  n_internal <- component_tree$Nnode

  ## ----------------------------------------------------------
  ## Map terminal nodes
  ## ----------------------------------------------------------

  original_tip_rows <- which(
    node_map$node_type == "tip"
  )

  tip_lookup <- setNames(
    original_tip_rows,
    node_map$original_label[
      original_tip_rows
    ]
  )

  tip_rows <- unname(
    tip_lookup[
      component_tree$tip.label
    ]
  )

  if (anyNA(tip_rows)) {

    missing_tips <- component_tree$tip.label[
      is.na(tip_rows)
    ]

    stop(
      "At least one component tip cannot be mapped to the original tree: ",
      paste(missing_tips, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  tip_map <- data.frame(
    local_node = seq_len(n_tip),
    node_type = "tip",
    stable_id = node_map$stable_id[
      tip_rows
    ],
    original_node = node_map$original_node[
      tip_rows
    ],
    original_label = node_map$original_label[
      tip_rows
    ],
    stringsAsFactors = FALSE
  )

  ## ----------------------------------------------------------
  ## Map internal nodes
  ## ----------------------------------------------------------

  if (n_internal > 0L) {

    if (is.null(component_tree$node.label)) {
      stop(
        "Derived component has lost its stable internal-node identifiers.",
        call. = FALSE
      )
    }

    if (length(component_tree$node.label) != n_internal) {
      stop(
        "Derived component has an invalid internal-node label count.",
        call. = FALSE
      )
    }

    internal_rows <- match(
      component_tree$node.label,
      node_map$stable_id
    )

    if (anyNA(internal_rows)) {

      unknown_ids <- component_tree$node.label[
        is.na(internal_rows)
      ]

      stop(
        "Derived component contains unknown stable node identifiers: ",
        paste(unknown_ids, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (any(
      node_map$node_type[
        internal_rows
      ] != "internal"
    )) {
      stop(
        "A component internal node mapped to a non-internal original node.",
        call. = FALSE
      )
    }

    internal_map <- data.frame(
      local_node =
        n_tip + seq_len(n_internal),
      node_type = "internal",
      stable_id =
        component_tree$node.label,
      original_node =
        node_map$original_node[
          internal_rows
        ],
      original_label =
        node_map$original_label[
          internal_rows
        ],
      stringsAsFactors = FALSE
    )

  } else {

    internal_map <- data.frame(
      local_node = integer(0L),
      node_type = character(0L),
      stable_id = character(0L),
      original_node = integer(0L),
      original_label = character(0L),
      stringsAsFactors = FALSE
    )
  }

  result <- rbind(
    tip_map,
    internal_map
  )

  if (anyDuplicated(result$stable_id)) {
    stop(
      "Derived component contains duplicated stable node identities.",
      call. = FALSE
    )
  }

  result
}


#' Look up stable IDs from original-tree node numbers
#'
#' @param node_map Original-tree node map.
#' @param node One or more original-tree `ape` node numbers.
#'
#' @return A character vector of stable node identifiers.
#'
#' @noRd
stable_id_from_original_node <- function(
    node_map,
    node
) {

  if (!is.numeric(node) ||
      anyNA(node) ||
      any(!is.finite(node)) ||
      any(node != floor(node))) {
    stop(
      "`node` must contain finite integer original-tree node numbers.",
      call. = FALSE
    )
  }

  rows <- match(
    as.integer(node),
    node_map$original_node
  )

  if (anyNA(rows)) {

    unknown <- node[
      is.na(rows)
    ]

    stop(
      "Unknown original-tree node number(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  node_map$stable_id[
    rows
  ]
}


#' Look up original-tree node numbers from stable IDs
#'
#' @param node_map Original-tree node map.
#' @param stable_id One or more stable node identifiers.
#'
#' @return An integer vector of original-tree `ape` node numbers.
#'
#' @noRd
original_node_from_stable_id <- function(
    node_map,
    stable_id
) {

  if (!is.character(stable_id) ||
      anyNA(stable_id) ||
      any(stable_id == "")) {
    stop(
      "`stable_id` must contain non-empty character identifiers.",
      call. = FALSE
    )
  }

  rows <- match(
    stable_id,
    node_map$stable_id
  )

  if (anyNA(rows)) {

    unknown <- stable_id[
      is.na(rows)
    ]

    stop(
      "Unknown stable node identifier(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  as.integer(
    node_map$original_node[
      rows
    ]
  )
}

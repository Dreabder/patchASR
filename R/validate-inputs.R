# Input validation infrastructure for patchASR.


#' Identify the root node of a rooted phylogeny
#'
#' @param tree A validated rooted `phylo` object.
#'
#' @return The original-tree ape node number of the root.
#'
#' @noRd
tree_root_node <- function(tree) {

  root_candidates <- setdiff(
    unique(tree$edge[, 1L]),
    tree$edge[, 2L]
  )

  if (length(root_candidates) != 1L) {
    stop(
      "Could not identify exactly one root node.",
      call. = FALSE
    )
  }

  as.integer(root_candidates)
}


#' Return direct children of one node
#'
#' @param tree A `phylo` object.
#' @param node An ape node number.
#'
#' @return An integer vector of direct child-node numbers.
#'
#' @noRd
tree_children <- function(tree, node) {

  as.integer(
    tree$edge[
      tree$edge[, 1L] == as.integer(node),
      2L
    ]
  )
}


#' Collect descendants including the focal node
#'
#' @param tree A rooted `phylo` object.
#' @param node An ape node number.
#'
#' @return An integer vector containing `node` and all descendants.
#'
#' @noRd
descendant_nodes_including_self <- function(
    tree,
    node
) {

  node <- as.integer(node)
  n_tip <- ape::Ntip(tree)

  stack <- node
  result <- integer(0L)

  while (length(stack) > 0L) {

    current <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    result <- c(
      result,
      current
    )

    if (current > n_tip) {
      stack <- c(
        stack,
        tree_children(
          tree,
          current
        )
      )
    }
  }

  sort(
    unique(result)
  )
}


#' Validate the input phylogeny
#'
#' @param tree Candidate phylogenetic tree.
#'
#' @return The validated tree, invisibly.
#'
#' @noRd
validate_phylo_tree <- function(tree) {

  if (!inherits(tree, "phylo")) {
    stop(
      "`tree` must inherit from class `phylo`.",
      call. = FALSE
    )
  }

  if (!ape::is.rooted(tree)) {
    stop(
      "`tree` must be rooted.",
      call. = FALSE
    )
  }

  if (!ape::is.binary(tree)) {
    stop(
      "`tree` must be fully bifurcating.",
      call. = FALSE
    )
  }

  if (is.null(tree$tip.label) ||
      length(tree$tip.label) < 2L) {
    stop(
      "`tree` must contain at least two terminal taxa.",
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

  if (is.null(tree$edge.length)) {
    stop(
      "`tree` must contain branch lengths.",
      call. = FALSE
    )
  }

  if (length(tree$edge.length) !=
      nrow(tree$edge)) {
    stop(
      "`tree$edge.length` has an invalid length.",
      call. = FALSE
    )
  }

  if (any(!is.finite(tree$edge.length))) {
    stop(
      "All branch lengths must be finite.",
      call. = FALSE
    )
  }

  if (any(tree$edge.length <= 0)) {
    stop(
      "All branch lengths must be strictly positive.",
      call. = FALSE
    )
  }

  ## Confirm that the root can be identified uniquely.
  tree_root_node(tree)

  invisible(tree)
}


#' Validate observed terminal states
#'
#' @param tree A validated `phylo` object.
#' @param states Candidate terminal-state vector.
#'
#' @return `states` reordered to `tree$tip.label`.
#'
#' @noRd
validate_tip_states <- function(
    tree,
    states
) {

  if (!is.numeric(states) ||
      is.null(names(states))) {
    stop(
      "`states` must be a named numeric vector.",
      call. = FALSE
    )
  }

  if (anyNA(names(states)) ||
      any(names(states) == "") ||
      anyDuplicated(names(states))) {
    stop(
      "`states` must have unique, non-empty tip names.",
      call. = FALSE
    )
  }

  if (any(!is.finite(states))) {
    stop(
      "All terminal trait values must be finite.",
      call. = FALSE
    )
  }

  missing_tips <- setdiff(
    tree$tip.label,
    names(states)
  )

  unknown_tips <- setdiff(
    names(states),
    tree$tip.label
  )

  if (length(missing_tips) > 0L) {
    stop(
      "Missing terminal trait value(s) for: ",
      paste(
        missing_tips,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (length(unknown_tips) > 0L) {
    stop(
      "Unknown terminal label(s) in `states`: ",
      paste(
        unknown_tips,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (length(states) !=
      length(tree$tip.label)) {
    stop(
      "`states` must contain exactly one value for every terminal taxon.",
      call. = FALSE
    )
  }

  states <- states[
    tree$tip.label
  ]

  if (!identical(
    names(states),
    tree$tip.label
  )) {
    stop(
      "Internal error while ordering terminal states.",
      call. = FALSE
    )
  }

  states
}


#' Validate shift-edge specification
#'
#' @param tree A validated rooted, bifurcating `phylo` object.
#' @param shifts A data frame with `parent` and `child` columns containing
#'   original-tree ape node numbers.
#'
#' @return A normalized data frame with integer `parent` and `child` columns.
#'
#' @noRd
validate_shifts <- function(
    tree,
    shifts
) {

  if (!is.data.frame(shifts)) {
    stop(
      "`shifts` must be a data frame.",
      call. = FALSE
    )
  }

  if (!all(
    c("parent", "child") %in%
    names(shifts)
  )) {
    stop(
      "`shifts` must contain `parent` and `child` columns.",
      call. = FALSE
    )
  }

  if (nrow(shifts) < 1L) {
    stop(
      "`shifts` must contain at least one requested shift edge.",
      call. = FALSE
    )
  }

  validate_node_reference <- function(
    x,
    field
  ) {

    if (!is.numeric(x) ||
        anyNA(x) ||
        any(!is.finite(x)) ||
        any(x != floor(x))) {
      stop(
        "`",
        field,
        "` must contain finite integer original-tree node numbers.",
        call. = FALSE
      )
    }

    as.integer(x)
  }

  parent <- validate_node_reference(
    shifts$parent,
    "parent"
  )

  child <- validate_node_reference(
    shifts$child,
    "child"
  )

  n_tip <- ape::Ntip(tree)
  n_total <- n_tip + tree$Nnode

  if (any(
    parent < 1L |
    parent > n_total
  )) {
    stop(
      "At least one `parent` node number is outside the input tree.",
      call. = FALSE
    )
  }

  if (any(
    child < 1L |
    child > n_total
  )) {
    stop(
      "At least one `child` node number is outside the input tree.",
      call. = FALSE
    )
  }

  if (any(parent <= n_tip)) {
    stop(
      "Every shift parent must be an internal node.",
      call. = FALSE
    )
  }

  root <- tree_root_node(tree)

  if (any(parent == root)) {
    stop(
      "The root cannot serve as a shift parent in the current implementation.",
      call. = FALSE
    )
  }

  ## Every parent-child pair must be one actual directed edge.
  for (i in seq_along(parent)) {

    edge_rows <- which(
      tree$edge[, 1L] == parent[[i]] &
        tree$edge[, 2L] == child[[i]]
    )

    if (length(edge_rows) != 1L) {
      stop(
        "Requested shift is not an existing directed edge: ",
        parent[[i]],
        " -> ",
        child[[i]],
        ".",
        call. = FALSE
      )
    }
  }

  edge_keys <- paste(
    parent,
    child,
    sep = "->"
  )

  if (anyDuplicated(edge_keys)) {
    stop(
      "Duplicated shift edges are not allowed.",
      call. = FALSE
    )
  }

  if (anyDuplicated(parent)) {
    stop(
      "Multiple requested shifts cannot share the same parent in the current implementation.",
      call. = FALSE
    )
  }

  ## Descendant patch clades must be pairwise disjoint.
  descendant_sets <- lapply(
    child,
    function(x) {
      descendant_nodes_including_self(
        tree,
        x
      )
    }
  )

  if (length(descendant_sets) > 1L) {

    for (i in seq_len(
      length(descendant_sets) - 1L
    )) {

      for (j in seq.int(
        i + 1L,
        length(descendant_sets)
      )) {

        overlap <- intersect(
          descendant_sets[[i]],
          descendant_sets[[j]]
        )

        if (length(overlap) > 0L) {
          stop(
            "Requested patch clades must be pairwise disjoint; ",
            "nested or overlapping patch clades are not supported ",
            "in the current implementation.",
            call. = FALSE
          )
        }
      }
    }
  }

  data.frame(
    parent = parent,
    child = child,
    stringsAsFactors = FALSE
  )
}


#' Validate all core patching inputs
#'
#' @param tree Input phylogeny.
#' @param states Named numeric terminal-state vector.
#' @param shifts Requested shift edges.
#' @param asr_fun Component-level ancestral-state reconstruction function.
#'
#' @return A normalized list containing validated inputs.
#'
#' @noRd
validate_patch_inputs <- function(
    tree,
    states,
    shifts,
    asr_fun
) {

  validate_phylo_tree(tree)

  validated_states <- validate_tip_states(
    tree,
    states
  )

  validated_shifts <- validate_shifts(
    tree,
    shifts
  )

  if (!is.function(asr_fun)) {
    stop(
      "`asr_fun` must be a function.",
      call. = FALSE
    )
  }

  list(
    tree = tree,
    states = validated_states,
    shifts = validated_shifts,
    asr_fun = asr_fun
  )
}

## ============================================================
## Step 2B — Partition / boundary-node audit
##
## Purpose:
##   Test whether, for currently supported patch configurations,
##   the original mother-side internal nodes removed by
##   drop.tip(..., collapse.singles = TRUE) are exactly the
##   specified shift parents (boundary nodes).
##
## This is a development audit script, not package code.
## ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for this audit.")
}

## ------------------------------------------------------------
## 1. Construct a small, fully labelled test tree
## ------------------------------------------------------------

tree_text <- paste0(
  "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
  "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
)

tree <- ape::read.tree(text = tree_text)

if (!ape::is.rooted(tree)) {
  stop("Audit tree is not rooted.")
}

if (!ape::is.binary(tree)) {
  stop("Audit tree is not fully bifurcating.")
}

if (is.null(tree$edge.length) ||
    any(!is.finite(tree$edge.length)) ||
    any(tree$edge.length <= 0)) {
  stop("Audit tree has invalid branch lengths.")
}


## ------------------------------------------------------------
## 2. Stable label <-> ape node helpers
## ------------------------------------------------------------

node_number_from_label <- function(tree, label) {

  tip_match <- match(label, tree$tip.label)

  if (!is.na(tip_match)) {
    return(as.integer(tip_match))
  }

  if (is.null(tree$node.label)) {
    stop("Internal node labels are missing.")
  }

  internal_match <- match(label, tree$node.label)

  if (is.na(internal_match)) {
    stop("Unknown node label: ", label)
  }

  as.integer(ape::Ntip(tree) + internal_match)
}


node_label_from_number <- function(tree, node) {

  node <- as.integer(node)
  n_tip <- ape::Ntip(tree)

  if (node <= n_tip) {
    return(tree$tip.label[[node]])
  }

  index <- node - n_tip

  if (index < 1L || index > tree$Nnode) {
    stop("Invalid node number: ", node)
  }

  tree$node.label[[index]]
}


children_of <- function(tree, node) {

  as.integer(
    tree$edge[
      tree$edge[, 1L] == as.integer(node),
      2L
    ]
  )
}


parent_of <- function(tree, node) {

  rows <- which(
    tree$edge[, 2L] == as.integer(node)
  )

  if (length(rows) == 0L) {
    return(NA_integer_)
  }

  if (length(rows) != 1L) {
    stop("Node has more than one parent: ", node)
  }

  as.integer(tree$edge[rows, 1L])
}


## ------------------------------------------------------------
## 3. Descendant-node helper
## ------------------------------------------------------------

descendants_including_self <- function(tree, node) {

  node <- as.integer(node)
  n_tip <- ape::Ntip(tree)

  stack <- node
  all_nodes <- integer(0L)

  while (length(stack) > 0L) {

    current <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    all_nodes <- c(all_nodes, current)

    if (current > n_tip) {
      stack <- c(
        stack,
        children_of(tree, current)
      )
    }
  }

  all_nodes <- sort(unique(all_nodes))

  list(
    all = all_nodes,
    tips = all_nodes[all_nodes <= n_tip],
    internal = all_nodes[all_nodes > n_tip]
  )
}


## ------------------------------------------------------------
## 4. Validate one shift configuration
## ------------------------------------------------------------

validate_shift_configuration <- function(tree, shifts) {

  if (!all(c("parent", "child") %in% names(shifts))) {
    stop("shifts must contain 'parent' and 'child'.")
  }

  if (nrow(shifts) < 1L) {
    stop("At least one shift is required.")
  }

  parent_nodes <- vapply(
    shifts$parent,
    function(x) node_number_from_label(tree, x),
    integer(1)
  )

  child_nodes <- vapply(
    shifts$child,
    function(x) node_number_from_label(tree, x),
    integer(1)
  )

  ## Parent must be internal.
  if (any(parent_nodes <= ape::Ntip(tree))) {
    stop("Every shift parent must be an internal node.")
  }

  ## Root cannot be a shift parent in v0.1.
  if (any(vapply(
    parent_nodes,
    function(x) is.na(parent_of(tree, x)),
    logical(1)
  ))) {
    stop("The root cannot serve as a shift parent in v0.1.")
  }

  ## Exact directed edge must exist.
  for (i in seq_len(nrow(shifts))) {

    rows <- which(
      tree$edge[, 1L] == parent_nodes[[i]] &
        tree$edge[, 2L] == child_nodes[[i]]
    )

    if (length(rows) != 1L) {
      stop(
        "Not an existing directed edge: ",
        shifts$parent[[i]], " -> ", shifts$child[[i]]
      )
    }
  }

  ## No duplicated shift edges.
  edge_keys <- paste(parent_nodes, child_nodes, sep = "->")

  if (anyDuplicated(edge_keys)) {
    stop("Duplicated shift edge.")
  }

  ## v0.1: no shared shift parent.
  if (anyDuplicated(parent_nodes)) {
    stop("Multiple shifts sharing the same parent are not supported.")
  }

  ## Patch descendant sets must be pairwise disjoint.
  descendant_sets <- lapply(
    child_nodes,
    function(x) descendants_including_self(tree, x)$all
  )

  if (length(descendant_sets) > 1L) {

    for (i in seq_len(length(descendant_sets) - 1L)) {

      for (j in (i + 1L):length(descendant_sets)) {

        if (length(intersect(
          descendant_sets[[i]],
          descendant_sets[[j]]
        )) > 0L) {

          stop(
            "Nested or overlapping patch clades detected."
          )
        }
      }
    }
  }

  invisible(TRUE)
}


## ------------------------------------------------------------
## 5. Audit one valid configuration
## ------------------------------------------------------------

audit_one_case <- function(tree, shifts, case_name) {

  validate_shift_configuration(tree, shifts)

  n_tip <- ape::Ntip(tree)

  parent_nodes <- vapply(
    shifts$parent,
    function(x) node_number_from_label(tree, x),
    integer(1)
  )

  child_nodes <- vapply(
    shifts$child,
    function(x) node_number_from_label(tree, x),
    integer(1)
  )

  descendant_info <- lapply(
    child_nodes,
    function(x) descendants_including_self(tree, x)
  )

  ## All tips assigned to patches.
  patch_tip_nodes <- sort(unique(unlist(
    lapply(descendant_info, `[[`, "tips"),
    use.names = FALSE
  )))

  patch_tip_labels <- tree$tip.label[patch_tip_nodes]

  ## All original internal nodes belonging to patch components.
  patch_internal_nodes <- sort(unique(unlist(
    lapply(descendant_info, `[[`, "internal"),
    use.names = FALSE
  )))

  patch_internal_labels <- if (length(patch_internal_nodes) == 0L) {
    character(0L)
  } else {
    vapply(
      patch_internal_nodes,
      function(x) node_label_from_number(tree, x),
      character(1)
    )
  }

  ## Construct shared mother tree.
  mother_tree <- ape::drop.tip(
    tree,
    tip = patch_tip_labels,
    collapse.singles = TRUE
  )

  original_internal_labels <- tree$node.label

  mother_internal_labels <- if (is.null(mother_tree$node.label)) {
    character(0L)
  } else {
    mother_tree$node.label
  }

  ## Internal nodes absent from the mother tree include:
  ##   1. patch internal nodes;
  ##   2. pruning-induced collapsed mother-side nodes.
  absent_from_mother <- setdiff(
    original_internal_labels,
    mother_internal_labels
  )

  collapsed_mother_side <- setdiff(
    absent_from_mother,
    patch_internal_labels
  )

  expected_boundary_labels <- vapply(
    parent_nodes,
    function(x) node_label_from_number(tree, x),
    character(1)
  )

  pass <- setequal(
    collapsed_mother_side,
    expected_boundary_labels
  )

  data.frame(
    case = case_name,
    patch_tips = paste(
      sort(patch_tip_labels),
      collapse = ","
    ),
    patch_internal = if (length(patch_internal_labels) == 0L) {
      "<none>"
    } else {
      paste(sort(patch_internal_labels), collapse = ",")
    },
    expected_boundary = paste(
      sort(expected_boundary_labels),
      collapse = ","
    ),
    collapsed_mother_side = if (length(collapsed_mother_side) == 0L) {
      "<none>"
    } else {
      paste(sort(collapsed_mother_side), collapse = ",")
    },
    pass = pass,
    stringsAsFactors = FALSE
  )
}


## ------------------------------------------------------------
## 6. Define audit cases
## ------------------------------------------------------------

cases <- list(

  ## Single internal-child shift:
  ## I1 -> J1; patch contains A/B and internal node J1.
  single_internal = data.frame(
    parent = "I1",
    child = "J1",
    stringsAsFactors = FALSE
  ),

  ## Single terminal-child shift:
  ## I1 -> C; patch consists only of observed tip C.
  single_terminal = data.frame(
    parent = "I1",
    child = "C",
    stringsAsFactors = FALSE
  ),

  ## Two distant internal patches.
  two_distant_internal = data.frame(
    parent = c("I1", "I2"),
    child = c("J1", "J2"),
    stringsAsFactors = FALSE
  ),

  ## Two adjacent but disjoint patches.
  ##
  ## I1 -> J1
  ## H1 -> D
  ##
  ## Boundaries I1 and H1 are adjacent in the original tree.
  adjacent_disjoint = data.frame(
    parent = c("I1", "H1"),
    child = c("J1", "D"),
    stringsAsFactors = FALSE
  ),

  ## Mixed internal + terminal patches.
  ##
  ## I1 -> J1  (internal patch)
  ## J2 -> E   (terminal patch)
  mixed_internal_terminal = data.frame(
    parent = c("I1", "J2"),
    child = c("J1", "E"),
    stringsAsFactors = FALSE
  )
)


## ------------------------------------------------------------
## 7. Run all valid audits
## ------------------------------------------------------------

audit_results <- do.call(
  rbind,
  lapply(
    names(cases),
    function(case_name) {
      audit_one_case(
        tree = tree,
        shifts = cases[[case_name]],
        case_name = case_name
      )
    }
  )
)

rownames(audit_results) <- NULL

cat("\n")
cat("============================================\n")
cat("Step 2B partition / boundary-node audit\n")
cat("============================================\n\n")

print(
  audit_results,
  row.names = FALSE
)

cat("\n")

if (all(audit_results$pass)) {

  cat(
    "PASS: In all audited valid configurations, ",
    "the mother-side internal nodes removed by pruning-induced ",
    "collapse are exactly the specified shift parents.\n",
    sep = ""
  )

} else {

  cat(
    "ATTENTION: At least one valid configuration caused ",
    "additional mother-side internal nodes to disappear.\n"
  )
}


## ------------------------------------------------------------
## 8. Confirm that a nested configuration is rejected
## ------------------------------------------------------------

nested_case <- data.frame(
  parent = c("H1", "I1"),
  child = c("I1", "J1"),
  stringsAsFactors = FALSE
)

nested_rejected <- FALSE

tryCatch(
  {
    validate_shift_configuration(
      tree,
      nested_case
    )
  },
  error = function(e) {

    nested_rejected <<- TRUE

    cat("\n")
    cat("Nested-shift validation check:\n")
    cat("  Correctly rejected nested configuration.\n")
    cat("  Message: ", conditionMessage(e), "\n", sep = "")
  }
)

if (!nested_rejected) {
  stop(
    "Audit failure: nested shift configuration was not rejected."
  )
}

cat("\nAudit script completed.\n")

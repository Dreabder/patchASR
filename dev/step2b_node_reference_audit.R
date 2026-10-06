## ============================================================
## Step 2B.2 — Node-reference API audit
##
## Purpose:
##   Verify the proposed v0.1 rule that shifts are supplied
##   using ape node numbers from the exact original input tree.
##
## This is a development audit script, not package code.
## ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for this audit.")
}


## ------------------------------------------------------------
## 1. Construct the same small audit tree
## ------------------------------------------------------------

tree_text <- paste0(
  "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
  "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
)

tree <- ape::read.tree(text = tree_text)


## ------------------------------------------------------------
## 2. Create an original-tree node map
## ------------------------------------------------------------

make_node_map <- function(tree) {

  n_tip <- ape::Ntip(tree)
  n_node <- tree$Nnode
  all_nodes <- seq_len(n_tip + n_node)

  node_type <- ifelse(
    all_nodes <= n_tip,
    "tip",
    "internal"
  )

  labels <- character(length(all_nodes))

  labels[all_nodes <= n_tip] <- tree$tip.label

  internal_indices <- all_nodes[all_nodes > n_tip] - n_tip

  if (!is.null(tree$node.label)) {
    labels[all_nodes > n_tip] <- tree$node.label[internal_indices]
  } else {
    labels[all_nodes > n_tip] <- NA_character_
  }

  stable_id <- ifelse(
    node_type == "tip",
    paste0("tip_", all_nodes),
    paste0("node_", all_nodes)
  )

  data.frame(
    original_node = all_nodes,
    node_type = node_type,
    original_label = labels,
    stable_id = stable_id,
    stringsAsFactors = FALSE
  )
}


node_map <- make_node_map(tree)

cat("\nOriginal-tree node map:\n\n")
print(node_map, row.names = FALSE)


## ------------------------------------------------------------
## 3. Helper for finding audit-tree nodes by label
##    This helper is used only to construct our test input.
## ------------------------------------------------------------

node_number_from_label <- function(tree, label) {

  tip_match <- match(label, tree$tip.label)

  if (!is.na(tip_match)) {
    return(as.integer(tip_match))
  }

  internal_match <- match(label, tree$node.label)

  if (is.na(internal_match)) {
    stop("Unknown label: ", label)
  }

  as.integer(ape::Ntip(tree) + internal_match)
}


## ------------------------------------------------------------
## 4. Proposed v0.1 shift-input validator
## ------------------------------------------------------------

validate_node_reference_shifts <- function(tree, shifts) {

  if (!is.data.frame(shifts)) {
    stop("`shifts` must be a data.frame.")
  }

  if (!all(c("parent", "child") %in% names(shifts))) {
    stop("`shifts` must contain `parent` and `child` columns.")
  }

  if (nrow(shifts) < 1L) {
    stop("`shifts` must contain at least one row.")
  }

  validate_integer_reference <- function(x, field) {

    if (!is.numeric(x) ||
        anyNA(x) ||
        any(!is.finite(x)) ||
        any(x != floor(x))) {
      stop(
        "`", field,
        "` must contain finite integer original-tree node numbers."
      )
    }

    as.integer(x)
  }

  parent <- validate_integer_reference(
    shifts$parent,
    "parent"
  )

  child <- validate_integer_reference(
    shifts$child,
    "child"
  )

  n_tip <- ape::Ntip(tree)
  n_total <- n_tip + tree$Nnode

  if (any(parent < 1L | parent > n_total)) {
    stop("At least one `parent` node number is outside the input tree.")
  }

  if (any(child < 1L | child > n_total)) {
    stop("At least one `child` node number is outside the input tree.")
  }

  if (any(parent <= n_tip)) {
    stop("Every shift parent must be an internal node.")
  }

  ## Root detection.
  child_nodes_in_tree <- tree$edge[, 2L]
  root_candidates <- setdiff(
    unique(tree$edge[, 1L]),
    child_nodes_in_tree
  )

  if (length(root_candidates) != 1L) {
    stop("Could not identify a unique root.")
  }

  root <- as.integer(root_candidates)

  if (any(parent == root)) {
    stop("The root cannot serve as a shift parent in version 0.1.")
  }

  ## Every requested pair must be an actual directed edge.
  for (i in seq_along(parent)) {

    edge_rows <- which(
      tree$edge[, 1L] == parent[[i]] &
        tree$edge[, 2L] == child[[i]]
    )

    if (length(edge_rows) != 1L) {
      stop(
        "Requested pair is not an existing directed edge: ",
        parent[[i]], " -> ", child[[i]]
      )
    }
  }

  data.frame(
    parent = parent,
    child = child,
    stringsAsFactors = FALSE
  )
}


## ------------------------------------------------------------
## 5. Construct valid shift inputs using original ape numbers
## ------------------------------------------------------------

I1 <- node_number_from_label(tree, "I1")
J1 <- node_number_from_label(tree, "J1")
C_tip <- node_number_from_label(tree, "C")

valid_internal_shift <- data.frame(
  parent = I1,
  child = J1
)

valid_terminal_shift <- data.frame(
  parent = I1,
  child = C_tip
)


## ------------------------------------------------------------
## 6. Test valid cases
## ------------------------------------------------------------

cat("\nValid internal-child shift:\n")
print(
  validate_node_reference_shifts(
    tree,
    valid_internal_shift
  ),
  row.names = FALSE
)

cat("\nValid terminal-child shift:\n")
print(
  validate_node_reference_shifts(
    tree,
    valid_terminal_shift
  ),
  row.names = FALSE
)


## ------------------------------------------------------------
## 7. Test expected failures
## ------------------------------------------------------------

expect_error <- function(expr, label) {

  rejected <- FALSE

  tryCatch(
    {
      force(expr)
    },
    error = function(e) {

      rejected <<- TRUE

      cat("\nPASS — correctly rejected: ", label, "\n", sep = "")
      cat("Message: ", conditionMessage(e), "\n", sep = "")
    }
  )

  if (!rejected) {
    stop(
      "Audit failure: expected rejection for ",
      label
    )
  }
}


expect_error(
  validate_node_reference_shifts(
    tree,
    data.frame(
      parent = I1 + 0.5,
      child = J1
    )
  ),
  "non-integer parent reference"
)


expect_error(
  validate_node_reference_shifts(
    tree,
    data.frame(
      parent = 9999,
      child = J1
    )
  ),
  "out-of-range parent reference"
)


A_tip <- node_number_from_label(tree, "A")

expect_error(
  validate_node_reference_shifts(
    tree,
    data.frame(
      parent = A_tip,
      child = J1
    )
  ),
  "tip used as shift parent"
)


H1 <- node_number_from_label(tree, "H1")

expect_error(
  validate_node_reference_shifts(
    tree,
    data.frame(
      parent = H1,
      child = J1
    )
  ),
  "parent-child pair that is not a direct edge"
)


ROOT <- node_number_from_label(tree, "ROOT")

expect_error(
  validate_node_reference_shifts(
    tree,
    data.frame(
      parent = ROOT,
      child = H1
    )
  ),
  "root used as shift parent"
)


## ------------------------------------------------------------
## 8. Demonstrate that edge-row index is not the node API
## ------------------------------------------------------------

cat("\nEdge matrix of original tree:\n\n")
print(tree$edge)

tree_reordered <- ape::reorder.phylo(
  tree,
  order = "postorder"
)

cat("\nEdge matrix after reorder.phylo(order = 'postorder'):\n\n")
print(tree_reordered$edge)

same_edge_row_order <- identical(
  tree$edge,
  tree_reordered$edge
)

cat(
  "\nIdentical edge-matrix row order after reordering? ",
  same_edge_row_order,
  "\n",
  sep = ""
)

cat(
  "\nThis illustrates why the core API identifies a shift by ",
  "parent and child node references, not by an edge-row index.\n",
  sep = ""
)

cat("\nStep 2B.2 node-reference audit completed.\n")

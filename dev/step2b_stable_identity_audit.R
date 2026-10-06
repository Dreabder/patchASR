## ============================================================
## Step 2B.3 — Stable node-identity audit
##
## Purpose:
##   Verify that original node identities can be recovered after
##   reorder.phylo(), extract.clade(), and drop.tip(), even when
##   the user's original tree has no internal node labels.
##
## This is a development audit script, not package code.
## ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for this audit.")
}


## ------------------------------------------------------------
## 1. Construct a labelled reference tree only for the audit
## ------------------------------------------------------------

tree_text <- paste0(
  "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
  "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
)

reference_tree <- ape::read.tree(text = tree_text)


## ------------------------------------------------------------
## 2. Helper used only to locate audit nodes in the reference tree
## ------------------------------------------------------------

node_number_from_label <- function(tree, label) {

  tip_match <- match(label, tree$tip.label)

  if (!is.na(tip_match)) {
    return(as.integer(tip_match))
  }

  if (is.null(tree$node.label)) {
    stop("Reference tree has no internal node labels.")
  }

  internal_match <- match(label, tree$node.label)

  if (is.na(internal_match)) {
    stop("Unknown reference-tree label: ", label)
  }

  as.integer(
    ape::Ntip(tree) + internal_match
  )
}


## Record important original node numbers before removing labels.

ROOT_node <- node_number_from_label(reference_tree, "ROOT")
H1_node   <- node_number_from_label(reference_tree, "H1")
I1_node   <- node_number_from_label(reference_tree, "I1")
J1_node   <- node_number_from_label(reference_tree, "J1")
I2_node   <- node_number_from_label(reference_tree, "I2")
J2_node   <- node_number_from_label(reference_tree, "J2")
S2_node   <- node_number_from_label(reference_tree, "S2")


## ------------------------------------------------------------
## 3. Simulate a real user tree with no internal node labels
## ------------------------------------------------------------

user_tree <- reference_tree

## Remove all user-supplied internal labels.
user_tree$node.label <- NULL

if (!is.null(user_tree$node.label)) {
  stop("Audit setup failed: internal node labels were not removed.")
}

original_edge_matrix <- user_tree$edge
original_tip_labels <- user_tree$tip.label
original_edge_lengths <- user_tree$edge.length


## ------------------------------------------------------------
## 4. Freeze original node identity
## ------------------------------------------------------------

freeze_node_identity <- function(tree) {

  n_tip <- ape::Ntip(tree)
  n_internal <- tree$Nnode
  n_total <- n_tip + n_internal

  original_nodes <- seq_len(n_total)

  node_type <- ifelse(
    original_nodes <= n_tip,
    "tip",
    "internal"
  )

  original_label <- rep(
    NA_character_,
    n_total
  )

  original_label[seq_len(n_tip)] <- tree$tip.label

  if (!is.null(tree$node.label)) {

    original_label[
      (n_tip + 1L):n_total
    ] <- tree$node.label
  }

  stable_id <- ifelse(
    node_type == "tip",
    paste0("tip_", original_nodes),
    paste0("node_", original_nodes)
  )

  node_map <- data.frame(
    original_node = original_nodes,
    node_type = node_type,
    original_label = original_label,
    stable_id = stable_id,
    stringsAsFactors = FALSE
  )

  ## Work on a copy, not the user's original tree.
  working_tree <- tree

  ## Internal stable IDs are temporarily carried in node.label.
  working_tree$node.label <- node_map$stable_id[
    node_map$node_type == "internal"
  ]

  list(
    tree = working_tree,
    node_map = node_map
  )
}


frozen <- freeze_node_identity(user_tree)

identity_tree <- frozen$tree
node_map <- frozen$node_map


cat("\n")
cat("============================================\n")
cat("Step 2B.3 stable node-identity audit\n")
cat("============================================\n\n")

cat("Original user tree has internal node labels? ")
cat(!is.null(user_tree$node.label), "\n")

cat("Internal working tree has stable node labels? ")
cat(!is.null(identity_tree$node.label), "\n\n")

cat("Original node map:\n\n")
print(node_map, row.names = FALSE)


## ------------------------------------------------------------
## 5. Confirm that freezing identity was non-destructive
## ------------------------------------------------------------

if (!is.null(user_tree$node.label)) {
  stop(
    "Identity freezing modified the user's original internal labels."
  )
}

if (!identical(user_tree$edge, original_edge_matrix)) {
  stop(
    "Identity freezing modified the user's original edge matrix."
  )
}

if (!identical(user_tree$tip.label, original_tip_labels)) {
  stop(
    "Identity freezing modified the user's original tip labels."
  )
}

if (!identical(user_tree$edge.length, original_edge_lengths)) {
  stop(
    "Identity freezing modified the user's original branch lengths."
  )
}

cat(
  "\nPASS — freezing stable identity did not modify the original user tree.\n"
)


## ------------------------------------------------------------
## 6. Map a derived component back to original stable identities
## ------------------------------------------------------------

map_component_nodes <- function(component_tree, node_map) {

  n_tip_component <- ape::Ntip(component_tree)
  n_internal_component <- component_tree$Nnode

  ## ---- Tips ----

  tip_rows <- match(
    component_tree$tip.label,
    node_map$original_label
  )

  if (anyNA(tip_rows)) {
    stop(
      "At least one component tip cannot be mapped to the original tree."
    )
  }

  if (any(node_map$node_type[tip_rows] != "tip")) {
    stop(
      "A component tip mapped to a non-tip original node."
    )
  }

  tip_map <- data.frame(
    local_node = seq_len(n_tip_component),
    node_type = "tip",
    current_label = component_tree$tip.label,
    stable_id = node_map$stable_id[tip_rows],
    original_node = node_map$original_node[tip_rows],
    stringsAsFactors = FALSE
  )

  ## ---- Internal nodes ----

  if (n_internal_component > 0L) {

    if (is.null(component_tree$node.label)) {
      stop(
        "Derived component lost all internal stable labels."
      )
    }

    if (length(component_tree$node.label) != n_internal_component) {
      stop(
        "Internal label count does not match component Nnode."
      )
    }

    internal_rows <- match(
      component_tree$node.label,
      node_map$stable_id
    )

    if (anyNA(internal_rows)) {
      stop(
        "At least one component internal node has an unknown stable ID."
      )
    }

    if (any(node_map$node_type[internal_rows] != "internal")) {
      stop(
        "A component internal node mapped to a non-internal original node."
      )
    }

    internal_map <- data.frame(
      local_node =
        n_tip_component + seq_len(n_internal_component),
      node_type = "internal",
      current_label = component_tree$node.label,
      stable_id = component_tree$node.label,
      original_node = node_map$original_node[internal_rows],
      stringsAsFactors = FALSE
    )

  } else {

    internal_map <- data.frame(
      local_node = integer(0L),
      node_type = character(0L),
      current_label = character(0L),
      stable_id = character(0L),
      original_node = integer(0L),
      stringsAsFactors = FALSE
    )
  }

  result <- rbind(
    tip_map,
    internal_map
  )

  if (anyDuplicated(result$stable_id)) {
    stop(
      "A derived component contains duplicated stable node IDs."
    )
  }

  result
}


## ------------------------------------------------------------
## 7. Helper for checking an expected component node set
## ------------------------------------------------------------

stable_id_for_original_node <- function(node_map, node) {

  row <- match(
    as.integer(node),
    node_map$original_node
  )

  if (is.na(row)) {
    stop("Unknown original node: ", node)
  }

  node_map$stable_id[[row]]
}


tip_stable_id <- function(node_map, tip_label) {

  rows <- which(
    node_map$node_type == "tip" &
      node_map$original_label == tip_label
  )

  if (length(rows) != 1L) {
    stop(
      "Could not uniquely identify original tip: ",
      tip_label
    )
  }

  node_map$stable_id[[rows]]
}


check_component <- function(
    component_tree,
    node_map,
    expected_stable_ids,
    case_name
) {

  component_map <- map_component_nodes(
    component_tree,
    node_map
  )

  observed <- sort(component_map$stable_id)
  expected <- sort(expected_stable_ids)

  pass <- identical(
    observed,
    expected
  )

  cat("\n--------------------------------------------\n")
  cat("Case: ", case_name, "\n", sep = "")
  cat("--------------------------------------------\n\n")

  print(
    component_map,
    row.names = FALSE
  )

  cat("\nExpected stable IDs:\n")
  print(expected)

  cat("\nObserved stable IDs:\n")
  print(observed)

  cat("\nPASS? ", pass, "\n", sep = "")

  if (!pass) {
    stop(
      "Stable-identity audit failed for case: ",
      case_name
    )
  }

  invisible(component_map)
}


## ------------------------------------------------------------
## 8. Original stable IDs used in expected results
## ------------------------------------------------------------

SID_ROOT <- stable_id_for_original_node(node_map, ROOT_node)
SID_H1   <- stable_id_for_original_node(node_map, H1_node)
SID_I1   <- stable_id_for_original_node(node_map, I1_node)
SID_J1   <- stable_id_for_original_node(node_map, J1_node)
SID_I2   <- stable_id_for_original_node(node_map, I2_node)
SID_J2   <- stable_id_for_original_node(node_map, J2_node)
SID_S2   <- stable_id_for_original_node(node_map, S2_node)

SID_A <- tip_stable_id(node_map, "A")
SID_B <- tip_stable_id(node_map, "B")
SID_C <- tip_stable_id(node_map, "C")
SID_D <- tip_stable_id(node_map, "D")
SID_E <- tip_stable_id(node_map, "E")
SID_F <- tip_stable_id(node_map, "F")
SID_G <- tip_stable_id(node_map, "G")
SID_H <- tip_stable_id(node_map, "H")


## ------------------------------------------------------------
## 9. Audit reorder.phylo()
## ------------------------------------------------------------

reordered_tree <- ape::reorder.phylo(
  identity_tree,
  order = "postorder"
)

check_component(
  component_tree = reordered_tree,
  node_map = node_map,
  expected_stable_ids = node_map$stable_id,
  case_name = "reorder.phylo"
)


## ------------------------------------------------------------
## 10. Audit extract.clade()
##
## Extract patch rooted at original J1.
## Expected nodes:
##   tips A, B
##   internal J1
## ------------------------------------------------------------

patch_J1 <- ape::extract.clade(
  identity_tree,
  node = J1_node
)

patch_map <- check_component(
  component_tree = patch_J1,
  node_map = node_map,
  expected_stable_ids = c(
    SID_A,
    SID_B,
    SID_J1
  ),
  case_name = "extract.clade at J1"
)


## ------------------------------------------------------------
## 11. Audit mother tree after one internal patch
##
## Remove descendants of J1: A and B.
##
## J1 belongs to the patch.
## I1 is collapsed as the boundary.
##
## Mother should retain:
##   tips C,D,E,F,G,H
##   internal ROOT,H1,I2,J2,S2
## ------------------------------------------------------------

mother_internal_patch <- ape::drop.tip(
  identity_tree,
  tip = c("A", "B"),
  collapse.singles = TRUE
)

mother_internal_map <- check_component(
  component_tree = mother_internal_patch,
  node_map = node_map,
  expected_stable_ids = c(
    SID_C,
    SID_D,
    SID_E,
    SID_F,
    SID_G,
    SID_H,
    SID_ROOT,
    SID_H1,
    SID_I2,
    SID_J2,
    SID_S2
  ),
  case_name = "drop.tip internal patch J1"
)


## ------------------------------------------------------------
## 12. Audit mother tree after one terminal patch
##
## Remove terminal shift child C.
##
## I1 collapses.
## J1 remains in the mother component.
## ------------------------------------------------------------

mother_terminal_patch <- ape::drop.tip(
  identity_tree,
  tip = "C",
  collapse.singles = TRUE
)

mother_terminal_map <- check_component(
  component_tree = mother_terminal_patch,
  node_map = node_map,
  expected_stable_ids = c(
    SID_A,
    SID_B,
    SID_D,
    SID_E,
    SID_F,
    SID_G,
    SID_H,
    SID_ROOT,
    SID_H1,
    SID_J1,
    SID_I2,
    SID_J2,
    SID_S2
  ),
  case_name = "drop.tip terminal patch C"
)


## ------------------------------------------------------------
## 13. Audit mixed internal + terminal patching
##
## Internal patch:
##   I1 -> J1
##   remove A,B; J1 belongs to patch; I1 collapses.
##
## Terminal patch:
##   J2 -> E
##   remove E; J2 collapses.
##
## Mother should retain:
##   tips C,D,F,G,H
##   internal ROOT,H1,I2,S2
## ------------------------------------------------------------

mother_mixed <- ape::drop.tip(
  identity_tree,
  tip = c("A", "B", "E"),
  collapse.singles = TRUE
)

mother_mixed_map <- check_component(
  component_tree = mother_mixed,
  node_map = node_map,
  expected_stable_ids = c(
    SID_C,
    SID_D,
    SID_F,
    SID_G,
    SID_H,
    SID_ROOT,
    SID_H1,
    SID_I2,
    SID_S2
  ),
  case_name = "mixed internal + terminal patches"
)


## ------------------------------------------------------------
## 14. Demonstrate that local ape numbers can differ
## ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("Local versus original node-number example\n")
cat("============================================\n\n")

cat("Patch rooted at original J1:\n\n")

print(
  patch_map[
    patch_map$node_type == "internal",
    c(
      "local_node",
      "stable_id",
      "original_node"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\nThe local ape node number may differ from the original node number,\n",
  "but the stable ID preserves the biological identity.\n",
  sep = ""
)


## ------------------------------------------------------------
## 15. Final audit message
## ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("Stable identity audit result\n")
cat("============================================\n\n")

cat(
  "PASS: stable node identities were recovered correctly after ",
  "reorder.phylo(), extract.clade(), and drop.tip() in all audited ",
  "configurations.\n",
  sep = ""
)

cat(
  "PASS: the user's original tree remained unchanged.\n"
)

cat(
  "PASS: internal user-supplied node labels were not required.\n"
)

cat("\nStep 2B.3 stable node-identity audit completed.\n")

## ============================================================
## Step 2B.4 — ASR interface contract audit
##
## Purpose:
##   Validate the proposed component-level ASR interface:
##
##       asr_fun(tree, states)
##
##   The audit checks:
##     1. component states are selected and reordered correctly;
##     2. valid ASR results are accepted and standardized;
##     3. malformed ASR results are rejected;
##     4. return order does not determine node identity;
##     5. terminal patches do not call asr_fun.
##
## This is a development audit script, not package code.
## ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for this audit.")
}


## ------------------------------------------------------------
## 1. Construct the same reference tree
## ------------------------------------------------------------

tree_text <- paste0(
  "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
  "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
)

reference_tree <- ape::read.tree(text = tree_text)


## ------------------------------------------------------------
## 2. Helper used only to locate audit nodes
## ------------------------------------------------------------

node_number_from_label <- function(tree, label) {

  tip_match <- match(label, tree$tip.label)

  if (!is.na(tip_match)) {
    return(as.integer(tip_match))
  }

  internal_match <- match(
    label,
    tree$node.label
  )

  if (is.na(internal_match)) {
    stop("Unknown reference-tree label: ", label)
  }

  as.integer(
    ape::Ntip(tree) + internal_match
  )
}


J1_node <- node_number_from_label(
  reference_tree,
  "J1"
)


## ------------------------------------------------------------
## 3. Simulate a user tree without internal node labels
## ------------------------------------------------------------

user_tree <- reference_tree
user_tree$node.label <- NULL


## ------------------------------------------------------------
## 4. Freeze stable node identities
## ------------------------------------------------------------

freeze_node_identity <- function(tree) {

  n_tip <- ape::Ntip(tree)
  n_internal <- tree$Nnode
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

  stable_id <- ifelse(
    node_type == "tip",
    paste0("tip_", original_node),
    paste0("node_", original_node)
  )

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


frozen <- freeze_node_identity(
  user_tree
)

identity_tree <- frozen$tree
node_map <- frozen$node_map


## ------------------------------------------------------------
## 5. Observed tip states
##
## Intentionally use a scrambled order.
## ------------------------------------------------------------

states <- c(
  H = 8,
  B = 2,
  F = 6,
  A = 1,
  D = 4,
  G = 7,
  C = 3,
  E = 5
)

cat("\n")
cat("============================================\n")
cat("Step 2B.4 ASR interface contract audit\n")
cat("============================================\n\n")

cat("Original state-vector order:\n")
print(names(states))


## ------------------------------------------------------------
## 6. Prepare component-specific states
## ------------------------------------------------------------

prepare_component_states <- function(
    component_tree,
    states
) {

  if (!inherits(component_tree, "phylo")) {
    stop("`component_tree` must inherit from class `phylo`.")
  }

  if (!is.numeric(states) ||
      is.null(names(states))) {
    stop("`states` must be a named numeric vector.")
  }

  if (anyNA(names(states)) ||
      any(names(states) == "") ||
      anyDuplicated(names(states))) {
    stop("`states` must have unique non-empty names.")
  }

  if (any(!is.finite(states))) {
    stop("All observed state values must be finite.")
  }

  component_tips <- component_tree$tip.label

  missing_tips <- setdiff(
    component_tips,
    names(states)
  )

  if (length(missing_tips) > 0L) {
    stop(
      "Missing observed state for component tip(s): ",
      paste(missing_tips, collapse = ", ")
    )
  }

  component_states <- states[
    component_tips
  ]

  if (!identical(
    names(component_states),
    component_tree$tip.label
  )) {
    stop(
      "Internal error: component states were not reordered correctly."
    )
  }

  component_states
}


## ------------------------------------------------------------
## 7. Validate and standardize an ASR result
## ------------------------------------------------------------

validate_asr_result <- function(
    component_tree,
    result
) {

  n_internal <- component_tree$Nnode

  if (n_internal < 1L) {
    stop(
      "A component with no internal nodes must not be sent to asr_fun."
    )
  }

  expected_ids <- component_tree$node.label

  if (is.null(expected_ids) ||
      length(expected_ids) != n_internal ||
      anyNA(expected_ids) ||
      any(expected_ids == "") ||
      anyDuplicated(expected_ids)) {

    stop(
      "Component internal stable IDs are missing or invalid."
    )
  }

  if (!is.numeric(result)) {
    stop(
      "`asr_fun` must return a numeric vector."
    )
  }

  if (is.null(names(result))) {
    stop(
      "`asr_fun` must return a named numeric vector."
    )
  }

  if (length(result) != n_internal) {
    stop(
      "`asr_fun` returned ",
      length(result),
      " estimate(s), but ",
      n_internal,
      " internal node(s) were expected."
    )
  }

  if (anyNA(names(result)) ||
      any(names(result) == "")) {
    stop(
      "Every ASR estimate must have a non-empty node identifier."
    )
  }

  if (anyDuplicated(names(result))) {
    stop(
      "`asr_fun` returned duplicated node identifiers."
    )
  }

  if (any(!is.finite(result))) {
    stop(
      "`asr_fun` returned non-finite ancestral-state estimates."
    )
  }

  missing_ids <- setdiff(
    expected_ids,
    names(result)
  )

  unexpected_ids <- setdiff(
    names(result),
    expected_ids
  )

  if (length(missing_ids) > 0L ||
      length(unexpected_ids) > 0L) {

    stop(
      "`asr_fun` returned an incorrect internal-node set. ",
      "Missing: ",
      if (length(missing_ids) == 0L) {
        "<none>"
      } else {
        paste(missing_ids, collapse = ", ")
      },
      "; unexpected: ",
      if (length(unexpected_ids) == 0L) {
        "<none>"
      } else {
        paste(unexpected_ids, collapse = ", ")
      },
      "."
    )
  }

  ## Result order is standardized using stable node identity.
  standardized <- result[
    expected_ids
  ]

  if (!identical(
    names(standardized),
    expected_ids
  )) {
    stop(
      "Internal error while standardizing ASR result order."
    )
  }

  standardized
}


## ------------------------------------------------------------
## 8. Component-level ASR runner
## ------------------------------------------------------------

run_component_asr <- function(
    component_tree = NULL,
    states,
    asr_fun,
    requires_asr = TRUE
) {

  if (!isTRUE(requires_asr)) {

    ## Terminal patch:
    ## no ancestral node exists and asr_fun must not be called.
    return(
      setNames(
        numeric(0L),
        character(0L)
      )
    )
  }

  if (!inherits(component_tree, "phylo")) {
    stop(
      "An ASR-requiring component must contain a `phylo` tree."
    )
  }

  if (component_tree$Nnode < 1L) {
    stop(
      "An ASR-requiring component must contain at least one internal node."
    )
  }

  if (!is.function(asr_fun)) {
    stop(
      "`asr_fun` must be a function."
    )
  }

  component_states <- prepare_component_states(
    component_tree,
    states
  )

  raw_result <- asr_fun(
    component_tree,
    component_states
  )

  validate_asr_result(
    component_tree,
    raw_result
  )
}


## ------------------------------------------------------------
## 9. Construct an internal patch component
## ------------------------------------------------------------

patch_J1 <- ape::extract.clade(
  identity_tree,
  node = J1_node
)

cat("\nPatch J1 tips:\n")
print(patch_J1$tip.label)

cat("\nPatch J1 stable internal IDs:\n")
print(patch_J1$node.label)


## ------------------------------------------------------------
## 10. Valid mock ASR backend
##
## Important:
##   This backend explicitly checks that tip states arrive in
##   tree$tip.label order.
##
## It deliberately returns internal estimates in reverse order,
## demonstrating that biological node identity comes from names,
## not from return-vector position.
## ------------------------------------------------------------

valid_asr_fun <- function(
    tree,
    states
) {

  if (!identical(
    names(states),
    tree$tip.label
  )) {
    stop(
      "Mock ASR backend received states in the wrong order."
    )
  }

  ids <- tree$node.label

  values <- seq_along(ids) + mean(states)

  result <- setNames(
    values,
    ids
  )

  ## Deliberately reverse output order.
  result[
    rev(seq_along(result))
  ]
}


cat("\n")
cat("--------------------------------------------\n")
cat("Valid ASR backend\n")
cat("--------------------------------------------\n\n")

valid_result <- run_component_asr(
  component_tree = patch_J1,
  states = states,
  asr_fun = valid_asr_fun,
  requires_asr = TRUE
)

print(valid_result)

if (!identical(
  names(valid_result),
  patch_J1$node.label
)) {
  stop(
    "Audit failure: valid ASR result was not standardized correctly."
  )
}

cat(
  "\nPASS — component states were correctly selected and ordered.\n"
)

cat(
  "PASS — ASR results were restored to stable-node order by name.\n"
)


## ------------------------------------------------------------
## 11. Error-testing helper
## ------------------------------------------------------------

expect_error <- function(
    expr,
    label
) {

  rejected <- FALSE

  tryCatch(
    {
      force(expr)
    },
    error = function(e) {

      rejected <<- TRUE

      cat(
        "\nPASS — correctly rejected: ",
        label,
        "\n",
        sep = ""
      )

      cat(
        "Message: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )
    }
  )

  if (!rejected) {
    stop(
      "Audit failure: expected rejection for ",
      label
    )
  }
}


## ------------------------------------------------------------
## 12. Missing internal-node estimate
## ------------------------------------------------------------

missing_node_asr_fun <- function(
    tree,
    states
) {

  ids <- tree$node.label

  if (length(ids) == 1L) {
    return(
      setNames(
        numeric(0L),
        character(0L)
      )
    )
  }

  setNames(
    seq_len(length(ids) - 1L),
    ids[-length(ids)]
  )
}


expect_error(
  run_component_asr(
    component_tree = patch_J1,
    states = states,
    asr_fun = missing_node_asr_fun
  ),
  "missing internal-node estimate"
)


## ------------------------------------------------------------
## 13. Unnamed result
## ------------------------------------------------------------

unnamed_asr_fun <- function(
    tree,
    states
) {

  rep(
    1,
    tree$Nnode
  )
}


expect_error(
  run_component_asr(
    component_tree = patch_J1,
    states = states,
    asr_fun = unnamed_asr_fun
  ),
  "unnamed ASR result"
)


## ------------------------------------------------------------
## 14. Duplicated node identifier
##
## Use the mother component because it has several internal nodes.
## ------------------------------------------------------------

mother_tree <- ape::drop.tip(
  identity_tree,
  tip = c("A", "B"),
  collapse.singles = TRUE
)

## ------------------------------------------------------------
## 14A. Valid multi-node result returned in reverse order
##
## This is the actual order-invariance test.
## ------------------------------------------------------------

cat("\n")
cat("--------------------------------------------\n")
cat("Valid multi-node ASR backend with reversed output order\n")
cat("--------------------------------------------\n\n")

valid_mother_result <- run_component_asr(
  component_tree = mother_tree,
  states = states,
  asr_fun = valid_asr_fun,
  requires_asr = TRUE
)

mother_states <- prepare_component_states(
  mother_tree,
  states
)

expected_mother_result <- setNames(
  seq_along(mother_tree$node.label) +
    mean(mother_states),
  mother_tree$node.label
)

cat("Expected stable-node order:\n")
print(names(expected_mother_result))

cat("\nStandardized returned order:\n")
print(names(valid_mother_result))

cat("\nExpected values:\n")
print(expected_mother_result)

cat("\nStandardized values:\n")
print(valid_mother_result)

if (!identical(
  names(valid_mother_result),
  names(expected_mother_result)
)) {
  stop(
    "Audit failure: multi-node ASR result was not restored ",
    "to stable-node order."
  )
}

if (!isTRUE(all.equal(
  unname(valid_mother_result),
  unname(expected_mother_result),
  tolerance = 1e-12
))) {
  stop(
    "Audit failure: node-associated ASR values changed during ",
    "order standardization."
  )
}

cat(
  "\nPASS — multi-node ASR output order was standardized ",
  "by stable node identity.\n",
  sep = ""
)

cat(
  "PASS — node-associated values remained attached to the ",
  "correct stable node IDs.\n",
  sep = ""
)

duplicated_id_asr_fun <- function(
    tree,
    states
) {

  ids <- tree$node.label

  bad_names <- ids

  if (length(bad_names) < 2L) {
    stop(
      "Audit setup requires at least two internal nodes."
    )
  }

  bad_names[[2L]] <- bad_names[[1L]]

  setNames(
    seq_along(ids),
    bad_names
  )
}


expect_error(
  run_component_asr(
    component_tree = mother_tree,
    states = states,
    asr_fun = duplicated_id_asr_fun
  ),
  "duplicated internal-node identifier"
)


## ------------------------------------------------------------
## 15. Unknown internal-node identifier
## ------------------------------------------------------------

unknown_id_asr_fun <- function(
    tree,
    states
) {

  ids <- tree$node.label

  ids[[1L]] <- "unknown_node"

  setNames(
    seq_along(ids),
    ids
  )
}


expect_error(
  run_component_asr(
    component_tree = mother_tree,
    states = states,
    asr_fun = unknown_id_asr_fun
  ),
  "unknown internal-node identifier"
)


## ------------------------------------------------------------
## 16. Non-finite estimate
## ------------------------------------------------------------

nonfinite_asr_fun <- function(
    tree,
    states
) {

  ids <- tree$node.label

  values <- rep(
    1,
    length(ids)
  )

  values[[1L]] <- NA_real_

  setNames(
    values,
    ids
  )
}


expect_error(
  run_component_asr(
    component_tree = mother_tree,
    states = states,
    asr_fun = nonfinite_asr_fun
  ),
  "non-finite ancestral-state estimate"
)


## ------------------------------------------------------------
## 17. Missing observed tip state
## ------------------------------------------------------------

states_missing_A <- states[
  names(states) != "A"
]

expect_error(
  run_component_asr(
    component_tree = patch_J1,
    states = states_missing_A,
    asr_fun = valid_asr_fun
  ),
  "missing observed component-tip state"
)


## ------------------------------------------------------------
## 18. Terminal patch must not call asr_fun
## ------------------------------------------------------------

terminal_call_count <- 0L

must_not_be_called <- function(
    tree,
    states
) {

  terminal_call_count <<-
    terminal_call_count + 1L

  stop(
    "This ASR function must not be called for a terminal patch."
  )
}


terminal_result <- run_component_asr(
  component_tree = NULL,
  states = states,
  asr_fun = must_not_be_called,
  requires_asr = FALSE
)

if (terminal_call_count != 0L) {
  stop(
    "Audit failure: asr_fun was called for a terminal patch."
  )
}

if (length(terminal_result) != 0L) {
  stop(
    "Audit failure: terminal patch returned ancestral estimates."
  )
}

cat(
  "\nPASS — terminal patch did not call asr_fun.\n"
)

cat(
  "PASS — terminal patch returned an empty ancestral-state result.\n"
)


## ------------------------------------------------------------
## 19. Final result
## ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("ASR interface audit result\n")
cat("============================================\n\n")

cat(
  "PASS: component tip states were subset and reordered correctly.\n"
)

cat(
  "PASS: valid ASR estimates were matched by stable node identity, ",
  "not by vector position.\n",
  sep = ""
)

cat(
  "PASS: malformed ASR outputs were rejected.\n"
)

cat(
  "PASS: missing observed component states were rejected.\n"
)

cat(
  "PASS: terminal patches bypassed the ASR backend.\n"
)

cat("\nStep 2B.4 ASR interface contract audit completed.\n")

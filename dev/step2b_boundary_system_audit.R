## ============================================================
## Step 2B.5 — Boundary-system audit
##
## Purpose:
##   Validate the branch-length interpolation and the joint
##   multi-boundary linear system:
##
##       A x = b
##
##   The audit checks:
##     1. a single boundary against direct interpolation;
##     2. two adjacent coupled boundaries;
##     3. mixed internal + terminal shifts;
##     4. shift-row order invariance;
##     5. missing known states are rejected;
##     6. unresolved boundaries cannot be treated as known states.
##
## This is development audit code, not package code.
## ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required for this audit.")
}


## ------------------------------------------------------------
## 1. Construct a test tree with unequal branch lengths
## ------------------------------------------------------------

tree_text <- paste0(
  "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
  "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
)

reference_tree <- ape::read.tree(
  text = tree_text
)


## ------------------------------------------------------------
## 2. Locate audit nodes before removing user labels
## ------------------------------------------------------------

node_number_from_label <- function(
    tree,
    label
) {

  tip_match <- match(
    label,
    tree$tip.label
  )

  if (!is.na(tip_match)) {
    return(as.integer(tip_match))
  }

  internal_match <- match(
    label,
    tree$node.label
  )

  if (is.na(internal_match)) {
    stop(
      "Unknown node label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) + internal_match
  )
}


ROOT <- node_number_from_label(reference_tree, "ROOT")
H1   <- node_number_from_label(reference_tree, "H1")
I1   <- node_number_from_label(reference_tree, "I1")
J1   <- node_number_from_label(reference_tree, "J1")
I2   <- node_number_from_label(reference_tree, "I2")
J2   <- node_number_from_label(reference_tree, "J2")
S2   <- node_number_from_label(reference_tree, "S2")

C_tip <- node_number_from_label(reference_tree, "C")
D_tip <- node_number_from_label(reference_tree, "D")
E_tip <- node_number_from_label(reference_tree, "E")
F_tip <- node_number_from_label(reference_tree, "F")


## ------------------------------------------------------------
## 3. Simulate an unlabeled user tree and freeze stable IDs
## ------------------------------------------------------------

user_tree <- reference_tree
user_tree$node.label <- NULL


freeze_node_identity <- function(tree) {

  n_tip <- ape::Ntip(tree)
  n_total <- n_tip + tree$Nnode

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

tree <- frozen$tree
node_map <- frozen$node_map


stable_id_for_node <- function(
    node_map,
    node
) {

  row <- match(
    as.integer(node),
    node_map$original_node
  )

  if (is.na(row)) {
    stop(
      "Unknown original node: ",
      node
    )
  }

  node_map$stable_id[[row]]
}


SID_ROOT <- stable_id_for_node(node_map, ROOT)
SID_H1   <- stable_id_for_node(node_map, H1)
SID_I1   <- stable_id_for_node(node_map, I1)
SID_J1   <- stable_id_for_node(node_map, J1)
SID_I2   <- stable_id_for_node(node_map, I2)
SID_J2   <- stable_id_for_node(node_map, J2)
SID_S2   <- stable_id_for_node(node_map, S2)


## ------------------------------------------------------------
## 4. Basic topology helpers
## ------------------------------------------------------------

parent_of <- function(
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
      "Node does not have exactly one parent."
    )
  }

  as.integer(
    tree$edge[rows, 1L]
  )
}


children_of <- function(
    tree,
    node
) {

  as.integer(
    tree$edge[
      tree$edge[, 1L] ==
        as.integer(node),
      2L
    ]
  )
}


edge_length_between <- function(
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
      "Requested parent-child edge was not found uniquely."
    )
  }

  value <- tree$edge.length[[rows]]

  if (!is.finite(value) ||
      value <= 0) {
    stop(
      "Boundary branch length must be finite and positive."
    )
  }

  as.numeric(value)
}


## ------------------------------------------------------------
## 5. Construct a synthetic unified known-state pool
##
## Tips represent observed states.
## Internal values represent mock component-level ASR estimates.
## ------------------------------------------------------------

observed_tip_states <- c(
  A = 1,
  B = 2,
  C = 3,
  D = 4,
  E = 5,
  F = 6,
  G = 7,
  H = 8
)


tip_rows <- which(
  node_map$node_type == "tip"
)

tip_values <- observed_tip_states[
  node_map$original_label[tip_rows]
]

tip_known_states <- setNames(
  as.numeric(tip_values),
  node_map$stable_id[tip_rows]
)


internal_known_states <- setNames(
  c(
    10,  ## ROOT
    20,  ## H1
    30,  ## I1
    40,  ## J1
    50,  ## I2
    60,  ## J2
    70   ## S2
  ),
  c(
    SID_ROOT,
    SID_H1,
    SID_I1,
    SID_J1,
    SID_I2,
    SID_J2,
    SID_S2
  )
)


all_known_states <- c(
  tip_known_states,
  internal_known_states
)


build_known_states <- function(
    boundary_nodes
) {

  boundary_ids <- vapply(
    boundary_nodes,
    function(x) {
      stable_id_for_node(
        node_map,
        x
      )
    },
    character(1)
  )

  all_known_states[
    setdiff(
      names(all_known_states),
      boundary_ids
    )
  ]
}


## ------------------------------------------------------------
## 6. Build one complete boundary linear system
## ------------------------------------------------------------

build_boundary_system <- function(
    tree,
    node_map,
    shifts,
    known_states
) {

  if (!is.data.frame(shifts) ||
      !all(c("parent", "child") %in% names(shifts))) {
    stop(
      "`shifts` must contain `parent` and `child` columns."
    )
  }

  if (nrow(shifts) < 1L) {
    stop(
      "At least one boundary equation is required."
    )
  }

  boundary_nodes <- as.integer(
    shifts$parent
  )

  if (anyDuplicated(boundary_nodes)) {
    stop(
      "Boundary nodes must be unique."
    )
  }

  boundary_ids <- vapply(
    boundary_nodes,
    function(x) {
      stable_id_for_node(
        node_map,
        x
      )
    },
    character(1)
  )

  if (!is.numeric(known_states) ||
      is.null(names(known_states)) ||
      anyDuplicated(names(known_states)) ||
      any(!is.finite(known_states))) {
    stop(
      "`known_states` must be a finite named numeric vector ",
      "with unique names."
    )
  }

  if (any(
    boundary_ids %in%
    names(known_states)
  )) {
    stop(
      "An unresolved boundary node must not be present ",
      "in the known-state pool."
    )
  }

  m <- length(boundary_nodes)

  A <- diag(m)

  dimnames(A) <- list(
    boundary_ids,
    boundary_ids
  )

  b <- setNames(
    numeric(m),
    boundary_ids
  )

  details <- vector(
    "list",
    m
  )

  for (i in seq_len(m)) {

    I <- boundary_nodes[[i]]
    J <- as.integer(
      shifts$child[[i]]
    )

    H <- parent_of(
      tree,
      I
    )

    if (is.na(H)) {
      stop(
        "A boundary parent cannot be the root."
      )
    }

    children <- children_of(
      tree,
      I
    )

    if (length(children) != 2L) {
      stop(
        "Boundary parent must have exactly two children."
      )
    }

    if (!(J %in% children)) {
      stop(
        "Shift child is not a direct child of its boundary parent."
      )
    }

    sibling <- setdiff(
      children,
      J
    )

    if (length(sibling) != 1L) {
      stop(
        "Could not identify a unique sibling-side node."
      )
    }

    d_HI <- edge_length_between(
      tree,
      H,
      I
    )

    d_IS <- edge_length_between(
      tree,
      I,
      sibling
    )

    w_H <- d_IS /
      (d_HI + d_IS)

    w_S <- d_HI /
      (d_HI + d_IS)

    if (!isTRUE(all.equal(
      w_H + w_S,
      1,
      tolerance = 1e-12
    ))) {
      stop(
        "Boundary interpolation weights do not sum to one."
      )
    }

    adjacent_nodes <- c(
      H,
      sibling
    )

    adjacent_weights <- c(
      w_H,
      w_S
    )

    for (q in seq_along(adjacent_nodes)) {

      adjacent_node <-
        adjacent_nodes[[q]]

      weight <-
        adjacent_weights[[q]]

      adjacent_id <-
        stable_id_for_node(
          node_map,
          adjacent_node
        )

      if (adjacent_node %in%
          boundary_nodes) {

        column <- match(
          adjacent_node,
          boundary_nodes
        )

        A[i, column] <-
          A[i, column] - weight

      } else {

        if (!(adjacent_id %in%
              names(known_states))) {
          stop(
            "Required known state is missing for node ",
            adjacent_id,
            "."
          )
        }

        b[[i]] <-
          b[[i]] +
          weight *
          known_states[[adjacent_id]]
      }
    }

    details[[i]] <- data.frame(
      boundary_id =
        boundary_ids[[i]],
      parent_side_id =
        stable_id_for_node(
          node_map,
          H
        ),
      sibling_side_id =
        stable_id_for_node(
          node_map,
          sibling
        ),
      d_parent_boundary =
        d_HI,
      d_boundary_sibling =
        d_IS,
      weight_parent =
        w_H,
      weight_sibling =
        w_S,
      stringsAsFactors = FALSE
    )
  }

  list(
    A = A,
    b = b,
    boundary_nodes = boundary_nodes,
    boundary_ids = boundary_ids,
    details = do.call(
      rbind,
      details
    )
  )
}


## ------------------------------------------------------------
## 7. Solve and validate the boundary system
## ------------------------------------------------------------

solve_boundary_system <- function(
    tree,
    node_map,
    shifts,
    known_states
) {

  system <- build_boundary_system(
    tree = tree,
    node_map = node_map,
    shifts = shifts,
    known_states = known_states
  )

  if (qr(system$A)$rank <
      nrow(system$A)) {
    stop(
      "Boundary coefficient matrix is singular."
    )
  }

  raw_solution <- solve(
    system$A,
    system$b
  )

  solution <- setNames(
    as.numeric(raw_solution),
    system$boundary_ids
  )

  if (any(!is.finite(solution))) {
    stop(
      "Boundary system produced non-finite estimates."
    )
  }

  residual <- max(
    abs(
      as.numeric(
        system$A %*% solution
      ) -
        as.numeric(system$b)
    )
  )

  if (!is.finite(residual) ||
      residual > 1e-10) {
    stop(
      "Boundary-system residual is too large: ",
      residual
    )
  }

  list(
    solution = solution,
    system = system,
    residual = residual
  )
}


order_by_name <- function(x) {

  x[
    order(names(x))
  ]
}


## ------------------------------------------------------------
## 8. Single-boundary audit
##
## Shift:
##   I1 -> J1
##
## Original geometry:
##
## H1 --3--> I1 --1--> C
##
## Therefore:
##
## X_I1 = 0.25 X_H1 + 0.75 X_C
##
## with H1 = 20 and C = 3:
##
## X_I1 = 7.25
## ------------------------------------------------------------

single_shift <- data.frame(
  parent = I1,
  child = J1
)

single_known <- build_known_states(
  boundary_nodes = I1
)

single_result <- solve_boundary_system(
  tree = tree,
  node_map = node_map,
  shifts = single_shift,
  known_states = single_known
)

cat("\n")
cat("============================================\n")
cat("Single-boundary audit\n")
cat("============================================\n\n")

cat("Boundary geometry:\n")
print(
  single_result$system$details,
  row.names = FALSE
)

cat("\nA matrix:\n")
print(single_result$system$A)

cat("\nb vector:\n")
print(single_result$system$b)

cat("\nBoundary solution:\n")
print(single_result$solution)

cat(
  "\nResidual: ",
  single_result$residual,
  "\n",
  sep = ""
)

expected_single <- setNames(
  0.25 * 20 +
    0.75 * 3,
  SID_I1
)

if (!isTRUE(all.equal(
  single_result$solution,
  expected_single,
  tolerance = 1e-12,
  check.attributes = TRUE
))) {
  stop(
    "Single-boundary interpolation does not match ",
    "the direct expected value."
  )
}

cat(
  "\nPASS — single-boundary solution matches direct ",
  "branch-length interpolation.\n",
  sep = ""
)


## ------------------------------------------------------------
## 9. Adjacent coupled-boundary audit
##
## Shifts:
##
##   I1 -> J1
##   H1 -> D
##
## Boundary equations:
##
## I1 = 0.25 H1 + 0.75 C
## H1 = 0.60 ROOT + 0.40 I1
##
## with C = 3 and ROOT = 10:
##
## I1 = 25/6
## H1 = 23/3
## ------------------------------------------------------------

adjacent_shifts <- data.frame(
  parent = c(
    I1,
    H1
  ),
  child = c(
    J1,
    D_tip
  )
)

adjacent_known <- build_known_states(
  boundary_nodes =
    adjacent_shifts$parent
)

adjacent_result <- solve_boundary_system(
  tree = tree,
  node_map = node_map,
  shifts = adjacent_shifts,
  known_states = adjacent_known
)

cat("\n")
cat("============================================\n")
cat("Adjacent coupled-boundary audit\n")
cat("============================================\n\n")

cat("Boundary geometry:\n")
print(
  adjacent_result$system$details,
  row.names = FALSE
)

cat("\nA matrix:\n")
print(adjacent_result$system$A)

cat("\nb vector:\n")
print(adjacent_result$system$b)

cat("\nBoundary solution:\n")
print(adjacent_result$solution)

cat(
  "\nResidual: ",
  adjacent_result$residual,
  "\n",
  sep = ""
)

expected_adjacent <- setNames(
  c(
    25 / 6,
    23 / 3
  ),
  c(
    SID_I1,
    SID_H1
  )
)

if (!isTRUE(all.equal(
  adjacent_result$solution,
  expected_adjacent,
  tolerance = 1e-12,
  check.attributes = TRUE
))) {
  stop(
    "Adjacent coupled-boundary solution does not match ",
    "the analytically expected result."
  )
}

cat(
  "\nPASS — adjacent coupled boundaries match the ",
  "analytical solution.\n",
  sep = ""
)


## ------------------------------------------------------------
## 10. Shift-row order invariance
## ------------------------------------------------------------

adjacent_shifts_reversed <-
  adjacent_shifts[
    nrow(adjacent_shifts):1L,
    ,
    drop = FALSE
  ]

adjacent_known_reversed <-
  build_known_states(
    boundary_nodes =
      adjacent_shifts_reversed$parent
  )

adjacent_reversed_result <-
  solve_boundary_system(
    tree = tree,
    node_map = node_map,
    shifts =
      adjacent_shifts_reversed,
    known_states =
      adjacent_known_reversed
  )

if (!isTRUE(all.equal(
  order_by_name(
    adjacent_result$solution
  ),
  order_by_name(
    adjacent_reversed_result$solution
  ),
  tolerance = 1e-12,
  check.attributes = TRUE
))) {
  stop(
    "Boundary solution changed when shift rows were reordered."
  )
}

cat(
  "\nPASS — reversing shift-row order did not change ",
  "the biological boundary solution.\n",
  sep = ""
)


## ------------------------------------------------------------
## 11. Mixed internal + terminal shift audit
##
## Shifts:
##
##   I1 -> J1   internal-child patch
##   J2 -> E    terminal-child patch
##
## Expected:
##
## I1 = 0.25*20 + 0.75*3 = 7.25
##
## J2:
##   I2 --2--> J2 --3--> F
##
## X_J2 = 0.60*50 + 0.40*6 = 32.4
## ------------------------------------------------------------

mixed_shifts <- data.frame(
  parent = c(
    I1,
    J2
  ),
  child = c(
    J1,
    E_tip
  )
)

mixed_known <- build_known_states(
  boundary_nodes =
    mixed_shifts$parent
)

mixed_result <- solve_boundary_system(
  tree = tree,
  node_map = node_map,
  shifts = mixed_shifts,
  known_states = mixed_known
)

cat("\n")
cat("============================================\n")
cat("Mixed internal + terminal boundary audit\n")
cat("============================================\n\n")

cat("Boundary geometry:\n")
print(
  mixed_result$system$details,
  row.names = FALSE
)

cat("\nBoundary solution:\n")
print(mixed_result$solution)

expected_mixed <- setNames(
  c(
    7.25,
    32.4
  ),
  c(
    SID_I1,
    SID_J2
  )
)

if (!isTRUE(all.equal(
  mixed_result$solution,
  expected_mixed,
  tolerance = 1e-12,
  check.attributes = TRUE
))) {
  stop(
    "Mixed internal + terminal boundary solution ",
    "does not match the expected result."
  )
}

cat(
  "\nPASS — mixed internal + terminal shifts were ",
  "reconstructed correctly.\n",
  sep = ""
)


## ------------------------------------------------------------
## 12. Error-testing helper
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
## 13. Missing required known state
## ------------------------------------------------------------

missing_C_known <-
  single_known[
    names(single_known) !=
      stable_id_for_node(
        node_map,
        C_tip
      )
  ]

expect_error(
  solve_boundary_system(
    tree = tree,
    node_map = node_map,
    shifts = single_shift,
    known_states = missing_C_known
  ),
  "missing required adjacent known state"
)


## ------------------------------------------------------------
## 14. Unresolved boundary incorrectly inserted as known
## ------------------------------------------------------------

bad_boundary_known <- c(
  single_known,
  setNames(
    999,
    SID_I1
  )
)

expect_error(
  solve_boundary_system(
    tree = tree,
    node_map = node_map,
    shifts = single_shift,
    known_states = bad_boundary_known
  ),
  "unresolved boundary present in known-state pool"
)


## ------------------------------------------------------------
## 15. Final result
## ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("Boundary-system audit result\n")
cat("============================================\n\n")

cat(
  "PASS: single-boundary interpolation matched the direct formula.\n"
)

cat(
  "PASS: adjacent boundary nodes were solved jointly and correctly.\n"
)

cat(
  "PASS: mixed internal and terminal shift boundaries were handled.\n"
)

cat(
  "PASS: the biological solution was invariant to shift-row order.\n"
)

cat(
  "PASS: missing adjacent known states were rejected.\n"
)

cat(
  "PASS: unresolved boundaries could not be treated as known states.\n"
)

cat("\nStep 2B.5 boundary-system audit completed.\n")

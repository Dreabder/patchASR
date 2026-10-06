make_validation_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
      "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
    )
  )
}


validation_node <- function(
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
      "Unknown test-tree label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


make_valid_states <- function() {

  c(
    H = 8,
    B = 2,
    F = 6,
    A = 1,
    D = 4,
    G = 7,
    C = 3,
    E = 5
  )
}


test_that(
  "valid core inputs are accepted and normalized",
  {

    tree <- make_validation_test_tree()

    shifts <- data.frame(
      parent = validation_node(
        tree,
        "I1"
      ),
      child = validation_node(
        tree,
        "J1"
      )
    )

    result <- validate_patch_inputs(
      tree = tree,
      states = make_valid_states(),
      shifts = shifts,
      asr_fun = function(tree, states) {
        numeric(0L)
      }
    )

    expect_identical(
      names(result$states),
      tree$tip.label
    )

    expect_type(
      result$shifts$parent,
      "integer"
    )

    expect_type(
      result$shifts$child,
      "integer"
    )

    expect_true(
      is.function(result$asr_fun)
    )
  }
)


test_that(
  "valid internal and terminal child shifts are accepted",
  {

    tree <- make_validation_test_tree()

    internal_shift <- data.frame(
      parent = validation_node(
        tree,
        "I1"
      ),
      child = validation_node(
        tree,
        "J1"
      )
    )

    terminal_shift <- data.frame(
      parent = validation_node(
        tree,
        "I1"
      ),
      child = validation_node(
        tree,
        "C"
      )
    )

    expect_no_error(
      validate_shifts(
        tree,
        internal_shift
      )
    )

    expect_no_error(
      validate_shifts(
        tree,
        terminal_shift
      )
    )
  }
)


test_that(
  "valid pairwise-disjoint multiple shifts are accepted",
  {

    tree <- make_validation_test_tree()

    shifts <- data.frame(
      parent = c(
        validation_node(tree, "I1"),
        validation_node(tree, "I2")
      ),
      child = c(
        validation_node(tree, "J1"),
        validation_node(tree, "J2")
      )
    )

    expect_no_error(
      validate_shifts(
        tree,
        shifts
      )
    )
  }
)


test_that(
  "invalid phylogenetic trees are rejected",
  {

    expect_error(
      validate_phylo_tree(
        list()
      ),
      "inherit"
    )

    unrooted <- ape::unroot(
      make_validation_test_tree()
    )

    expect_error(
      validate_phylo_tree(
        unrooted
      ),
      "rooted"
    )

    polytomy <- ape::read.tree(
      text = "((A:1,B:1,C:1):1,D:1);"
    )

    expect_true(
      ape::is.rooted(polytomy)
    )

    expect_false(
      ape::is.binary(polytomy)
    )

    expect_error(
      validate_phylo_tree(
        polytomy
      ),
      "bifurcating"
    )

    no_lengths <-
      make_validation_test_tree()

    no_lengths$edge.length <- NULL

    expect_error(
      validate_phylo_tree(
        no_lengths
      ),
      "branch lengths"
    )

    zero_length <-
      make_validation_test_tree()

    zero_length$edge.length[[1L]] <- 0

    expect_error(
      validate_phylo_tree(
        zero_length
      ),
      "strictly positive"
    )

    negative_length <-
      make_validation_test_tree()

    negative_length$edge.length[[1L]] <- -1

    expect_error(
      validate_phylo_tree(
        negative_length
      ),
      "strictly positive"
    )

    nonfinite_length <-
      make_validation_test_tree()

    nonfinite_length$edge.length[[1L]] <-
      Inf

    expect_error(
      validate_phylo_tree(
        nonfinite_length
      ),
      "finite"
    )

    duplicated_tip <-
      make_validation_test_tree()

    duplicated_tip$tip.label[[2L]] <-
      duplicated_tip$tip.label[[1L]]

    expect_error(
      validate_phylo_tree(
        duplicated_tip
      ),
      "unique"
    )
  }
)


test_that(
  "invalid terminal-state vectors are rejected",
  {

    tree <- make_validation_test_tree()

    expect_error(
      validate_tip_states(
        tree,
        1:8
      ),
      "named numeric"
    )

    duplicated_names <-
      make_valid_states()

    names(duplicated_names)[[2L]] <-
      names(duplicated_names)[[1L]]

    expect_error(
      validate_tip_states(
        tree,
        duplicated_names
      ),
      "unique"
    )

    missing_state <-
      make_valid_states()

    missing_state <-
      missing_state[
        names(missing_state) != "A"
      ]

    expect_error(
      validate_tip_states(
        tree,
        missing_state
      ),
      "Missing"
    )

    unknown_state <-
      make_valid_states()

    names(unknown_state)[
      names(unknown_state) == "A"
    ] <- "UNKNOWN"

    expect_error(
      validate_tip_states(
        tree,
        unknown_state
      ),
      "Missing|Unknown"
    )

    nonfinite_state <-
      make_valid_states()

    nonfinite_state[["A"]] <-
      NA_real_

    expect_error(
      validate_tip_states(
        tree,
        nonfinite_state
      ),
      "finite"
    )
  }
)


test_that(
  "malformed shift tables are rejected",
  {

    tree <- make_validation_test_tree()

    I1 <- validation_node(
      tree,
      "I1"
    )

    J1 <- validation_node(
      tree,
      "J1"
    )

    expect_error(
      validate_shifts(
        tree,
        c(I1, J1)
      ),
      "data frame"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = I1
        )
      ),
      "parent.*child"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = numeric(0),
          child = numeric(0)
        )
      ),
      "at least one"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = I1 + 0.5,
          child = J1
        )
      ),
      "finite integer"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = 9999,
          child = J1
        )
      ),
      "outside"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = I1,
          child = 9999
        )
      ),
      "outside"
    )
  }
)


test_that(
  "topologically invalid shift edges are rejected",
  {

    tree <- make_validation_test_tree()

    ROOT <- validation_node(
      tree,
      "ROOT"
    )

    H1 <- validation_node(
      tree,
      "H1"
    )

    I1 <- validation_node(
      tree,
      "I1"
    )

    J1 <- validation_node(
      tree,
      "J1"
    )

    A <- validation_node(
      tree,
      "A"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = A,
          child = J1
        )
      ),
      "internal"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = ROOT,
          child = H1
        )
      ),
      "root"
    )

    expect_error(
      validate_shifts(
        tree,
        data.frame(
          parent = H1,
          child = J1
        )
      ),
      "existing directed edge"
    )
  }
)


test_that(
  "duplicated and shared-parent shifts are rejected",
  {

    tree <- make_validation_test_tree()

    I1 <- validation_node(
      tree,
      "I1"
    )

    J1 <- validation_node(
      tree,
      "J1"
    )

    C <- validation_node(
      tree,
      "C"
    )

    duplicate_edges <- data.frame(
      parent = c(I1, I1),
      child = c(J1, J1)
    )

    expect_error(
      validate_shifts(
        tree,
        duplicate_edges
      ),
      "Duplicated"
    )

    shared_parent <- data.frame(
      parent = c(I1, I1),
      child = c(J1, C)
    )

    expect_error(
      validate_shifts(
        tree,
        shared_parent
      ),
      "same parent"
    )
  }
)


test_that(
  "nested or overlapping patch clades are rejected",
  {

    tree <- make_validation_test_tree()

    nested <- data.frame(
      parent = c(
        validation_node(tree, "H1"),
        validation_node(tree, "I1")
      ),
      child = c(
        validation_node(tree, "I1"),
        validation_node(tree, "J1")
      )
    )

    expect_error(
      validate_shifts(
        tree,
        nested
      ),
      "pairwise disjoint"
    )
  }
)


test_that(
  "non-function ASR backends are rejected",
  {

    tree <- make_validation_test_tree()

    shifts <- data.frame(
      parent = validation_node(
        tree,
        "I1"
      ),
      child = validation_node(
        tree,
        "J1"
      )
    )

    expect_error(
      validate_patch_inputs(
        tree = tree,
        states = make_valid_states(),
        shifts = shifts,
        asr_fun = "not_a_function"
      ),
      "must be a function"
    )
  }
)

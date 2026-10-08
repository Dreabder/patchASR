make_patch_asr_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
      "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
    )
  )
}


patch_asr_test_node <- function(
    tree,
    label
) {

  tip_match <- match(
    label,
    tree$tip.label
  )

  if (!is.na(
    tip_match
  )) {
    return(
      as.integer(
        tip_match
      )
    )
  }

  internal_match <- match(
    label,
    tree$node.label
  )

  if (is.na(
    internal_match
  )) {
    stop(
      "Unknown patch_asr test label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


make_patch_asr_states <- function() {

  ## Deliberately scrambled to confirm that patch_asr()
  ## does not depend on the input order of terminal states.
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


fixed_patch_asr_backend <- function(
    tree,
    states
) {

  if (!identical(
    names(states),
    tree$tip.label
  )) {
    stop(
      "Backend received terminal states in the wrong order."
    )
  }

  values <- c(
    node_9 = 10,
    node_10 = 20,
    node_11 = 30,
    node_12 = 40,
    node_13 = 50,
    node_14 = 60,
    node_15 = 70
  )

  ids <- tree$node.label

  if (!all(
    ids %in%
    names(values)
  )) {
    stop(
      "Backend encountered an unknown stable node identifier."
    )
  }

  result <-
    values[
      ids
    ]

  ## Deliberately reverse output order.
  result[
    rev(
      seq_along(result)
    )
  ]
}


test_that(
  "patch_asr performs a complete single-shift reconstruction",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <-
      reference_tree

    tree$node.label <-
      NULL

    shifts <- data.frame(
      parent =
        patch_asr_test_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_asr_test_node(
          reference_tree,
          "J1"
        )
    )

    fit <- patch_asr(
      tree = tree,
      states = make_patch_asr_states(),
      shifts = shifts,
      asr_fun = fixed_patch_asr_backend
    )

    expect_s3_class(
      fit,
      "patch_asr"
    )

    expect_equal(
      nrow(
        fit$ancestral_states
      ),
      tree$Nnode
    )

    expect_identical(
      fit$ancestral_states$node,
      paste0(
        "node_",
        ape::Ntip(tree) +
          seq_len(
            tree$Nnode
          )
      )
    )

    expect_false(
      anyDuplicated(
        fit$ancestral_states$node
      ) > 0L
    )

    expect_true(
      all(
        is.finite(
          fit$ancestral_states$
            estimate
        )
      )
    )

    boundary_row <-
      fit$ancestral_states[
        fit$ancestral_states$node ==
          "node_11",
        ,
        drop = FALSE
      ]

    expect_equal(
      boundary_row$estimate,
      7.25,
      tolerance =
        1e-12
    )

    expect_identical(
      boundary_row$source,
      "boundary_reconstruction"
    )

    expect_identical(
      boundary_row$component,
      "boundary"
    )
  }
)


test_that(
  "patch_asr returns the specified audit-supporting fields",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    fit <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        data.frame(
          parent =
            patch_asr_test_node(
              reference_tree,
              "I1"
            ),
          child =
            patch_asr_test_node(
              reference_tree,
              "J1"
            )
        ),
      asr_fun =
        fixed_patch_asr_backend
    )

    expect_true(
      all(
        c(
          "ancestral_states",
          "patch_plan",
          "boundary_states",
          "component_states",
          "node_map",
          "call",
          "diagnostics"
        ) %in%
          names(fit)
      )
    )

    expect_s3_class(
      fit$patch_plan,
      "patchASR_patch_plan"
    )

    expect_s3_class(
      fit$diagnostics$partition,
      "patchASR_partition"
    )

    expect_s3_class(
      fit$diagnostics$component_asr,
      "patchASR_component_asr"
    )

    expect_s3_class(
      fit$diagnostics$boundary_result,
      "patchASR_boundary_result"
    )

    expect_lt(
      fit$diagnostics$
        boundary_residual,
      1e-10
    )
  }
)


test_that(
  "patch_asr supports terminal-child shifts end to end",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    fit <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        data.frame(
          parent =
            patch_asr_test_node(
              reference_tree,
              "I1"
            ),
          child =
            patch_asr_test_node(
              reference_tree,
              "C"
            )
        ),
      asr_fun =
        fixed_patch_asr_backend
    )

    expect_equal(
      nrow(
        fit$ancestral_states
      ),
      tree$Nnode
    )

    expect_equal(
      fit$patch_plan$
        patches[[1L]]$
        child_type,
      "terminal"
    )

    expect_false(
      fit$patch_plan$
        patches[[1L]]$
        requires_asr
    )

    expect_length(
      fit$diagnostics$
        component_asr$
        patch_states[[1L]],
      0L
    )

    expect_equal(
      nrow(
        fit$boundary_states
      ),
      1L
    )
  }
)


test_that(
  "patch_asr supports adjacent coupled boundaries end to end",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    fit <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        data.frame(
          parent = c(
            patch_asr_test_node(
              reference_tree,
              "I1"
            ),
            patch_asr_test_node(
              reference_tree,
              "H1"
            )
          ),
          child = c(
            patch_asr_test_node(
              reference_tree,
              "J1"
            ),
            patch_asr_test_node(
              reference_tree,
              "D"
            )
          )
        ),
      asr_fun =
        fixed_patch_asr_backend
    )

    expected_boundary_states <- c(
      node_11 = 25 / 6,
      node_10 = 23 / 3
    )

    observed <-
      fit$diagnostics$
      boundary_result$
      states

    expect_equal(
      observed,
      expected_boundary_states,
      tolerance =
        1e-12
    )

    expect_equal(
      nrow(
        fit$ancestral_states
      ),
      tree$Nnode
    )
  }
)


test_that(
  "final numerical reconstruction is invariant to compatible shift-row order",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    shifts_1 <- data.frame(
      parent = c(
        patch_asr_test_node(
          reference_tree,
          "I1"
        ),
        patch_asr_test_node(
          reference_tree,
          "H1"
        )
      ),
      child = c(
        patch_asr_test_node(
          reference_tree,
          "J1"
        ),
        patch_asr_test_node(
          reference_tree,
          "D"
        )
      )
    )

    shifts_2 <-
      shifts_1[
        2:1,
        ,
        drop = FALSE
      ]

    fit_1 <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        shifts_1,
      asr_fun =
        fixed_patch_asr_backend
    )

    fit_2 <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        shifts_2,
      asr_fun =
        fixed_patch_asr_backend
    )

    compare_1 <-
      fit_1$ancestral_states[
        ,
        c(
          "node",
          "estimate",
          "source"
        ),
        drop = FALSE
      ]

    compare_2 <-
      fit_2$ancestral_states[
        ,
        c(
          "node",
          "estimate",
          "source"
        ),
        drop = FALSE
      ]

    expect_equal(
      compare_1,
      compare_2,
      tolerance =
        1e-12
    )
  }
)


test_that(
  "patch_asr does not modify the original user inputs",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    states <-
      make_patch_asr_states()

    shifts <- data.frame(
      parent =
        patch_asr_test_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_asr_test_node(
          reference_tree,
          "J1"
        )
    )

    tree_before <- tree
    states_before <- states
    shifts_before <- shifts

    invisible(
      patch_asr(
        tree = tree,
        states = states,
        shifts = shifts,
        asr_fun =
          fixed_patch_asr_backend
      )
    )

    expect_identical(
      tree,
      tree_before
    )

    expect_identical(
      states,
      states_before
    )

    expect_identical(
      shifts,
      shifts_before
    )
  }
)


test_that(
  "patch_asr rejects invalid residual tolerances",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    shifts <- data.frame(
      parent =
        patch_asr_test_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_asr_test_node(
          reference_tree,
          "J1"
        )
    )

    expect_error(
      patch_asr(
        tree = tree,
        states =
          make_patch_asr_states(),
        shifts =
          shifts,
        asr_fun =
          fixed_patch_asr_backend,
        residual_tolerance =
          0
      ),
      "finite positive"
    )

    expect_error(
      patch_asr(
        tree = tree,
        states =
          make_patch_asr_states(),
        shifts =
          shifts,
        asr_fun =
          fixed_patch_asr_backend,
        residual_tolerance =
          NA_real_
      ),
      "finite positive"
    )
  }
)


test_that(
  "patch_asr public interface rejects invalid shift configurations",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    nested_shifts <- data.frame(
      parent = c(
        patch_asr_test_node(
          reference_tree,
          "H1"
        ),
        patch_asr_test_node(
          reference_tree,
          "I1"
        )
      ),
      child = c(
        patch_asr_test_node(
          reference_tree,
          "I1"
        ),
        patch_asr_test_node(
          reference_tree,
          "J1"
        )
      )
    )

    expect_error(
      patch_asr(
        tree = tree,
        states =
          make_patch_asr_states(),
        shifts =
          nested_shifts,
        asr_fun =
          fixed_patch_asr_backend
      ),
      "pairwise disjoint"
    )
  }
)


test_that(
  "print method returns patch_asr object invisibly",
  {

    reference_tree <-
      make_patch_asr_test_tree()

    tree <- reference_tree
    tree$node.label <- NULL

    fit <- patch_asr(
      tree = tree,
      states =
        make_patch_asr_states(),
      shifts =
        data.frame(
          parent =
            patch_asr_test_node(
              reference_tree,
              "I1"
            ),
          child =
            patch_asr_test_node(
              reference_tree,
              "J1"
            )
        ),
      asr_fun =
        fixed_patch_asr_backend
    )

    expect_invisible(
      print(fit)
    )
  }
)

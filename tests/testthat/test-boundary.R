make_boundary_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
      "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
    )
  )
}


boundary_node <- function(
    tree,
    label
) {

  tip_match <- match(
    label,
    tree$tip.label
  )

  if (!is.na(tip_match)) {
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
      "Unknown boundary-test label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


make_boundary_states <- function() {

  c(
    A = 1,
    B = 2,
    C = 3,
    D = 4,
    E = 5,
    F = 6,
    G = 7,
    H = 8
  )
}


fixed_boundary_asr <- function(
    tree,
    states
) {

  values <- c(
    node_9 = 10,
    node_10 = 20,
    node_11 = 30,
    node_12 = 40,
    node_13 = 50,
    node_14 = 60,
    node_15 = 70
  )

  ids <-
    tree$node.label

  if (!all(
    ids %in%
    names(values)
  )) {
    stop(
      "Mock boundary backend encountered an unknown stable ID."
    )
  }

  result <-
    values[
      ids
    ]

  ## Return in reverse order so the runner must standardize it.
  result[
    rev(
      seq_along(result)
    )
  ]
}


prepare_boundary_fixture <- function(
    shifts
) {

  reference_tree <-
    make_boundary_test_tree()

  user_tree <-
    reference_tree

  user_tree$node.label <-
    NULL

  frozen <-
    freeze_node_identity(
      user_tree
    )

  shift_table <-
    shifts(
      reference_tree
    )

  plan <-
    build_patch_plan(
      tree =
        user_tree,
      node_map =
        frozen$node_map,
      shifts =
        shift_table
    )

  partition <-
    partition_phylogeny(
      tree =
        user_tree,
      node_map =
        frozen$node_map,
      patch_plan =
        plan
    )

  component_asr <-
    run_partition_asr(
      partition =
        partition,
      states =
        make_boundary_states(),
      asr_fun =
        fixed_boundary_asr
    )

  list(
    reference_tree =
      reference_tree,
    user_tree =
      user_tree,
    node_map =
      frozen$node_map,
    plan =
      plan,
    partition =
      partition,
    component_asr =
      component_asr
  )
}


test_that(
  "unified known-state pool contains tips and non-boundary internal states",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              boundary_node(
                tree,
                "I1"
              ),
            child =
              boundary_node(
                tree,
                "J1"
              )
          )
        }
      )

    known <-
      build_unified_known_states(
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr,
        boundary_ids =
          fixture$plan$boundary_ids
      )

    expect_false(
      any(
        fixture$plan$boundary_ids %in%
          names(known)
      )
    )

    expect_equal(
      length(known),
      nrow(
        fixture$node_map
      ) -
        length(
          fixture$plan$boundary_ids
        )
    )

    expect_true(
      all(
        fixture$component_asr$
          all_internal_states %in%
          known
      )
    )
  }
)


test_that(
  "single boundary matches direct branch-length interpolation",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              boundary_node(
                tree,
                "I1"
              ),
            child =
              boundary_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <-
      run_boundary_reconstruction(
        patch_plan =
          fixture$plan,
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr
      )

    expect_s3_class(
      result,
      "phyloPatch_boundary_result"
    )

    ## H1 = 20
    ## C = 3
    ##
    ## H1 --3--> I1 --1--> C
    ##
    ## I1 = 0.25*20 + 0.75*3 = 7.25

    expect_equal(
      unname(
        result$states
      ),
      7.25,
      tolerance =
        1e-12
    )

    expect_lt(
      result$residual,
      1e-10
    )
  }
)


test_that(
  "adjacent boundary nodes are solved jointly",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent = c(
              boundary_node(
                tree,
                "I1"
              ),
              boundary_node(
                tree,
                "H1"
              )
            ),
            child = c(
              boundary_node(
                tree,
                "J1"
              ),
              boundary_node(
                tree,
                "D"
              )
            )
          )
        }
      )

    result <-
      run_boundary_reconstruction(
        patch_plan =
          fixture$plan,
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr
      )

    expected <- c(
      node_11 = 25 / 6,
      node_10 = 23 / 3
    )

    expect_equal(
      result$states,
      expected,
      tolerance =
        1e-12
    )

    expect_equal(
      result$system$A,
      matrix(
        c(
          1,
          -0.4,
          -0.25,
          1
        ),
        nrow = 2L,
        byrow = FALSE,
        dimnames = list(
          c(
            "node_11",
            "node_10"
          ),
          c(
            "node_11",
            "node_10"
          )
        )
      ),
      tolerance =
        1e-12
    )

    expect_lt(
      result$residual,
      1e-10
    )
  }
)


test_that(
  "mixed internal and terminal shift boundaries are reconstructed correctly",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent = c(
              boundary_node(
                tree,
                "I1"
              ),
              boundary_node(
                tree,
                "J2"
              )
            ),
            child = c(
              boundary_node(
                tree,
                "J1"
              ),
              boundary_node(
                tree,
                "E"
              )
            )
          )
        }
      )

    result <-
      run_boundary_reconstruction(
        patch_plan =
          fixture$plan,
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr
      )

    expected <- c(
      node_11 = 7.25,
      node_14 = 32.4
    )

    expect_equal(
      result$states,
      expected,
      tolerance =
        1e-12
    )
  }
)


test_that(
  "biological boundary solution is invariant to shift-row order",
  {

    make_shifts_1 <- function(
    tree
    ) {

      data.frame(
        parent = c(
          boundary_node(
            tree,
            "I1"
          ),
          boundary_node(
            tree,
            "H1"
          )
        ),
        child = c(
          boundary_node(
            tree,
            "J1"
          ),
          boundary_node(
            tree,
            "D"
          )
        )
      )
    }

    make_shifts_2 <- function(
    tree
    ) {

      make_shifts_1(
        tree
      )[
        2:1,
        ,
        drop = FALSE
      ]
    }

    fixture_1 <-
      prepare_boundary_fixture(
        make_shifts_1
      )

    fixture_2 <-
      prepare_boundary_fixture(
        make_shifts_2
      )

    result_1 <-
      run_boundary_reconstruction(
        patch_plan =
          fixture_1$plan,
        node_map =
          fixture_1$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture_1$component_asr
      )

    result_2 <-
      run_boundary_reconstruction(
        patch_plan =
          fixture_2$plan,
        node_map =
          fixture_2$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture_2$component_asr
      )

    order_states <- function(
    x
    ) {

      x[
        order(
          names(x)
        )
      ]
    }

    expect_equal(
      order_states(
        result_1$states
      ),
      order_states(
        result_2$states
      ),
      tolerance =
        1e-12
    )
  }
)


test_that(
  "missing required known states are rejected",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              boundary_node(
                tree,
                "I1"
              ),
            child =
              boundary_node(
                tree,
                "J1"
              )
          )
        }
      )

    known <-
      build_unified_known_states(
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr,
        boundary_ids =
          fixture$plan$boundary_ids
      )

    ## tip_3 is C, required as the sibling-side state.
    known_without_C <-
      known[
        names(known) !=
          "tip_3"
      ]

    expect_error(
      build_boundary_system(
        patch_plan =
          fixture$plan,
        known_states =
          known_without_C
      ),
      "Required known state"
    )
  }
)


test_that(
  "unresolved boundaries cannot be inserted into known-state pool",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              boundary_node(
                tree,
                "I1"
              ),
            child =
              boundary_node(
                tree,
                "J1"
              )
          )
        }
      )

    known <-
      build_unified_known_states(
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr,
        boundary_ids =
          fixture$plan$boundary_ids
      )

    bad_known <- c(
      known,
      setNames(
        999,
        fixture$plan$
          boundary_ids[[1L]]
      )
    )

    expect_error(
      build_boundary_system(
        patch_plan =
          fixture$plan,
        known_states =
          bad_known
      ),
      "unresolved boundary"
    )
  }
)


test_that(
  "singular boundary systems are rejected",
  {

    system <- structure(
      list(
        A = matrix(
          c(
            1,
            -1,
            -1,
            1
          ),
          nrow = 2L,
          byrow = TRUE,
          dimnames = list(
            c(
              "node_a",
              "node_b"
            ),
            c(
              "node_a",
              "node_b"
            )
          )
        ),

        b = c(
          node_a = 1,
          node_b = 1
        ),

        boundary_ids = c(
          "node_a",
          "node_b"
        ),

        details =
          data.frame()
      ),
      class =
        "phyloPatch_boundary_system"
    )

    expect_error(
      solve_boundary_system(
        system
      ),
      "singular"
    )
  }
)


test_that(
  "boundary estimate table records provenance",
  {

    fixture <-
      prepare_boundary_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              boundary_node(
                tree,
                "I1"
              ),
            child =
              boundary_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <-
      run_boundary_reconstruction(
        patch_plan =
          fixture$plan,
        node_map =
          fixture$node_map,
        states =
          make_boundary_states(),
        component_asr =
          fixture$component_asr
      )

    expect_true(
      all(
        result$estimate_table$
          source ==
          "boundary_reconstruction"
      )
    )

    expect_true(
      all(
        result$estimate_table$
          component ==
          "boundary"
      )
    )

    expect_equal(
      result$estimate_table$node,
      names(
        result$states
      )
    )
  }
)

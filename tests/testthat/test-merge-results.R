make_merge_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
      "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
    )
  )
}


merge_node <- function(
    tree,
    label
) {

  tip_match <- match(
    label,
    tree$tip.label
  )

  if (!is.na(tip_match)) {
    return(
      as.integer(tip_match)
    )
  }

  internal_match <- match(
    label,
    tree$node.label
  )

  if (is.na(internal_match)) {
    stop(
      "Unknown merge-test label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


make_merge_states <- function() {

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


fixed_merge_asr <- function(
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

  ids <- tree$node.label

  values[
    ids
  ]
}


prepare_merge_fixture <- function(
    shifts
) {

  reference_tree <-
    make_merge_test_tree()

  user_tree <-
    reference_tree

  user_tree$node.label <-
    NULL

  frozen <-
    freeze_node_identity(
      user_tree
    )

  plan <-
    build_patch_plan(
      tree =
        user_tree,
      node_map =
        frozen$node_map,
      shifts =
        shifts(
          reference_tree
        )
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
        make_merge_states(),
      asr_fun =
        fixed_merge_asr
    )

  boundary_result <-
    run_boundary_reconstruction(
      patch_plan =
        plan,
      node_map =
        frozen$node_map,
      states =
        make_merge_states(),
      component_asr =
        component_asr
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
      component_asr,
    boundary_result =
      boundary_result
  )
}


test_that(
  "single-shift merge contains every original internal node exactly once",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <-
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          fixture$component_asr,
        boundary_result =
          fixture$boundary_result
      )

    expected_ids <-
      fixture$node_map$stable_id[
        fixture$node_map$node_type ==
          "internal"
      ]

    expect_equal(
      nrow(result),
      length(
        expected_ids
      )
    )

    expect_identical(
      result$node,
      expected_ids
    )

    expect_false(
      anyDuplicated(
        result$node
      ) > 0L
    )

    expect_true(
      all(
        is.finite(
          result$estimate
        )
      )
    )
  }
)


test_that(
  "single-shift merge preserves component and boundary provenance",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <-
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          fixture$component_asr,
        boundary_result =
          fixture$boundary_result
      )

    boundary_row <-
      result[
        result$node ==
          "node_11",
        ,
        drop = FALSE
      ]

    expect_equal(
      nrow(
        boundary_row
      ),
      1L
    )

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

    patch_row <-
      result[
        result$node ==
          "node_12",
        ,
        drop = FALSE
      ]

    expect_equal(
      patch_row$estimate,
      40
    )

    expect_identical(
      patch_row$source,
      "component_asr"
    )

    expect_identical(
      patch_row$component,
      "patch_1"
    )
  }
)


test_that(
  "mixed shifts still merge to the complete original internal-node set",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent = c(
              merge_node(
                tree,
                "I1"
              ),
              merge_node(
                tree,
                "J2"
              )
            ),
            child = c(
              merge_node(
                tree,
                "J1"
              ),
              merge_node(
                tree,
                "E"
              )
            )
          )
        }
      )

    result <-
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          fixture$component_asr,
        boundary_result =
          fixture$boundary_result
      )

    expected_ids <-
      fixture$node_map$stable_id[
        fixture$node_map$node_type ==
          "internal"
      ]

    expect_identical(
      result$node,
      expected_ids
    )

    boundary_rows <-
      result[
        result$source ==
          "boundary_reconstruction",
        ,
        drop = FALSE
      ]

    expect_equal(
      sort(
        boundary_rows$node
      ),
      sort(
        c(
          "node_11",
          "node_14"
        )
      )
    )

    expect_equal(
      boundary_rows$estimate[
        boundary_rows$node ==
          "node_11"
      ],
      7.25,
      tolerance =
        1e-12
    )

    expect_equal(
      boundary_rows$estimate[
        boundary_rows$node ==
          "node_14"
      ],
      32.4,
      tolerance =
        1e-12
    )
  }
)


test_that(
  "missing component estimate-table rows are rejected",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    bad_component <-
      fixture$component_asr

    bad_component$estimate_table <-
      bad_component$estimate_table[
        -1L,
        ,
        drop = FALSE
      ]

    expect_error(
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          bad_component,
        boundary_result =
          fixture$boundary_result
      ),
      "does not match"
    )
  }
)


test_that(
  "duplicated merged node identities are rejected",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    bad_boundary <-
      fixture$boundary_result

    bad_boundary$estimate_table$node <-
      fixture$component_asr$
      estimate_table$node[[1L]]

    names(
      bad_boundary$states
    ) <-
      fixture$component_asr$
      estimate_table$node[[1L]]

    expect_error(
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          fixture$component_asr,
        boundary_result =
          bad_boundary
      ),
      "occurs more than once|incorrect internal-node set"
    )
  }
)


test_that(
  "non-finite merged estimates are rejected",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    bad_component <-
      fixture$component_asr

    bad_component$estimate_table$
      estimate[[1L]] <-
      NA_real_

    expect_error(
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          bad_component,
        boundary_result =
          fixture$boundary_result
      ),
      "non-finite"
    )
  }
)


test_that(
  "inconsistent estimate-table values are rejected",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    bad_component <-
      fixture$component_asr

    bad_component$estimate_table$
      estimate[[1L]] <-
      bad_component$estimate_table$
      estimate[[1L]] +
      100

    expect_error(
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          bad_component,
        boundary_result =
          fixture$boundary_result
      ),
      "inconsistent"
    )
  }
)


test_that(
  "invalid provenance is rejected",
  {

    fixture <-
      prepare_merge_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              merge_node(
                tree,
                "I1"
              ),
            child =
              merge_node(
                tree,
                "J1"
              )
          )
        }
      )

    bad_boundary <-
      fixture$boundary_result

    bad_boundary$estimate_table$
      source <-
      "component_asr"

    expect_error(
      merge_ancestral_results(
        node_map =
          fixture$node_map,
        component_asr =
          fixture$component_asr,
        boundary_result =
          bad_boundary
      ),
      "provenance"
    )
  }
)

make_partition_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
      "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
    )
  )
}


partition_node <- function(
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
      "Unknown partition-test label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


prepare_partition_fixture <- function() {

  reference_tree <-
    make_partition_test_tree()

  user_tree <-
    reference_tree

  user_tree$node.label <-
    NULL

  frozen <-
    freeze_node_identity(
      user_tree
    )

  list(
    reference_tree =
      reference_tree,
    user_tree =
      user_tree,
    node_map =
      frozen$node_map
  )
}


partition_sid <- function(
    fixture,
    label
) {

  node <- partition_node(
    fixture$reference_tree,
    label
  )

  stable_id_from_original_node(
    fixture$node_map,
    node
  )
}


build_partition_for_shifts <- function(
    fixture,
    shifts
) {

  plan <- build_patch_plan(
    tree =
      fixture$user_tree,
    node_map =
      fixture$node_map,
    shifts =
      shifts
  )

  partition_phylogeny(
    tree =
      fixture$user_tree,
    node_map =
      fixture$node_map,
    patch_plan =
      plan
  )
}


test_that(
  "single internal patch partitions mother, patch, and boundary correctly",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent =
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
      child =
        partition_node(
          fixture$reference_tree,
          "J1"
        )
    )

    partition <-
      build_partition_for_shifts(
        fixture,
        shifts
      )

    expect_s3_class(
      partition,
      "patchASR_partition"
    )

    expect_equal(
      sort(
        partition$mother_tree$tip.label
      ),
      c(
        "C",
        "D",
        "E",
        "F",
        "G",
        "H"
      )
    )

    expect_equal(
      sort(
        partition$boundary_ids
      ),
      partition_sid(
        fixture,
        "I1"
      )
    )

    expect_equal(
      sort(
        partition$collapsed_mother_side_ids
      ),
      sort(
        partition$boundary_ids
      )
    )

    expected_mother_internal <- c(
      partition_sid(fixture, "ROOT"),
      partition_sid(fixture, "H1"),
      partition_sid(fixture, "I2"),
      partition_sid(fixture, "J2"),
      partition_sid(fixture, "S2")
    )

    expect_equal(
      sort(
        partition$mother_internal_ids
      ),
      sort(
        expected_mother_internal
      )
    )

    expect_equal(
      partition$patch_internal_ids,
      partition_sid(
        fixture,
        "J1"
      )
    )

    patch <-
      partition$patch_components[[1L]]

    expect_true(
      patch$requires_asr
    )

    expect_s3_class(
      patch$tree,
      "phylo"
    )

    expect_equal(
      sort(
        patch$tree$tip.label
      ),
      c("A", "B")
    )
  }
)


test_that(
  "single terminal patch contains no patch internal nodes",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent =
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
      child =
        partition_node(
          fixture$reference_tree,
          "C"
        )
    )

    partition <-
      build_partition_for_shifts(
        fixture,
        shifts
      )

    expect_equal(
      sort(
        partition$mother_tree$tip.label
      ),
      c(
        "A",
        "B",
        "D",
        "E",
        "F",
        "G",
        "H"
      )
    )

    expect_length(
      partition$patch_internal_ids,
      0L
    )

    expect_equal(
      partition$boundary_ids,
      partition_sid(
        fixture,
        "I1"
      )
    )

    expect_equal(
      partition$collapsed_mother_side_ids,
      partition$boundary_ids
    )

    ## J1 remains in the mother tree for the terminal shift I1 -> C.
    expect_true(
      partition_sid(
        fixture,
        "J1"
      ) %in%
        partition$mother_internal_ids
    )

    patch <-
      partition$patch_components[[1L]]

    expect_false(
      patch$requires_asr
    )

    expect_null(
      patch$tree
    )

    expect_equal(
      nrow(
        patch$component_map
      ),
      1L
    )

    expect_identical(
      patch$component_map$node_type,
      "tip"
    )

    expect_identical(
      patch$component_map$original_label,
      "C"
    )
  }
)


test_that(
  "two distant internal patches partition correctly",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent = c(
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
        partition_node(
          fixture$reference_tree,
          "I2"
        )
      ),
      child = c(
        partition_node(
          fixture$reference_tree,
          "J1"
        ),
        partition_node(
          fixture$reference_tree,
          "J2"
        )
      )
    )

    partition <-
      build_partition_for_shifts(
        fixture,
        shifts
      )

    expect_equal(
      sort(
        partition$mother_tree$tip.label
      ),
      c(
        "C",
        "D",
        "G",
        "H"
      )
    )

    expect_equal(
      sort(
        partition$boundary_ids
      ),
      sort(
        c(
          partition_sid(
            fixture,
            "I1"
          ),
          partition_sid(
            fixture,
            "I2"
          )
        )
      )
    )

    expect_equal(
      sort(
        partition$patch_internal_ids
      ),
      sort(
        c(
          partition_sid(
            fixture,
            "J1"
          ),
          partition_sid(
            fixture,
            "J2"
          )
        )
      )
    )

    expect_equal(
      sort(
        partition$collapsed_mother_side_ids
      ),
      sort(
        partition$boundary_ids
      )
    )

    expected_mother_internal <- c(
      partition_sid(
        fixture,
        "ROOT"
      ),
      partition_sid(
        fixture,
        "H1"
      ),
      partition_sid(
        fixture,
        "S2"
      )
    )

    expect_equal(
      sort(
        partition$mother_internal_ids
      ),
      sort(
        expected_mother_internal
      )
    )
  }
)


test_that(
  "adjacent disjoint shifts collapse exactly their boundary parents",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent = c(
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
        partition_node(
          fixture$reference_tree,
          "H1"
        )
      ),
      child = c(
        partition_node(
          fixture$reference_tree,
          "J1"
        ),
        partition_node(
          fixture$reference_tree,
          "D"
        )
      )
    )

    partition <-
      build_partition_for_shifts(
        fixture,
        shifts
      )

    expected_boundaries <- c(
      partition_sid(
        fixture,
        "I1"
      ),
      partition_sid(
        fixture,
        "H1"
      )
    )

    expect_equal(
      sort(
        partition$boundary_ids
      ),
      sort(
        expected_boundaries
      )
    )

    expect_equal(
      sort(
        partition$collapsed_mother_side_ids
      ),
      sort(
        expected_boundaries
      )
    )

    expect_equal(
      sort(
        partition$mother_tree$tip.label
      ),
      c(
        "C",
        "E",
        "F",
        "G",
        "H"
      )
    )
  }
)


test_that(
  "mixed internal and terminal patches partition correctly",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent = c(
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
        partition_node(
          fixture$reference_tree,
          "J2"
        )
      ),
      child = c(
        partition_node(
          fixture$reference_tree,
          "J1"
        ),
        partition_node(
          fixture$reference_tree,
          "E"
        )
      )
    )

    partition <-
      build_partition_for_shifts(
        fixture,
        shifts
      )

    expect_equal(
      sort(
        partition$mother_tree$tip.label
      ),
      c(
        "C",
        "D",
        "F",
        "G",
        "H"
      )
    )

    expected_boundaries <- c(
      partition_sid(
        fixture,
        "I1"
      ),
      partition_sid(
        fixture,
        "J2"
      )
    )

    expect_equal(
      sort(
        partition$boundary_ids
      ),
      sort(
        expected_boundaries
      )
    )

    expect_equal(
      sort(
        partition$collapsed_mother_side_ids
      ),
      sort(
        expected_boundaries
      )
    )

    expect_equal(
      partition$patch_components[[1L]]$
        requires_asr,
      TRUE
    )

    expect_equal(
      partition$patch_components[[2L]]$
        requires_asr,
      FALSE
    )
  }
)


test_that(
  "partitioning does not modify original inputs",
  {

    fixture <-
      prepare_partition_fixture()

    tree_before <-
      fixture$user_tree

    map_before <-
      fixture$node_map

    shifts <- data.frame(
      parent =
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
      child =
        partition_node(
          fixture$reference_tree,
          "J1"
        )
    )

    plan <- build_patch_plan(
      tree =
        fixture$user_tree,
      node_map =
        fixture$node_map,
      shifts =
        shifts
    )

    plan_before <-
      plan

    invisible(
      partition_phylogeny(
        tree =
          fixture$user_tree,
        node_map =
          fixture$node_map,
        patch_plan =
          plan
      )
    )

    expect_identical(
      fixture$user_tree,
      tree_before
    )

    expect_identical(
      fixture$node_map,
      map_before
    )

    expect_identical(
      plan,
      plan_before
    )
  }
)


test_that(
  "biological partition is invariant to shift-row order",
  {

    fixture <-
      prepare_partition_fixture()

    shifts_1 <- data.frame(
      parent = c(
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
        partition_node(
          fixture$reference_tree,
          "I2"
        )
      ),
      child = c(
        partition_node(
          fixture$reference_tree,
          "J1"
        ),
        partition_node(
          fixture$reference_tree,
          "J2"
        )
      )
    )

    shifts_2 <-
      shifts_1[
        2:1,
        ,
        drop = FALSE
      ]

    partition_1 <-
      build_partition_for_shifts(
        fixture,
        shifts_1
      )

    partition_2 <-
      build_partition_for_shifts(
        fixture,
        shifts_2
      )

    expect_equal(
      sort(
        partition_1$mother_internal_ids
      ),
      sort(
        partition_2$mother_internal_ids
      )
    )

    expect_equal(
      sort(
        partition_1$patch_internal_ids
      ),
      sort(
        partition_2$patch_internal_ids
      )
    )

    expect_equal(
      sort(
        partition_1$boundary_ids
      ),
      sort(
        partition_2$boundary_ids
      )
    )

    expect_equal(
      sort(
        partition_1$mother_tip_ids
      ),
      sort(
        partition_2$mother_tip_ids
      )
    )

    expect_equal(
      sort(
        partition_1$patch_tip_ids
      ),
      sort(
        partition_2$patch_tip_ids
      )
    )
  }
)


test_that(
  "tampered boundary assignments are rejected",
  {

    fixture <-
      prepare_partition_fixture()

    shifts <- data.frame(
      parent =
        partition_node(
          fixture$reference_tree,
          "I1"
        ),
      child =
        partition_node(
          fixture$reference_tree,
          "J1"
        )
    )

    plan <- build_patch_plan(
      tree =
        fixture$user_tree,
      node_map =
        fixture$node_map,
      shifts =
        shifts
    )

    ## Deliberately replace the true boundary I1 with H1.
    plan$boundary_ids <-
      partition_sid(
        fixture,
        "H1"
      )

    expect_error(
      partition_phylogeny(
        tree =
          fixture$user_tree,
        node_map =
          fixture$node_map,
        patch_plan =
          plan
      ),
      "pruning-induced collapse"
    )
  }
)

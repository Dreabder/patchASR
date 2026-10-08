make_patch_plan_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
      "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
    )
  )
}


patch_plan_node <- function(
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
      "Unknown test-tree label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


prepare_patch_plan_fixture <- function() {

  reference_tree <-
    make_patch_plan_test_tree()

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


test_that(
  "internal-child patch plans contain the expected topology",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    shifts <- data.frame(
      parent =
        patch_plan_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_plan_node(
          reference_tree,
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

    expect_s3_class(
      plan,
      "patchASR_patch_plan"
    )

    expect_equal(
      plan$n_patches,
      1L
    )

    patch <- plan$patches[[1L]]

    expect_equal(
      patch$parent_node,
      patch_plan_node(
        reference_tree,
        "I1"
      )
    )

    expect_equal(
      patch$child_node,
      patch_plan_node(
        reference_tree,
        "J1"
      )
    )

    expect_equal(
      patch$upstream_node,
      patch_plan_node(
        reference_tree,
        "H1"
      )
    )

    expect_equal(
      patch$sibling_node,
      patch_plan_node(
        reference_tree,
        "C"
      )
    )

    expect_identical(
      patch$child_type,
      "internal"
    )

    expect_true(
      patch$requires_asr
    )

    expect_equal(
      sort(
        patch$descendant_tip_labels
      ),
      c("A", "B")
    )

    expect_equal(
      patch$descendant_internal_nodes,
      patch_plan_node(
        reference_tree,
        "J1"
      )
    )
  }
)


test_that(
  "terminal-child patch plans contain no internal patch nodes",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    shifts <- data.frame(
      parent =
        patch_plan_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_plan_node(
          reference_tree,
          "C"
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

    patch <- plan$patches[[1L]]

    expect_identical(
      patch$child_type,
      "terminal"
    )

    expect_false(
      patch$requires_asr
    )

    expect_equal(
      patch$descendant_tip_labels,
      "C"
    )

    expect_length(
      patch$descendant_internal_nodes,
      0L
    )

    expect_length(
      patch$descendant_internal_ids,
      0L
    )

    expect_equal(
      patch$sibling_node,
      patch_plan_node(
        reference_tree,
        "J1"
      )
    )
  }
)


test_that(
  "multiple disjoint shifts produce separate patch entries",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    shifts <- data.frame(
      parent = c(
        patch_plan_node(
          reference_tree,
          "I1"
        ),
        patch_plan_node(
          reference_tree,
          "I2"
        )
      ),
      child = c(
        patch_plan_node(
          reference_tree,
          "J1"
        ),
        patch_plan_node(
          reference_tree,
          "J2"
        )
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

    expect_equal(
      plan$n_patches,
      2L
    )

    expect_length(
      plan$boundary_ids,
      2L
    )

    expect_false(
      anyDuplicated(
        plan$boundary_ids
      ) > 0L
    )

    expect_equal(
      sort(
        plan$patch_tip_labels
      ),
      c(
        "A",
        "B",
        "E",
        "F"
      )
    )

    expect_false(
      anyDuplicated(
        c(
          plan$patches[[1L]]$
            descendant_tip_ids,
          plan$patches[[2L]]$
            descendant_tip_ids,
          plan$patches[[1L]]$
            descendant_internal_ids,
          plan$patches[[2L]]$
            descendant_internal_ids
        )
      ) > 0L
    )
  }
)


test_that(
  "patch plans retain original boundary branch lengths",
  {

    reference_tree <- ape::read.tree(
      text = paste0(
        "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
        "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
      )
    )

    I1 <- patch_plan_node(
      reference_tree,
      "I1"
    )

    J1 <- patch_plan_node(
      reference_tree,
      "J1"
    )

    user_tree <- reference_tree
    user_tree$node.label <- NULL

    frozen <- freeze_node_identity(
      user_tree
    )

    plan <- build_patch_plan(
      tree = user_tree,
      node_map = frozen$node_map,
      shifts = data.frame(
        parent = I1,
        child = J1
      )
    )

    patch <- plan$patches[[1L]]

    expect_equal(
      patch$branch_length_parent_boundary,
      3
    )

    expect_equal(
      patch$branch_length_boundary_sibling,
      1
    )

    expect_equal(
      patch$branch_length_shift_edge,
      4
    )
  }
)


test_that(
  "terminal shifts retain the correct sibling-side branch length",
  {

    reference_tree <- ape::read.tree(
      text = paste0(
        "((((A:1,B:1)J1:4,C:1)I1:3,D:5)H1:2,",
        "((E:2,F:3)J2:2,(G:1,H:1)S2:1)I2:2)ROOT;"
      )
    )

    I1 <- patch_plan_node(
      reference_tree,
      "I1"
    )

    C_tip <- patch_plan_node(
      reference_tree,
      "C"
    )

    user_tree <- reference_tree
    user_tree$node.label <- NULL

    frozen <- freeze_node_identity(
      user_tree
    )

    plan <- build_patch_plan(
      tree = user_tree,
      node_map = frozen$node_map,
      shifts = data.frame(
        parent = I1,
        child = C_tip
      )
    )

    patch <- plan$patches[[1L]]

    ## For I1 -> C, the sibling side is J1.
    ## The original I1 -> J1 branch length is 4.
    expect_equal(
      patch$branch_length_parent_boundary,
      3
    )

    expect_equal(
      patch$branch_length_boundary_sibling,
      4
    )

    expect_equal(
      patch$branch_length_shift_edge,
      1
    )
  }
)


test_that(
  "patch-plan construction does not modify its inputs",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    tree_before <-
      fixture$user_tree

    node_map_before <-
      fixture$node_map

    shifts <- data.frame(
      parent =
        patch_plan_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_plan_node(
          reference_tree,
          "J1"
        )
    )

    shifts_before <- shifts

    invisible(
      build_patch_plan(
        tree =
          fixture$user_tree,
        node_map =
          fixture$node_map,
        shifts =
          shifts
      )
    )

    expect_identical(
      fixture$user_tree,
      tree_before
    )

    expect_identical(
      fixture$node_map,
      node_map_before
    )

    expect_identical(
      shifts,
      shifts_before
    )
  }
)


test_that(
  "biological patch content is invariant to shift-row order",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    shifts_1 <- data.frame(
      parent = c(
        patch_plan_node(
          reference_tree,
          "I1"
        ),
        patch_plan_node(
          reference_tree,
          "I2"
        )
      ),
      child = c(
        patch_plan_node(
          reference_tree,
          "J1"
        ),
        patch_plan_node(
          reference_tree,
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

    plan_1 <- build_patch_plan(
      tree =
        fixture$user_tree,
      node_map =
        fixture$node_map,
      shifts =
        shifts_1
    )

    plan_2 <- build_patch_plan(
      tree =
        fixture$user_tree,
      node_map =
        fixture$node_map,
      shifts =
        shifts_2
    )

    patch_signature <- function(
    plan
    ) {

      signatures <- vapply(
        plan$patches,
        function(x) {

          paste(
            x$parent_id,
            x$child_id,
            x$child_type,
            paste(
              sort(
                x$descendant_tip_ids
              ),
              collapse = ","
            ),
            paste(
              sort(
                x$descendant_internal_ids
              ),
              collapse = ","
            ),
            sep = "|"
          )
        },
        character(1)
      )

      sort(signatures)
    }

    expect_identical(
      patch_signature(plan_1),
      patch_signature(plan_2)
    )

    expect_identical(
      sort(plan_1$boundary_ids),
      sort(plan_2$boundary_ids)
    )

    expect_identical(
      sort(plan_1$patch_tip_ids),
      sort(plan_2$patch_tip_ids)
    )
  }
)


test_that(
  "invalid node maps are rejected",
  {

    fixture <-
      prepare_patch_plan_fixture()

    reference_tree <-
      fixture$reference_tree

    shifts <- data.frame(
      parent =
        patch_plan_node(
          reference_tree,
          "I1"
        ),
      child =
        patch_plan_node(
          reference_tree,
          "J1"
        )
    )

    broken_map <-
      fixture$node_map[-1L, ]

    expect_error(
      build_patch_plan(
        tree =
          fixture$user_tree,
        node_map =
          broken_map,
        shifts =
          shifts
      ),
      "does not correspond"
    )
  }
)

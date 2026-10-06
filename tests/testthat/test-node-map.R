make_node_map_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
      "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
    )
  )
}


test_that(
  "stable identity freezing does not modify the user tree",
  {

    tree <- make_node_map_test_tree()

    tree$node.label <- NULL

    original_tree <- tree

    frozen <- freeze_node_identity(
      tree
    )

    expect_identical(
      tree,
      original_tree
    )

    expect_null(
      tree$node.label
    )

    expect_equal(
      nrow(frozen$node_map),
      ape::Ntip(tree) + tree$Nnode
    )

    expect_false(
      anyDuplicated(
        frozen$node_map$stable_id
      ) > 0L
    )

    expect_equal(
      frozen$tree$node.label,
      frozen$node_map$stable_id[
        frozen$node_map$node_type ==
          "internal"
      ]
    )
  }
)


test_that(
  "existing internal labels are preserved in the node map",
  {

    tree <- make_node_map_test_tree()

    original_labels <-
      tree$node.label

    frozen <- freeze_node_identity(
      tree
    )

    internal_map <-
      frozen$node_map[
        frozen$node_map$node_type ==
          "internal",
        ,
        drop = FALSE
      ]

    expect_equal(
      internal_map$original_label,
      original_labels
    )

    expect_equal(
      frozen$tree$node.label,
      internal_map$stable_id
    )

    expect_equal(
      tree$node.label,
      original_labels
    )
  }
)


test_that(
  "node identity survives phylo reordering",
  {

    tree <- make_node_map_test_tree()
    tree$node.label <- NULL

    frozen <- freeze_node_identity(
      tree
    )

    reordered <- ape::reorder.phylo(
      frozen$tree,
      order = "postorder"
    )

    mapped <- map_component_nodes(
      reordered,
      frozen$node_map
    )

    expect_equal(
      sort(mapped$stable_id),
      sort(frozen$node_map$stable_id)
    )
  }
)


test_that(
  "node identity survives patch extraction",
  {

    reference_tree <-
      make_node_map_test_tree()

    j1_node <-
      ape::Ntip(reference_tree) +
      match(
        "J1",
        reference_tree$node.label
      )

    user_tree <- reference_tree
    user_tree$node.label <- NULL

    frozen <- freeze_node_identity(
      user_tree
    )

    patch <- ape::extract.clade(
      frozen$tree,
      node = j1_node
    )

    mapped <- map_component_nodes(
      patch,
      frozen$node_map
    )

    expected_original_nodes <- c(
      match("A", reference_tree$tip.label),
      match("B", reference_tree$tip.label),
      j1_node
    )

    expected_stable_ids <-
      stable_id_from_original_node(
        frozen$node_map,
        expected_original_nodes
      )

    expect_equal(
      sort(mapped$stable_id),
      sort(expected_stable_ids)
    )

    internal_row <-
      mapped[
        mapped$node_type == "internal",
        ,
        drop = FALSE
      ]

    expect_equal(
      internal_row$original_node,
      j1_node
    )

    expect_equal(
      internal_row$local_node,
      3L
    )
  }
)


test_that(
  "node identity survives pruning an internal patch",
  {

    reference_tree <-
      make_node_map_test_tree()

    n_tip <-
      ape::Ntip(reference_tree)

    i1_node <-
      n_tip +
      match(
        "I1",
        reference_tree$node.label
      )

    j1_node <-
      n_tip +
      match(
        "J1",
        reference_tree$node.label
      )

    user_tree <- reference_tree
    user_tree$node.label <- NULL

    frozen <- freeze_node_identity(
      user_tree
    )

    mother <- ape::drop.tip(
      frozen$tree,
      tip = c("A", "B"),
      collapse.singles = TRUE
    )

    mapped <- map_component_nodes(
      mother,
      frozen$node_map
    )

    removed_original_nodes <- c(
      match("A", reference_tree$tip.label),
      match("B", reference_tree$tip.label),
      i1_node,
      j1_node
    )

    expected_ids <-
      frozen$node_map$stable_id[
        !frozen$node_map$original_node %in%
          removed_original_nodes
      ]

    expect_equal(
      sort(mapped$stable_id),
      sort(expected_ids)
    )
  }
)


test_that(
  "node identity survives pruning a terminal patch",
  {

    reference_tree <-
      make_node_map_test_tree()

    n_tip <-
      ape::Ntip(reference_tree)

    i1_node <-
      n_tip +
      match(
        "I1",
        reference_tree$node.label
      )

    c_tip <-
      match(
        "C",
        reference_tree$tip.label
      )

    user_tree <- reference_tree
    user_tree$node.label <- NULL

    frozen <- freeze_node_identity(
      user_tree
    )

    mother <- ape::drop.tip(
      frozen$tree,
      tip = "C",
      collapse.singles = TRUE
    )

    mapped <- map_component_nodes(
      mother,
      frozen$node_map
    )

    removed_original_nodes <- c(
      c_tip,
      i1_node
    )

    expected_ids <-
      frozen$node_map$stable_id[
        !frozen$node_map$original_node %in%
          removed_original_nodes
      ]

    expect_equal(
      sort(mapped$stable_id),
      sort(expected_ids)
    )
  }
)


test_that(
  "original and stable node identifiers round-trip correctly",
  {

    tree <- make_node_map_test_tree()
    tree$node.label <- NULL

    frozen <- freeze_node_identity(
      tree
    )

    original_nodes <- c(
      1L,
      3L,
      ape::Ntip(tree) + 1L,
      ape::Ntip(tree) + tree$Nnode
    )

    stable_ids <-
      stable_id_from_original_node(
        frozen$node_map,
        original_nodes
      )

    recovered_nodes <-
      original_node_from_stable_id(
        frozen$node_map,
        stable_ids
      )

    expect_identical(
      recovered_nodes,
      original_nodes
    )
  }
)


test_that(
  "invalid tree and identity inputs are rejected",
  {

    expect_error(
      freeze_node_identity(
        list()
      ),
      "must inherit"
    )

    tree <- make_node_map_test_tree()

    tree$tip.label[2] <-
      tree$tip.label[1]

    expect_error(
      freeze_node_identity(
        tree
      ),
      "unique"
    )

    good_tree <-
      make_node_map_test_tree()

    good_tree$node.label <- NULL

    frozen <-
      freeze_node_identity(
        good_tree
      )

    expect_error(
      stable_id_from_original_node(
        frozen$node_map,
        9999
      ),
      "Unknown"
    )

    expect_error(
      original_node_from_stable_id(
        frozen$node_map,
        "unknown_node"
      ),
      "Unknown"
    )
  }
)


test_that(
  "derived components without stable internal IDs are rejected",
  {

    tree <- make_node_map_test_tree()
    tree$node.label <- NULL

    frozen <- freeze_node_identity(
      tree
    )

    broken_component <-
      frozen$tree

    broken_component$node.label <-
      NULL

    expect_error(
      map_component_nodes(
        broken_component,
        frozen$node_map
      ),
      "lost"
    )
  }
)

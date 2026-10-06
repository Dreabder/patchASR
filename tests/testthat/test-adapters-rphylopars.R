make_rphylopars_adapter_tree <- function() {

  ape::read.tree(
    text = paste0(
      "(((((A:1,B:1)J1:1,(C:1,D:1)J2:1)J:1,",
      "E:1)I:1,F:1)H0:1,(G:1,H:1)K:1)ROOT;"
    )
  )
}


prepare_rphylopars_adapter_tree <- function() {

  tree <-
    make_rphylopars_adapter_tree()

  tree$node.label <-
    paste0(
      "node_",
      ape::Ntip(tree) +
        seq_len(
          tree$Nnode
        )
    )

  tree
}


make_rphylopars_adapter_states <- function(
    tree
) {

  stats::setNames(
    seq_len(
      ape::Ntip(tree)
    ),
    tree$tip.label
  )
}


test_that(
  "Rphylopars BM adapter returns stable finite internal-node estimates",
  {

    skip_if_not_installed(
      "Rphylopars"
    )

    tree <-
      prepare_rphylopars_adapter_tree()

    states <-
      make_rphylopars_adapter_states(
        tree
      )

    backend <-
      make_asr_rphylopars_bm(
        REML = TRUE
      )

    result <-
      backend(
        tree,
        states
      )

    expect_type(
      result,
      "double"
    )

    expect_length(
      result,
      tree$Nnode
    )

    expect_identical(
      names(result),
      tree$node.label
    )

    expect_true(
      all(
        is.finite(
          result
        )
      )
    )
  }
)


test_that(
  "Rphylopars BM adapter supports the minimum two-tip component",
  {

    skip_if_not_installed(
      "Rphylopars"
    )

    tree <-
      ape::read.tree(
        text = "(A:1,B:2);"
      )

    tree$node.label <-
      "node_3"

    states <-
      c(
        A = 1,
        B = 3
      )

    backend <-
      make_asr_rphylopars_bm()

    result <-
      backend(
        tree,
        states
      )

    expect_identical(
      names(result),
      "node_3"
    )

    expect_length(
      result,
      1L
    )

    expect_true(
      is.finite(
        result[[1]]
      )
    )
  }
)


test_that(
  "Rphylopars BM adapter rejects invalid REML arguments",
  {

    expect_error(
      make_asr_rphylopars_bm(
        REML = NA
      ),
      "REML"
    )

    expect_error(
      make_asr_rphylopars_bm(
        REML = 1
      ),
      "REML"
    )

    expect_error(
      make_asr_rphylopars_bm(
        REML = c(
          TRUE,
          FALSE
        )
      ),
      "REML"
    )
  }
)


test_that(
  "Rphylopars BM adapter integrates with patch_asr end to end",
  {

    skip_if_not_installed(
      "Rphylopars"
    )

    reference_tree <-
      make_rphylopars_adapter_tree()

    I_node <-
      ape::Ntip(
        reference_tree
      ) +
      match(
        "I",
        reference_tree$node.label
      )

    J_node <-
      ape::Ntip(
        reference_tree
      ) +
      match(
        "J",
        reference_tree$node.label
      )

    tree <-
      reference_tree

    tree$node.label <-
      NULL

    states <-
      stats::setNames(
        seq_len(
          ape::Ntip(tree)
        ),
        tree$tip.label
      )

    shifts <-
      data.frame(
        parent =
          I_node,
        child =
          J_node
      )

    backend <-
      make_asr_rphylopars_bm(
        REML = TRUE
      )

    fit <-
      patch_asr(
        tree = tree,
        states = states,
        shifts = shifts,
        asr_fun = backend
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

    expect_equal(
      nrow(
        fit$boundary_states
      ),
      1L
    )
  }
)

test_that(
  "Rphylopars BM adapter handles numeric tip labels that collide with local node IDs",
  {

    skip_if_not_installed(
      "Rphylopars"
    )

    ## Six tips imply local internal ape IDs 7:11.
    ## Tip labels 7, 8, 9, and 10 therefore deliberately
    ## collide with local internal-node numbers.
    tree <-
      ape::read.tree(
        text = paste0(
          "(((1:1,2:1):1,(7:1,8:1):1):1,",
          "(9:1,10:1):2);"
        )
      )


    tree$node.label <-
      paste0(
        "node_",
        101L +
          seq_len(
            tree$Nnode
          )
      )


    states <-
      stats::setNames(
        seq_len(
          ape::Ntip(tree)
        ),
        tree$tip.label
      )


    local_internal_ids <-
      as.character(
        ape::Ntip(tree) +
          seq_len(
            tree$Nnode
          )
      )


    ## Confirm that this test really contains the collision
    ## that caused the real 128-tip failure.
    expect_true(
      length(
        intersect(
          tree$tip.label,
          local_internal_ids
        )
      ) > 0L
    )


    backend <-
      make_asr_rphylopars_bm(
        REML = TRUE
      )


    result <-
      backend(
        tree,
        states
      )


    expect_type(
      result,
      "double"
    )


    expect_length(
      result,
      tree$Nnode
    )


    expect_identical(
      names(result),
      tree$node.label
    )


    expect_true(
      all(
        is.finite(
          result
        )
      )
    )
  }
)

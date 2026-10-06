make_ape_adapter_tree <- function() {

  ape::read.tree(
    text = paste0(
      "(((((A:1,B:1)J1:1,(C:1,D:1)J2:1)J:1,",
      "E:1)I:1,F:1)H0:1,(G:1,H:1)K:1)ROOT;"
    )
  )
}


prepare_ape_adapter_tree <- function() {

  tree <-
    make_ape_adapter_tree()

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


make_ape_adapter_states <- function(
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
  "ape adapters return exactly one finite estimate per stable internal node",
  {

    tree <-
      prepare_ape_adapter_tree()

    states <-
      make_ape_adapter_states(
        tree
      )

    backends <- list(
      PIC =
        asr_ape_pic,
      ML_BM =
        asr_ape_ml_bm,
      GLS_BM =
        asr_ape_gls_bm
    )

    for (backend_name in
         names(backends)) {

      result <-
        backends[[backend_name]](
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
  }
)


test_that(
  "ape adapters reject terminal states that are not in exact tree-tip order",
  {

    tree <-
      prepare_ape_adapter_tree()

    states <-
      make_ape_adapter_states(
        tree
      )

    bad_states <-
      rev(
        states
      )

    expect_error(
      asr_ape_pic(
        tree,
        bad_states
      ),
      "exact tree-tip order"
    )

    expect_error(
      asr_ape_ml_bm(
        tree,
        bad_states
      ),
      "exact tree-tip order"
    )

    expect_error(
      asr_ape_gls_bm(
        tree,
        bad_states
      ),
      "exact tree-tip order"
    )
  }
)


test_that(
  "ape result standardization rejects an unrecognized node identity map",
  {

    tree <-
      prepare_ape_adapter_tree()

    bad_fit <- list(
      ace =
        stats::setNames(
          seq_len(
            tree$Nnode
          ),
          paste0(
            "wrong_",
            seq_len(
              tree$Nnode
            )
          )
        )
    )

    expect_error(
      standardize_ape_ace_result(
        bad_fit,
        tree,
        "test method"
      ),
      "cannot be mapped safely"
    )
  }
)


test_that(
  "ape adapters integrate with patch_asr end to end",
  {

    reference_tree <-
      make_ape_adapter_tree()

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

    backends <- list(
      PIC =
        asr_ape_pic,
      ML_BM =
        asr_ape_ml_bm,
      GLS_BM =
        asr_ape_gls_bm
    )

    for (backend_name in
         names(backends)) {

      fit <-
        patch_asr(
          tree = tree,
          states = states,
          shifts = shifts,
          asr_fun =
            backends[[backend_name]]
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
  }
)

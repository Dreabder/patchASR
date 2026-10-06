make_bayes_adapter_tree <- function() {

  ape::read.tree(
    text = paste0(
      "(((((A:1,B:1)J1:1,(C:1,D:1)J2:1)J:1,",
      "E:1)I:1,F:1)H0:1,(G:1,H:1)K:1)ROOT;"
    )
  )
}


prepare_bayes_adapter_tree <- function() {

  tree <-
    make_bayes_adapter_tree()

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


make_bayes_adapter_states <- function(
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
  "phytools Bayesian adapter returns stable finite estimates",
  {

    skip_if_not_installed(
      "phytools"
    )

    tree <-
      prepare_bayes_adapter_tree()

    states <-
      make_bayes_adapter_states(
        tree
      )

    backend <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
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
  "phytools Bayesian adapter is reproducible for the same base seed",
  {

    skip_if_not_installed(
      "phytools"
    )

    tree <-
      prepare_bayes_adapter_tree()

    states <-
      make_bayes_adapter_states(
        tree
      )

    backend <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
      )

    result_1 <-
      backend(
        tree,
        states
      )

    result_2 <-
      backend(
        tree,
        states
      )

    expect_identical(
      result_1,
      result_2
    )
  }
)


test_that(
  "different Bayesian base seeds produce different stochastic results",
  {

    skip_if_not_installed(
      "phytools"
    )

    tree <-
      prepare_bayes_adapter_tree()

    states <-
      make_bayes_adapter_states(
        tree
      )

    backend_1 <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
      )

    backend_2 <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 124L
      )

    result_1 <-
      backend_1(
        tree,
        states
      )

    result_2 <-
      backend_2(
        tree,
        states
      )

    expect_false(
      identical(
        result_1,
        result_2
      )
    )
  }
)


test_that(
  "Bayesian adapter preserves the caller RNG state",
  {

    skip_if_not_installed(
      "phytools"
    )

    tree <-
      prepare_bayes_adapter_tree()

    states <-
      make_bayes_adapter_states(
        tree
      )

    backend <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
      )

    set.seed(
      999L
    )

    rng_before <-
      .Random.seed

    invisible(
      backend(
        tree,
        states
      )
    )

    rng_after <-
      .Random.seed

    expect_identical(
      rng_after,
      rng_before
    )
  }
)


test_that(
  "Bayesian adapter rejects a one-internal-node component explicitly",
  {

    skip_if_not_installed(
      "phytools"
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
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
      )

    expect_error(
      backend(
        tree,
        states
      ),
      "one internal node"
    )
  }
)


test_that(
  "Bayesian adapter validates factory parameters",
  {

    expect_error(
      make_asr_phytools_bayes(
        ngen = 0
      ),
      "ngen"
    )

    expect_error(
      make_asr_phytools_bayes(
        sample_freq = 0
      ),
      "sample_freq"
    )

    expect_error(
      make_asr_phytools_bayes(
        burnin_frac = 1
      ),
      "burnin_frac"
    )

    expect_error(
      make_asr_phytools_bayes(
        seed = NA
      ),
      "seed"
    )

    expect_error(
      make_asr_phytools_bayes(
        ngen = 100,
        sample_freq = 200
      ),
      "cannot exceed"
    )
  }
)


test_that(
  "Bayesian adapter integrates with patch_asr end to end",
  {

    skip_if_not_installed(
      "phytools"
    )

    reference_tree <-
      make_bayes_adapter_tree()

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
        parent = I_node,
        child = J_node
      )

    backend <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
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

    expect_true(
      all(
        is.finite(
          fit$ancestral_states$estimate
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

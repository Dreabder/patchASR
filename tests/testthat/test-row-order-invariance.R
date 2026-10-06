make_row_order_invariance_case <- function() {

  tree_text <- paste0(
    "(",

    "(",
    "(",
    "((A:1,B:1):1,(C:1,D:1):1)JA:1,",
    "(E:1,F:1)SA:1",
    ")IA:1,",
    "(G:1,H:1)GA:1",
    ")HA:1,",

    "(",
    "(",
    "((I:1,J:1):1,(K:1,L:1):1)JB:1,",
    "(M:1,N:1)SB:1",
    ")IB:1,",
    "(O:1,P:1)GB:1",
    ")HB:1",

    ")ROOT;"
  )


  reference_tree <-
    ape::read.tree(
      text = tree_text
    )


  node_number_from_label <- function(
    tree,
    label
  ) {

    index <- match(
      label,
      tree$node.label
    )

    if (is.na(index)) {

      stop(
        "Could not find internal-node label: ",
        label
      )
    }

    ape::Ntip(tree) +
      index
  }


  IA <-
    node_number_from_label(
      reference_tree,
      "IA"
    )

  JA <-
    node_number_from_label(
      reference_tree,
      "JA"
    )

  IB <-
    node_number_from_label(
      reference_tree,
      "IB"
    )

  JB <-
    node_number_from_label(
      reference_tree,
      "JB"
    )


  tree <-
    reference_tree

  tree$node.label <-
    NULL


  states <-
    stats::setNames(
      c(
        1.2,
        2.1,
        3.7,
        4.4,
        2.8,
        5.6,
        4.9,
        7.3,
        6.1,
        8.8,
        7.6,
        10.2,
        9.4,
        11.7,
        10.9,
        13.5
      ),
      tree$tip.label
    )


  shifts_AB <-
    data.frame(
      parent = c(
        IA,
        IB
      ),
      child = c(
        JA,
        JB
      )
    )


  shifts_BA <-
    data.frame(
      parent = c(
        IB,
        IA
      ),
      child = c(
        JB,
        JA
      )
    )


  list(
    tree = tree,
    states = states,
    shifts_AB = shifts_AB,
    shifts_BA = shifts_BA
  )
}


extract_order_invariance_estimates <- function(
    fit
) {

  x <-
    fit$ancestral_states[
      ,
      c(
        "original_node",
        "estimate"
      ),
      drop = FALSE
    ]


  x$original_node <-
    as.integer(
      x$original_node
    )


  x <-
    x[
      order(
        x$original_node
      ),
      ,
      drop = FALSE
    ]


  stats::setNames(
    as.numeric(
      x$estimate
    ),
    as.character(
      x$original_node
    )
  )
}


deterministic_identity_backend <- function(
    tree,
    states
) {

  stable_ids <-
    tree$node.label


  values <-
    as.numeric(
      sub(
        "^node_",
        "",
        stable_ids
      )
    )


  stats::setNames(
    values,
    stable_ids
  )
}


test_that(
  "patch_asr is invariant to compatible shift-row order",
  {

    case <-
      make_row_order_invariance_case()


    fit_AB <-
      patch_asr(
        tree = case$tree,
        states = case$states,
        shifts = case$shifts_AB,
        asr_fun =
          deterministic_identity_backend
      )


    fit_BA <-
      patch_asr(
        tree = case$tree,
        states = case$states,
        shifts = case$shifts_BA,
        asr_fun =
          deterministic_identity_backend
      )


    estimates_AB <-
      extract_order_invariance_estimates(
        fit_AB
      )


    estimates_BA <-
      extract_order_invariance_estimates(
        fit_BA
      )


    expect_identical(
      estimates_AB,
      estimates_BA
    )


    expect_equal(
      nrow(
        fit_AB$boundary_states
      ),
      2L
    )


    expect_equal(
      nrow(
        fit_BA$boundary_states
      ),
      2L
    )


    expect_true(
      all(
        is.finite(
          estimates_AB
        )
      )
    )
  }
)


test_that(
  "Bayesian component-specific seeding is invariant to shift-row order",
  {

    skip_if_not_installed(
      "phytools"
    )


    case <-
      make_row_order_invariance_case()


    backend <-
      make_asr_phytools_bayes(
        ngen = 400L,
        sample_freq = 100L,
        burnin_frac = 0.20,
        seed = 123L
      )


    fit_AB <-
      patch_asr(
        tree = case$tree,
        states = case$states,
        shifts = case$shifts_AB,
        asr_fun = backend
      )


    fit_BA <-
      patch_asr(
        tree = case$tree,
        states = case$states,
        shifts = case$shifts_BA,
        asr_fun = backend
      )


    estimates_AB <-
      extract_order_invariance_estimates(
        fit_AB
      )


    estimates_BA <-
      extract_order_invariance_estimates(
        fit_BA
      )


    expect_identical(
      estimates_AB,
      estimates_BA
    )


    expect_equal(
      nrow(
        fit_AB$boundary_states
      ),
      2L
    )


    expect_equal(
      nrow(
        fit_BA$boundary_states
      ),
      2L
    )
  }
)

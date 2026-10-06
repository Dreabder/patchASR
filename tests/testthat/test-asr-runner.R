make_asr_runner_test_tree <- function() {

  ape::read.tree(
    text = paste0(
      "((((A:1,B:1)J1:1,C:1)I1:1,D:1)H1:1,",
      "((E:1,F:1)J2:1,(G:1,H:1)S2:1)I2:1)ROOT;"
    )
  )
}


asr_runner_node <- function(
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
      "Unknown ASR-runner test label: ",
      label
    )
  }

  as.integer(
    ape::Ntip(tree) +
      internal_match
  )
}


make_asr_runner_states <- function() {

  ## Deliberately scrambled.
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


prepare_asr_runner_fixture <- function(
    shifts
) {

  reference_tree <-
    make_asr_runner_test_tree()

  user_tree <-
    reference_tree

  user_tree$node.label <-
    NULL

  frozen <-
    freeze_node_identity(
      user_tree
    )

  plan <- build_patch_plan(
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
      partition
  )
}


valid_mock_asr <- function(
    tree,
    states
) {

  ## The runner must supply states in exact component-tip order.
  if (!identical(
    names(states),
    tree$tip.label
  )) {
    stop(
      "Mock backend received component states in the wrong order."
    )
  }

  ids <-
    tree$node.label

  values <-
    seq_along(ids) +
    mean(states)

  result <- setNames(
    values,
    ids
  )

  ## Deliberately return the estimates in reverse order.
  result[
    rev(
      seq_along(result)
    )
  ]
}


test_that(
  "component tip states are selected and reordered correctly",
  {

    tree <- ape::read.tree(
      text =
        "((A:1,B:1):1,C:1);"
    )

    tree$node.label <-
      c(
        "node_4",
        "node_5"
      )

    states <- c(
      C = 3,
      B = 2,
      A = 1
    )

    result <-
      prepare_component_states(
        tree,
        states
      )

    expect_identical(
      names(result),
      tree$tip.label
    )

    expect_equal(
      unname(result),
      c(1, 2, 3)
    )
  }
)


test_that(
  "valid ASR output is standardized by stable node identity",
  {

    tree <- ape::read.tree(
      text =
        "(((A:1,B:1):1,C:1):1,D:1);"
    )

    tree$node.label <-
      c(
        "node_5",
        "node_6",
        "node_7"
      )

    states <- c(
      D = 4,
      B = 2,
      A = 1,
      C = 3
    )

    result <- run_component_asr(
      component_tree =
        tree,
      states =
        states,
      asr_fun =
        valid_mock_asr,
      requires_asr =
        TRUE
    )

    expect_identical(
      names(result),
      tree$node.label
    )

    expected <- setNames(
      seq_along(
        tree$node.label
      ) +
        mean(
          states[
            tree$tip.label
          ]
        ),
      tree$node.label
    )

    expect_equal(
      result,
      expected
    )
  }
)


test_that(
  "missing component tip states are rejected",
  {

    tree <- ape::read.tree(
      text =
        "((A:1,B:1):1,C:1);"
    )

    tree$node.label <-
      c(
        "node_4",
        "node_5"
      )

    states <- c(
      A = 1,
      C = 3
    )

    expect_error(
      prepare_component_states(
        tree,
        states
      ),
      "Missing observed state"
    )
  }
)


test_that(
  "malformed ASR outputs are rejected",
  {

    tree <- ape::read.tree(
      text =
        "(((A:1,B:1):1,C:1):1,D:1);"
    )

    tree$node.label <-
      c(
        "node_5",
        "node_6",
        "node_7"
      )

    states <- c(
      A = 1,
      B = 2,
      C = 3,
      D = 4
    )

    unnamed_backend <- function(
    tree,
    states
    ) {
      rep(
        1,
        tree$Nnode
      )
    }

    expect_error(
      run_component_asr(
        tree,
        states,
        unnamed_backend
      ),
      "named numeric vector"
    )

    missing_backend <- function(
    tree,
    states
    ) {

      ids <-
        tree$node.label

      setNames(
        rep(
          1,
          length(ids) - 1L
        ),
        ids[
          -length(ids)
        ]
      )
    }

    expect_error(
      run_component_asr(
        tree,
        states,
        missing_backend
      ),
      "estimate"
    )

    duplicated_backend <- function(
    tree,
    states
    ) {

      ids <-
        tree$node.label

      bad_names <-
        ids

      bad_names[[2L]] <-
        bad_names[[1L]]

      setNames(
        seq_along(ids),
        bad_names
      )
    }

    expect_error(
      run_component_asr(
        tree,
        states,
        duplicated_backend
      ),
      "duplicated"
    )

    unknown_backend <- function(
    tree,
    states
    ) {

      ids <-
        tree$node.label

      ids[[1L]] <-
        "unknown_node"

      setNames(
        seq_along(ids),
        ids
      )
    }

    expect_error(
      run_component_asr(
        tree,
        states,
        unknown_backend
      ),
      "incorrect internal-node set"
    )

    nonfinite_backend <- function(
    tree,
    states
    ) {

      ids <-
        tree$node.label

      values <-
        rep(
          1,
          length(ids)
        )

      values[[1L]] <-
        NA_real_

      setNames(
        values,
        ids
      )
    }

    expect_error(
      run_component_asr(
        tree,
        states,
        nonfinite_backend
      ),
      "non-finite"
    )
  }
)


test_that(
  "terminal patches bypass the ASR backend",
  {

    call_count <- 0L

    must_not_run <- function(
    tree,
    states
    ) {

      call_count <<-
        call_count + 1L

      stop(
        "Terminal patch incorrectly called ASR."
      )
    }

    result <- run_component_asr(
      component_tree =
        NULL,
      states =
        make_asr_runner_states(),
      asr_fun =
        must_not_run,
      requires_asr =
        FALSE
    )

    expect_equal(
      call_count,
      0L
    )

    expect_length(
      result,
      0L
    )
  }
)


test_that(
  "partition ASR runs mother and internal patch components",
  {

    fixture <-
      prepare_asr_runner_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              asr_runner_node(
                tree,
                "I1"
              ),
            child =
              asr_runner_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <- run_partition_asr(
      partition =
        fixture$partition,
      states =
        make_asr_runner_states(),
      asr_fun =
        valid_mock_asr
    )

    expect_s3_class(
      result,
      "phyloPatch_component_asr"
    )

    expect_equal(
      sort(
        names(
          result$mother_states
        )
      ),
      sort(
        fixture$partition$
          mother_internal_ids
      )
    )

    expect_equal(
      sort(
        names(
          result$patch_states[[1L]]
        )
      ),
      sort(
        fixture$partition$
          patch_components[[1L]]$internal_ids
      )
    )

    expect_equal(
      sort(
        names(
          result$all_internal_states
        )
      ),
      sort(
        c(
          fixture$partition$
            mother_internal_ids,
          fixture$partition$
            patch_internal_ids
        )
      )
    )

    expect_false(
      anyDuplicated(
        names(
          result$all_internal_states
        )
      ) > 0L
    )
  }
)


test_that(
  "mixed internal and terminal patches call ASR only where required",
  {

    fixture <-
      prepare_asr_runner_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent = c(
              asr_runner_node(
                tree,
                "I1"
              ),
              asr_runner_node(
                tree,
                "J2"
              )
            ),
            child = c(
              asr_runner_node(
                tree,
                "J1"
              ),
              asr_runner_node(
                tree,
                "E"
              )
            )
          )
        }
      )

    call_log <-
      character(0L)

    logging_backend <- function(
    tree,
    states
    ) {

      call_log <<-
        c(
          call_log,
          paste(
            tree$tip.label,
            collapse = ","
          )
        )

      ids <-
        tree$node.label

      setNames(
        seq_along(ids) +
          mean(states),
        ids
      )
    }

    result <- run_partition_asr(
      partition =
        fixture$partition,
      states =
        make_asr_runner_states(),
      asr_fun =
        logging_backend
    )

    ## Mother + one internal patch.
    ## The terminal patch must not invoke the backend.
    expect_equal(
      length(call_log),
      2L
    )

    expect_length(
      result$patch_states[[2L]],
      0L
    )

    expect_true(
      length(
        result$patch_states[[1L]]
      ) > 0L
    )
  }
)


test_that(
  "component estimate table records source and component provenance",
  {

    fixture <-
      prepare_asr_runner_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              asr_runner_node(
                tree,
                "I1"
              ),
            child =
              asr_runner_node(
                tree,
                "J1"
              )
          )
        }
      )

    result <- run_partition_asr(
      partition =
        fixture$partition,
      states =
        make_asr_runner_states(),
      asr_fun =
        valid_mock_asr
    )

    expect_true(
      all(
        result$estimate_table$
          source ==
          "component_asr"
      )
    )

    expect_true(
      "mother" %in%
        result$estimate_table$
        component
    )

    expect_true(
      "patch_1" %in%
        result$estimate_table$
        component
    )

    expect_equal(
      sort(
        result$estimate_table$
          node
      ),
      sort(
        names(
          result$all_internal_states
        )
      )
    )
  }
)


test_that(
  "invalid partition and ASR backend inputs are rejected",
  {

    expect_error(
      run_partition_asr(
        partition =
          list(),
        states =
          make_asr_runner_states(),
        asr_fun =
          valid_mock_asr
      ),
      "partition"
    )

    fixture <-
      prepare_asr_runner_fixture(
        shifts = function(
    tree
        ) {

          data.frame(
            parent =
              asr_runner_node(
                tree,
                "I1"
              ),
            child =
              asr_runner_node(
                tree,
                "J1"
              )
          )
        }
      )

    expect_error(
      run_partition_asr(
        partition =
          fixture$partition,
        states =
          make_asr_runner_states(),
        asr_fun =
          "not_a_function"
      ),
      "must be a function"
    )
  }
)

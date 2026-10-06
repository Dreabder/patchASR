#!/usr/bin/env Rscript

## ============================================================
## Audit phytools::anc.Bayes() output
##
## Goals:
##   1. Inspect MCMC column names used for internal nodes.
##   2. Determine whether tree$node.label is retained.
##   3. Test regular, small, and minimum binary components.
##   4. Test reproducibility under an explicit random seed.
##
## Development audit only.
## ============================================================


if (!requireNamespace(
  "ape",
  quietly = TRUE
)) {
  stop("Package 'ape' is required.")
}

if (!requireNamespace(
  "phytools",
  quietly = TRUE
)) {
  stop("Package 'phytools' is required.")
}


audit_one_bayes_case <- function(
    case_name,
    tree,
    states,
    seed = 123L,
    ngen = 2000L,
    sample_freq = 100L
) {

  cat(
    "\n\n============================================================\n",
    "CASE: ", case_name, "\n",
    "============================================================\n",
    sep = ""
  )

  stable_ids <- paste0(
    "node_",
    ape::Ntip(tree) +
      seq_len(tree$Nnode)
  )

  tree$node.label <- stable_ids

  states <- states[
    tree$tip.label
  ]

  stopifnot(
    identical(
      names(states),
      tree$tip.label
    )
  )

  local_internal_ids <- as.character(
    ape::Ntip(tree) +
      seq_len(tree$Nnode)
  )

  cat("\nNtip:\n")
  print(
    ape::Ntip(tree)
  )

  cat("\nNnode:\n")
  print(
    tree$Nnode
  )

  cat("\ntree$tip.label:\n")
  print(
    tree$tip.label
  )

  cat("\ntree$node.label:\n")
  print(
    tree$node.label
  )

  cat("\nExpected local ape internal IDs:\n")
  print(
    local_internal_ids
  )

  cat("\nSeed:\n")
  print(
    seed
  )

  ## ----------------------------------------------------------
  ## First run
  ## ----------------------------------------------------------

  set.seed(
    seed
  )

  fit1 <- tryCatch(
    phytools::anc.Bayes(
      tree = tree,
      x = states,
      ngen = as.integer(ngen),
      control = list(
        sample = as.integer(sample_freq)
      )
    ),
    error = function(e) {
      e
    }
  )

  if (inherits(
    fit1,
    "error"
  )) {

    cat(
      "\nanc.Bayes FAILED:\n"
    )

    cat(
      conditionMessage(fit1),
      "\n"
    )

    return(
      invisible(NULL)
    )
  }

  cat("\nanc.Bayes succeeded.\n")

  cat("\nTop-level fit object names:\n")
  print(
    names(fit1)
  )

  if (is.null(
    fit1$mcmc
  )) {
    cat(
      "\nERROR: fit$mcmc is NULL.\n"
    )

    return(
      invisible(NULL)
    )
  }

  mcmc1 <- as.data.frame(
    fit1$mcmc
  )

  cat("\nclass(fit$mcmc):\n")
  print(
    class(
      fit1$mcmc
    )
  )

  cat("\ndim(mcmc):\n")
  print(
    dim(mcmc1)
  )

  cat("\ncolnames(mcmc):\n")
  print(
    colnames(mcmc1)
  )

  cat("\nFirst rows of mcmc:\n")
  print(
    utils::head(
      mcmc1,
      5L
    )
  )

  column_names <- colnames(
    mcmc1
  )

  cat(
    "\nAll stable node labels occur in MCMC columns: ",
    all(
      stable_ids %in%
        column_names
    ),
    "\n",
    sep = ""
  )

  cat(
    "All local ape internal IDs occur in MCMC columns: ",
    all(
      local_internal_ids %in%
        column_names
    ),
    "\n",
    sep = ""
  )

  cat(
    "Number of stable node columns found: ",
    sum(
      stable_ids %in%
        column_names
    ),
    " / ",
    length(stable_ids),
    "\n",
    sep = ""
  )

  cat(
    "Number of local ape ID columns found: ",
    sum(
      local_internal_ids %in%
        column_names
    ),
    " / ",
    length(local_internal_ids),
    "\n",
    sep = ""
  )

  ## ----------------------------------------------------------
  ## Same seed reproducibility check
  ## ----------------------------------------------------------

  set.seed(
    seed
  )

  fit2 <- tryCatch(
    phytools::anc.Bayes(
      tree = tree,
      x = states,
      ngen = as.integer(ngen),
      control = list(
        sample = as.integer(sample_freq)
      )
    ),
    error = function(e) {
      e
    }
  )

  if (inherits(
    fit2,
    "error"
  )) {

    cat(
      "\nSecond anc.Bayes run FAILED:\n"
    )

    cat(
      conditionMessage(fit2),
      "\n"
    )

  } else {

    mcmc2 <- as.data.frame(
      fit2$mcmc
    )

    cat(
      "\nSame seed gives identical MCMC output: ",
      identical(
        mcmc1,
        mcmc2
      ),
      "\n",
      sep = ""
    )
  }

  ## ----------------------------------------------------------
  ## Different seed diagnostic
  ## ----------------------------------------------------------

  set.seed(
    seed + 1L
  )

  fit3 <- tryCatch(
    phytools::anc.Bayes(
      tree = tree,
      x = states,
      ngen = as.integer(ngen),
      control = list(
        sample = as.integer(sample_freq)
      )
    ),
    error = function(e) {
      e
    }
  )

  if (!inherits(
    fit3,
    "error"
  )) {

    mcmc3 <- as.data.frame(
      fit3$mcmc
    )

    cat(
      "Different seed gives different MCMC output: ",
      !identical(
        mcmc1,
        mcmc3
      ),
      "\n",
      sep = ""
    )
  }

  cat(
    "\nCASE COMPLETE: ",
    case_name,
    "\n",
    sep = ""
  )

  invisible(
    list(
      tree = tree,
      states = states,
      fit = fit1,
      mcmc = mcmc1
    )
  )
}


## ============================================================
## Case 1: regular six-tip component
## ============================================================

tree_regular <- ape::read.tree(
  text = paste0(
    "((A:1,B:1):1,",
    "((C:1,D:1):1,(E:1,F:1):1):1);"
  )
)

states_regular <- stats::setNames(
  1:6,
  tree_regular$tip.label
)

audit_one_bayes_case(
  case_name = "regular_six_tip_tree",
  tree = tree_regular,
  states = states_regular
)


## ============================================================
## Case 2: three-tip component
## ============================================================

tree_three <- ape::read.tree(
  text = "((A:1,B:1):1,C:1);"
)

states_three <- stats::setNames(
  1:3,
  tree_three$tip.label
)

audit_one_bayes_case(
  case_name = "three_tip_tree",
  tree = tree_three,
  states = states_three
)


## ============================================================
## Case 3: two-tip / one-internal-node component
## ============================================================

tree_two <- ape::read.tree(
  text = "(A:1,B:2);"
)

states_two <- c(
  A = 1,
  B = 3
)

audit_one_bayes_case(
  case_name = "two_tip_one_internal_node_tree",
  tree = tree_two,
  states = states_two
)


cat(
  "\n\n============================================================\n",
  "PHYTOOLS anc.Bayes OUTPUT AUDIT FINISHED\n",
  "============================================================\n"
)

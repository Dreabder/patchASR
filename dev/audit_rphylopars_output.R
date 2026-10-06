#!/usr/bin/env Rscript

## ============================================================
## Audit Rphylopars ancestral reconstruction output
##
## Purpose:
##   Determine how fit$anc_recon identifies internal nodes
##   under the Rphylopars version installed on this machine.
##
## This is a development audit only.
## It does NOT modify phyloPatch.
## ============================================================


if (!requireNamespace(
  "ape",
  quietly = TRUE
)) {
  stop("Package 'ape' is required.")
}

if (!requireNamespace(
  "Rphylopars",
  quietly = TRUE
)) {
  stop("Package 'Rphylopars' is required.")
}


audit_one_rphylopars_case <- function(
    case_name,
    tree,
    states,
    REML = TRUE
) {

  cat(
    "\n\n",
    "============================================================\n",
    "CASE: ", case_name, "\n",
    "============================================================\n",
    sep = ""
  )

  ## ----------------------------------------------------------
  ## Assign stable labels similar to those used by phyloPatch
  ## ----------------------------------------------------------

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

  ## ----------------------------------------------------------
  ## Show component-tree identity information
  ## ----------------------------------------------------------

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

  local_internal_ids <- as.character(
    ape::Ntip(tree) +
      seq_len(tree$Nnode)
  )

  cat("\nExpected local ape internal IDs:\n")
  print(
    local_internal_ids
  )

  ## ----------------------------------------------------------
  ## Construct Rphylopars input
  ## ----------------------------------------------------------

  trait_data <- data.frame(
    species = names(states),
    trait = as.numeric(states),
    stringsAsFactors = FALSE
  )

  cat("\ntrait_data:\n")
  print(
    trait_data
  )

  ## ----------------------------------------------------------
  ## Run Rphylopars
  ## ----------------------------------------------------------

  fit <- tryCatch(
    Rphylopars::phylopars(
      trait_data = trait_data,
      tree = tree,
      model = "BM",
      pheno_error = FALSE,
      REML = REML
    ),
    error = function(e) {
      e
    }
  )

  if (inherits(
    fit,
    "error"
  )) {

    cat(
      "\nRphylopars FAILED for this case:\n"
    )

    cat(
      conditionMessage(fit),
      "\n"
    )

    return(
      invisible(NULL)
    )
  }

  cat("\nRphylopars fit succeeded.\n")

  cat("\nTop-level fit object names:\n")
  print(
    names(fit)
  )

  ## ----------------------------------------------------------
  ## Audit anc_recon
  ## ----------------------------------------------------------

  if (is.null(
    fit$anc_recon
  )) {
    cat(
      "\nERROR: fit$anc_recon is NULL.\n"
    )

    return(
      invisible(NULL)
    )
  }

  anc <- as.data.frame(
    fit$anc_recon
  )

  cat("\nclass(fit$anc_recon):\n")
  print(
    class(
      fit$anc_recon
    )
  )

  cat("\ndim(anc_recon):\n")
  print(
    dim(anc)
  )

  cat("\ncolnames(anc_recon):\n")
  print(
    colnames(anc)
  )

  cat("\nrownames(anc_recon):\n")
  print(
    rownames(anc)
  )

  cat("\nFirst rows of anc_recon:\n")
  print(
    utils::head(
      anc,
      10L
    )
  )

  cat("\nLast rows of anc_recon:\n")
  print(
    utils::tail(
      anc,
      10L
    )
  )

  ## ----------------------------------------------------------
  ## Explicit mapping checks
  ## ----------------------------------------------------------

  rn <- rownames(anc)

  expected_total_rows <-
    ape::Ntip(tree) +
    tree$Nnode

  cat(
    "\nExpected Ntip + Nnode rows: ",
    expected_total_rows,
    "\n",
    sep = ""
  )

  cat(
    "Observed rows: ",
    nrow(anc),
    "\n",
    sep = ""
  )

  cat(
    "Row count equals Ntip + Nnode: ",
    identical(
      nrow(anc),
      expected_total_rows
    ),
    "\n",
    sep = ""
  )

  if (is.null(rn)) {

    cat(
      "anc_recon has no row names.\n"
    )

  } else {

    cat(
      "All stable node labels occur in row names: ",
      all(
        stable_ids %in% rn
      ),
      "\n",
      sep = ""
    )

    cat(
      "All local ape internal IDs occur in row names: ",
      all(
        local_internal_ids %in% rn
      ),
      "\n",
      sep = ""
    )

    cat(
      "All tip labels occur in row names: ",
      all(
        tree$tip.label %in% rn
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
      fit = fit,
      anc = anc
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
  c(
    1,
    2,
    3,
    4,
    5,
    6
  ),
  tree_regular$tip.label
)

audit_one_rphylopars_case(
  case_name = "regular_six_tip_tree",
  tree = tree_regular,
  states = states_regular
)


## ============================================================
## Case 2: small three-tip component
## ============================================================

tree_three <- ape::read.tree(
  text = "((A:1,B:1):1,C:1);"
)

states_three <- stats::setNames(
  c(
    1,
    2,
    3
  ),
  tree_three$tip.label
)

audit_one_rphylopars_case(
  case_name = "three_tip_tree",
  tree = tree_three,
  states = states_three
)


## ============================================================
## Case 3: minimum two-tip / one-internal-node component
## ============================================================

tree_two <- ape::read.tree(
  text = "(A:1,B:2);"
)

states_two <- stats::setNames(
  c(
    1,
    3
  ),
  tree_two$tip.label
)

audit_one_rphylopars_case(
  case_name = "two_tip_one_internal_node_tree",
  tree = tree_two,
  states = states_two
)


cat(
  "\n\n============================================================\n",
  "RPHYLOPARS OUTPUT AUDIT FINISHED\n",
  "============================================================\n"
)

#!/usr/bin/env Rscript

## ============================================================
## Multiple-SBLC shift-row-order invariance validation
##
## Goal:
##
## Compare
##
##   shifts = A, B
##
## against
##
##   shifts = B, A
##
## while keeping the biological shift edges identical.
##
## Methods:
##   PIC
##   ML-BM
##   GLS-BM
##   Rphylopars-BM
##   phytools Bayesian
##
## Expected:
##   Final ancestral-state estimates must not depend on the
##   input row order of the shift table.
## ============================================================


## ------------------------------------------------------------
## Dependencies
## ------------------------------------------------------------

if (!requireNamespace(
  "ape",
  quietly = TRUE
)) {
  stop("Package 'ape' is required.")
}

if (!requireNamespace(
  "pkgload",
  quietly = TRUE
)) {
  stop("Package 'pkgload' is required.")
}


pkgload::load_all(
  ".",
  quiet = TRUE
)


## ============================================================
## Construct a controlled fully bifurcating 16-tip tree
##
## Shift A:
##   IA -> JA
##
## Shift B:
##   IB -> JB
##
## Each patch contains four tips.
## The two patch descendant sets are disjoint.
## ============================================================

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


reference_tree <- ape::read.tree(
  text = tree_text
)


if (!ape::is.rooted(
  reference_tree
)) {
  stop(
    "Validation tree is not rooted."
  )
}


if (!ape::is.binary.tree(
  reference_tree
)) {
  stop(
    "Validation tree is not fully bifurcating."
  )
}


## ------------------------------------------------------------
## Resolve original ape node numbers from temporary node labels
## ------------------------------------------------------------

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


IA <- node_number_from_label(
  reference_tree,
  "IA"
)

JA <- node_number_from_label(
  reference_tree,
  "JA"
)

IB <- node_number_from_label(
  reference_tree,
  "IB"
)

JB <- node_number_from_label(
  reference_tree,
  "JB"
)


cat(
  "Ntip:  ",
  ape::Ntip(reference_tree),
  "\n",
  sep = ""
)

cat(
  "Nnode: ",
  reference_tree$Nnode,
  "\n",
  sep = ""
)

cat(
  "Shift A: ",
  IA,
  " -> ",
  JA,
  "\n",
  sep = ""
)

cat(
  "Shift B: ",
  IB,
  " -> ",
  JB,
  "\n",
  sep = ""
)


## ------------------------------------------------------------
## Confirm both shift edges really exist
## ------------------------------------------------------------

check_edge <- function(
    tree,
    parent,
    child
) {

  sum(
    tree$edge[, 1] == parent &
      tree$edge[, 2] == child
  ) == 1L
}


if (!check_edge(
  reference_tree,
  IA,
  JA
)) {

  stop(
    "Shift A is not an exact edge."
  )
}


if (!check_edge(
  reference_tree,
  IB,
  JB
)) {

  stop(
    "Shift B is not an exact edge."
  )
}


## ------------------------------------------------------------
## Remove the temporary biological labels.
##
## patchASR should work from original ape node numbers and
## create its own stable working identity internally.
## ------------------------------------------------------------

tree <- reference_tree

tree$node.label <- NULL


## ------------------------------------------------------------
## Continuous tip states
##
## Deliberately non-linear values avoid an overly symmetric
## numerical test case.
## ------------------------------------------------------------

states <- stats::setNames(
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


if (!identical(
  names(states),
  tree$tip.label
)) {

  stop(
    "Tip states are not in exact tree tip order."
  )
}


## ============================================================
## Same biological shifts, opposite row order
## ============================================================

shifts_AB <- data.frame(
  parent = c(
    IA,
    IB
  ),
  child = c(
    JA,
    JB
  )
)


shifts_BA <- data.frame(
  parent = c(
    IB,
    IA
  ),
  child = c(
    JB,
    JA
  )
)


cat(
  "\nshifts_AB:\n"
)

print(
  shifts_AB
)


cat(
  "\nshifts_BA:\n"
)

print(
  shifts_BA
)


## ============================================================
## ASR backends
## ============================================================

backends <- list(
  PIC =
    asr_ape_pic,

  ML_BM =
    asr_ape_ml_bm,

  GLS_BM =
    asr_ape_gls_bm,

  Rphylopars_BM =
    make_asr_rphylopars_bm(
      REML = TRUE
    ),

  Bayesian =
    make_asr_phytools_bayes(
      ngen = 1000L,
      sample_freq = 100L,
      burnin_frac = 0.20,
      seed = 123L
    )
)


## ============================================================
## Helpers
## ============================================================

extract_final_estimates <- function(
    fit
) {

  x <- fit$ancestral_states[
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


  x <- x[
    order(
      x$original_node
    ),
    ,
    drop = FALSE
  ]


  rownames(x) <- NULL

  x
}


## Deterministic methods allow a conservative floating-point
## tolerance. Bayesian is expected to be bit-for-bit identical
## because each biological component receives the same
## component-specific seed regardless of shift-row order.

deterministic_tolerance <- 1e-8


summary_list <- list()

nodewise_list <- list()


## ============================================================
## Run each backend twice:
##
##   AB
##   BA
## ============================================================

for (method_name in
     names(backends)) {

  cat(
    "\n\n============================================================\n",
    "METHOD: ",
    method_name,
    "\n",
    "============================================================\n",
    sep = ""
  )


  backend <-backends[[method_name]]


  ## ----------------------------------------------------------
  ## AB
  ## ----------------------------------------------------------

  fit_AB <-
    patch_asr(
      tree = tree,
      states = states,
      shifts = shifts_AB,
      asr_fun = backend
    )


  ## ----------------------------------------------------------
  ## BA
  ## ----------------------------------------------------------

  fit_BA <-
    patch_asr(
      tree = tree,
      states = states,
      shifts = shifts_BA,
      asr_fun = backend
    )


  ## ----------------------------------------------------------
  ## Final estimates
  ## ----------------------------------------------------------

  result_AB <-
    extract_final_estimates(
      fit_AB
    )


  result_BA <-
    extract_final_estimates(
      fit_BA
    )


  if (!identical(
    result_AB$original_node,
    result_BA$original_node
  )) {

    stop(
      method_name,
      ": AB and BA runs contain different node-ID sets."
    )
  }


  ## ----------------------------------------------------------
  ## Basic completeness
  ## ----------------------------------------------------------

  complete_AB <-
    nrow(result_AB) ==
    tree$Nnode

  complete_BA <-
    nrow(result_BA) ==
    tree$Nnode


  finite_AB <-
    all(
      is.finite(
        result_AB$estimate
      )
    )

  finite_BA <-
    all(
      is.finite(
        result_BA$estimate
      )
    )


  boundary_count_AB <-
    nrow(
      fit_AB$boundary_states
    )

  boundary_count_BA <-
    nrow(
      fit_BA$boundary_states
    )


  ## ----------------------------------------------------------
  ## Numerical comparison
  ## ----------------------------------------------------------

  difference <-
    result_BA$estimate -
    result_AB$estimate


  absolute_difference <-
    abs(
      difference
    )


  max_abs_diff <-
    max(
      absolute_difference
    )


  mean_abs_diff <-
    mean(
      absolute_difference
    )


  rmse <-
    sqrt(
      mean(
        difference^2
      )
    )


  exact_identical <-
    identical(
      result_AB$estimate,
      result_BA$estimate
    )


  n_over_tolerance <-
    sum(
      absolute_difference >
        deterministic_tolerance
    )


  ## ----------------------------------------------------------
  ## PASS criterion
  ## ----------------------------------------------------------

  if (identical(
    method_name,
    "Bayesian"
  )) {

    numerical_pass <-
      exact_identical &&
      max_abs_diff == 0

  } else {

    numerical_pass <-
      all(
        absolute_difference <=
          deterministic_tolerance
      )
  }


  pass <-
    complete_AB &&
    complete_BA &&
    finite_AB &&
    finite_BA &&
    boundary_count_AB == 2L &&
    boundary_count_BA == 2L &&
    numerical_pass


  ## ----------------------------------------------------------
  ## Print
  ## ----------------------------------------------------------

  cat(
    "AB internal-node count: ",
    nrow(result_AB),
    "\n",
    sep = ""
  )


  cat(
    "BA internal-node count: ",
    nrow(result_BA),
    "\n",
    sep = ""
  )


  cat(
    "AB boundary count: ",
    boundary_count_AB,
    "\n",
    sep = ""
  )


  cat(
    "BA boundary count: ",
    boundary_count_BA,
    "\n",
    sep = ""
  )


  cat(
    "Exact estimate-vector identity: ",
    exact_identical,
    "\n",
    sep = ""
  )


  cat(
    "Maximum absolute difference: ",
    format(
      max_abs_diff,
      scientific = TRUE,
      digits = 16
    ),
    "\n",
    sep = ""
  )


  cat(
    "Method status: ",
    if (pass) {
      "PASS"
    } else {
      "FAIL"
    },
    "\n",
    sep = ""
  )


  ## ----------------------------------------------------------
  ## Save in memory
  ## ----------------------------------------------------------

  summary_list[[method_name]] <-
    data.frame(
      method =
        method_name,

      n_internal_nodes =
        nrow(result_AB),

      boundary_count_AB =
        boundary_count_AB,

      boundary_count_BA =
        boundary_count_BA,

      exact_identical =
        exact_identical,

      max_abs_diff =
        max_abs_diff,

      mean_abs_diff =
        mean_abs_diff,

      rmse =
        rmse,

      n_over_1e_8 =
        n_over_tolerance,

      status =
        if (pass) {
          "PASS"
        } else {
          "FAIL"
        },

      stringsAsFactors = FALSE
    )


  nodewise_list[[method_name]] <-
    data.frame(
      method =
        method_name,

      node_id =
        result_AB$original_node,

      estimate_AB =
        result_AB$estimate,

      estimate_BA =
        result_BA$estimate,

      difference =
        difference,

      abs_difference =
        absolute_difference,

      stringsAsFactors = FALSE
    )
}


## ============================================================
## Combine reports
## ============================================================

summary_table <-
  do.call(
    rbind,
    summary_list
  )

rownames(
  summary_table
) <- NULL


nodewise_table <-
  do.call(
    rbind,
    nodewise_list
  )

rownames(
  nodewise_table
) <- NULL


## ============================================================
## Overall validation status
## ============================================================

method_pass_vector <-
  summary_table$status ==
  "PASS"


cat(
  "\nPer-method PASS vector:\n"
)

print(
  stats::setNames(
    method_pass_vector,
    summary_table$method
  )
)


overall_pass <-
  length(
    method_pass_vector
  ) ==
  length(
    backends
  ) &&
  all(
    method_pass_vector
  )


## ============================================================
## Save validation evidence
## ============================================================

summary_file <- file.path(
  "dev",
  "regression",
  "shift_row_order_invariance_summary.csv"
)


nodewise_file <- file.path(
  "dev",
  "regression",
  "shift_row_order_invariance_nodewise.csv"
)


write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)


write.csv(
  nodewise_table,
  nodewise_file,
  row.names = FALSE
)


## ============================================================
## Final report
## ============================================================

cat(
  "\n\n============================================================\n",
  "MULTIPLE-SBLC SHIFT-ROW-ORDER INVARIANCE SUMMARY\n",
  "============================================================\n",
  sep = ""
)


print(
  summary_table,
  row.names = FALSE,
  digits = 16
)


cat(
  "\nOverall shift-row-order invariance status: ",
  if (overall_pass) {
    "PASS"
   } else {
    "FAIL"
  },
  "\n",
  sep = ""
)


cat(
  "\nSummary file:\n",
  summary_file,
  "\n",
  sep = ""
)


cat(
  "\nNode-wise file:\n",
  nodewise_file,
  "\n",
  sep = ""
)

#!/usr/bin/env Rscript

## ============================================================
## Bayesian reproducibility validation for phyloPatch
##
## Real dataset:
##   BM/balanced/jump_16/rep_001
##
## This script validates software-level stochastic behaviour:
##
##   1. Same input + same base seed -> identical result.
##   2. Result is independent of the caller's RNG state.
##   3. The caller's RNG state is restored after patch_asr().
##   4. Different base seeds -> different Bayesian result.
##
## IMPORTANT:
## Short MCMC settings are used here for software validation,
## not for scientific posterior inference or convergence claims.
## ============================================================


## ------------------------------------------------------------
## Required packages
## ------------------------------------------------------------

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

if (!requireNamespace(
  "pkgload",
  quietly = TRUE
)) {
  stop("Package 'pkgload' is required.")
}


## ------------------------------------------------------------
## Load current development version of phyloPatch
## ------------------------------------------------------------

pkgload::load_all(
  ".",
  quiet = TRUE
)


## ------------------------------------------------------------
## Files
## ------------------------------------------------------------

case_dir <- paste0(
  "/home/chengyq/work/ASE/data/jump_simu/simu/",
  "BM/balanced/jump_16/rep_001"
)

tree_file <- file.path(
  case_dir,
  "tree.nwk"
)

trait_file <- file.path(
  case_dir,
  "node_traits.csv"
)

run1_file <- file.path(
  case_dir,
  "bayesian_validation_same_seed_run1.csv"
)

run2_file <- file.path(
  case_dir,
  "bayesian_validation_same_seed_run2.csv"
)

different_seed_file <- file.path(
  case_dir,
  "bayesian_validation_different_seed.csv"
)

summary_file <- file.path(
  case_dir,
  "bayesian_validation_summary.csv"
)


## ------------------------------------------------------------
## Read inputs
## ------------------------------------------------------------

tree <- ape::read.tree(
  tree_file
)

dat <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


get_one <- function(
    x,
    name
) {

  values <- unique(
    x[
      !is.na(x)
    ]
  )

  if (length(values) != 1L) {

    stop(
      name,
      " must contain exactly one unique non-missing value."
    )
  }

  values[[1]]
}


I <- as.integer(
  get_one(
    dat$I,
    "I"
  )
)

J <- as.integer(
  get_one(
    dat$J,
    "J"
  )
)


## ------------------------------------------------------------
## Confirm shift
## ------------------------------------------------------------

edge_rows <- which(
  tree$edge[, 1] == I &
    tree$edge[, 2] == J
)

if (length(edge_rows) != 1L) {

  stop(
    "I -> J does not identify exactly one edge."
  )
}


shifts <- data.frame(
  parent = I,
  child = J
)


## ------------------------------------------------------------
## Tip states
## ------------------------------------------------------------

trait_lookup <- stats::setNames(
  as.numeric(
    dat$trait_value
  ),
  as.character(
    dat$node_id
  )
)

states <- trait_lookup[
  tree$tip.label
]


if (length(states) !=
    ape::Ntip(tree) ||
    any(!is.finite(states)) ||
    !identical(
      names(states),
      tree$tip.label
    )) {

  stop(
    "Tip-state extraction failed."
  )
}


cat(
  "Tree tips:      ",
  ape::Ntip(tree),
  "\n",
  sep = ""
)

cat(
  "Internal nodes: ",
  tree$Nnode,
  "\n",
  sep = ""
)

cat(
  "Shift:          ",
  I,
  " -> ",
  J,
  "\n",
  sep = ""
)


## ============================================================
## Validation settings
##
## These settings are deliberately short because this is a
## reproducibility/software test, not an MCMC convergence study.
## ============================================================

validation_ngen <- 1000L
validation_sample_freq <- 100L
validation_burnin_frac <- 0.20

base_seed_same <- 123L
base_seed_different <- 124L


cat(
  "\nValidation MCMC settings:\n"
)

cat(
  "ngen        = ",
  validation_ngen,
  "\n",
  sep = ""
)

cat(
  "sample_freq = ",
  validation_sample_freq,
  "\n",
  sep = ""
)

cat(
  "burnin_frac = ",
  validation_burnin_frac,
  "\n",
  sep = ""
)

cat(
  "base seed 1 = ",
  base_seed_same,
  "\n",
  sep = ""
)

cat(
  "base seed 2 = ",
  base_seed_different,
  "\n",
  sep = ""
)


## ------------------------------------------------------------
## Backend with fixed base seed
## ------------------------------------------------------------

backend_same_seed <-
  make_asr_phytools_bayes(
    ngen = validation_ngen,
    sample_freq =
      validation_sample_freq,
    burnin_frac =
      validation_burnin_frac,
    seed = base_seed_same
  )


## ============================================================
## RUN 1
##
## Caller RNG state = 111
## ============================================================

cat(
  "\n============================================================\n",
  "BAYESIAN VALIDATION RUN 1\n",
  "============================================================\n",
  sep = ""
)


set.seed(
  111L
)

rng_before_run1 <-
  .Random.seed


fit1 <-
  patch_asr(
    tree = tree,
    states = states,
    shifts = shifts,
    asr_fun = backend_same_seed
  )


rng_after_run1 <-
  .Random.seed


rng_preserved_run1 <-
  identical(
    rng_before_run1,
    rng_after_run1
  )


cat(
  "Caller RNG preserved in run 1: ",
  rng_preserved_run1,
  "\n",
  sep = ""
)


## ============================================================
## RUN 2
##
## Same Bayesian base seed, but deliberately different caller
## RNG state = 222.
## ============================================================

cat(
  "\n============================================================\n",
  "BAYESIAN VALIDATION RUN 2\n",
  "============================================================\n",
  sep = ""
)


set.seed(
  222L
)

rng_before_run2 <-
  .Random.seed


fit2 <-
  patch_asr(
    tree = tree,
    states = states,
    shifts = shifts,
    asr_fun = backend_same_seed
  )


rng_after_run2 <-
  .Random.seed


rng_preserved_run2 <-
  identical(
    rng_before_run2,
    rng_after_run2
  )


cat(
  "Caller RNG preserved in run 2: ",
  rng_preserved_run2,
  "\n",
  sep = ""
)


## ============================================================
## Prepare comparable tables
## ============================================================

extract_estimates <- function(
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


result1 <-
  extract_estimates(
    fit1
  )

result2 <-
  extract_estimates(
    fit2
  )


if (!identical(
  result1$original_node,
  result2$original_node
)) {

  stop(
    "Run 1 and run 2 do not contain identical node IDs."
  )
}


same_seed_diff <-
  result2$estimate -
  result1$estimate


same_seed_max_abs_diff <-
  max(
    abs(
      same_seed_diff
    )
  )


same_seed_identical <-
  identical(
    result1$estimate,
    result2$estimate
  )


cat(
  "\nSame base seed despite different caller RNG:\n"
)

cat(
  "Exact estimate vector identical: ",
  same_seed_identical,
  "\n",
  sep = ""
)

cat(
  "Maximum absolute difference: ",
  format(
    same_seed_max_abs_diff,
    scientific = TRUE,
    digits = 16
  ),
  "\n",
  sep = ""
)


## ============================================================
## RUN 3
##
## Different Bayesian base seed.
## ============================================================

cat(
  "\n============================================================\n",
  "BAYESIAN VALIDATION RUN 3: DIFFERENT BASE SEED\n",
  "============================================================\n",
  sep = ""
)


backend_different_seed <-
  make_asr_phytools_bayes(
    ngen = validation_ngen,
    sample_freq =
      validation_sample_freq,
    burnin_frac =
      validation_burnin_frac,
    seed = base_seed_different
  )


set.seed(
  333L
)

rng_before_run3 <-
  .Random.seed


fit3 <-
  patch_asr(
    tree = tree,
    states = states,
    shifts = shifts,
    asr_fun = backend_different_seed
  )


rng_after_run3 <-
  .Random.seed


rng_preserved_run3 <-
  identical(
    rng_before_run3,
    rng_after_run3
  )


result3 <-
  extract_estimates(
    fit3
  )


if (!identical(
  result1$original_node,
  result3$original_node
)) {

  stop(
    "Different-seed run does not contain identical node IDs."
  )
}


different_seed_diff <-
  result3$estimate -
  result1$estimate


different_seed_max_abs_diff <-
  max(
    abs(
      different_seed_diff
    )
  )


different_seed_mean_abs_diff <-
  mean(
    abs(
      different_seed_diff
    )
  )


different_seed_rmse <-
  sqrt(
    mean(
      different_seed_diff^2
    )
  )


different_seed_changes_result <-
  any(
    result1$estimate !=
      result3$estimate
  )


cat(
  "Caller RNG preserved in run 3: ",
  rng_preserved_run3,
  "\n",
  sep = ""
)

cat(
  "Different base seed changes result: ",
  different_seed_changes_result,
  "\n",
  sep = ""
)

cat(
  "Different-seed maximum absolute difference: ",
  format(
    different_seed_max_abs_diff,
    scientific = TRUE,
    digits = 10
  ),
  "\n",
  sep = ""
)


## ============================================================
## General completeness checks
## ============================================================

expected_n_nodes <-
  tree$Nnode


complete_node_count <-
  nrow(result1) ==
  expected_n_nodes &&
  nrow(result2) ==
  expected_n_nodes &&
  nrow(result3) ==
  expected_n_nodes


all_finite <-
  all(
    is.finite(
      result1$estimate
    )
  ) &&
  all(
    is.finite(
      result2$estimate
    )
  ) &&
  all(
    is.finite(
      result3$estimate
    )
  )


unique_nodes <-
  !anyDuplicated(
    result1$original_node
  ) &&
  !anyDuplicated(
    result2$original_node
  ) &&
  !anyDuplicated(
    result3$original_node
  )


## ============================================================
## Save individual results
## ============================================================

names(result1) <-
  c(
    "node_id",
    "estimate"
  )

names(result2) <-
  c(
    "node_id",
    "estimate"
  )

names(result3) <-
  c(
    "node_id",
    "estimate"
  )


write.csv(
  result1,
  run1_file,
  row.names = FALSE
)

write.csv(
  result2,
  run2_file,
  row.names = FALSE
)

write.csv(
  result3,
  different_seed_file,
  row.names = FALSE
)


## ============================================================
## Summary
## ============================================================

summary_table <- data.frame(
  check = c(
    "complete_internal_node_count",
    "all_estimates_finite",
    "internal_node_ids_unique",
    "caller_rng_preserved_run1",
    "caller_rng_preserved_run2",
    "caller_rng_preserved_run3",
    "same_seed_exact_reproducibility",
    "different_seed_changes_result"
  ),

  value = c(
    complete_node_count,
    all_finite,
    unique_nodes,
    rng_preserved_run1,
    rng_preserved_run2,
    rng_preserved_run3,
    same_seed_identical,
    different_seed_changes_result
  ),

  expected = rep(
    TRUE,
    8L
  ),

  status = ifelse(
    c(
      complete_node_count,
      all_finite,
      unique_nodes,
      rng_preserved_run1,
      rng_preserved_run2,
      rng_preserved_run3,
      same_seed_identical,
      different_seed_changes_result
    ),
    "PASS",
    "FAIL"
  ),

  stringsAsFactors = FALSE
)


write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)


cat(
  "\n============================================================\n",
  "BAYESIAN VALIDATION SUMMARY\n",
  "============================================================\n",
  sep = ""
)


print(
  summary_table,
  row.names = FALSE
)


cat(
  "\nSame-seed max abs difference = ",
  format(
    same_seed_max_abs_diff,
    scientific = TRUE,
    digits = 16
  ),
  "\n",
  sep = ""
)


cat(
  "Different-seed max abs difference = ",
  format(
    different_seed_max_abs_diff,
    scientific = TRUE,
    digits = 10
  ),
  "\n",
  sep = ""
)


cat(
  "Different-seed mean abs difference = ",
  format(
    different_seed_mean_abs_diff,
    scientific = TRUE,
    digits = 10
  ),
  "\n",
  sep = ""
)


cat(
  "Different-seed RMSE = ",
  format(
    different_seed_rmse,
    scientific = TRUE,
    digits = 10
  ),
  "\n",
  sep = ""
)


overall_pass <-
  all(
    summary_table$status ==
      "PASS"
  )


cat(
  "\nOverall Bayesian validation status: ",
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

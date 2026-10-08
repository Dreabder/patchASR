#!/usr/bin/env Rscript

## ============================================================
## Direct validation:
## patchASR Bayesian adapter
## versus
## raw phytools::anc.Bayes posterior means
##
## Dataset:
##   BM/balanced/jump_16/rep_001
##
## Purpose:
##   For every ASR component created by patch_asr():
##
##   1. run make_asr_phytools_bayes()
##   2. independently run phytools::anc.Bayes()
##   3. manually discard burn-in
##   4. manually extract local internal-node columns
##   5. manually calculate posterior means
##   6. compare the two results exactly
##
## This validates the adapter implementation itself.
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
## Load current patchASR development version
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

component_summary_file <- file.path(
  case_dir,
  "bayesian_adapter_vs_raw_component_summary.csv"
)

nodewise_file <- file.path(
  case_dir,
  "bayesian_adapter_vs_raw_nodewise.csv"
)


## ------------------------------------------------------------
## Validation settings
##
## These intentionally match Step 4.6.2A.
## ------------------------------------------------------------

validation_ngen <- 1000L

validation_sample_freq <- 100L

validation_burnin_frac <- 0.20

base_seed <- 123L


## ------------------------------------------------------------
## Read tree and traits
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
## Verify shift
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
## Obtain the SAME component-seed helper used by the adapter.
##
## This is an internal helper and is accessed here only because
## this is a development validation script, not public user code.
## ============================================================

stable_component_seed_fun <-
  getFromNamespace(
    "stable_component_seed",
    "patchASR"
  )


## ------------------------------------------------------------
## Official patchASR Bayesian backend
## ------------------------------------------------------------

official_backend <-
  make_asr_phytools_bayes(
    ngen = validation_ngen,
    sample_freq =
      validation_sample_freq,
    burnin_frac =
      validation_burnin_frac,
    seed = base_seed
  )


## ============================================================
## Independent raw anc.Bayes calculation
## ============================================================

run_raw_bayes <- function(
    tree,
    states
) {

  if (tree$Nnode == 1L) {

    stop(
      "Raw Bayesian validation does not support ",
      "one-internal-node components."
    )
  }


  local_internal_ids <-
    as.character(
      ape::Ntip(tree) +
        seq_len(
          tree$Nnode
        )
    )


  component_seed <-
    stable_component_seed_fun(
      base_seed,
      tree
    )


  ## ----------------------------------------------------------
  ## Preserve caller RNG state during the independent raw call
  ## ----------------------------------------------------------

  had_random_seed <-
    exists(
      ".Random.seed",
      envir = globalenv(),
      inherits = FALSE
    )


  if (had_random_seed) {

    old_random_seed <-
      get(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
      )
  }


  on.exit(
    {

      if (had_random_seed) {

        assign(
          ".Random.seed",
          old_random_seed,
          envir = globalenv()
        )

      } else if (exists(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
      )) {

        rm(
          ".Random.seed",
          envir = globalenv()
        )
      }
    },
    add = TRUE
  )


  set.seed(
    component_seed
  )


  ## ----------------------------------------------------------
  ## Direct phytools call
  ## ----------------------------------------------------------

  fit <-
    phytools::anc.Bayes(
      tree = tree,
      x = states,
      ngen = validation_ngen,
      control = list(
        sample =
          validation_sample_freq
      )
    )


  if (is.null(
    fit$mcmc
  )) {

    stop(
      "Raw anc.Bayes() returned NULL mcmc."
    )
  }


  mcmc <-
    as.data.frame(
      fit$mcmc
    )


  ## ----------------------------------------------------------
  ## Manual burn-in
  ## ----------------------------------------------------------

  burnin_n <-
    floor(
      nrow(mcmc) *
        validation_burnin_frac
    )


  if (burnin_n >=
      nrow(mcmc)) {

    stop(
      "Burn-in removes every MCMC row."
    )
  }


  post_mcmc <-
    mcmc[
      (burnin_n + 1L):nrow(mcmc),
      ,
      drop = FALSE
    ]


  ## ----------------------------------------------------------
  ## Explicit local-node mapping
  ## ----------------------------------------------------------

  if (!all(
    local_internal_ids %in%
    colnames(post_mcmc)
  )) {

    stop(
      "Raw anc.Bayes output does not contain every ",
      "local internal-node column."
    )
  }


  values <-
    colMeans(
      post_mcmc[
        ,
        local_internal_ids,
        drop = FALSE
      ],
      na.rm = TRUE
    )


  values <-
    stats::setNames(
      as.numeric(values),
      tree$node.label
    )


  if (any(!is.finite(
    values
  ))) {

    stop(
      "Raw posterior means contain NA/Inf."
    )
  }


  list(
    estimates = values,
    component_seed = component_seed,
    n_mcmc_rows = nrow(mcmc),
    burnin_rows = burnin_n,
    retained_rows =
      nrow(post_mcmc)
  )
}


## ============================================================
## Validation backend
##
## patch_asr() constructs the REAL mother and patch components.
## Each component is independently evaluated by:
##
##   official patchASR adapter
##   versus
##   direct raw anc.Bayes calculation
## ============================================================

component_counter <- 0L

component_summaries <- list()

nodewise_results <- list()


validation_backend <- function(
    tree,
    states
) {

  component_counter <<-
    component_counter + 1L


  component_name <-
    paste0(
      "component_",
      component_counter
    )


  cat(
    "\n\n============================================================\n",
    "VALIDATING ",
    component_name,
    "\n",
    "============================================================\n",
    sep = ""
  )


  cat(
    "Ntip  = ",
    ape::Ntip(tree),
    "\n",
    sep = ""
  )


  cat(
    "Nnode = ",
    tree$Nnode,
    "\n",
    sep = ""
  )


  ## ----------------------------------------------------------
  ## Official adapter
  ## ----------------------------------------------------------

  adapter_values <-
    official_backend(
      tree,
      states
    )


  ## ----------------------------------------------------------
  ## Independent raw implementation
  ## ----------------------------------------------------------

  raw <-
    run_raw_bayes(
      tree,
      states
    )


  raw_values <-
    raw$estimates


  ## ----------------------------------------------------------
  ## Identity checks
  ## ----------------------------------------------------------

  if (!identical(
    names(adapter_values),
    names(raw_values)
  )) {

    stop(
      component_name,
      ": adapter and raw node identities differ."
    )
  }


  difference <-
    adapter_values -
    raw_values


  absolute_difference <-
    abs(
      difference
    )


  exact_identical <-
    identical(
      adapter_values,
      raw_values
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


  pass <-
    exact_identical &&
    max_abs_diff == 0


  cat(
    "Component seed = ",
    raw$component_seed,
    "\n",
    sep = ""
  )


  cat(
    "MCMC rows = ",
    raw$n_mcmc_rows,
    "\n",
    sep = ""
  )


  cat(
    "Burn-in rows = ",
    raw$burnin_rows,
    "\n",
    sep = ""
  )


  cat(
    "Retained rows = ",
    raw$retained_rows,
    "\n",
    sep = ""
  )


  cat(
    "Exact adapter/raw identity = ",
    exact_identical,
    "\n",
    sep = ""
  )


  cat(
    "Maximum absolute difference = ",
    format(
      max_abs_diff,
      scientific = TRUE,
      digits = 16
    ),
    "\n",
    sep = ""
  )


  cat(
    "Component status = ",
    if (pass) {
      "PASS"
    } else {
      "FAIL"
    },
    "\n",
    sep = ""
  )


  component_summaries[[component_counter]] <<-
    data.frame(
      component =
        component_name,

      n_tip =
        ape::Ntip(tree),

      n_node =
        tree$Nnode,

      component_seed =
        raw$component_seed,

      n_mcmc_rows =
        raw$n_mcmc_rows,

      burnin_rows =
        raw$burnin_rows,

      retained_rows =
        raw$retained_rows,

      exact_identical =
        exact_identical,

      max_abs_diff =
        max_abs_diff,

      mean_abs_diff =
        mean_abs_diff,

      rmse =
        rmse,

      status =
        if (pass) {
          "PASS"
        } else {
          "FAIL"
        },

      stringsAsFactors = FALSE
    )


  nodewise_results[[component_counter]] <<-
    data.frame(
      component =
        component_name,

      stable_node =
        names(
          adapter_values
        ),

      adapter_estimate =
        as.numeric(
          adapter_values
        ),

      raw_estimate =
        as.numeric(
          raw_values
        ),

      difference =
        as.numeric(
          difference
        ),

      abs_difference =
        as.numeric(
          absolute_difference
        ),

      stringsAsFactors = FALSE
    )


  ## Return the official adapter values to patch_asr().
  adapter_values
}


## ============================================================
## Run validation through REAL patch_asr component construction
## ============================================================

cat(
  "\nStarting raw-vs-adapter Bayesian validation...\n"
)


validation_fit <-
  patch_asr(
    tree = tree,
    states = states,
    shifts = shifts,
    asr_fun =
      validation_backend
  )


## ============================================================
## Combine reports
## ============================================================

component_summary <-
  do.call(
    rbind,
    component_summaries
  )

rownames(
  component_summary
) <- NULL


nodewise <-
  do.call(
    rbind,
    nodewise_results
  )

rownames(
  nodewise
) <- NULL


## ------------------------------------------------------------
## Final checks
## ------------------------------------------------------------

expected_component_calls <- 2L


component_count_ok <-
  component_counter ==
  expected_component_calls


all_component_pass <-
  all(
    component_summary$status ==
      "PASS"
  )


all_nodewise_zero <-
  all(
    nodewise$abs_difference ==
      0
  )


all_component_seeds_unique <-
  length(
    unique(
      component_summary$component_seed
    )
  ) ==
  nrow(
    component_summary
  )


overall_pass <-
  component_count_ok &&
  all_component_pass &&
  all_nodewise_zero &&
  all_component_seeds_unique


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

write.csv(
  component_summary,
  component_summary_file,
  row.names = FALSE
)


write.csv(
  nodewise,
  nodewise_file,
  row.names = FALSE
)


## ------------------------------------------------------------
## Print report
## ------------------------------------------------------------

cat(
  "\n\n============================================================\n",
  "BAYESIAN ADAPTER vs RAW anc.Bayes SUMMARY\n",
  "============================================================\n",
  sep = ""
)


print(
  component_summary,
  row.names = FALSE,
  digits = 16
)


cat(
  "\nExpected component calls: ",
  expected_component_calls,
  "\n",
  sep = ""
)


cat(
  "Observed component calls: ",
  component_counter,
  "\n",
  sep = ""
)


cat(
  "Component count correct: ",
  component_count_ok,
  "\n",
  sep = ""
)


cat(
  "All component seeds unique: ",
  all_component_seeds_unique,
  "\n",
  sep = ""
)


cat(
  "All node-wise differences exactly zero: ",
  all_nodewise_zero,
  "\n",
  sep = ""
)


cat(
  "\nOverall adapter-vs-raw status: ",
  if (overall_pass) {
    "PASS"
  } else {
    "FAIL"
  },
  "\n",
  sep = ""
)


cat(
  "\nComponent summary file:\n",
  component_summary_file,
  "\n",
  sep = ""
)


cat(
  "\nNode-wise comparison file:\n",
  nodewise_file,
  "\n",
  sep = ""
)

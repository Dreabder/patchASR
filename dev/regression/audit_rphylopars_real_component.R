#!/usr/bin/env Rscript

## ============================================================
## Audit Rphylopars node identities on REAL components passed
## by patchASR::patch_asr()
##
## IMPORTANT:
## This script is diagnostic only.
## The dummy return values from the audit backend must NOT be
## interpreted as ancestral-state estimates.
## ============================================================


## ------------------------------------------------------------
## Required packages
## ------------------------------------------------------------

if (!requireNamespace(
  "ape",
  quietly = TRUE
)) {
  stop(
    "Package 'ape' is required."
  )
}

if (!requireNamespace(
  "Rphylopars",
  quietly = TRUE
)) {
  stop(
    "Package 'Rphylopars' is required."
  )
}

if (!requireNamespace(
  "pkgload",
  quietly = TRUE
)) {
  stop(
    "Package 'pkgload' is required."
  )
}


## ------------------------------------------------------------
## Load current development version of patchASR
## ------------------------------------------------------------

pkgload::load_all(
  ".",
  quiet = TRUE
)


## ------------------------------------------------------------
## Input files
## ------------------------------------------------------------

source(
  file.path(
    "dev",
    "regression",
    "regression_paths.R"
  )
)

case_dir <- fixture_dir

tree_file <- file.path(
  case_dir,
  "tree.nwk"
)

trait_file <- file.path(
  case_dir,
  "node_traits.csv"
)


tree <- ape::read.tree(
  tree_file
)

dat <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## ------------------------------------------------------------
## Helper: extract one unique metadata value
## ------------------------------------------------------------

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


## ------------------------------------------------------------
## Shift edge
## ------------------------------------------------------------

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


cat(
  "Original tree Ntip  = ",
  ape::Ntip(tree),
  "\n",
  sep = ""
)

cat(
  "Original tree Nnode = ",
  tree$Nnode,
  "\n",
  sep = ""
)

cat(
  "Shift edge          = ",
  I,
  " -> ",
  J,
  "\n",
  sep = ""
)


## Verify that I -> J is exactly one edge
edge_match <- which(
  tree$edge[, 1] == I &
    tree$edge[, 2] == J
)

if (length(edge_match) != 1L) {

  stop(
    "I -> J does not identify exactly one edge ",
    "in the input tree."
  )
}


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
    ape::Ntip(tree)) {

  stop(
    "Incorrect number of tip states."
  )
}

if (any(!is.finite(
  states
))) {

  stop(
    "Tip states contain NA/Inf."
  )
}

if (!identical(
  names(states),
  tree$tip.label
)) {

  stop(
    "Tip-state names are not in exact tree tip order."
  )
}


shifts <- data.frame(
  parent = I,
  child = J
)


## ============================================================
## Audit backend
##
## patch_asr() itself constructs mother/patch trees and then
## passes each real component into this function.
##
## We run Rphylopars only to inspect anc_recon identity.
##
## We then return finite dummy values carrying the correct
## stable node names so that patch_asr() can continue.
## ============================================================

audit_call_number <- 0L


audit_rphylopars_backend <- function(
    tree,
    states
) {

  audit_call_number <<-
    audit_call_number + 1L


  cat(
    "\n\n",
    "============================================================\n",
    "ASR COMPONENT CALL: ",
    audit_call_number,
    "\n",
    "============================================================\n",
    sep = ""
  )


  ## ----------------------------------------------------------
  ## Basic component information
  ## ----------------------------------------------------------

  n_tip <- ape::Ntip(
    tree
  )

  n_node <- tree$Nnode


  cat(
    "Ntip  = ",
    n_tip,
    "\n",
    sep = ""
  )

  cat(
    "Nnode = ",
    n_node,
    "\n",
    sep = ""
  )


  cat(
    "\nFirst component tip labels:\n"
  )

  print(
    utils::head(
      tree$tip.label,
      20L
    )
  )


  cat(
    "\nFirst component stable node labels:\n"
  )

  print(
    utils::head(
      tree$node.label,
      30L
    )
  )


  cat(
    "\nLast component stable node labels:\n"
  )

  print(
    utils::tail(
      tree$node.label,
      30L
    )
  )


  ## ----------------------------------------------------------
  ## Three candidate identity systems
  ## ----------------------------------------------------------

  stable_ids <-
    tree$node.label


  stripped_stable_ids <- sub(
    "^node_",
    "",
    stable_ids
  )


  local_ape_ids <- as.character(
    n_tip +
      seq_len(
        n_node
      )
  )


  cat(
    "\nFirst local ape internal IDs:\n"
  )

  print(
    utils::head(
      local_ape_ids,
      30L
    )
  )


  ## ----------------------------------------------------------
  ## Rphylopars input
  ## ----------------------------------------------------------

  trait_data <- data.frame(
    species = names(states),
    trait = as.numeric(states),
    stringsAsFactors = FALSE
  )


  ## ----------------------------------------------------------
  ## Run Rphylopars
  ## ----------------------------------------------------------

  fit <- Rphylopars::phylopars(
    trait_data = trait_data,
    tree = tree,
    model = "BM",
    pheno_error = FALSE,
    REML = TRUE
  )


  if (is.null(
    fit$anc_recon
  )) {

    stop(
      "Rphylopars returned NULL anc_recon."
    )
  }


  anc <- as.data.frame(
    fit$anc_recon
  )

  rn <- rownames(
    anc
  )


  cat(
    "\nanc_recon dimensions:\n"
  )

  print(
    dim(anc)
  )


  cat(
    "\nFirst anc_recon row names:\n"
  )

  print(
    utils::head(
      rn,
      40L
    )
  )


  cat(
    "\nLast anc_recon row names:\n"
  )

  print(
    utils::tail(
      rn,
      40L
    )
  )


  ## ----------------------------------------------------------
  ## Identity checks
  ## ----------------------------------------------------------

  cat(
    "\n===== identity checks =====\n"
  )


  exact_stable_ok <-
    all(
      stable_ids %in% rn
    )


  stripped_stable_ok <-
    all(
      stripped_stable_ids %in% rn
    )


  local_ape_ok <-
    all(
      local_ape_ids %in% rn
    )


  cat(
    "All exact stable IDs present: ",
    exact_stable_ok,
    "\n",
    sep = ""
  )


  cat(
    "All stripped stable IDs present: ",
    stripped_stable_ok,
    "\n",
    sep = ""
  )


  cat(
    "All local ape IDs present: ",
    local_ape_ok,
    "\n",
    sep = ""
  )


  cat(
    "Exact stable IDs matched: ",
    sum(
      stable_ids %in% rn
    ),
    " / ",
    length(stable_ids),
    "\n",
    sep = ""
  )


  cat(
    "Stripped stable IDs matched: ",
    sum(
      stripped_stable_ids %in% rn
    ),
    " / ",
    length(stripped_stable_ids),
    "\n",
    sep = ""
  )


  cat(
    "Local ape IDs matched: ",
    sum(
      local_ape_ids %in% rn
    ),
    " / ",
    length(local_ape_ids),
    "\n",
    sep = ""
  )


  ## ----------------------------------------------------------
  ## Example identity table
  ## ----------------------------------------------------------

  cat(
    "\nExample mapping table:\n"
  )


  mapping <- data.frame(
    stable_id = stable_ids,
    stripped_stable_id =
      stripped_stable_ids,
    local_ape_id =
      local_ape_ids,
    stringsAsFactors = FALSE
  )


  print(
    utils::head(
      mapping,
      30L
    )
  )


  ## ----------------------------------------------------------
  ## IMPORTANT:
  ## Return diagnostic placeholder values.
  ##
  ## These are NOT ancestral reconstructions.
  ## They only allow patch_asr() to finish its control flow.
  ## ----------------------------------------------------------

  if (n_node == 0L) {

    return(
      stats::setNames(
        numeric(0),
        character(0)
      )
    )
  }


  stats::setNames(
    rep(
      0,
      n_node
    ),
    stable_ids
  )
}


## ============================================================
## Let patch_asr() construct the REAL components
## ============================================================

cat(
  "\n\nStarting patch_asr component audit...\n"
)


audit_result <- patchASR::patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = audit_rphylopars_backend
)


cat(
  "\n\n============================================================\n",
  "REAL-COMPONENT RPHYLOPARS AUDIT FINISHED\n",
  "Number of ASR component calls: ",
  audit_call_number,
  "\n",
  "============================================================\n",
  sep = ""
)

#!/usr/bin/env Rscript

## ============================================================
## phyloPatch deterministic regression
##
## Same exact dataset and shift as the legacy reference.
## ============================================================


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
  stop("Package 'pkgload' is required for this development regression.")
}


## ------------------------------------------------------------
## Load current development version of phyloPatch
## ------------------------------------------------------------

pkgload::load_all(
  ".",
  quiet = TRUE
)


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

output_file <- file.path(
  case_dir,
  "phylopatch_patch_4methods.csv"
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


## ------------------------------------------------------------
## Verify node identity
## ------------------------------------------------------------

n_tip <- ape::Ntip(
  tree
)

expected_internal_ids <- as.character(
  n_tip +
    seq_len(
      tree$Nnode
    )
)

if (!identical(
  tree$node.label,
  expected_internal_ids
)) {

  stop(
    "Internal Newick labels do not exactly match current ape ",
    "internal-node numbers. Stop regression rather than assuming ",
    "node identity."
  )
}


## ------------------------------------------------------------
## Helper for one metadata value
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
      " must contain exactly one unique value."
    )
  }

  values[[1]]
}


## ------------------------------------------------------------
## Derive the shift directly from node_traits.csv
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

edge_rows <- which(
  tree$edge[, 1] == I &
    tree$edge[, 2] == J
)

if (length(edge_rows) != 1L) {

  stop(
    "Metadata I -> J does not identify exactly one tree edge."
  )
}

shifts <- data.frame(
  parent = I,
  child = J
)


## ------------------------------------------------------------
## Extract observed tip states
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
  "Tree tips:",
  ape::Ntip(tree),
  "\n"
)

cat(
  "Internal nodes:",
  tree$Nnode,
  "\n"
)

cat(
  "Shift:",
  I,
  "->",
  J,
  "\n"
)


## ------------------------------------------------------------
## Run phyloPatch
## ------------------------------------------------------------

fit_pic <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_pic
)

fit_ml <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_ml_bm
)

fit_gls <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = asr_ape_gls_bm
)

rph_backend <- make_asr_rphylopars_bm(
  REML = TRUE
)

fit_rph <- patch_asr(
  tree = tree,
  states = states,
  shifts = shifts,
  asr_fun = rph_backend
)


## ------------------------------------------------------------
## Extract estimates
## ------------------------------------------------------------

extract_one <- function(
    fit,
    column_name
) {

  x <- fit$ancestral_states[
    ,
    c(
      "original_node",
      "estimate"
    ),
    drop = FALSE
  ]

  names(x) <- c(
    "node_id",
    column_name
  )

  x$node_id <- as.integer(
    x$node_id
  )

  x
}


pic_df <- extract_one(
  fit_pic,
  "prediction_PIC_style"
)

ml_df <- extract_one(
  fit_ml,
  "prediction_ML"
)

gls_df <- extract_one(
  fit_gls,
  "prediction_GLS"
)

rph_df <- extract_one(
  fit_rph,
  "prediction_Rphylopars"
)


result <- Reduce(
  function(
    x,
    y
  ) {

    merge(
      x,
      y,
      by = "node_id",
      all = TRUE,
      sort = FALSE
    )
  },
  list(
    pic_df,
    ml_df,
    gls_df,
    rph_df
  )
)

result <- result[
  order(
    result$node_id
  ),
  ,
  drop = FALSE
]

rownames(result) <- NULL


## ------------------------------------------------------------
## Safety checks
## ------------------------------------------------------------

if (nrow(result) !=
    tree$Nnode) {

  stop(
    "phyloPatch output contains an unexpected number of nodes."
  )
}

prediction_columns <- setdiff(
  colnames(result),
  "node_id"
)

if (any(
  !is.finite(
    as.matrix(
      result[
        ,
        prediction_columns,
        drop = FALSE
      ]
    )
  )
)) {

  stop(
    "phyloPatch regression output contains NA/Inf."
  )
}


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

write.csv(
  result,
  output_file,
  row.names = FALSE
)

cat(
  "phyloPatch deterministic regression completed.\n"
)

cat(
  "Rows:",
  nrow(result),
  "\n"
)

cat(
  "Output:",
  output_file,
  "\n"
)

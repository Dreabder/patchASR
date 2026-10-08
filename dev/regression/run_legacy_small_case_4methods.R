#!/usr/bin/env Rscript

## ============================================================
## Legacy deterministic patch-ASR reference
##
## Dataset:
## BM/balanced/jump_16/rep_001
##
## Methods:
##   PIC
##   ML-BM
##   GLS-BM
##   Rphylopars-BM
##
## Bayesian is intentionally excluded from this regression step.
## ============================================================


legacy_script <- Sys.getenv(
  "PATCHASR_LEGACY_SCRIPT",
  unset = ""
)

if (!nzchar(legacy_script)) {
  stop(
    paste0(
      "To regenerate the frozen legacy reference, set the ",
      "PATCHASR_LEGACY_SCRIPT environment variable to the path ",
      "of the original research implementation."
    )
  )
}

if (!file.exists(legacy_script)) {
  stop(
    "Legacy implementation file does not exist: ",
    legacy_script
  )
}

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

output_file <- file.path(
  results_dir,
  "legacy_patch_4methods_regenerated.csv"
)


## ------------------------------------------------------------
## Load the actual legacy implementation
## ------------------------------------------------------------

source(
  legacy_script
)


## ------------------------------------------------------------
## Read the exact regression dataset
## ------------------------------------------------------------

dat <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

my_tree <- ape::read.tree(
  tree_file
)


## ------------------------------------------------------------
## Prepare and partition using the legacy implementation
## ------------------------------------------------------------

obj <- prepare_patch_ase_input(
  shifted_df = dat,
  my_tree = my_tree
)

pruned <- execute_patch_prune(
  obj
)


## ------------------------------------------------------------
## PIC
## ------------------------------------------------------------

res_pic <- run_ase_pic_local_patch_single(
  obj,
  pruned
)

res_pic <- build_patch_final_output(
  res_pic,
  obj
)


## ------------------------------------------------------------
## ML-BM
## ------------------------------------------------------------

res_ml <- run_ase_ml_global_patch_single(
  obj,
  pruned
)

res_ml <- build_patch_final_output(
  res_ml,
  obj
)


## ------------------------------------------------------------
## GLS-BM
## ------------------------------------------------------------

res_gls <- run_ase_gls_global_patch_single(
  obj,
  pruned
)

res_gls <- build_patch_final_output(
  res_gls,
  obj
)


## ------------------------------------------------------------
## Rphylopars-BM
## ------------------------------------------------------------

res_rph <- run_ase_rphylopars_global_patch_single(
  obj,
  pruned,
  REML = TRUE
)

res_rph <- build_patch_final_output(
  res_rph,
  obj
)


## ------------------------------------------------------------
## Extract one prediction column per method
## ------------------------------------------------------------

pic_df <- extract_patch_method_df(
  res_pic,
  "PIC_style_local"
)

ml_df <- extract_patch_method_df(
  res_ml,
  "ML_global"
)

gls_df <- extract_patch_method_df(
  res_gls,
  "GLS_global"
)

rph_df <- extract_patch_method_df(
  res_rph,
  "Rphylopars_like_global"
)


## ------------------------------------------------------------
## Merge by original node ID
## ------------------------------------------------------------

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

result$node_id <- as.integer(
  result$node_id
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

expected_n <- my_tree$Nnode

if (nrow(result) !=
    expected_n) {

  stop(
    "Legacy result contains ",
    nrow(result),
    " internal nodes; expected ",
    expected_n,
    "."
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
    "Legacy regression result contains NA/Inf."
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
  "Legacy deterministic regression completed.\n"
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

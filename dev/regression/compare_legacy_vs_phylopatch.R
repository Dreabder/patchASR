#!/usr/bin/env Rscript

## ============================================================
## Numerical regression:
## legacy patch implementation vs phyloPatch
##
## Dataset:
##   BM/balanced/jump_16/rep_001
##
## Deterministic methods:
##   PIC
##   ML-BM
##   GLS-BM
##   Rphylopars-BM
##
## Bayesian is intentionally excluded here because it requires
## a separate stochastic/RNG regression analysis.
## ============================================================


## ------------------------------------------------------------
## Settings
## ------------------------------------------------------------

case_dir <- paste0(
  "/home/chengyq/work/ASE/data/jump_simu/simu/",
  "BM/balanced/jump_16/rep_001"
)

legacy_file <- file.path(
  case_dir,
  "legacy_patch_4methods.csv"
)

new_file <- file.path(
  case_dir,
  "phylopatch_patch_4methods.csv"
)

trait_file <- file.path(
  case_dir,
  "node_traits.csv"
)

summary_file <- file.path(
  case_dir,
  "legacy_vs_phylopatch_4methods_summary.csv"
)

nodewise_file <- file.path(
  case_dir,
  "legacy_vs_phylopatch_4methods_nodewise.csv"
)


## Initial numerical regression tolerance.
tolerance <- 1e-8


## ------------------------------------------------------------
## Input existence checks
## ------------------------------------------------------------

required_files <- c(
  legacy_file,
  new_file,
  trait_file
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0L) {

  stop(
    "Required regression file(s) missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


## ------------------------------------------------------------
## Read results
## ------------------------------------------------------------

legacy <- read.csv(
  legacy_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

new <- read.csv(
  new_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

traits <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Legacy rows:     ",
  nrow(legacy),
  "\n",
  sep = ""
)

cat(
  "phyloPatch rows: ",
  nrow(new),
  "\n",
  sep = ""
)


## ------------------------------------------------------------
## Basic node-ID checks
## ------------------------------------------------------------

if (!"node_id" %in%
    names(legacy)) {

  stop(
    "Legacy result has no `node_id` column."
  )
}

if (!"node_id" %in%
    names(new)) {

  stop(
    "phyloPatch result has no `node_id` column."
  )
}


legacy$node_id <- as.integer(
  legacy$node_id
)

new$node_id <- as.integer(
  new$node_id
)


if (anyNA(
  legacy$node_id
) ||
anyNA(
  new$node_id
)) {

  stop(
    "At least one node_id could not be converted to integer."
  )
}


if (anyDuplicated(
  legacy$node_id
)) {

  stop(
    "Legacy result contains duplicated node_id values."
  )
}

if (anyDuplicated(
  new$node_id
)) {

  stop(
    "phyloPatch result contains duplicated node_id values."
  )
}


if (!setequal(
  legacy$node_id,
  new$node_id
)) {

  missing_from_new <- setdiff(
    legacy$node_id,
    new$node_id
  )

  missing_from_legacy <- setdiff(
    new$node_id,
    legacy$node_id
  )

  stop(
    "Legacy and phyloPatch node sets differ.\n",
    "Missing from phyloPatch: ",
    paste(
      missing_from_new,
      collapse = ", "
    ),
    "\n",
    "Missing from legacy: ",
    paste(
      missing_from_legacy,
      collapse = ", "
    )
  )
}


cat(
  "Node-ID sets match exactly: TRUE\n"
)


## ------------------------------------------------------------
## Find the actual legacy column names safely
## ------------------------------------------------------------

find_one_column <- function(
    df,
    candidates,
    method_name
) {

  hits <- candidates[
    candidates %in%
      names(df)
  ]

  if (length(hits) == 0L) {

    stop(
      "Could not identify the legacy ",
      method_name,
      " column.\n",
      "Candidate names were: ",
      paste(
        candidates,
        collapse = ", "
      ),
      "\nAvailable columns are: ",
      paste(
        names(df),
        collapse = ", "
      )
    )
  }

  if (length(hits) > 1L) {

    stop(
      "Multiple candidate columns were found for ",
      method_name,
      ": ",
      paste(
        hits,
        collapse = ", "
      ),
      ". Resolve this ambiguity before comparison."
    )
  }

  hits[[1]]
}


legacy_pic_col <- find_one_column(
  legacy,
  c(
    "PIC_style_local",
    "prediction_PIC_style",
    "prediction_PIC_style_local",
    "PIC_style"
  ),
  "PIC"
)

legacy_ml_col <- find_one_column(
  legacy,
  c(
    "ML_global",
    "prediction_ML",
    "prediction_ML_global",
    "ML"
  ),
  "ML"
)

legacy_gls_col <- find_one_column(
  legacy,
  c(
    "GLS_global",
    "prediction_GLS",
    "prediction_GLS_global",
    "GLS"
  ),
  "GLS"
)

legacy_rph_col <- find_one_column(
  legacy,
  c(
    "Rphylopars_like_global",
    "prediction_Rphylopars",
    "prediction_Rphylopars_like_global",
    "Rphylopars_global",
    "Rphylopars"
  ),
  "Rphylopars"
)


## ------------------------------------------------------------
## Required phyloPatch columns
## ------------------------------------------------------------

new_columns <- c(
  PIC = "prediction_PIC_style",
  ML = "prediction_ML",
  GLS = "prediction_GLS",
  Rphylopars = "prediction_Rphylopars"
)

missing_new_columns <- setdiff(
  unname(new_columns),
  names(new)
)

if (length(missing_new_columns) > 0L) {

  stop(
    "phyloPatch output is missing expected column(s): ",
    paste(
      missing_new_columns,
      collapse = ", "
    )
  )
}


legacy_columns <- c(
  PIC = legacy_pic_col,
  ML = legacy_ml_col,
  GLS = legacy_gls_col,
  Rphylopars = legacy_rph_col
)


cat(
  "\nResolved legacy columns:\n"
)

print(
  legacy_columns
)


## ------------------------------------------------------------
## Sort both tables by original node ID
## ------------------------------------------------------------

legacy <- legacy[
  order(
    legacy$node_id
  ),
  ,
  drop = FALSE
]

new <- new[
  order(
    new$node_id
  ),
  ,
  drop = FALSE
]


if (!identical(
  legacy$node_id,
  new$node_id
)) {

  stop(
    "Node IDs do not align after sorting."
  )
}


## ------------------------------------------------------------
## Identify the boundary node I
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


boundary_I <- as.integer(
  get_one(
    traits$I,
    "I"
  )
)


cat(
  "\nBoundary node I: ",
  boundary_I,
  "\n",
  sep = ""
)


## ------------------------------------------------------------
## Build node-wise comparison table
## ------------------------------------------------------------

nodewise <- data.frame(
  node_id = legacy$node_id,
  is_boundary_I =
    legacy$node_id == boundary_I,
  stringsAsFactors = FALSE
)


summary_list <- list()


for (method in
     names(legacy_columns)) {

  legacy_values <- as.numeric(
    legacy[[legacy_columns[[method]]]]
  )

  new_values <- as.numeric(
    new[[new_columns[[method]]]]
  )


  if (any(!is.finite(
    legacy_values
  ))) {

    stop(
      method,
      ": legacy result contains non-finite values."
    )
  }


  if (any(!is.finite(
    new_values
  ))) {

    stop(
      method,
      ": phyloPatch result contains non-finite values."
    )
  }


  difference <-
    new_values -
    legacy_values

  absolute_difference <-
    abs(
      difference
    )


  max_index <- which.max(
    absolute_difference
  )


  boundary_index <- match(
    boundary_I,
    legacy$node_id
  )


  summary_list[[method]] <- data.frame(
    method = method,
    n_nodes =
      length(
        difference
      ),
    max_abs_diff =
      max(
        absolute_difference
      ),
    node_with_max_diff =
      legacy$node_id[
        max_index
      ],
    mean_abs_diff =
      mean(
        absolute_difference
      ),
    rmse =
      sqrt(
        mean(
          difference^2
        )
      ),
    boundary_I_abs_diff =
      absolute_difference[
        boundary_index
      ],
    n_over_tolerance =
      sum(
        absolute_difference >
          tolerance
      ),
    tolerance =
      tolerance,
    status =
      if (
        all(
          absolute_difference <=
          tolerance
        )
      ) {
        "PASS"
      } else {
        "FAIL"
      },
    stringsAsFactors = FALSE
  )


  prefix <- tolower(
    method
  )


  nodewise[[
      paste0(
        "legacy_",
        prefix
      )]] <- legacy_values


  nodewise[[
      paste0(
        "phylopatch_",
        prefix
      )]] <- new_values


  nodewise[[paste0(
        "diff_",
        prefix
      )]] <- difference


  nodewise[[
      paste0(
        "abs_diff_",
        prefix
      )]] <- absolute_difference
}


summary_table <- do.call(
  rbind,
  summary_list
)

rownames(
  summary_table
) <- NULL


## ------------------------------------------------------------
## Save results
## ------------------------------------------------------------

write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)

write.csv(
  nodewise,
  nodewise_file,
  row.names = FALSE
)


## ------------------------------------------------------------
## Print concise regression report
## ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "LEGACY vs phyloPatch NUMERICAL REGRESSION\n",
  "============================================================\n",
  sep = ""
)

print(
  summary_table,
  row.names = FALSE,
  digits = 12
)


cat(
  "\nOverall status: ",
  if (
    all(
      summary_table$status ==
      "PASS"
    )
  ) {
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

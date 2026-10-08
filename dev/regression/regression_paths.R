## ============================================================
## Shared paths for self-contained regression validation
##
## Run regression scripts from the patchASR repository root.
## ============================================================

repo_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

description_file <- file.path(
  repo_root,
  "DESCRIPTION"
)

if (!file.exists(description_file)) {
  stop(
    "Regression scripts must be run from the patchASR repository root."
  )
}

description_text <- readLines(
  description_file,
  warn = FALSE
)

if (!any(description_text == "Package: patchASR")) {
  stop(
    "Current working directory is not the patchASR repository root."
  )
}

fixture_dir <- file.path(
  repo_root,
  "dev",
  "regression",
  "fixtures",
  "BM_balanced_jump16_rep001"
)

results_dir <- file.path(
  repo_root,
  "dev",
  "regression",
  "results"
)

if (!dir.exists(fixture_dir)) {
  stop(
    "Regression fixture directory is missing: ",
    fixture_dir
  )
}

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

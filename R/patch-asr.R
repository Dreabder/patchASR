# User-facing orchestration for patchASR.


#' Patch-aware ancestral-state reconstruction
#'
#' Perform ancestral-state reconstruction after partitioning a rooted
#' phylogeny across one or more designated shift edges.
#'
#' `patch_asr()` validates the input data, freezes stable node identities,
#' constructs a patch plan, partitions the phylogeny, performs
#' component-wise ancestral-state reconstruction, reconstructs boundary
#' nodes, and finally merges all estimates back to the internal nodes of
#' the original phylogeny.
#'
#' @param tree A rooted, fully bifurcating phylogenetic tree of class
#'   `phylo` with finite, strictly positive branch lengths and unique
#'   terminal labels.
#'
#' @param states A named numeric vector containing one finite observed
#'   continuous trait value for every terminal taxon in `tree`.
#'
#' @param shifts A data frame containing `parent` and `child` columns.
#'   Each row identifies one directed shift edge using the `ape` node
#'   numbers of the original input `tree`.
#'
#' @param asr_fun A function with interface `asr_fun(tree, states)` that
#'   returns one finite ancestral-state estimate for every internal node
#'   of the supplied component. The returned value must be a named numeric
#'   vector whose names are the stable internal-node identifiers carried
#'   by the supplied component tree.
#'
#' @param residual_tolerance A finite positive numeric scalar specifying
#'   the maximum allowed absolute residual for the boundary linear system.
#'   The default is `1e-10`.
#'
#' @return An object of class `patch_asr`. The object contains the final
#'   ancestral-state table, patch plan, boundary estimates, component-level
#'   estimates, stable node map, original call, and diagnostic information.
#'
#' @importFrom stats setNames
#' @export
patch_asr <- function(
    tree,
    states,
    shifts,
    asr_fun,
    residual_tolerance = 1e-10
) {

  ## ----------------------------------------------------------
  ## 1. Validate user inputs
  ## ----------------------------------------------------------

  validated <- validate_patch_inputs(
    tree = tree,
    states = states,
    shifts = shifts,
    asr_fun = asr_fun
  )

  if (!is.numeric(
    residual_tolerance
  ) ||
  length(
    residual_tolerance
  ) != 1L ||
  is.na(
    residual_tolerance
  ) ||
  !is.finite(
    residual_tolerance
  ) ||
  residual_tolerance <= 0) {
    stop(
      "`residual_tolerance` must be one finite positive number.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## 2. Freeze original-tree node identity
  ## ----------------------------------------------------------

  frozen <- freeze_node_identity(
    validated$tree
  )

  ## ----------------------------------------------------------
  ## 3. Construct patch plan
  ## ----------------------------------------------------------

  patch_plan <- build_patch_plan(
    tree = validated$tree,
    node_map = frozen$node_map,
    shifts = validated$shifts
  )

  ## ----------------------------------------------------------
  ## 4. Partition the phylogeny
  ## ----------------------------------------------------------

  partition <- partition_phylogeny(
    tree = validated$tree,
    node_map = frozen$node_map,
    patch_plan = patch_plan
  )

  ## ----------------------------------------------------------
  ## 5. Perform component-level ASR
  ## ----------------------------------------------------------

  component_asr <- run_partition_asr(
    partition = partition,
    states = validated$states,
    asr_fun = validated$asr_fun
  )

  ## ----------------------------------------------------------
  ## 6. Reconstruct boundary nodes
  ## ----------------------------------------------------------

  boundary_result <- run_boundary_reconstruction(
    patch_plan = patch_plan,
    node_map = frozen$node_map,
    states = validated$states,
    component_asr = component_asr,
    residual_tolerance = residual_tolerance
  )

  ## ----------------------------------------------------------
  ## 7. Merge all original internal-node estimates
  ## ----------------------------------------------------------

  ancestral_states <- merge_ancestral_results(
    node_map = frozen$node_map,
    component_asr = component_asr,
    boundary_result = boundary_result
  )

  ## ----------------------------------------------------------
  ## 8. Construct final S3 object
  ## ----------------------------------------------------------

  result <- list(
    ancestral_states =
      ancestral_states,

    patch_plan =
      patch_plan,

    boundary_states =
      boundary_result$estimate_table,

    component_states =
      component_asr$estimate_table,

    node_map =
      frozen$node_map,

    call =
      match.call(),

    diagnostics =
      list(
        partition =
          partition,

        component_asr =
          component_asr,

        boundary_result =
          boundary_result,

        boundary_residual =
          boundary_result$residual
      )
  )

  class(result) <- "patch_asr"

  result
}


#' Print a patch-aware ancestral-state reconstruction
#'
#' @param x An object of class `patch_asr`.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.patch_asr <- function(
    x,
    ...
) {

  if (!inherits(
    x,
    "patch_asr"
  )) {
    stop(
      "`x` must inherit from class `patch_asr`.",
      call. = FALSE
    )
  }

  n_internal <-
    if (is.data.frame(
      x$ancestral_states
    )) {
      nrow(
        x$ancestral_states
      )
    } else {
      NA_integer_
    }

  n_patches <-
    if (inherits(
      x$patch_plan,
      "patchASR_patch_plan"
    )) {
      x$patch_plan$n_patches
    } else {
      NA_integer_
    }

  n_boundaries <-
    if (is.data.frame(
      x$boundary_states
    )) {
      nrow(
        x$boundary_states
      )
    } else {
      NA_integer_
    }

  cat(
    "<patch_asr>\n",
    sep = ""
  )

  cat(
    "  Internal nodes: ",
    n_internal,
    "\n",
    sep = ""
  )

  cat(
    "  Patches: ",
    n_patches,
    "\n",
    sep = ""
  )

  cat(
    "  Boundary nodes: ",
    n_boundaries,
    "\n",
    sep = ""
  )

  if (!is.null(
    x$diagnostics$
    boundary_residual
  )) {

    cat(
      "  Boundary residual: ",
      format(
        x$diagnostics$
          boundary_residual,
        scientific = TRUE
      ),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

# Boundary-state reconstruction infrastructure for patchASR.
#
# Boundary nodes are reconstructed jointly from original-tree branch
# geometry and the unified pool of already known component and terminal
# states.


#' Build the unified known-state pool
#'
#' @param node_map Stable original-tree node map.
#' @param states Named numeric observed states for all original tips.
#' @param component_asr Result produced by `run_partition_asr()`.
#' @param boundary_ids Stable identifiers of unresolved boundary nodes.
#'
#' @return A named numeric vector containing all known non-boundary states.
#'
#' @noRd
build_unified_known_states <- function(
    node_map,
    states,
    component_asr,
    boundary_ids
) {

  required_columns <- c(
    "original_node",
    "node_type",
    "original_label",
    "stable_id"
  )

  if (!is.data.frame(node_map) ||
      !all(
        required_columns %in%
        names(node_map)
      )) {
    stop(
      "`node_map` has an invalid structure.",
      call. = FALSE
    )
  }

  if (anyNA(node_map$stable_id) ||
      any(node_map$stable_id == "") ||
      anyDuplicated(node_map$stable_id)) {
    stop(
      "`node_map` contains invalid stable node identifiers.",
      call. = FALSE
    )
  }

  if (!is.numeric(states) ||
      is.null(names(states))) {
    stop(
      "`states` must be a named numeric vector.",
      call. = FALSE
    )
  }

  if (anyNA(names(states)) ||
      any(names(states) == "") ||
      anyDuplicated(names(states))) {
    stop(
      "`states` must have unique, non-empty names.",
      call. = FALSE
    )
  }

  if (any(!is.finite(states))) {
    stop(
      "All observed terminal states must be finite.",
      call. = FALSE
    )
  }

  if (!inherits(
    component_asr,
    "patchASR_component_asr"
  )) {
    stop(
      "`component_asr` must be a valid patchASR component-ASR object.",
      call. = FALSE
    )
  }

  if (!is.character(boundary_ids) ||
      anyNA(boundary_ids) ||
      any(boundary_ids == "") ||
      anyDuplicated(boundary_ids)) {
    stop(
      "`boundary_ids` must contain unique, non-empty stable identifiers.",
      call. = FALSE
    )
  }

  if (!all(
    boundary_ids %in%
    node_map$stable_id
  )) {
    stop(
      "At least one boundary identifier is absent from `node_map`.",
      call. = FALSE
    )
  }

  boundary_rows <- match(
    boundary_ids,
    node_map$stable_id
  )

  if (any(
    node_map$node_type[
      boundary_rows
    ] != "internal"
  )) {
    stop(
      "Every boundary identifier must refer to an original internal node.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Observed terminal states
  ## ----------------------------------------------------------

  tip_rows <- which(
    node_map$node_type ==
      "tip"
  )

  tip_labels <-
    node_map$original_label[
      tip_rows
    ]

  missing_tips <- setdiff(
    tip_labels,
    names(states)
  )

  unknown_tips <- setdiff(
    names(states),
    tip_labels
  )

  if (length(missing_tips) > 0L) {
    stop(
      "Missing observed terminal state(s): ",
      paste(
        missing_tips,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (length(unknown_tips) > 0L) {
    stop(
      "Unknown terminal state label(s): ",
      paste(
        unknown_tips,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (length(states) !=
      length(tip_rows)) {
    stop(
      "`states` must contain exactly one value for every original tip.",
      call. = FALSE
    )
  }

  tip_states <- setNames(
    as.numeric(
      states[
        tip_labels
      ]
    ),
    node_map$stable_id[
      tip_rows
    ]
  )

  ## ----------------------------------------------------------
  ## Component-derived internal states
  ## ----------------------------------------------------------

  internal_states <-
    component_asr$all_internal_states

  if (!is.numeric(internal_states) ||
      is.null(names(internal_states))) {
    stop(
      "`component_asr$all_internal_states` must be a named numeric vector.",
      call. = FALSE
    )
  }

  if (anyNA(names(internal_states)) ||
      any(names(internal_states) == "") ||
      anyDuplicated(names(internal_states))) {
    stop(
      "Component ASR states contain invalid node identifiers.",
      call. = FALSE
    )
  }

  if (any(!is.finite(
    internal_states
  ))) {
    stop(
      "Component ASR states must all be finite.",
      call. = FALSE
    )
  }

  if (any(
    names(internal_states) %in%
    boundary_ids
  )) {
    stop(
      "Unresolved boundary nodes must not occur in component ASR states.",
      call. = FALSE
    )
  }

  original_internal_ids <-
    node_map$stable_id[
      node_map$node_type ==
        "internal"
    ]

  expected_internal_ids <- setdiff(
    original_internal_ids,
    boundary_ids
  )

  if (!setequal(
    names(internal_states),
    expected_internal_ids
  )) {
    stop(
      "Component ASR states do not match the expected non-boundary ",
      "internal-node set.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Unified pool
  ## ----------------------------------------------------------

  known_states <- c(
    tip_states,
    internal_states
  )

  if (anyDuplicated(
    names(known_states)
  )) {
    stop(
      "Unified known-state pool contains duplicated node identities.",
      call. = FALSE
    )
  }

  if (any(
    boundary_ids %in%
    names(known_states)
  )) {
    stop(
      "Unresolved boundary nodes must not be inserted into the ",
      "known-state pool.",
      call. = FALSE
    )
  }

  expected_known_ids <- c(
    node_map$stable_id[
      tip_rows
    ],
    expected_internal_ids
  )

  if (!setequal(
    names(known_states),
    expected_known_ids
  )) {
    stop(
      "Unified known-state pool is incomplete or contains unknown nodes.",
      call. = FALSE
    )
  }

  known_states
}


#' Construct the joint boundary linear system
#'
#' @param patch_plan Patch plan produced by `build_patch_plan()`.
#' @param known_states Unified known-state vector.
#'
#' @return A list containing `A`, `b`, boundary identifiers, and equation
#'   details.
#'
#' @noRd
build_boundary_system <- function(
    patch_plan,
    known_states
) {

  if (!inherits(
    patch_plan,
    "patchASR_patch_plan"
  )) {
    stop(
      "`patch_plan` must be a valid patchASR patch-plan object.",
      call. = FALSE
    )
  }

  if (!is.numeric(known_states) ||
      is.null(names(known_states)) ||
      anyNA(names(known_states)) ||
      any(names(known_states) == "") ||
      anyDuplicated(names(known_states)) ||
      any(!is.finite(known_states))) {
    stop(
      "`known_states` must be a finite named numeric vector ",
      "with unique node identifiers.",
      call. = FALSE
    )
  }

  patches <-
    patch_plan$patches

  boundary_ids <-
    as.character(
      patch_plan$boundary_ids
    )

  if (!is.list(patches) ||
      length(patches) < 1L) {
    stop(
      "`patch_plan` contains no boundary equations.",
      call. = FALSE
    )
  }

  if (length(boundary_ids) !=
      length(patches)) {
    stop(
      "Patch-plan boundary count does not match the number of patches.",
      call. = FALSE
    )
  }

  if (anyNA(boundary_ids) ||
      any(boundary_ids == "") ||
      anyDuplicated(boundary_ids)) {
    stop(
      "Patch plan contains invalid boundary identifiers.",
      call. = FALSE
    )
  }

  if (any(
    boundary_ids %in%
    names(known_states)
  )) {
    stop(
      "An unresolved boundary node must not be present in the ",
      "known-state pool.",
      call. = FALSE
    )
  }

  n_boundary <-
    length(boundary_ids)

  A <- diag(
    n_boundary
  )

  dimnames(A) <- list(
    boundary_ids,
    boundary_ids
  )

  b <- setNames(
    numeric(
      n_boundary
    ),
    boundary_ids
  )

  details <- vector(
    "list",
    n_boundary
  )

  ## ----------------------------------------------------------
  ## Construct one interpolation equation per boundary
  ## ----------------------------------------------------------

  for (i in seq_len(
    n_boundary
  )) {

    patch <-
      patches[[i]]

    boundary_id <-
      patch$parent_id

    if (!identical(
      boundary_id,
      boundary_ids[[i]]
    )) {
      stop(
        "Patch-plan boundary identifiers are inconsistent with ",
        "their patch entries.",
        call. = FALSE
      )
    }

    upstream_id <-
      patch$upstream_id

    sibling_id <-
      patch$sibling_id

    d_parent_boundary <-
      patch$
      branch_length_parent_boundary

    d_boundary_sibling <-
      patch$
      branch_length_boundary_sibling

    if (length(
      d_parent_boundary
    ) != 1L ||
    length(
      d_boundary_sibling
    ) != 1L ||
    !is.finite(
      d_parent_boundary
    ) ||
    !is.finite(
      d_boundary_sibling
    ) ||
    d_parent_boundary <= 0 ||
    d_boundary_sibling <= 0) {
      stop(
        "Boundary interpolation requires finite, strictly positive ",
        "original-tree branch lengths.",
        call. = FALSE
      )
    }

    weight_parent <-
      d_boundary_sibling /
      (
        d_parent_boundary +
          d_boundary_sibling
      )

    weight_sibling <-
      d_parent_boundary /
      (
        d_parent_boundary +
          d_boundary_sibling
      )

    if (!isTRUE(all.equal(
      weight_parent +
      weight_sibling,
      1,
      tolerance = 1e-12
    ))) {
      stop(
        "Boundary interpolation weights do not sum to one.",
        call. = FALSE
      )
    }

    adjacent_ids <- c(
      upstream_id,
      sibling_id
    )

    adjacent_weights <- c(
      weight_parent,
      weight_sibling
    )

    for (j in seq_along(
      adjacent_ids
    )) {

      adjacent_id <-
        adjacent_ids[[j]]

      weight <-
        adjacent_weights[[j]]

      if (adjacent_id %in%
          boundary_ids) {

        column <- match(
          adjacent_id,
          boundary_ids
        )

        A[i, column] <-
          A[i, column] -
          weight

      } else {

        if (!(adjacent_id %in%
              names(known_states))) {
          stop(
            "Required known state is missing for node ",
            adjacent_id,
            ".",
            call. = FALSE
          )
        }

        b[[i]] <-
          b[[i]] +
          weight *
          known_states[[adjacent_id]]
      }
    }

    details[[i]] <- data.frame(
      boundary_id =
        boundary_id,

      upstream_id =
        upstream_id,

      sibling_id =
        sibling_id,

      branch_length_parent_boundary =
        d_parent_boundary,

      branch_length_boundary_sibling =
        d_boundary_sibling,

      weight_parent =
        weight_parent,

      weight_sibling =
        weight_sibling,

      stringsAsFactors =
        FALSE
    )
  }

  structure(
    list(
      A =
        A,

      b =
        b,

      boundary_ids =
        boundary_ids,

      details =
        do.call(
          rbind,
          details
        )
    ),
    class =
      "patchASR_boundary_system"
  )
}


#' Solve a joint boundary system
#'
#' @param boundary_system Result produced by `build_boundary_system()`.
#' @param residual_tolerance Maximum allowed absolute system residual.
#'
#' @return A list containing the named boundary solution and residual.
#'
#' @noRd
solve_boundary_system <- function(
    boundary_system,
    residual_tolerance = 1e-10
) {

  if (!inherits(
    boundary_system,
    "patchASR_boundary_system"
  )) {
    stop(
      "`boundary_system` must be a valid patchASR boundary system.",
      call. = FALSE
    )
  }

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

  A <-
    boundary_system$A

  b <-
    boundary_system$b

  boundary_ids <-
    boundary_system$
    boundary_ids

  if (!is.matrix(A) ||
      !is.numeric(A) ||
      nrow(A) < 1L ||
      nrow(A) != ncol(A) ||
      any(!is.finite(A))) {
    stop(
      "Boundary coefficient matrix is invalid.",
      call. = FALSE
    )
  }

  if (!is.numeric(b) ||
      length(b) != nrow(A) ||
      any(!is.finite(b))) {
    stop(
      "Boundary right-hand-side vector is invalid.",
      call. = FALSE
    )
  }

  if (length(boundary_ids) !=
      nrow(A) ||
      anyNA(boundary_ids) ||
      any(boundary_ids == "") ||
      anyDuplicated(boundary_ids)) {
    stop(
      "Boundary-system node identifiers are invalid.",
      call. = FALSE
    )
  }

  if (qr(A)$rank <
      nrow(A)) {
    stop(
      "Boundary coefficient matrix is singular.",
      call. = FALSE
    )
  }

  raw_solution <- solve(
    A,
    b
  )

  solution <- setNames(
    as.numeric(
      raw_solution
    ),
    boundary_ids
  )

  if (any(
    !is.finite(solution)
  )) {
    stop(
      "Boundary system produced non-finite estimates.",
      call. = FALSE
    )
  }

  residual_vector <-
    as.numeric(
      A %*% solution
    ) -
    as.numeric(
      b
    )

  max_abs_residual <-
    max(
      abs(
        residual_vector
      )
    )

  if (!is.finite(
    max_abs_residual
  ) ||
  max_abs_residual >
  residual_tolerance) {
    stop(
      "Boundary-system numerical residual exceeds the allowed tolerance.",
      call. = FALSE
    )
  }

  list(
    states =
      solution,

    residual =
      max_abs_residual,

    residual_vector =
      setNames(
        residual_vector,
        boundary_ids
      )
  )
}


#' Construct boundary provenance table
#'
#' @param boundary_states Named numeric boundary-state vector.
#'
#' @return A data frame containing boundary estimates and provenance.
#'
#' @noRd
boundary_estimate_table <- function(
    boundary_states
) {

  data.frame(
    node =
      names(
        boundary_states
      ),

    estimate =
      as.numeric(
        boundary_states
      ),

    source =
      rep(
        "boundary_reconstruction",
        length(
          boundary_states
        )
      ),

    component =
      rep(
        "boundary",
        length(
          boundary_states
        )
      ),

    stringsAsFactors =
      FALSE
  )
}


#' Reconstruct all boundary ancestral states
#'
#' @param patch_plan Patch plan produced by `build_patch_plan()`.
#' @param node_map Stable original-tree node map.
#' @param states Named numeric observed terminal states.
#' @param component_asr Component-level reconstruction result.
#' @param residual_tolerance Maximum allowed absolute residual.
#'
#' @return A `patchASR_boundary_result` object.
#'
#' @noRd
run_boundary_reconstruction <- function(
    patch_plan,
    node_map,
    states,
    component_asr,
    residual_tolerance = 1e-10
) {

  known_states <-
    build_unified_known_states(
      node_map =
        node_map,
      states =
        states,
      component_asr =
        component_asr,
      boundary_ids =
        patch_plan$boundary_ids
    )

  boundary_system <-
    build_boundary_system(
      patch_plan =
        patch_plan,
      known_states =
        known_states
    )

  solved <-
    solve_boundary_system(
      boundary_system =
        boundary_system,
      residual_tolerance =
        residual_tolerance
    )

  if (!setequal(
    names(
      solved$states
    ),
    patch_plan$boundary_ids
  )) {
    stop(
      "Solved boundary-state set does not match the patch plan.",
      call. = FALSE
    )
  }

  if (length(
    solved$states
  ) != length(
    patch_plan$boundary_ids
  )) {
    stop(
      "Boundary reconstruction did not produce exactly one state ",
      "per designated boundary node.",
      call. = FALSE
    )
  }

  estimate_table <-
    boundary_estimate_table(
      solved$states
    )

  structure(
    list(
      states =
        solved$states,

      estimate_table =
        estimate_table,

      known_states =
        known_states,

      system =
        boundary_system,

      residual =
        solved$residual,

      residual_vector =
        solved$residual_vector
    ),
    class =
      "patchASR_boundary_result"
  )
}

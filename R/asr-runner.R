# Component-wise ancestral-state reconstruction infrastructure
# for phyloPatch.


#' Prepare observed states for one phylogenetic component
#'
#' @param component_tree A derived `phylo` component.
#' @param states Named numeric observed terminal states from the
#'   original input phylogeny.
#'
#' @return A named numeric vector ordered exactly as
#'   `component_tree$tip.label`.
#'
#' @noRd
prepare_component_states <- function(
    component_tree,
    states
) {

  if (!inherits(
    component_tree,
    "phylo"
  )) {
    stop(
      "`component_tree` must inherit from class `phylo`.",
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
      "All observed state values must be finite.",
      call. = FALSE
    )
  }

  component_tips <-
    component_tree$tip.label

  missing_tips <- setdiff(
    component_tips,
    names(states)
  )

  if (length(
    missing_tips
  ) > 0L) {
    stop(
      "Missing observed state for component tip(s): ",
      paste(
        missing_tips,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  component_states <-
    states[
      component_tips
    ]

  if (!identical(
    names(component_states),
    component_tree$tip.label
  )) {
    stop(
      "Internal error: component states were not reordered correctly.",
      call. = FALSE
    )
  }

  component_states
}


#' Validate and standardize component-level ASR output
#'
#' @param component_tree A phylogenetic component containing one or more
#'   internal nodes carrying stable node identifiers.
#' @param result Raw result returned by `asr_fun`.
#'
#' @return A named numeric vector ordered by the stable internal-node
#'   identities in `component_tree`.
#'
#' @noRd
validate_asr_result <- function(
    component_tree,
    result
) {

  n_internal <-
    component_tree$Nnode

  if (n_internal < 1L) {
    stop(
      "A component with no internal nodes must not be sent to `asr_fun`.",
      call. = FALSE
    )
  }

  expected_ids <-
    component_tree$node.label

  if (is.null(expected_ids) ||
      length(expected_ids) !=
      n_internal ||
      anyNA(expected_ids) ||
      any(expected_ids == "") ||
      anyDuplicated(expected_ids)) {
    stop(
      "Component internal stable node identifiers are missing or invalid.",
      call. = FALSE
    )
  }

  if (!is.numeric(result)) {
    stop(
      "`asr_fun` must return a numeric vector.",
      call. = FALSE
    )
  }

  if (is.null(
    names(result)
  )) {
    stop(
      "`asr_fun` must return a named numeric vector.",
      call. = FALSE
    )
  }

  if (length(result) !=
      n_internal) {
    stop(
      "`asr_fun` returned ",
      length(result),
      " estimate(s), but ",
      n_internal,
      " internal node(s) were expected.",
      call. = FALSE
    )
  }

  if (anyNA(names(result)) ||
      any(names(result) == "")) {
    stop(
      "Every ASR estimate must have a non-empty stable node identifier.",
      call. = FALSE
    )
  }

  if (anyDuplicated(
    names(result)
  )) {
    stop(
      "`asr_fun` returned duplicated node identifiers.",
      call. = FALSE
    )
  }

  if (any(
    !is.finite(result)
  )) {
    stop(
      "`asr_fun` returned non-finite ancestral-state estimates.",
      call. = FALSE
    )
  }

  missing_ids <- setdiff(
    expected_ids,
    names(result)
  )

  unexpected_ids <- setdiff(
    names(result),
    expected_ids
  )

  if (length(missing_ids) > 0L ||
      length(unexpected_ids) > 0L) {

    stop(
      "`asr_fun` returned an incorrect internal-node set. ",
      "Missing: ",
      if (length(missing_ids) == 0L) {
        "<none>"
      } else {
        paste(
          missing_ids,
          collapse = ", "
        )
      },
      "; unexpected: ",
      if (length(unexpected_ids) == 0L) {
        "<none>"
      } else {
        paste(
          unexpected_ids,
          collapse = ", "
        )
      },
      ".",
      call. = FALSE
    )
  }

  standardized <-
    result[
      expected_ids
    ]

  if (!identical(
    names(standardized),
    expected_ids
  )) {
    stop(
      "Internal error while standardizing ASR-result node order.",
      call. = FALSE
    )
  }

  standardized
}


#' Run ASR for one component
#'
#' @param component_tree A phylogenetic component, or `NULL` for a
#'   terminal patch.
#' @param states Named numeric observed terminal states.
#' @param asr_fun Standardized component-level ASR function.
#' @param requires_asr Whether this component requires ASR.
#'
#' @return A named numeric vector containing standardized ancestral-state
#'   estimates. Terminal patches return an empty vector.
#'
#' @noRd
run_component_asr <- function(
    component_tree = NULL,
    states,
    asr_fun,
    requires_asr = TRUE
) {

  ## ----------------------------------------------------------
  ## Terminal patch
  ## ----------------------------------------------------------

  if (!isTRUE(
    requires_asr
  )) {

    return(
      setNames(
        numeric(0L),
        character(0L)
      )
    )
  }

  ## ----------------------------------------------------------
  ## Components requiring ASR
  ## ----------------------------------------------------------

  if (!inherits(
    component_tree,
    "phylo"
  )) {
    stop(
      "An ASR-requiring component must contain a `phylo` tree.",
      call. = FALSE
    )
  }

  if (component_tree$Nnode <
      1L) {
    stop(
      "An ASR-requiring component must contain at least one internal node.",
      call. = FALSE
    )
  }

  if (!is.function(
    asr_fun
  )) {
    stop(
      "`asr_fun` must be a function.",
      call. = FALSE
    )
  }

  component_states <-
    prepare_component_states(
      component_tree,
      states
    )

  raw_result <-
    asr_fun(
      component_tree,
      component_states
    )

  validate_asr_result(
    component_tree,
    raw_result
  )
}


#' Convert one component result to a provenance table
#'
#' @param estimates Named numeric internal-node estimates.
#' @param component Component identifier.
#'
#' @return A data frame containing component-level estimates.
#'
#' @noRd
component_estimate_table <- function(
    estimates,
    component
) {

  n <- length(
    estimates
  )

  data.frame(
    node =
      names(estimates),

    estimate =
      as.numeric(estimates),

    source =
      rep(
        "component_asr",
        n
      ),

    component =
      rep(
        component,
        n
      ),

    stringsAsFactors =
      FALSE
  )
}


#' Run component-wise ASR across an entire partition
#'
#' @param partition A `phyloPatch_partition` object.
#' @param states Named numeric observed terminal states from the
#'   original tree.
#' @param asr_fun Standardized component-level ASR backend.
#'
#' @return A `phyloPatch_component_asr` object.
#'
#' @noRd
run_partition_asr <- function(
    partition,
    states,
    asr_fun
) {

  if (!inherits(
    partition,
    "phyloPatch_partition"
  )) {
    stop(
      "`partition` must be a valid phyloPatch partition object.",
      call. = FALSE
    )
  }

  if (!is.function(
    asr_fun
  )) {
    stop(
      "`asr_fun` must be a function.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Mother component
  ## ----------------------------------------------------------

  mother_states <-
    run_component_asr(
      component_tree =
        partition$mother_tree,
      states =
        states,
      asr_fun =
        asr_fun,
      requires_asr =
        TRUE
    )

  if (!setequal(
    names(mother_states),
    partition$mother_internal_ids
  )) {
    stop(
      "Mother-component ASR result does not match the ",
      "partition's mother internal-node set.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Patch components
  ## ----------------------------------------------------------

  patch_ids <- vapply(
    partition$patch_components,
    function(x) {
      x$patch_id
    },
    character(1)
  )

  if (anyNA(patch_ids) ||
      any(patch_ids == "") ||
      anyDuplicated(patch_ids)) {
    stop(
      "Partition contains invalid patch identifiers.",
      call. = FALSE
    )
  }

  patch_states <- lapply(
    partition$patch_components,
    function(component) {

      result <-
        run_component_asr(
          component_tree =
            component$tree,
          states =
            states,
          asr_fun =
            asr_fun,
          requires_asr =
            component$requires_asr
        )

      if (!setequal(
        names(result),
        component$internal_ids
      )) {
        stop(
          "Patch-component ASR result does not match ",
          "the partition's patch internal-node set.",
          call. = FALSE
        )
      }

      result
    }
  )

  names(patch_states) <-
    patch_ids

  ## ----------------------------------------------------------
  ## Combine all component-derived ancestral states
  ## ----------------------------------------------------------

  patch_internal_states <-
    if (length(
      patch_states
    ) == 0L) {

      setNames(
        numeric(0L),
        character(0L)
      )

    } else {

      do.call(
        c,
        unname(
          patch_states
        )
      )
    }

  if (is.null(
    patch_internal_states
  )) {
    patch_internal_states <-
      setNames(
        numeric(0L),
        character(0L)
      )
  }

  all_internal_states <- c(
    mother_states,
    patch_internal_states
  )

  if (anyDuplicated(
    names(all_internal_states)
  )) {
    stop(
      "Component-wise ASR produced duplicated stable node identities.",
      call. = FALSE
    )
  }

  expected_component_ids <- c(
    partition$mother_internal_ids,
    partition$patch_internal_ids
  )

  if (!setequal(
    names(all_internal_states),
    expected_component_ids
  )) {
    stop(
      "Combined component-level ASR results do not match the ",
      "partition's non-boundary internal-node set.",
      call. = FALSE
    )
  }

  if (length(
    all_internal_states
  ) != length(
    expected_component_ids
  )) {
    stop(
      "Combined component-level ASR results do not contain exactly ",
      "one estimate per non-boundary internal node.",
      call. = FALSE
    )
  }

  if (any(
    !is.finite(
      all_internal_states
    )
  )) {
    stop(
      "Combined component-level ASR results contain non-finite estimates.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Build provenance table
  ## ----------------------------------------------------------

  tables <- list(
    component_estimate_table(
      mother_states,
      "mother"
    )
  )

  for (i in seq_along(
    patch_states
  )) {

    tables[[length(tables) + 1L]] <-
      component_estimate_table(
        patch_states[[i]],
        names(patch_states)[[i]]
      )
  }

  estimate_table <-
    do.call(
      rbind,
      tables
    )

  rownames(
    estimate_table
  ) <- NULL

  ## ----------------------------------------------------------
  ## Return result
  ## ----------------------------------------------------------

  structure(
    list(
      mother_states =
        mother_states,

      patch_states =
        patch_states,

      all_internal_states =
        all_internal_states,

      estimate_table =
        estimate_table
    ),
    class =
      "phyloPatch_component_asr"
  )
}

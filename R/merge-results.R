# Final ancestral-state result merging for phyloPatch.


#' Validate one ancestral-state estimate table
#'
#' @param x Estimate table.
#' @param object_name Name used in diagnostic messages.
#'
#' @return `x`, invisibly.
#'
#' @noRd
validate_estimate_table <- function(
    x,
    object_name
) {

  required_columns <- c(
    "node",
    "estimate",
    "source",
    "component"
  )

  if (!is.data.frame(x) ||
      !all(required_columns %in%
           names(x))) {
    stop(
      "`",
      object_name,
      "` has an invalid estimate-table structure.",
      call. = FALSE
    )
  }

  if (!is.character(x$node) ||
      anyNA(x$node) ||
      any(x$node == "")) {
    stop(
      "`",
      object_name,
      "` contains invalid node identifiers.",
      call. = FALSE
    )
  }

  if (anyDuplicated(x$node)) {
    stop(
      "`",
      object_name,
      "` contains duplicated node identifiers.",
      call. = FALSE
    )
  }

  if (!is.numeric(x$estimate) ||
      any(!is.finite(x$estimate))) {
    stop(
      "`",
      object_name,
      "` contains non-finite or non-numeric estimates.",
      call. = FALSE
    )
  }

  if (!is.character(x$source) ||
      anyNA(x$source) ||
      any(x$source == "")) {
    stop(
      "`",
      object_name,
      "` contains invalid estimate sources.",
      call. = FALSE
    )
  }

  if (!is.character(x$component) ||
      anyNA(x$component) ||
      any(x$component == "")) {
    stop(
      "`",
      object_name,
      "` contains invalid component identifiers.",
      call. = FALSE
    )
  }

  invisible(x)
}


#' Merge component and boundary ancestral-state estimates
#'
#' @param node_map Stable original-tree node map.
#' @param component_asr Result produced by `run_partition_asr()`.
#' @param boundary_result Result produced by
#'   `run_boundary_reconstruction()`.
#'
#' @return A data frame containing exactly one estimate for every
#'   original internal node.
#'
#' @noRd
merge_ancestral_results <- function(
    node_map,
    component_asr,
    boundary_result
) {

  ## ----------------------------------------------------------
  ## Validate node map
  ## ----------------------------------------------------------

  required_map_columns <- c(
    "original_node",
    "node_type",
    "original_label",
    "stable_id"
  )

  if (!is.data.frame(node_map) ||
      !all(required_map_columns %in%
           names(node_map))) {
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

  ## ----------------------------------------------------------
  ## Validate upstream result objects
  ## ----------------------------------------------------------

  if (!inherits(
    component_asr,
    "phyloPatch_component_asr"
  )) {
    stop(
      "`component_asr` must be a valid phyloPatch component-ASR object.",
      call. = FALSE
    )
  }

  if (!inherits(
    boundary_result,
    "phyloPatch_boundary_result"
  )) {
    stop(
      "`boundary_result` must be a valid phyloPatch boundary result.",
      call. = FALSE
    )
  }

  component_table <-
    component_asr$estimate_table

  boundary_table <-
    boundary_result$estimate_table

  validate_estimate_table(
    component_table,
    "component_asr$estimate_table"
  )

  validate_estimate_table(
    boundary_table,
    "boundary_result$estimate_table"
  )

  ## ----------------------------------------------------------
  ## Cross-check component table against its state vector
  ## ----------------------------------------------------------

  component_states <-
    component_asr$all_internal_states

  if (!is.numeric(component_states) ||
      is.null(names(component_states)) ||
      anyNA(names(component_states)) ||
      any(names(component_states) == "") ||
      anyDuplicated(names(component_states)) ||
      any(!is.finite(component_states))) {
    stop(
      "`component_asr$all_internal_states` is invalid.",
      call. = FALSE
    )
  }

  if (!setequal(
    component_table$node,
    names(component_states)
  )) {
    stop(
      "Component estimate table does not match the component state vector.",
      call. = FALSE
    )
  }

  expected_component_values <-
    as.numeric(
      component_states[
        component_table$node
      ]
    )

  if (!isTRUE(all.equal(
    component_table$estimate,
    expected_component_values,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop(
      "Component estimate table contains values inconsistent with ",
      "the component state vector.",
      call. = FALSE
    )
  }

  if (!all(
    component_table$source ==
    "component_asr"
  )) {
    stop(
      "Component estimate table contains invalid provenance sources.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Cross-check boundary table against its state vector
  ## ----------------------------------------------------------

  boundary_states <-
    boundary_result$states

  if (!is.numeric(boundary_states) ||
      is.null(names(boundary_states)) ||
      anyNA(names(boundary_states)) ||
      any(names(boundary_states) == "") ||
      anyDuplicated(names(boundary_states)) ||
      any(!is.finite(boundary_states))) {
    stop(
      "`boundary_result$states` is invalid.",
      call. = FALSE
    )
  }

  if (!setequal(
    boundary_table$node,
    names(boundary_states)
  )) {
    stop(
      "Boundary estimate table does not match the boundary state vector.",
      call. = FALSE
    )
  }

  expected_boundary_values <-
    as.numeric(
      boundary_states[
        boundary_table$node
      ]
    )

  if (!isTRUE(all.equal(
    boundary_table$estimate,
    expected_boundary_values,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop(
      "Boundary estimate table contains values inconsistent with ",
      "the boundary state vector.",
      call. = FALSE
    )
  }

  if (!all(
    boundary_table$source ==
    "boundary_reconstruction"
  )) {
    stop(
      "Boundary estimate table contains invalid provenance sources.",
      call. = FALSE
    )
  }

  if (!all(
    boundary_table$component ==
    "boundary"
  )) {
    stop(
      "Boundary estimate table contains invalid component provenance.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Combine both result sources
  ## ----------------------------------------------------------

  combined <- rbind(
    component_table[
      ,
      c(
        "node",
        "estimate",
        "source",
        "component"
      ),
      drop = FALSE
    ],
    boundary_table[
      ,
      c(
        "node",
        "estimate",
        "source",
        "component"
      ),
      drop = FALSE
    ]
  )

  rownames(combined) <- NULL

  if (anyDuplicated(
    combined$node
  )) {
    stop(
      "An original internal node occurs more than once in the merged result.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Require exactly the complete original internal-node set
  ## ----------------------------------------------------------

  internal_rows <- which(
    node_map$node_type ==
      "internal"
  )

  expected_internal_ids <-
    node_map$stable_id[
      internal_rows
    ]

  missing_ids <- setdiff(
    expected_internal_ids,
    combined$node
  )

  unexpected_ids <- setdiff(
    combined$node,
    expected_internal_ids
  )

  if (length(missing_ids) > 0L ||
      length(unexpected_ids) > 0L) {

    stop(
      "Merged ancestral-state result has an incorrect internal-node set. ",
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

  if (nrow(combined) !=
      length(expected_internal_ids)) {
    stop(
      "Merged ancestral-state result does not contain exactly one ",
      "estimate per original internal node.",
      call. = FALSE
    )
  }

  if (any(!is.finite(
    combined$estimate
  ))) {
    stop(
      "Merged ancestral-state result contains non-finite estimates.",
      call. = FALSE
    )
  }

  ## ----------------------------------------------------------
  ## Map final rows back to the original tree
  ## ----------------------------------------------------------

  map_rows <- match(
    combined$node,
    node_map$stable_id
  )

  if (anyNA(map_rows)) {
    stop(
      "Merged result contains a node that cannot be mapped to the ",
      "original tree.",
      call. = FALSE
    )
  }

  if (any(
    node_map$node_type[
      map_rows
    ] != "internal"
  )) {
    stop(
      "Merged result contains a node that is not an original ",
      "internal node.",
      call. = FALSE
    )
  }

  final <- data.frame(
    node =
      combined$node,

    original_node =
      as.integer(
        node_map$original_node[
          map_rows
        ]
      ),

    original_label =
      node_map$original_label[
        map_rows
      ],

    estimate =
      combined$estimate,

    source =
      combined$source,

    component =
      combined$component,

    stringsAsFactors =
      FALSE
  )

  ## ----------------------------------------------------------
  ## Standardize final order to original-tree node order
  ## ----------------------------------------------------------

  final <- final[
    order(
      final$original_node
    ),
    ,
    drop = FALSE
  ]

  rownames(final) <- NULL

  expected_order <-
    node_map$stable_id[
      internal_rows
    ]

  if (!identical(
    final$node,
    expected_order
  )) {
    stop(
      "Internal error while restoring final original-tree node order.",
      call. = FALSE
    )
  }

  final
}

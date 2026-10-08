# Built-in ASR adapters based on ape::ace().


#' Validate inputs supplied to a built-in ASR adapter
#'
#' @param tree Component phylogeny.
#' @param states Named terminal-state vector.
#' @param method_name Method name used in diagnostics.
#'
#' @return Stable internal-node identifiers, invisibly.
#'
#' @noRd
validate_builtin_asr_inputs <- function(
    tree,
    states,
    method_name
) {

  if (!inherits(
    tree,
    "phylo"
  )) {
    stop(
      method_name,
      ": `tree` must inherit from class `phylo`.",
      call. = FALSE
    )
  }

  if (!is.numeric(states) ||
      is.null(names(states)) ||
      anyNA(names(states)) ||
      any(names(states) == "") ||
      anyDuplicated(names(states))) {
    stop(
      method_name,
      ": `states` must be a uniquely named numeric vector.",
      call. = FALSE
    )
  }

  if (length(states) !=
      ape::Ntip(tree)) {
    stop(
      method_name,
      ": number of terminal states does not equal Ntip(tree).",
      call. = FALSE
    )
  }

  if (any(!is.finite(
    states
  ))) {
    stop(
      method_name,
      ": terminal states must all be finite.",
      call. = FALSE
    )
  }

  if (!identical(
    names(states),
    tree$tip.label
  )) {
    stop(
      method_name,
      ": terminal states are not in exact tree-tip order.",
      call. = FALSE
    )
  }

  if (tree$Nnode < 1L) {
    return(
      invisible(
        character(0)
      )
    )
  }

  stable_ids <-
    tree$node.label

  if (is.null(stable_ids) ||
      length(stable_ids) !=
      tree$Nnode ||
      anyNA(stable_ids) ||
      any(stable_ids == "") ||
      anyDuplicated(stable_ids)) {
    stop(
      method_name,
      ": component tree does not contain a complete unique ",
      "stable internal-node identity map.",
      call. = FALSE
    )
  }

  invisible(
    stable_ids
  )
}


#' Map an ape::ace result to stable component node identifiers
#'
#' @param fit Result returned by `ape::ace()`.
#' @param tree Component phylogeny.
#' @param method_name Method name used in diagnostics.
#'
#' @return Named numeric vector in stable-node order.
#'
#' @noRd
standardize_ape_ace_result <- function(
    fit,
    tree,
    method_name
) {

  if (is.null(
    fit$ace
  )) {
    stop(
      method_name,
      ": `ape::ace()` did not return an `ace` result.",
      call. = FALSE
    )
  }

  values <-
    fit$ace

  if (!is.numeric(values) ||
      length(values) !=
      tree$Nnode) {
    stop(
      method_name,
      ": `ape::ace()` returned an unexpected number ",
      "of internal-node estimates.",
      call. = FALSE
    )
  }

  stable_ids <-
    tree$node.label

  local_ids <-
    as.character(
      ape::Ntip(tree) +
        seq_len(
          tree$Nnode
        )
    )

  value_names <-
    names(values)

  if (is.null(value_names) ||
      anyNA(value_names) ||
      any(value_names == "") ||
      anyDuplicated(value_names)) {
    stop(
      method_name,
      ": `ape::ace()` returned estimates without a usable ",
      "internal-node identity map.",
      call. = FALSE
    )
  }

  ## ape versions normally identify ace estimates by the
  ## component tree's local ape node numbers.
  if (setequal(
    value_names,
    local_ids
  )) {

    values <-
      values[
        local_ids
      ]

    ## Retain support if an implementation returns the supplied
    ## internal labels directly.
  } else if (setequal(
    value_names,
    stable_ids
  )) {

    values <-
      values[
        stable_ids
      ]

  } else {

    stop(
      method_name,
      ": `ape::ace()` internal-node names cannot be mapped ",
      "safely to the component tree. Observed names: ",
      paste(
        utils::head(
          value_names,
          20L
        ),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  values <-
    as.numeric(
      values
    )

  if (any(!is.finite(
    values
  ))) {
    stop(
      method_name,
      ": `ape::ace()` produced non-finite ancestral estimates.",
      call. = FALSE
    )
  }

  stats::setNames(
    values,
    stable_ids
  )
}


#' PIC-style ancestral-state reconstruction using ape
#'
#' Runs continuous-trait ancestral-state reconstruction with
#' `ape::ace(method = "pic")` and maps the resulting component-tree
#' estimates to the stable internal-node identifiers used by
#' `patchASR`.
#'
#' This function satisfies the `asr_fun(tree, states)` interface required
#' by [patch_asr()].
#'
#' @param tree A component phylogeny supplied by `patchASR`.
#' @param states A named numeric vector of terminal states in exact
#'   `tree$tip.label` order.
#'
#' @return A named numeric vector containing one finite estimate for every
#'   internal node of `tree`.
#'
#' @export
asr_ape_pic <- function(
    tree,
    states
) {

  validate_builtin_asr_inputs(
    tree,
    states,
    "ape PIC-style ASR"
  )

  fit <- ape::ace(
    x = states,
    phy = tree,
    type = "continuous",
    method = "pic",
    CI = FALSE
  )

  standardize_ape_ace_result(
    fit,
    tree,
    "ape PIC-style ASR"
  )
}


#' Maximum-likelihood BM ancestral-state reconstruction using ape
#'
#' Runs continuous-trait ancestral-state reconstruction with
#' `ape::ace(method = "ML", model = "BM")`.
#'
#' This function satisfies the `asr_fun(tree, states)` interface required
#' by [patch_asr()].
#'
#' @inheritParams asr_ape_pic
#'
#' @return A named numeric vector containing one finite estimate for every
#'   internal node of `tree`.
#'
#' @export
asr_ape_ml_bm <- function(
    tree,
    states
) {

  validate_builtin_asr_inputs(
    tree,
    states,
    "ape ML-BM ASR"
  )

  fit <- ape::ace(
    x = states,
    phy = tree,
    type = "continuous",
    method = "ML",
    model = "BM",
    CI = FALSE
  )

  standardize_ape_ace_result(
    fit,
    tree,
    "ape ML-BM ASR"
  )
}


#' GLS BM ancestral-state reconstruction using ape
#'
#' Runs continuous-trait ancestral-state reconstruction with
#' `ape::ace(method = "GLS")` and a Brownian correlation structure
#' created with `ape::corBrownian()`.
#'
#' This function satisfies the `asr_fun(tree, states)` interface required
#' by [patch_asr()].
#'
#' @inheritParams asr_ape_pic
#'
#' @return A named numeric vector containing one finite estimate for every
#'   internal node of `tree`.
#'
#' @export
asr_ape_gls_bm <- function(
    tree,
    states
) {

  validate_builtin_asr_inputs(
    tree,
    states,
    "ape GLS-BM ASR"
  )

  cor_bm <- ape::corBrownian(
    value = 1,
    phy = tree
  )

  fit <- withCallingHandlers(
    ape::ace(
      x = states,
      phy = tree,
      type = "continuous",
      method = "GLS",
      corStruct = cor_bm,
      CI = FALSE
    ),
    warning = function(w) {

      message_text <-
        conditionMessage(w)

      nuisance_warning <-
        grepl(
          paste0(
            "No covariate specified, species will be taken ",
            "as ordered in the data frame"
          ),
          message_text,
          fixed = TRUE
        )

      if (nuisance_warning) {
        invokeRestart(
          "muffleWarning"
        )
      }
    }
  )

  standardize_ape_ace_result(
    fit,
    tree,
    "ape GLS-BM ASR"
  )
}

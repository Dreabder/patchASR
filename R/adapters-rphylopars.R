# Rphylopars ancestral-state reconstruction adapter.


#' Brownian-motion ancestral-state reconstruction using Rphylopars
#'
#' Creates an ancestral-state reconstruction backend based on
#' `Rphylopars::phylopars()` under a Brownian-motion model.
#'
#' The returned function satisfies the `asr_fun(tree, states)` interface
#' required by [patch_asr()].
#'
#' Internally, the adapter works on a temporary copy of the component tree
#' whose tip labels are replaced by collision-safe identifiers. This avoids
#' ambiguity when original numeric tip labels coincide with local `ape`
#' internal-node numbers after pruning. The original tree is not modified.
#'
#' The adapter verifies the complete row-identity system returned by
#' `Rphylopars::phylopars()` before mapping reconstructed internal nodes
#' back to the stable identifiers stored in `tree$node.label`.
#'
#' @param REML Logical. Whether restricted maximum likelihood should be used
#'   by `Rphylopars::phylopars()`. Defaults to `TRUE`.
#'
#' @return A function with arguments `tree` and `states` that returns one
#'   finite ancestral-state estimate for every internal node of `tree`,
#'   named with the stable node identifiers in `tree$node.label`.
#'
#' @export
make_asr_rphylopars_bm <- function(
    REML = TRUE
) {

  if (!is.logical(REML) ||
      length(REML) != 1L ||
      is.na(REML)) {

    stop(
      "`REML` must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  function(
    tree,
    states
  ) {

    if (!requireNamespace(
      "Rphylopars",
      quietly = TRUE
    )) {

      stop(
        "Package `Rphylopars` is required for ",
        "`make_asr_rphylopars_bm()` backends. ",
        "Install it before using this adapter.",
        call. = FALSE
      )
    }


    validate_builtin_asr_inputs(
      tree,
      states,
      "Rphylopars BM ASR"
    )


    if (tree$Nnode == 0L) {

      return(
        stats::setNames(
          numeric(0),
          character(0)
        )
      )
    }


    ## --------------------------------------------------------
    ## Construct a collision-safe temporary working tree.
    ##
    ## Numeric biological tip labels can collide with local ape
    ## internal-node numbers after pruning. Rphylopars combines
    ## tip and ancestral reconstructions in one row-identity
    ## system, so such collisions can make the returned rows
    ## ambiguous.
    ##
    ## Relabeling tips on an internal copy is safe because the
    ## topology, branch lengths, trait values, and their exact
    ## correspondence are preserved.
    ## --------------------------------------------------------

    rph_tree <-
      tree


    rph_tree$tip.label <-
      paste0(
        "phylopatch_tip_",
        seq_len(
          ape::Ntip(rph_tree)
        )
      )


    ## Rphylopars does not need phyloPatch stable internal labels.
    ## Remove them from the temporary working copy so that local
    ## ape node identity is the only internal-node convention
    ## relevant to this backend.
    rph_tree$node.label <-
      NULL


    rph_states <-
      stats::setNames(
        as.numeric(states),
        rph_tree$tip.label
      )


    trait_data <-
      data.frame(
        species = names(rph_states),
        trait = as.numeric(rph_states),
        stringsAsFactors = FALSE
      )


    ## --------------------------------------------------------
    ## Run Rphylopars
    ## --------------------------------------------------------

    fit <-
      Rphylopars::phylopars(
        trait_data = trait_data,
        tree = rph_tree,
        model = "BM",
        pheno_error = FALSE,
        REML = REML
      )


    if (is.null(
      fit$anc_recon
    )) {

      stop(
        "Rphylopars BM ASR: `phylopars()` did not return ",
        "`anc_recon`.",
        call. = FALSE
      )
    }


    anc <-
      as.data.frame(
        fit$anc_recon
      )


    expected_rows <-
      ape::Ntip(rph_tree) +
      rph_tree$Nnode


    if (nrow(anc) !=
        expected_rows) {

      stop(
        "Rphylopars BM ASR: `anc_recon` has ",
        nrow(anc),
        " rows; expected ",
        expected_rows,
        " (= Ntip + Nnode).",
        call. = FALSE
      )
    }


    if (!"trait" %in%
        colnames(anc)) {

      stop(
        "Rphylopars BM ASR: `anc_recon` does not contain ",
        "the expected `trait` column.",
        call. = FALSE
      )
    }


    row_ids <-
      rownames(anc)


    if (is.null(row_ids) ||
        anyNA(row_ids) ||
        any(row_ids == "") ||
        anyDuplicated(row_ids)) {

      stop(
        "Rphylopars BM ASR: `anc_recon` does not have ",
        "a complete unique row identity map.",
        call. = FALSE
      )
    }


    ## --------------------------------------------------------
    ## Expected identities on the temporary component tree
    ## --------------------------------------------------------

    local_internal_ids <-
      as.character(
        ape::Ntip(rph_tree) +
          seq_len(
            rph_tree$Nnode
          )
      )


    expected_raw_ids <-
      c(
        rph_tree$tip.label,
        local_internal_ids
      )


    ## Because temporary tip names are generated by this adapter,
    ## this should never happen. Keep the check as an invariant.
    if (anyDuplicated(
      expected_raw_ids
    )) {

      stop(
        "Rphylopars BM ASR internal error: temporary tip ",
        "identifiers collide with local node identifiers.",
        call. = FALSE
      )
    }


    expected_syntactic_ids <-
      make.names(
        expected_raw_ids,
        unique = TRUE
      )


    ## --------------------------------------------------------
    ## Resolve the complete Rphylopars row-identity convention.
    ##
    ## Accept either:
    ##
    ##   1. raw temporary tip labels + raw local ape node IDs;
    ##   2. the complete syntactically normalized form.
    ##
    ## No positional fallback is used.
    ## --------------------------------------------------------

    if (setequal(
      row_ids,
      expected_raw_ids
    )) {

      internal_row_ids <-
        local_internal_ids

    } else if (setequal(
      row_ids,
      expected_syntactic_ids
    )) {

      internal_positions <-
        match(
          local_internal_ids,
          expected_raw_ids
        )


      internal_row_ids <-
        expected_syntactic_ids[
          internal_positions
        ]

    } else {

      stop(
        "Rphylopars BM ASR: `anc_recon` row identities cannot ",
        "be mapped safely to the temporary component tree. ",
        "Neither the complete raw identity set nor its ",
        "syntactic `make.names()` form matches the returned ",
        "row-name set.",
        call. = FALSE
      )
    }


    if (length(internal_row_ids) !=
        tree$Nnode ||
        anyNA(internal_row_ids) ||
        !all(
          internal_row_ids %in%
          row_ids
        )) {

      stop(
        "Rphylopars BM ASR: resolved internal-node row ",
        "identities are incomplete.",
        call. = FALSE
      )
    }


    ## --------------------------------------------------------
    ## Extract internal-node estimates
    ## --------------------------------------------------------

    values <-
      anc[
        internal_row_ids,
        "trait",
        drop = TRUE
      ]


    values <-
      as.numeric(
        values
      )


    if (length(values) !=
        tree$Nnode) {

      stop(
        "Rphylopars BM ASR: extracted an unexpected number ",
        "of ancestral estimates.",
        call. = FALSE
      )
    }


    if (any(!is.finite(
      values
    ))) {

      stop(
        "Rphylopars BM ASR produced non-finite ancestral ",
        "estimates.",
        call. = FALSE
      )
    }


    ## --------------------------------------------------------
    ## Restore phyloPatch stable biological node identity
    ## --------------------------------------------------------

    stats::setNames(
      values,
      tree$node.label
    )
  }
}

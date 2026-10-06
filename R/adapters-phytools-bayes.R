# Bayesian ancestral-state reconstruction adapter based on
# phytools::anc.Bayes().


#' Validate a positive integer-like scalar
#'
#' @param x Value to validate.
#' @param name Argument name used in diagnostics.
#'
#' @return Integer value.
#'
#' @noRd
validate_positive_integer_scalar <- function(
    x,
    name
) {

  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x <= 0 ||
      x != floor(x) ||
      x > .Machine$integer.max) {

    stop(
      "`", name,
      "` must be a single positive integer.",
      call. = FALSE
    )
  }

  as.integer(x)
}


#' Derive a deterministic seed for one phylogenetic component
#'
#' The seed depends on the user-supplied base seed and the stable
#' biological identity of the component rather than on the order in
#' which components happen to be evaluated.
#'
#' @param base_seed Positive integer base seed.
#' @param tree Component phylogeny.
#'
#' @return A positive integer seed.
#'
#' @noRd
stable_component_seed <- function(
    base_seed,
    tree
) {

  key <- paste(
    c(
      "tips",
      sort(tree$tip.label),
      "nodes",
      sort(tree$node.label)
    ),
    collapse = "|"
  )

  modulus <- 1900000000

  value <-
    as.double(base_seed) %% modulus

  code_points <-
    utf8ToInt(
      enc2utf8(key)
    )

  for (code_point in
       code_points) {

    value <-
      (
        value * 131 +
          code_point
      ) %% modulus
  }

  as.integer(
    max(
      1,
      floor(value)
    )
  )
}


#' Bayesian ancestral-state reconstruction using phytools
#'
#' Creates a Bayesian ancestral-state reconstruction backend based on
#' `phytools::anc.Bayes()`.
#'
#' The returned function satisfies the `asr_fun(tree, states)` interface
#' required by [patch_asr()].
#'
#' Posterior ancestral-state estimates are calculated as the column means
#' of retained MCMC samples after discarding the requested fraction of
#' sampled rows as burn-in.
#'
#' A deterministic component-specific random seed is derived from the
#' supplied base seed and the stable identities of the taxa and internal
#' nodes in each component. This avoids making stochastic results depend
#' on the order in which compatible components are evaluated.
#'
#' The current adapter does not silently replace `anc.Bayes()` with another
#' ancestral-state reconstruction method when `anc.Bayes()` cannot handle
#' a component.
#'
#' @param ngen Positive integer number of MCMC generations. Defaults to
#'   `10000`.
#' @param sample_freq Positive integer MCMC sampling frequency. Defaults to
#'   `100`.
#' @param burnin_frac Numeric value satisfying `0 <= burnin_frac < 1`.
#'   This fraction of sampled MCMC rows is discarded as burn-in.
#' @param seed Positive integer base random seed. A deterministic
#'   component-specific seed is derived from this value. Defaults to `123`.
#'
#' @return A function with arguments `tree` and `states` that returns one
#'   finite posterior-mean ancestral-state estimate for every internal node
#'   of `tree`, named with the stable identifiers in `tree$node.label`.
#'
#' @export
make_asr_phytools_bayes <- function(
    ngen = 10000L,
    sample_freq = 100L,
    burnin_frac = 0.20,
    seed = 123L
) {

  ngen <-
    validate_positive_integer_scalar(
      ngen,
      "ngen"
    )

  sample_freq <-
    validate_positive_integer_scalar(
      sample_freq,
      "sample_freq"
    )

  seed <-
    validate_positive_integer_scalar(
      seed,
      "seed"
    )

  if (!is.numeric(burnin_frac) ||
      length(burnin_frac) != 1L ||
      !is.finite(burnin_frac) ||
      burnin_frac < 0 ||
      burnin_frac >= 1) {

    stop(
      "`burnin_frac` must be a single finite numeric value ",
      "satisfying 0 <= burnin_frac < 1.",
      call. = FALSE
    )
  }

  if (sample_freq > ngen) {

    stop(
      "`sample_freq` cannot exceed `ngen`.",
      call. = FALSE
    )
  }

  function(
    tree,
    states
  ) {

    if (!requireNamespace(
      "phytools",
      quietly = TRUE
    )) {

      stop(
        "Package `phytools` is required for ",
        "`make_asr_phytools_bayes()` backends. ",
        "Install it before using this adapter.",
        call. = FALSE
      )
    }

    validate_builtin_asr_inputs(
      tree,
      states,
      "phytools Bayesian ASR"
    )

    if (tree$Nnode == 0L) {

      return(
        stats::setNames(
          numeric(0),
          character(0)
        )
      )
    }

    ## Audited phytools 2.4.4 behaviour:
    ## anc.Bayes() fails on a binary two-tip tree with one
    ## internal node. Do not silently substitute another ASR method.
    if (tree$Nnode == 1L) {

      stop(
        "phytools Bayesian ASR cannot be applied safely to a ",
        "component containing only one internal node. ",
        "`phytools::anc.Bayes()` failed on this minimum binary ",
        "component during package validation. No fallback method ",
        "has been substituted.",
        call. = FALSE
      )
    }

    local_internal_ids <-
      as.character(
        ape::Ntip(tree) +
          seq_len(
            tree$Nnode
          )
      )

    component_seed <-
      stable_component_seed(
        seed,
        tree
      )

    ## Preserve the caller's global RNG state.
    had_random_seed <-
      exists(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
      )

    if (had_random_seed) {

      old_random_seed <-
        get(
          ".Random.seed",
          envir = globalenv(),
          inherits = FALSE
        )
    }

    on.exit(
      {

        if (had_random_seed) {

          assign(
            ".Random.seed",
            old_random_seed,
            envir = globalenv()
          )

        } else if (exists(
          ".Random.seed",
          envir = globalenv(),
          inherits = FALSE
        )) {

          rm(
            ".Random.seed",
            envir = globalenv()
          )
        }
      },
      add = TRUE
    )

    set.seed(
      component_seed
    )

    fit <- tryCatch(
      phytools::anc.Bayes(
        tree = tree,
        x = states,
        ngen = ngen,
        control = list(
          sample = sample_freq
        )
      ),
      error = function(e) {

        stop(
          "phytools Bayesian ASR failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (is.null(
      fit$mcmc
    )) {

      stop(
        "phytools Bayesian ASR did not return an MCMC table.",
        call. = FALSE
      )
    }

    mcmc <-
      as.data.frame(
        fit$mcmc
      )

    if (nrow(mcmc) < 2L) {

      stop(
        "phytools Bayesian ASR returned too few MCMC samples.",
        call. = FALSE
      )
    }

    burnin_n <-
      floor(
        nrow(mcmc) *
          burnin_frac
      )

    if (burnin_n >=
        nrow(mcmc)) {

      stop(
        "Bayesian burn-in removes all MCMC samples.",
        call. = FALSE
      )
    }

    post_mcmc <-
      mcmc[
        (burnin_n + 1L):nrow(mcmc),
        ,
        drop = FALSE
      ]

    mcmc_names <-
      colnames(
        post_mcmc
      )

    if (is.null(mcmc_names) ||
        anyNA(mcmc_names) ||
        any(mcmc_names == "") ||
        anyDuplicated(mcmc_names)) {

      stop(
        "phytools Bayesian ASR returned MCMC output without ",
        "a complete unique column identity map.",
        call. = FALSE
      )
    }

    if (!all(
      local_internal_ids %in%
      mcmc_names
    )) {

      stop(
        "phytools Bayesian ASR could not identify all local ",
        "`ape` internal-node columns in the MCMC output. Missing IDs: ",
        paste(
          setdiff(
            local_internal_ids,
            mcmc_names
          ),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    values <-
      colMeans(
        post_mcmc[
          ,
          local_internal_ids,
          drop = FALSE
        ],
        na.rm = TRUE
      )

    values <-
      as.numeric(
        values
      )

    if (length(values) !=
        tree$Nnode) {

      stop(
        "phytools Bayesian ASR extracted an unexpected number ",
        "of ancestral estimates.",
        call. = FALSE
      )
    }

    if (any(!is.finite(
      values
    ))) {

      stop(
        "phytools Bayesian ASR produced non-finite posterior means.",
        call. = FALSE
      )
    }

    stats::setNames(
      values,
      tree$node.label
    )
  }
}

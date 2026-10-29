# ============================================================================
# compute_custom_bic()
# BIC for hybrid (mixed continuous / discrete) DAGs without the
# Conditional-Gaussian (CG) restriction.
#
# bnlearn's score(..., type = "bic-cg") requires that discrete nodes have no
# continuous parents and stops with "arc X -> Y violates the assumptions of
# the model" otherwise. Networks learned with tpc() + mixCItest can contain
# such arcs, so the BIC is computed node by node instead:
#
#   continuous node : lm(node ~ parents)              (any parent types;
#                                                      factors are dummy-coded)
#   discrete node   : nnet::multinom(node ~ parents)  (any parent types)
#
# Log-likelihoods and parameter counts are summed over nodes:
#   BIC = -2 * logLik + k * log(n)        (LOWER is better)
# Note: this is the opposite sign convention to bnlearn::score(),
# where HIGHER is better.
#
# Arguments
#   dag     : bnlearn "bn" object
#   data    : data frame containing all nodes of dag
#   verbose : print per-node log-likelihood and parameter count
#
# Returns a list: loglik, k (number of parameters), BIC, per_node (data frame)
# ============================================================================

if (!requireNamespace("nnet", quietly = TRUE)) {
  install.packages("nnet")
}
library(bnlearn)
library(nnet)

compute_custom_bic <- function(dag, data, verbose = FALSE) {

  node_list    <- nodes(dag)
  n            <- nrow(data)
  total_loglik <- 0
  total_params <- 0
  per_node     <- data.frame(node = character(), loglik = numeric(),
                             k = numeric(), stringsAsFactors = FALSE)

  for (nd in node_list) {

    par <- bnlearn::parents(dag, nd)
    y   <- data[[nd]]

    # Model formula: intercept-only for root nodes
    formula_str <- if (length(par) == 0) {
      paste0(nd, " ~ 1")
    } else {
      paste0(nd, " ~ ", paste(par, collapse = " + "))
    }
    f <- as.formula(formula_str)

    if (is.numeric(y)) {
      # Continuous node: OLS regression on its parents
      fit <- lm(f, data = data)
      ll  <- as.numeric(stats::logLik(fit))
      k   <- attr(stats::logLik(fit), "df")

    } else {
      # Discrete node: multinomial logistic regression on its parents
      fit <- nnet::multinom(f, data = data, trace = FALSE)
      ll  <- as.numeric(stats::logLik(fit))
      k   <- attr(stats::logLik(fit), "df")
    }

    if (verbose) cat(sprintf("  %-22s parents=%s  logLik=%.2f  k=%d\n",
                             nd, paste(par, collapse = ","), ll, k))

    total_loglik <- total_loglik + ll
    total_params <- total_params + k
    per_node <- rbind(per_node, data.frame(node = nd, loglik = ll, k = k))
  }

  bic <- -2 * total_loglik + total_params * log(n)

  list(loglik   = total_loglik,
       k        = total_params,
       BIC      = bic,
       per_node = per_node)
}

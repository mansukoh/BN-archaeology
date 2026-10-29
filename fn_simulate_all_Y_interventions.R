# ============================================================================
# Monte Carlo simulation of interventional distributions
#
# Functions defined in this file:
#   simulate_all_Y_interventions() : run simulations for every intervention
#                                    variable Y at once
#   simulate_doY_region()          : simulate all nodes under one intervention
#                                    on Y
#   draw_Y_empirical()             : sample Y from the observed data within
#                                    the intervention region
# ============================================================================


# ============================================================================
# simulate_all_Y_interventions()
# For each intervention variable Y:
#   continuous Y  : simulate under do(Y <= ay) and do(Y > ay)
#   categorical Y : simulate under do(Y = k) for each level k
#
# Arguments
#   fitted_list : fitted local model for each node (see main.rmd)
#   avg_graph   : consensus DAG (bnlearn "bn" object)
#   dat         : data used to fit the models
#   Y_vars      : intervention variables
#   ay_list     : named list; threshold (continuous Y) or vector of levels
#                 (categorical Y) for each Y
#   n           : number of Monte Carlo samples per intervention
#
# Returns a named list (one element per Y):
#   continuous Y  : list(le = <data frame>, gt = <data frame>, ay, n)
#   categorical Y : list(simulations = list(<level> = <data frame>), levels, n)
# ============================================================================

simulate_all_Y_interventions <- function(fitted_list, avg_graph, dat,
                                         Y_vars,  # intervention variables
                                         ay_list,
                                         n = 1e6) {

  results <- list()

  for (Y in Y_vars) {
    ay <- ay_list[[Y]]

    if(is.numeric(dat[[Y]])) {
      # Continuous Y: do(Y <= ay) and do(Y > ay)
      sim_le <- simulate_doY_region(fitted_list, avg_graph, dat,
                                    Y, region = "le", ay = ay, n = n)

      sim_gt <- simulate_doY_region(fitted_list, avg_graph, dat,
                                    Y, region = "gt", ay = ay, n = n)

      results[[Y]] <- list(
        le = sim_le,
        gt = sim_gt,
        ay = ay,
        n = n
      )
    } else {
      # Categorical Y: do(Y = level) for each level
      sim_list <- list()
      for(level in ay) {  # here ay = levels(dat[[Y]])
        sim_list[[level]] <- simulate_doY_region(fitted_list, avg_graph, dat,
                                                 Y, region = "eq", ay = level, n = n)
      }

      results[[Y]] <- list(
        simulations = sim_list,
        levels = ay,
        n = n
      )
    }
  }

  return(results)
}


# ============================================================================
# simulate_doY_region()
# Simulate the joint distribution of all nodes under an intervention on Y,
# using the fitted local models in fitted_list.
#
# Simulation proceeds in two passes over the topological ordering of the DAG:
#   Pass 1: generate all non-descendants of Y from their fitted models;
#           set Y by sampling from the intervention region (draw_Y_empirical).
#   Pass 2: generate the descendants of Y in topological order, conditioning
#           on the simulated values of their parents (which now reflect the
#           intervention on Y).
#
# Arguments
#   region : "le" = do(Y <= ay), "gt" = do(Y > ay), "eq" = do(Y = ay)
#   ay     : threshold (continuous Y) or level (categorical Y)
#   n      : number of Monte Carlo samples
#
# Returns a data frame with n rows and one column per node.
# ============================================================================

simulate_doY_region <- function(fitted_list, avg_graph, dat,
                                Y, region = c("le", "gt", "eq"), ay, n = 1e6) {

  region <- match.arg(region)
  topo   <- bnlearn::node.ordering(avg_graph)   # topological ordering of all nodes

  discrete_vars <- c("SettlementType", "CommOrg", "HouseholdOrg")

  # Descendants of Y: the only nodes affected by the intervention
  descY <- bnlearn::descendants(avg_graph, Y)

  # Simulation container (one element per node)
  sim        <- vector("list", length(topo))
  names(sim) <- topo

  # ── Pass 1: non-descendants of Y; set Y by the intervention ────────────────
  for (nd in topo) {

    # Intervention node Y: sample from the specified region
    if (nd == Y) {
      sim[[Y]] <- draw_Y_empirical(dat, Y, region, ay, n)
      next
    }

    # Descendants of Y: generated in Pass 2
    if (nd %in% descY) next

    # All other nodes: generate from their fitted conditional model
    model <- fitted_list[[nd]]
    pa    <- bnlearn::parents(avg_graph, nd)

    if (length(pa) > 0) {
      par_df <- as.data.frame(sim[pa], stringsAsFactors = FALSE)
      # Restore factor levels for discrete parent nodes
      for (v in intersect(pa, discrete_vars)) {
        par_df[[v]] <- factor(par_df[[v]], levels = levels(dat[[v]]))
      }
    } else par_df <- NULL

    if (inherits(model, "lm")) {
      # Continuous node with parents: linear predictor + Gaussian residual noise
      mu        <- predict(model, newdata = par_df)
      sim[[nd]] <- mu + rnorm(n, 0, summary(model)$sigma)

    } else if (is.list(model) && all(c("coef", "sd") %in% names(model))) {
      # Continuous root node: unconditional Gaussian
      sim[[nd]] <- rnorm(n, mean = model$coef, sd = model$sd)

    } else if (inherits(model, "multinom")) {
      # Discrete node with parents: multinomial logistic regression
      probs <- predict(model, newdata = par_df, type = "probs")
      # Binary outcome: predict() returns a vector, so build the 2-column matrix
      if (is.null(dim(probs))) {
        lvls  <- levels(factor(dat[[nd]]))
        probs <- cbind(1 - probs, probs)
        colnames(probs) <- lvls
      }
      lvls      <- levels(factor(dat[[nd]]))
      probs     <- probs[, lvls, drop = FALSE]
      # Inverse-CDF sampling: one category per row
      cum       <- t(apply(probs, 1, cumsum))
      idx       <- max.col(runif(n) <= cum, ties.method = "first")
      sim[[nd]] <- factor(lvls[idx], levels = lvls)

    } else if (inherits(model, "table")) {
      # Discrete root node: sample from the empirical marginal distribution
      lvls      <- names(model)
      sim[[nd]] <- factor(sample(lvls, n, TRUE, as.numeric(model)), levels = lvls)

    } else {
      stop(sprintf("Node %s: unsupported model class '%s'",
                   nd, paste(class(model), collapse = ", ")))
    }
  }

  # ── Pass 2: descendants of Y in topological order ──────────────────────────
  # Each node is generated conditional on its parents, which now include
  # the intervened Y (directly or through other descendants)
  for (nd in topo) {
    if (!(nd %in% descY)) next

    model  <- fitted_list[[nd]]
    pa     <- bnlearn::parents(avg_graph, nd)
    par_df <- as.data.frame(sim[pa], stringsAsFactors = FALSE)

    # Restore factor levels for discrete parent nodes
    for (v in intersect(pa, discrete_vars)) {
      par_df[[v]] <- factor(par_df[[v]], levels = levels(dat[[v]]))
    }

    if (inherits(model, "lm")) {
      # Continuous descendant: linear predictor + Gaussian residual noise
      mu        <- predict(model, newdata = par_df)
      sim[[nd]] <- mu + rnorm(n, 0, summary(model)$sigma)

    } else if (inherits(model, "multinom")) {
      # Discrete descendant: multinomial logistic regression
      probs <- predict(model, newdata = par_df, type = "probs")
      if (is.null(dim(probs))) {
        lvls  <- levels(factor(dat[[nd]]))
        probs <- cbind(1 - probs, probs)
        colnames(probs) <- lvls
      }
      lvls      <- levels(factor(dat[[nd]]))
      probs     <- probs[, lvls, drop = FALSE]
      cum       <- t(apply(probs, 1, cumsum))
      idx       <- max.col(runif(n) <= cum, ties.method = "first")
      sim[[nd]] <- factor(lvls[idx], levels = lvls)

    } else {
      stop(sprintf(
        "Descendant node %s must be generated from a conditional model (lm or multinom).", nd
      ))
    }
  }

  as.data.frame(sim, stringsAsFactors = FALSE)
}


# ============================================================================
# draw_Y_empirical()
# Sample n values of Y (with replacement) from the observed data within the
# intervention region:
#   "le" : do(Y <= ay)  observed values at or below the threshold ay
#   "gt" : do(Y >  ay)  observed values above the threshold ay
#   "eq" : do(Y =  ay)  observed values equal to ay (categorical Y)
# Stops with an error if fewer than 10 observations are eligible.
# ============================================================================

draw_Y_empirical <- function(dat, Y, region = c("le", "gt", "eq"), ay, n) {
  region <- match.arg(region)
  ydat   <- dat[[Y]]

  if (region == "eq") {
    # Categorical Y: pool = observations with Y == ay
    ydat <- ydat[!is.na(ydat)]
    pool <- ydat[ydat == ay]
  } else {
    # Continuous Y: pool = observations satisfying the inequality
    ydat <- ydat[is.finite(ydat)]
    pool <- if (region == "le") ydat[ydat <= ay] else ydat[ydat > ay]
  }

  cat("Y=", Y, " region=", region, " ay=", ay,
      " pool=", length(pool), " / n=", n, "\n")

  if (length(pool) < 10) {
    stop("Y intervention pool too small. Try a different ay value or check Y data.")
  }

  sample(pool, n, replace = TRUE)
}

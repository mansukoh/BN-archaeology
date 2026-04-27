# ============================================================================
# Function to generate simulations for multiple Y interventions (all at once)
# ============================================================================

simulate_all_Y_interventions <- function(fitted_list, avg_graph, dat,
                                         Y_vars,  # e.g., c("Subsistence", "Ecology")
                                         ay_list, # named list: list(Subsistence=0.5, Ecology=0.3)
                                         n = 1e6) {
  
  results <- list()
  
  for (Y in Y_vars) {
    ay <- ay_list[[Y]]
    
    if(is.numeric(dat[[Y]])) {
      # 연속형 변수: do(Y <= ay)와 do(Y > ay)
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
      # 범주형 변수: 각 level마다 do(Y = level)
      sim_list <- list()
      for(level in ay) {  # ay는 levels(dat[[Y]])
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


#============================================================================================
# # Simulate the joint distribution of all nodes under the intervention do(Y = ay)
# using the fitted local models in fitted_list.
#
# Simulation proceeds in two passes over the topological ordering of the DAG:
#   Pass 1: generate all non-descendant nodes from their fitted models;
#           fix Y by forced sampling from the intervention region.
#   Pass 2: generate descendant nodes of Y in topological order,
#           conditioning on the already-simulated values of their parents
#           (which now reflect the intervention on Y).
#
# Returns a data frame with n rows and one column per node.
#===========================================================================================


simulate_doY_region <- function(fitted_list, avg_graph, dat,
                                Y, region = c("le", "gt", "eq"), ay, n = 1e6) {
  
  region <- match.arg(region)
  topo   <- bnlearn::node.ordering(avg_graph)   # topological ordering of all nodes
  
  discrete_vars <- c("SettlementType", "CommOrg", "HouseholdOrg")
  
  # Identify descendants of Y; only these need to be regenerated after the intervention
  descY <- bnlearn::descendants(avg_graph, Y)
  
  # Initialise simulation container (one element per node)
  sim        <- vector("list", length(topo))
  names(sim) <- topo
  
  # ── Pass 1: generate non-descendant nodes; fix Y via forced sampling ─────────
  for (nd in topo) {
    
    # Intervention node Y: force-sample from the specified region
    if (nd == Y) {
      sim[[Y]] <- draw_Y_empirical(dat, Y, region, ay, n)
      next
    }
    
    # Descendants of Y: skip for now; generated in Pass 2 using the forced Y value
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
      # Continuous node with parents: Gaussian linear model + residual noise
      mu        <- predict(model, newdata = par_df)
      sim[[nd]] <- mu + rnorm(n, 0, summary(model)$sigma)
      
    } else if (is.list(model) && all(c("coef", "sd") %in% names(model))) {
      # Continuous root node (no parents): unconditional Gaussian
      sim[[nd]] <- rnorm(n, mean = model$coef, sd = model$sd)
      
    } else if (inherits(model, "multinom")) {
      # Discrete node: multinomial logistic regression
      probs <- predict(model, newdata = par_df, type = "probs")
      # Handle binary outcome returned as a vector instead of a matrix
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
      
    } else if (inherits(model, "table")) {
      # Discrete root node: sample from empirical marginal distribution
      lvls      <- names(model)
      sim[[nd]] <- factor(sample(lvls, n, TRUE, as.numeric(model)), levels = lvls)
      
    } else {
      stop(sprintf("Node %s: unsupported model class '%s'",
                   nd, paste(class(model), collapse = ", ")))
    }
  }
  
  # ── Pass 2: generate descendants of Y in topological order ──────────────────
  # These nodes are conditioned on their parents, which now include the forced Y value
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
      # Continuous descendant: Gaussian linear model + residual noise
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


# =============================================================================================
# Sample n values of Y from the observed data under a specified intervention region:
#   "le" : do(Y <= ay)  — sample from observed values at or below threshold ay
#   "gt" : do(Y >  ay)  — sample from observed values above threshold ay
#   "eq" : do(Y =  ay)  — sample from observed values equal to ay (categorical)
# Stops with an error if the eligible pool contains fewer than 10 observations.
# =============================================================================================

draw_Y_empirical <- function(dat, Y, region = c("le", "gt", "eq"), ay, n) {
  region <- match.arg(region)
  ydat   <- dat[[Y]]
  
  if (region == "eq") {
    # Categorical Y: pool = all observations with Y == ay
    ydat <- ydat[!is.na(ydat)]
    pool <- ydat[ydat == ay]
  } else {
    # Continuous Y: pool = observations satisfying the inequality constraint
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
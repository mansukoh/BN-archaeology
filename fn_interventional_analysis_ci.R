
#============================================================================
# Interventional probability computation with 95% confidence intervals
# Handles all combinations of continuous/categorical X and Y variables
# ============================================================================

compute_interventional_analysis_ci <- function(sim_results, avg_graph, dat,
                                               X_vars, ax_list) {
  # Storage for results across all X variables
  all_results <- list()
  
  for (X in X_vars) {
    
    # 1. Find ancestor nodes of X in the DAG (strict ancestors only)
    anc <- bnlearn::ancestors(avg_graph, X)
    anc <- setdiff(anc, X)
    
    if (length(anc) == 0) {
      cat("  No ancestors found for", X, "\n")
      next
    }
    
    # 2. Determine whether X is categorical or continuous
    is_categorical_X <- is.factor(dat[[X]]) || is.character(dat[[X]])
    
    if (is_categorical_X) {
      # Categorical X: compute P(X = level | do(Y)) for each level
      levels_X <- if (is.factor(dat[[X]])) levels(dat[[X]]) else unique(dat[[X]])
      cat(X, "  Variable type: Categorical\n")
      cat("  Levels:", paste(levels_X, collapse = ", "), "\n\n")
    } else {
      # Continuous X: use specified threshold or default to median
      ax <- if (is.null(ax_list) || is.null(ax_list[[X]])) {
        median(dat[[X]], na.rm = TRUE)
      } else {
        ax_list[[X]]
      }
      cat(X, "  Variable type: Continuous\n")
    }
    
    # 3. Compute interventional probabilities for each ancestor Y of X
    X_results <- list()
    
    for (Y in anc) {
      if (!(Y %in% names(sim_results))) {
        cat("    Skipping", Y, "- no simulation results available\n")
        next
      }
      
      # Determine whether Y is categorical or continuous
      is_categorical_Y <- is.factor(dat[[Y]]) || is.character(dat[[Y]])
      
      Y_results <- list()
      
      # ==================================================================
      # Case 1: Y is continuous — compare do(Y <= ay) vs do(Y > ay)
      # ==================================================================
      if (!is_categorical_Y) {
        
        if (is_categorical_X) {
          # X categorical, Y continuous: P(X = level | do(Y <= ay)) and P(X = level | do(Y > ay))
          for (level in levels_X) {
            cat(X, "\n  Level:", level, "\n")
            
            cat(X, Y, "    ")
            prob_le <- compute_interventional_prob_with_ci(
              sim_results, Y, X,
              region_Y    = "le",
              condition_X = "eq",
              ax          = level
            )
            
            cat("    ")
            prob_gt <- compute_interventional_prob_with_ci(
              sim_results, Y, X,
              region_Y    = "gt",
              condition_X = "eq",
              ax          = level
            )
            
            ratio_result <- compute_ratio_ci(prob_le, prob_gt)
            
            Y_results[[level]] <- list(
              prob_le       = prob_le$prob,
              prob_le_lower = prob_le$ci_lower,
              prob_le_upper = prob_le$ci_upper,
              prob_gt       = prob_gt$prob,
              prob_gt_lower = prob_gt$ci_lower,
              prob_gt_upper = prob_gt$ci_upper,
              diff          = prob_gt$prob - prob_le$prob,
              ratio         = ratio_result$ratio,
              ratio_lower   = ratio_result$ratio_lower,
              ratio_upper   = ratio_result$ratio_upper,
              Y_type        = "continuous"
            )
          }
          
        } else {
          # X continuous, Y continuous: P(X > ax | do(Y <= ay)) and P(X > ax | do(Y > ay))
          cat("\n    ")
          prob_le <- compute_interventional_prob_with_ci(
            sim_results, Y, X,
            region_Y    = "le",
            condition_X = "gt",
            ax          = ax
          )
          
          cat("    ")
          prob_gt <- compute_interventional_prob_with_ci(
            sim_results, Y, X,
            region_Y    = "gt",
            condition_X = "gt",
            ax          = ax
          )
          
          ratio_result <- compute_ratio_ci(prob_le, prob_gt)
          
          Y_results[["continuous"]] <- list(
            prob_le       = prob_le$prob,
            prob_le_lower = prob_le$ci_lower,
            prob_le_upper = prob_le$ci_upper,
            prob_gt       = prob_gt$prob,
            prob_gt_lower = prob_gt$ci_lower,
            prob_gt_upper = prob_gt$ci_upper,
            diff          = prob_gt$prob - prob_le$prob,
            ax            = ax,
            ratio         = ratio_result$ratio,
            ratio_lower   = ratio_result$ratio_lower,
            ratio_upper   = ratio_result$ratio_upper,
            Y_type        = "continuous"
          )
        }
        
        # ==================================================================
        # Case 2: Y is categorical — compute P(X | do(Y = k)) for each level k
        # ==================================================================
      } else {
        
        levels_Y <- if (is.factor(dat[[Y]])) levels(dat[[Y]]) else unique(dat[[Y]])
        cat(Y, "\n  Y is categorical with levels:", paste(levels_Y, collapse = ", "), "\n")
        
        if (is_categorical_X) {
          # X categorical, Y categorical: P(X = level_X | do(Y = level_Y))
          for (level_X in levels_X) {
            cat(X, "\n    X level:", level_X, "\n")
            level_results <- list()
            
            for (level_Y in levels_Y) {
              cat(Y, "      Y level:", level_Y, "\n      ")
              prob_Y <- compute_interventional_prob_with_ci_categorical_Y(
                sim_results, Y, X,
                Y_level     = level_Y,
                condition_X = "eq",
                ax          = level_X
              )
              level_results[[level_Y]] <- list(
                prob     = prob_Y$prob,
                ci_lower = prob_Y$ci_lower,
                ci_upper = prob_Y$ci_upper
              )
            }
            
            Y_results[[level_X]] <- list(
              probs_by_Y_level = level_results,
              Y_type           = "categorical"
            )
          }
          
        } else {
          # X continuous, Y categorical: P(X > ax | do(Y = level_Y)) for each level_Y
          level_results <- list()
          
          for (level_Y in levels_Y) {
            cat("\n    Y level:", level_Y, "\n    ")
            prob_Y <- compute_interventional_prob_with_ci_categorical_Y(
              sim_results, Y, X,
              Y_level     = level_Y,
              condition_X = "gt",
              ax          = ax
            )
            level_results[[level_Y]] <- list(
              prob     = prob_Y$prob,
              ci_lower = prob_Y$ci_lower,
              ci_upper = prob_Y$ci_upper
            )
          }
          
          Y_results[["continuous"]] <- list(
            probs_by_Y_level = level_results,
            ax               = ax,
            Y_type           = "categorical"
          )
        }
      }
      
      X_results[[Y]] <- Y_results
    }
    
    all_results[[X]] <- X_results
  }
  
  return(invisible(all_results))
}

# ============================================================================
# Compute P(X | do(Y)) and 95% CI when Y is continuous
# Uses Wilson score interval; falls back to rule-of-three for extreme probabilities
# region_Y:    "le" = do(Y <= ay),  "gt" = do(Y > ay)
# condition_X: "gt" = P(X > ax),   "eq" = P(X = ax),  "le" = P(X <= ax)
# ============================================================================
compute_interventional_prob_with_ci <- function(sim_results, Y, X,
                                                region_Y    = c("le", "gt"),
                                                condition_X = c("gt", "eq", "le"),
                                                ax, alpha = 0.05) {
  
  region_Y    <- match.arg(region_Y)
  condition_X <- match.arg(condition_X)
  
  # Retrieve simulated X values under the specified intervention on Y
  sim_data <- sim_results[[Y]][[region_Y]]
  X_values <- sim_data[[X]]
  n        <- sum(is.finite(X_values))
  
  # Count successes according to the condition on X
  if (condition_X == "gt") {
    successes <- sum(X_values > ax, na.rm = TRUE)
  } else if (condition_X == "eq") {
    if (is.factor(X_values) || is.character(X_values)) {
      successes <- sum(X_values == ax, na.rm = TRUE)
    } else {
      successes <- sum(abs(X_values - ax) < 1e-6, na.rm = TRUE)  # numeric equality with tolerance
    }
  } else {
    successes <- sum(X_values <= ax, na.rm = TRUE)
  }
  
  prob <- successes / n
  
  # 95% CI via Wilson score interval (preferred for proportions near 0 or 1)
  z <- qnorm(1 - alpha / 2)
  
  if (n > 0 && prob > 0 && prob < 1) {
    denominator <- 1 + z^2 / n
    center      <- (prob + z^2 / (2 * n)) / denominator
    margin      <- z * sqrt(prob * (1 - prob) / n + z^2 / (4 * n^2)) / denominator
    ci_lower    <- max(0, center - margin)
    ci_upper    <- min(1, center + margin)
  } else {
    # Boundary cases (prob = 0 or 1): use rule of three
    if (prob == 0) {
      ci_lower <- 0
      ci_upper <- 1 - (alpha / 2)^(1 / n)
    } else if (prob == 1) {
      ci_lower <- (alpha / 2)^(1 / n)
      ci_upper <- 1
    } else {
      ci_lower <- ci_upper <- prob
    }
  }
  
  return(list(
    prob        = prob,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    n           = n,
    Y           = Y,
    X           = X,
    region_Y    = region_Y,
    condition_X = condition_X,
    ax          = ax
  ))
}

# ============================================================================
# Compute P(X | do(Y = k)) and 95% CI when Y is categorical
# Retrieves simulation results for the specific Y level k
# ============================================================================
compute_interventional_prob_with_ci_categorical_Y <- function(sim_results, Y, X,
                                                              Y_level,
                                                              condition_X = c("gt", "eq", "le"),
                                                              ax, alpha = 0.05) {
  
  condition_X <- match.arg(condition_X)
  
  # Retrieve simulated X values under do(Y = Y_level)
  sim_data <- sim_results[[Y]]$simulations[[Y_level]]
  X_values <- sim_data[[X]]
  n        <- sum(is.finite(X_values))
  
  # Count successes according to the condition on X
  if (condition_X == "gt") {
    successes <- sum(X_values > ax, na.rm = TRUE)
  } else if (condition_X == "eq") {
    if (is.factor(X_values) || is.character(X_values)) {
      successes <- sum(X_values == ax, na.rm = TRUE)
    } else {
      successes <- sum(abs(X_values - ax) < 1e-6, na.rm = TRUE)
    }
  } else {
    successes <- sum(X_values <= ax, na.rm = TRUE)
  }
  
  prob <- successes / n
  
  # 95% CI via Wilson score interval
  z <- qnorm(1 - alpha / 2)
  
  if (n > 0 && prob > 0 && prob < 1) {
    denominator <- 1 + z^2 / n
    center      <- (prob + z^2 / (2 * n)) / denominator
    margin      <- z * sqrt(prob * (1 - prob) / n + z^2 / (4 * n^2)) / denominator
    ci_lower    <- max(0, center - margin)
    ci_upper    <- min(1, center + margin)
  } else {
    # Boundary cases: use rule of three
    if (prob == 0) {
      ci_lower <- 0
      ci_upper <- 1 - (alpha / 2)^(1 / n)
    } else if (prob == 1) {
      ci_lower <- (alpha / 2)^(1 / n)
      ci_upper <- 1
    } else {
      ci_lower <- ci_upper <- prob
    }
  }
  
  return(list(
    prob        = prob,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    n           = n,
    Y           = Y,
    X           = X,
    Y_level     = Y_level,
    condition_X = condition_X,
    ax          = ax
  ))
}

# ============================================================================
# Compute ratio P(X > ax | do(Y > ay)) / P(X > ax | do(Y <= ay)) and 95% CI
# Uses log-transformation method for CI; returns NA for degenerate cases
# ============================================================================
compute_ratio_ci <- function(prob_le_result, prob_gt_result, alpha = 0.05) {
  
  prob_le <- prob_le_result$prob
  prob_gt <- prob_gt_result$prob
  n_le    <- prob_le_result$n
  n_gt    <- prob_gt_result$n
  
  # Return NA if inputs are invalid or denominator is effectively zero
  if (is.null(prob_le) || is.null(prob_gt) ||
      !is.finite(prob_le) || !is.finite(prob_gt) ||
      prob_le < 1e-10) {
    return(list(ratio = NA_real_, ratio_lower = NA_real_, ratio_upper = NA_real_))
  }
  
  ratio <- prob_gt / prob_le
  
  # CI via delta method on log scale (valid when both probabilities are interior)
  if (prob_le > 0 && prob_gt > 0 && prob_le < 1 && prob_gt < 1) {
    se_log_ratio <- sqrt(
      (1 - prob_le) / (n_le * prob_le) +
        (1 - prob_gt) / (n_gt * prob_gt)
    )
    z             <- qnorm(1 - alpha / 2)
    log_ratio     <- log(ratio)
    ratio_lower   <- exp(log_ratio - z * se_log_ratio)
    ratio_upper   <- exp(log_ratio + z * se_log_ratio)
  } else {
    # Boundary probabilities: CI not computable via delta method
    ratio_lower <- NA_real_
    ratio_upper <- NA_real_
  }
  
  return(list(
    ratio       = ratio,
    ratio_lower = ratio_lower,
    ratio_upper = ratio_upper
  ))
}

# ============================================================================
# Summarise all interventional results into two data frames:
#   continuous_Y : results where Y is continuous (includes ratio and 95% CI)
#   categorical_Y: results where Y is categorical (includes P(X | do(Y = k)))
# ============================================================================
create_results_table <- function(results, digits = 3) {
  
  result_list_cont_Y <- list()   # Y continuous
  result_list_cat_Y  <- list()   # Y categorical
  
  for (X in names(results)) {
    for (Y in names(results[[X]])) {
      
      Y_data    <- results[[X]][[Y]]
      first_key <- names(Y_data)[1]
      
      # Determine Y type from stored metadata (default to continuous)
      Y_type <- if (!is.null(Y_data[[first_key]]$Y_type)) {
        Y_data[[first_key]]$Y_type
      } else {
        "continuous"
      }
      
      # ================================================================
      # Y continuous: extract prob_le, prob_gt, ratio, and their CIs
      # ================================================================
      if (Y_type == "continuous") {
        
        if ("continuous" %in% names(Y_data)) {
          # X continuous
          result_list_cont_Y[[length(result_list_cont_Y) + 1]] <- data.frame(
            X             = X,
            Y             = Y,
            Y_type        = "continuous",
            Y_value       = NA,
            X_type        = "continuous",
            X_value       = as.character(Y_data$continuous$ax),
            prob_Y_le     = round(Y_data$continuous$prob_le,       digits),
            prob_Y_le_lower = round(Y_data$continuous$prob_le_lower, digits),
            prob_Y_le_upper = round(Y_data$continuous$prob_le_upper, digits),
            prob_Y_gt     = round(Y_data$continuous$prob_gt,       digits),
            prob_Y_gt_lower = round(Y_data$continuous$prob_gt_lower, digits),
            prob_Y_gt_upper = round(Y_data$continuous$prob_gt_upper, digits),
            diff          = round(Y_data$continuous$diff,          digits),
            ratio         = round(Y_data$continuous$ratio,         digits),
            ratio_lower   = round(Y_data$continuous$ratio_lower,   digits),
            ratio_upper   = round(Y_data$continuous$ratio_upper,   digits),
            stringsAsFactors = FALSE
          )
        } else {
          # X categorical: one row per X level
          for (level in names(Y_data)) {
            result_list_cont_Y[[length(result_list_cont_Y) + 1]] <- data.frame(
              X             = X,
              Y             = Y,
              Y_type        = "continuous",
              Y_value       = NA,
              X_type        = "categorical",
              X_value       = as.character(level),
              prob_Y_le     = round(Y_data[[level]]$prob_le,       digits),
              prob_Y_le_lower = round(Y_data[[level]]$prob_le_lower, digits),
              prob_Y_le_upper = round(Y_data[[level]]$prob_le_upper, digits),
              prob_Y_gt     = round(Y_data[[level]]$prob_gt,       digits),
              prob_Y_gt_lower = round(Y_data[[level]]$prob_gt_lower, digits),
              prob_Y_gt_upper = round(Y_data[[level]]$prob_gt_upper, digits),
              diff          = round(Y_data[[level]]$diff,          digits),
              ratio         = round(Y_data[[level]]$ratio,         digits),
              ratio_lower   = round(Y_data[[level]]$ratio_lower,   digits),
              ratio_upper   = round(Y_data[[level]]$ratio_upper,   digits),
              stringsAsFactors = FALSE
            )
          }
        }
        
        # ================================================================
        # Y categorical: extract P(X | do(Y = k)) for each level k
        # ================================================================
      } else {
        
        if ("continuous" %in% names(Y_data)) {
          # X continuous: one row per Y level
          probs_by_Y <- Y_data$continuous$probs_by_Y_level
          for (Y_level in names(probs_by_Y)) {
            result_list_cat_Y[[length(result_list_cat_Y) + 1]] <- data.frame(
              X          = X,
              Y          = Y,
              Y_type     = "categorical",
              Y_value    = Y_level,
              X_type     = "continuous",
              X_value    = as.character(Y_data$continuous$ax),
              prob       = round(probs_by_Y[[Y_level]]$prob,     digits),
              prob_lower = round(probs_by_Y[[Y_level]]$ci_lower, digits),
              prob_upper = round(probs_by_Y[[Y_level]]$ci_upper, digits),
              stringsAsFactors = FALSE
            )
          }
        } else {
          # X categorical: one row per (X level, Y level) combination
          for (X_level in names(Y_data)) {
            probs_by_Y <- Y_data[[X_level]]$probs_by_Y_level
            for (Y_level in names(probs_by_Y)) {
              result_list_cat_Y[[length(result_list_cat_Y) + 1]] <- data.frame(
                X          = X,
                Y          = Y,
                Y_type     = "categorical",
                Y_value    = Y_level,
                X_type     = "categorical",
                X_value    = X_level,
                prob       = round(probs_by_Y[[Y_level]]$prob,     digits),
                prob_lower = round(probs_by_Y[[Y_level]]$ci_lower, digits),
                prob_upper = round(probs_by_Y[[Y_level]]$ci_upper, digits),
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  
  # Combine rows and add ratio direction label for continuous Y results
  df_cont_Y <- if (length(result_list_cont_Y) > 0) {
    df <- do.call(rbind, result_list_cont_Y)
    df$ratio_ind <- ifelse(is.na(df$ratio), "NA",
                           ifelse(df$ratio > 1, "Increase",
                                  ifelse(df$ratio < 1, "Decrease", "No change")))
    df
  } else NULL
  
  df_cat_Y <- if (length(result_list_cat_Y) > 0) {
    do.call(rbind, result_list_cat_Y)
  } else NULL
  
  return(list(
    continuous_Y  = df_cont_Y,
    categorical_Y = df_cat_Y
  ))
}



#========================================================================================
# Plot interventional probability ratios as grouped bar charts 
# Produces three figures saved to figures/:
#   barplot_subsistence.pdf/png  : one panel per Subsistence outcome variable
#   barplot_socialorg.pdf/png    : CommSize, HouseholdOrg, CommOrg side by side
#   barplot_settlementtype.pdf/png: one panel per SettlementType category
#========================================================================================

plot_ratio <- function(results_tables, Subistence_in_network, SocialOrg_in_network) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  
  dat <- as.data.frame(results_tables)
  
  # Bar fill colours: pink = ratio > 1 (increase), steelblue = ratio < 1 (decrease)
  color_map <- c("Increase" = "pink", "Decrease" = "steelblue")
  
  # ── Global font size constants ────────────────────────────────────────────
  BASE_SIZE     <- 14   # base font size for theme_bw()
  AXIS_SIZE     <- 13   # axis tick label size
  LABEL_SIZE    <- 14   # axis title size
  SUBTITLE_SIZE <- 14   # panel subtitle size
  
  # ── Shared bold theme ─────────────────────────────────────────────────────
  bold_theme <- function() {
    theme_bw(base_size = BASE_SIZE) +
      theme(
        axis.text.x   = element_text(angle = 30, hjust = 1,
                                     size = AXIS_SIZE, face = "bold"),
        axis.text.y   = element_text(size = AXIS_SIZE, face = "bold"),
        axis.title.x  = element_text(size = LABEL_SIZE, face = "bold"),
        axis.title.y  = element_text(size = LABEL_SIZE, face = "bold"),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, face = "bold"),
        legend.text   = element_text(size = AXIS_SIZE, face = "bold"),
        legend.title  = element_text(size = LABEL_SIZE, face = "bold"),
        plot.title    = element_blank(),
        legend.position = "none"
      )
  }
  
  # ── Helper: single bar chart for a subset of the results ─────────────────
  # Bars are sorted in descending order of ratio; dashed line marks ratio = 1
  make_bar <- function(sub) {
    sub <- sub %>% arrange(-ratio)
    sub$Y <- factor(sub$Y, levels = sub$Y)
    
    ggplot(sub, aes(x = Y, y = ratio, fill = ratio_ind)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
      scale_fill_manual(values = color_map, name = NULL) +
      labs(x = "Intervention variable",
           y = "Interventional probability ratio") +
      bold_theme()
  }
  
  # ── Helper: extract shared legend from a dummy plot ───────────────────────
  legend_plot <- ggplot(data.frame(ratio_ind = names(color_map)),
                        aes(x = ratio_ind, fill = ratio_ind)) +
    geom_bar() +
    scale_fill_manual(values = color_map, name = NULL) +
    guides(fill = guide_legend(nrow = 1))
  
  get_legend <- function(p) {
    g   <- ggplotGrob(p + theme(
      legend.position = "bottom",
      legend.text  = element_text(size = AXIS_SIZE, face = "bold"),
      legend.title = element_text(size = LABEL_SIZE, face = "bold")
    ))
    idx <- which(grepl("guide-box", sapply(g$grobs, function(x) x$name)))[1]
    g$grobs[[idx]]
  }
  shared_legend <- get_legend(legend_plot)
  
  # ── Helper: attach legend and save figure as PDF and PNG ─────────────────
  save_fig <- function(top_grob, legend, fname, width, height) {
    combined <- arrangeGrob(top_grob, legend, nrow = 2, heights = c(10, 1))
    ggsave(paste0("figures/", fname, ".pdf"), combined,
           width = width, height = height)
    ggsave(paste0("figures/", fname, ".png"), combined,
           width = width, height = height, dpi = 300)
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # Figure 1: Subsistence outcomes
  # One panel per Subsistence variable; panels arranged in 2 rows
  # ══════════════════════════════════════════════════════════════════════════
  vars1  <- Subistence_in_network
  plots1 <- lapply(vars1, function(v) {
    sub <- dplyr::filter(dat, X == v)
    make_bar(sub) + labs(subtitle = paste("Outcome:", v))
  })
  
  top1 <- arrangeGrob(grobs = plots1, nrow = 2)
  save_fig(top1, shared_legend, "barplot_subsistence", width = 14, height = 10)
  
  # ══════════════════════════════════════════════════════════════════════════
  # Figure 2: Social organisation outcomes
  # CommSize always included; HouseholdOrg and CommOrg added if present in DAG
  # Categorical outcomes use dodged bars, one bar per category level
  # ══════════════════════════════════════════════════════════════════════════
  p_comm <- dplyr::filter(dat, X == "CommSize") %>%
    arrange(-ratio) %>%
    mutate(Y = factor(Y, levels = unique(as.character(Y)))) %>%
    ggplot(aes(x = Y, y = ratio, fill = ratio_ind)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(values = color_map, name = NULL) +
    labs(x = "Intervention variable",
         y = "Interventional probability ratio",
         subtitle = "Outcome: CommSize") +
    bold_theme()
  
  plot_list  <- list(p_comm)
  width_list <- c(1)
  
  if ("HouseholdOrg" %in% SocialOrg_in_network) {
    
    # Colour each household organisation type distinctly
    ho_colors <- c("Large extended" = "#E74C3C",
                   "Nuclear"        = "#3498DB",
                   "Small extended" = "#2ECC71")
    
    # Sort Y by mean ratio across household types for a consistent ordering
    ho_sub <- dplyr::filter(dat, X == "HouseholdOrg") %>%
      group_by(Y) %>%
      mutate(mean_ratio = mean(ratio)) %>%
      ungroup() %>%
      arrange(mean_ratio) %>%
      mutate(Y = factor(Y, levels = unique(as.character(Y))))
    
    p_ho <- ggplot(ho_sub, aes(x = Y, y = ratio, fill = X_value, group = X_value)) +
      geom_bar(stat = "identity", position = position_dodge(0.7),
               color = "black", linewidth = 0.3, width = 0.65) +
      # geom_errorbar(aes(ymin = ratio_lower, ymax = ratio_upper),
      #               position = position_dodge(0.7), width = 0.2, linewidth = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
      scale_fill_manual(values = ho_colors, name = "Household type") +
      labs(x = "Intervention variable",
           y = "Interventional probability ratio",
           subtitle = "Outcome: HouseholdOrg") +
      bold_theme() +
      theme(legend.position = "bottom",
            legend.text  = element_text(size = AXIS_SIZE, face = "bold"),
            legend.title = element_text(size = LABEL_SIZE, face = "bold"))
    
    plot_list  <- c(plot_list, list(p_ho))
    width_list <- c(width_list, 1.5)
  }
  
  if ("CommOrg" %in% SocialOrg_in_network) {
    
    # Colour each community organisation type distinctly
    co_colors <- c("Clans"              = "#E74C3C",
                   "No exogamous clans" = "#3498DB")
    
    # Sort Y by mean ratio across community types for a consistent ordering
    co_sub <- dplyr::filter(dat, X == "CommOrg") %>%
      group_by(Y) %>%
      mutate(mean_ratio = mean(ratio)) %>%
      ungroup() %>%
      arrange(mean_ratio) %>%
      mutate(Y = factor(Y, levels = unique(as.character(Y))))
    
    p_co <- ggplot(co_sub, aes(x = Y, y = ratio, fill = X_value, group = X_value)) +
      geom_bar(stat = "identity", position = position_dodge(0.7),
               color = "black", linewidth = 0.3, width = 0.65) +
      # geom_errorbar(aes(ymin = ratio_lower, ymax = ratio_upper),
      #               position = position_dodge(0.7), width = 0.2, linewidth = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
      scale_fill_manual(values = co_colors, name = "Community type") +
      labs(x = "Intervention variable",
           y = "Interventional probability ratio",
           subtitle = "Outcome: CommOrg") +
      bold_theme() +
      theme(legend.position = "bottom",
            legend.text  = element_text(size = AXIS_SIZE, face = "bold"),
            legend.title = element_text(size = LABEL_SIZE, face = "bold"))
    
    plot_list  <- c(plot_list, list(p_co))
    width_list <- c(width_list, 1.5)
  }
  
  fig2 <- arrangeGrob(grobs = plot_list, nrow = 2,
                      layout_matrix = rbind(c(1, 2), c(3, NA)))
  ggsave("figures/barplot_socialorg.pdf", fig2, width = 9,  height = 9)
  ggsave("figures/barplot_socialorg.png", fig2, width = 9,  height = 9, dpi = 300)
  
  
  # ══════════════════════════════════════════════════════════════════════════
  # Figure 3: SettlementType outcomes
  # One panel per settlement category (Village, Hamlet, Homesteads, Camp),
  # arranged in a 2x2 grid
  # ══════════════════════════════════════════════════════════════════════════
  st_vals <- c("Village", "Hamlet", "Homesteads", "Camp")
  plots3  <- lapply(st_vals, function(v) {
    sub <- dplyr::filter(dat, X == "SettlementType", X_value == v) %>%
      arrange(-ratio) %>%
      mutate(Y = factor(Y, levels = unique(as.character(Y))))
    
    ggplot(sub, aes(x = Y, y = ratio, fill = ratio_ind)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6) +
      # geom_errorbar(aes(ymin = ratio_lower, ymax = ratio_upper),
      #               width = 0.2, linewidth = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
      scale_fill_manual(values = color_map, name = NULL) +
      labs(x = "Intervention variable",
           y = "Interventional probability ratio",
           subtitle = paste("Outcome: SettlementType =", v)) +
      bold_theme()
  })
  
  top3 <- arrangeGrob(grobs = plots3, nrow = 2, ncol = 2)
  save_fig(top3, shared_legend, "barplot_settlementtype", width = 13, height = 10)
  
  cat("All figures saved to figures/\n")
}

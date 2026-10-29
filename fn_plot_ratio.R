# ============================================================================
# plot_ratio()
# Bar charts of interventional results (input: output of create_results_table()).
#
# Figures saved to output/ (PDF and PNG):
#   barplot_<X>_<direction> : one figure per outcome X, continuous interventions Y.
#                             Bars = interventional probability ratio
#                             P(. | do(Y > ay)) / P(. | do(Y <= ay)), sorted in
#                             descending order; dashed line at ratio = 1.
#                             Categorical X: bars coloured by level of X
#                             (SettlementType is faceted by level).
#   barplot_<X>_HouseholdOrg: categorical intervention Y = HouseholdOrg;
#                             bars = P(X | do(HouseholdOrg = k)) for each level k.
#
# Arguments
#   results_tables         : list(continuous_Y, categorical_Y)
#   Subistence_in_network,
#   SocialOrg_in_network   : kept for compatibility (not used in the body)
#   direction              : blacklist scenario label used in file names
#
# Requires an existing output/ folder.
# ============================================================================

plot_ratio <- function(results_tables, Subistence_in_network, SocialOrg_in_network, direction) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  
  # results_tables is the list returned by create_results_table():
  #   results_tables$continuous_Y  : Y continuous  (do(Y<=ay) vs do(Y>ay), has `ratio`)
  #   results_tables$categorical_Y : Y categorical (P(X | do(Y=level)) per level, no `ratio`)
  dat     <- as.data.frame(results_tables$continuous_Y)
  dat_cat <- if (!is.null(results_tables$categorical_Y)) {
    as.data.frame(results_tables$categorical_Y)
  } else NULL
  
  # Bar fill colours for increase (ratio > 1) and decrease (ratio < 1);
  # currently the same colour for both
  color_map <- c("Increase" = "steelblue", "Decrease" = "steelblue")
  
  # ── Global font size constants ────────────────────────────────────────────
  BASE_SIZE     <- 13   # base font size for theme_bw()
  AXIS_SIZE     <- 10   # axis tick label size
  LABEL_SIZE    <- 12   # axis title size
  SUBTITLE_SIZE <- 13   # panel subtitle size
  
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
    sub <- sub %>% dplyr::arrange(dplyr::desc(ratio))
    sub$Y <- factor(as.character(sub$Y),
                    levels = unique(as.character(sub$Y)))

    ggplot(sub, aes(x = Y, y = ratio)) +
      geom_col(fill = "steelblue", color = "black",
               linewidth = 0.3, width = 0.6) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
      labs(x = "Intervention variable",
           y = "Interventional probability ratio") +
      bold_theme()
  }

  # ── Helper: extract shared legend from a dummy plot (currently unused) ────
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
  
  # ── Helper: attach legend and save figure as PDF and PNG (currently unused)
  save_fig <- function(top_grob, legend, fname, width, height) {
    combined <- arrangeGrob(top_grob, legend, nrow = 2, heights = c(10, 1))
    ggsave(paste0("output/", fname, ".pdf"), combined,
           width = width, height = height)
    ggsave(paste0("output/", fname, ".png"), combined,
           width = width, height = height, dpi = 300)
  }
  
  # ── Helpers: save one (outcome X, intervention Y) plot as its own file ─────
  # Files are saved as output/barplot_X_Y.pdf and output/barplot_X_Y.png
  safe_name <- function(x) {
    gsub("[^A-Za-z0-9_]+", "_", as.character(x))
  }
  
  save_individual_plot <- function(p, X, Y, width = 6, height = 5) {
    fname <- paste0("barplot_", safe_name(X), "_", safe_name(Y))
    ggsave(paste0("output/", fname, ".pdf"), p, width = width, height = height)
    ggsave(paste0("output/", fname, ".png"), p, width = width, height = height, dpi = 300)
  }
  
  
  # ══════════════════════════════════════════════════════════════════════════
  # Part 1: continuous intervention variables Y
  # One figure per outcome X (all intervention variables Y in one plot).
  # Categorical outcomes: levels shown by colour (dodged bars);
  # SettlementType is faceted by level to keep the bars readable.
  # Files: output/barplot_<X>_<direction>.pdf / .png
  # ══════════════════════════════════════════════════════════════════════════
  if (!is.null(dat) && nrow(dat) > 0) {
    for (x_name in unique(as.character(dat$X))) {
      sub <- dplyr::filter(dat, X == x_name) %>%
        mutate(Y = as.character(Y))
      if (nrow(sub) == 0) next

      x_is_cat <- any(sub$X_type == "categorical", na.rm = TRUE)

      # Arrange intervention variables by descending mean ratio.
      y_order <- sub %>%
        group_by(Y) %>%
        summarise(mean_ratio = mean(ratio, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_ratio)) %>%
        pull(Y)
      sub$Y <- factor(sub$Y, levels = y_order)

      if (x_is_cat) {
        if (x_name == "SettlementType") {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = c("Village", "Homesteads", "Hamlet"))
        } else if (x_name == "HouseholdOrg") {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = c("Large extended", "Small extended", "Nuclear"))
        } else {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = unique(as.character(sub$X_value)))
        }

        p_x <- ggplot(sub, aes(x = Y, y = ratio,
                               fill = X_value, group = X_value)) +
          geom_col(position = position_dodge(0.7),
                   color = "black", linewidth = 0.3, width = 0.65) +
          geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
          labs(x = "Intervention variable",
               y = "Interventional probability ratio",
               subtitle = paste("Outcome:", x_name), fill = x_name) +
          bold_theme() +
          theme(legend.position = "bottom")

        if (x_name == "SettlementType") {
          st_colors <- c("Village" = "#F8766D",
                         "Homesteads" = "#00BA38",
                         "Hamlet" = "#619CFF")
          p_x <- p_x +
            scale_fill_manual(values = st_colors, limits = names(st_colors)) +
            facet_wrap(~ X_value, nrow = 1, drop = FALSE) +
            theme(legend.position = "none")
        } else if (x_name == "HouseholdOrg") {
          ho_colors <- c("Large extended" = "#3498DB",
                         "Small extended" = "#2ECC71",
                         "Nuclear" = "#E74C3C")
          p_x <- p_x +
            scale_fill_manual(values = ho_colors,
                              breaks = names(ho_colors),
                              limits = names(ho_colors),
                              name = "Household type")
        }
      } else {
        p_x <- make_bar(sub) + labs(subtitle = paste("Outcome:", x_name))
      }

      output_name <- paste0("output/barplot_", safe_name(x_name), "_",
                            safe_name(direction))
      n_y <- length(unique(sub$Y))
      plot_width <- max(7, min(15, 1.0 * n_y + 3))
      if (x_name == "SettlementType") plot_width <- max(12, plot_width)
      ggsave(paste0(output_name, ".pdf"), p_x,
             width = plot_width, height = 5.5)
      ggsave(paste0(output_name, ".png"), p_x,
             width = plot_width, height = 5.5, dpi = 300)
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Part 2: categorical intervention variables Y
  # No ratio is defined (there is no le/gt split for a categorical
  # intervention), so each panel shows P(X | do(Y = level)), one bar per
  # level of Y. When X is also categorical, bars are grouped by X level and
  # coloured by Y level.
  # A panel is built for every (X, Y) combination, but only panels with
  # Y = HouseholdOrg are saved: output/barplot_<X>_HouseholdOrg.pdf / .png
  # ══════════════════════════════════════════════════════════════════════════
  if (!is.null(dat_cat) && nrow(dat_cat) > 0) {
    
    # Fixed palette for Y levels (the categorical intervention's categories);
    # extended automatically if more levels are present than pre-defined here.
    y_level_palette <- c("#3498DB", "#E74C3C", "#2ECC71", "#F1C40F", "#9B59B6", "#1ABC9C")
    
    make_cat_panel <- function(sub, x_name, y_name) {
      # Set intervention levels and colours explicitly when Y is HouseholdOrg.
      if (y_name == "HouseholdOrg") {
        y_levels <- c("Large extended", "Small extended", "Nuclear")
        pal <- c("Large extended" = "#3498DB",
                 "Small extended" = "#2ECC71",
                 "Nuclear" = "#E74C3C")
      } else {
        y_levels <- unique(as.character(sub$Y_value))
        pal <- setNames(rep(y_level_palette,
                            length.out = length(y_levels)), y_levels)
      }
      sub$Y_value <- factor(as.character(sub$Y_value), levels = y_levels)
      x_is_cat <- any(sub$X_type == "categorical", na.rm = TRUE)

      if (x_is_cat) {
        if (x_name == "SettlementType") {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = c("Village", "Homesteads", "Hamlet"))
        } else if (x_name == "HouseholdOrg") {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = c("Large extended", "Small extended", "Nuclear"))
        } else {
          sub$X_value <- factor(as.character(sub$X_value),
                                levels = unique(as.character(sub$X_value)))
        }

        p <- ggplot(sub, aes(x = X_value, y = prob,
                             fill = Y_value, group = Y_value)) +
          geom_col(position = position_dodge(0.7),
                   color = "black", linewidth = 0.3, width = 0.65) +
          scale_fill_manual(values = pal, breaks = y_levels,
                            name = y_name) +
          labs(x = y_name, y = "Interventional probability",
               subtitle = if (x_name == "SettlementType")
                 "Outcome: Settlement Type" else paste0("Outcome: ", x_name)) +
          bold_theme() +
          theme(legend.position = "bottom",
                legend.text = element_text(size = AXIS_SIZE, face = "bold"),
                legend.title = element_text(size = LABEL_SIZE, face = "bold"))

        if (x_name == "SettlementType") {
          label_data <- sub %>%
            group_by(X_value) %>%
            summarise(label_y = max(prob_upper, na.rm = TRUE) + 0.025,
                      .groups = "drop")
          p <- p +
            geom_text(data = label_data,
                      aes(x = X_value, y = label_y, label = X_value),
                      inherit.aes = FALSE, size = 4, fontface = "bold") +
            scale_x_discrete(labels = c("Village" = "",
                                        "Homesteads" = "", "Hamlet" = "")) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
        }
        return(p)
      }

      # Numeric outcome under a categorical intervention: one bar per Y level.
      # For HouseholdOrg the bar order and category-specific colours are fixed.
      ggplot(sub, aes(x = Y_value, y = prob, fill = Y_value)) +
        geom_col(color = "black", linewidth = 0.3, width = 0.6) +
        scale_x_discrete(limits = y_levels) +
        scale_fill_manual(values = pal, limits = y_levels, name = NULL) +
        labs(x = paste("Intervention:", y_name),
             y = "Interventional probability",
             subtitle = paste0("Outcome: ", x_name)) +
        bold_theme() +
        theme(legend.position = "none")
    }

    combos <- unique(dat_cat[, c("X", "Y")])
    plots4  <- lapply(seq_len(nrow(combos)), function(i) {
      sub <- dplyr::filter(dat_cat, X == combos$X[i], Y == combos$Y[i])
      make_cat_panel(sub, combos$X[i], combos$Y[i])
    })
    
    # Only HouseholdOrg interventions (categorical Y) get separate files.
    for (i in seq_len(nrow(combos))) {
      if (as.character(combos$Y[i]) == "HouseholdOrg") {
        save_individual_plot(plots4[[i]],
                             as.character(combos$X[i]),
                             as.character(combos$Y[i]))
      }
    }
    
    # Optional: combine all categorical-Y panels into one figure
    # n_plots <- length(plots4)
    # ncol4   <- min(2, n_plots)
    # nrow4   <- ceiling(n_plots / ncol4)
    
    # top4 <- arrangeGrob(grobs = plots4, ncol = ncol4, nrow = nrow4)
    # ggsave(paste0("output/barplot_categoricalY_", direction, ".pdf"), top4,
    #       width = 7 * ncol4, height = 5 * nrow4)
    # ggsave(paste0("output/barplot_categoricalY_", direction, ".png"), top4,
    #        width = 7 * ncol4, height = 5 * nrow4, dpi = 300)
    
  } else {
    cat("  No categorical_Y results to plot (results_tables$categorical_Y is NULL/empty)\n")
  }
  
  cat("All figures saved to output/\n")
}


#  Sensitivity to arc orientation (models M0, M1, M2, M3)
#
# Procedure
#   Step 1. Run Sections 1-5 to obtain the consensus DAG (M0).
#   Step 2. Identify the arcs of M0 with uncertain direction
#           (direction probability < 0.6).
#   Step 3. For each combination of orientations of these arcs:
#           3-1. Reverse the arcs in M0 and confirm the result is a DAG.
#           3-2. Run Sections 6-8 with this network.
#
# With threshold = 0.3 the uncertain arcs are
#   Gather_Agri -> CommSize  and  CommSize -> Hunt_Agri, giving
#   M0 : consensus DAG
#   M1 : Gather_Agri -> CommSize reversed
#   M2 : CommSize -> Hunt_Agri reversed
#   M3 : both arcs reversed

# ── Step 1: consensus DAG (M0) ───────────────────────────────────────────────
# Taken from the copy saved in Section 5 (e.g. avg_graph_none), so it stays
# the original DAG even after avg_graph is replaced below
avg_graph_original <- get(paste0("avg_graph_", direction))

# ── Step 2: arcs with uncertain direction ────────────────────────────────────
uncertain_arcs <- selected_arcs[selected_arcs$direction < 0.6, ]
uncertain_arcs

# ── Step 3-1: network for one combination ────────────────────────────────────
sens_model <- "M1"     # "M0", "M1", "M2" or "M3"

dag_sens <- avg_graph_original
if (sens_model %in% c("M1", "M3"))
  dag_sens <- reverse.arc(dag_sens, from = "Gather_Agri", to = "CommSize",
                          
                          check.cycles = TRUE)
if (sens_model %in% c("M2", "M3"))
  dag_sens <- reverse.arc(dag_sens, from = "CommSize", to = "Hunt_Agri",
                          check.cycles = TRUE)

# Confirm the result is still a DAG
# (reverse.arc() also stops with an error if a cycle is created)
acyclic(dag_sens, directed = TRUE)

avg_graph <- dag_sens

# Label for output files of Sections 7-8 (e.g. "none_M1"), so results of
# M1-M3 do not overwrite those of M0
out_label <- if (sens_model == "M0") direction else paste0(direction, "_", sens_model)

# ── Step 3-2: run Sections 6-8 with avg_graph ────────────────────────────────
# Repeat Step 3 for M0, M1, M2 and M3, then compare the BIC files
# output/BIC_<out_label>_thr3.csv (lower is better).

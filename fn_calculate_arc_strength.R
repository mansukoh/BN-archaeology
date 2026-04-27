# Calculate arc strength and direction probability from bootstrap networks
calculate_arc_strength <- function(boot_nets, nodes) {
  
  # Initialise results data frame: one row per ordered arc (from -> to)
  result <- data.frame(
    from      = character(),
    to        = character(),
    strength  = numeric(),   # proportion of bootstrap nets containing this edge (either direction)
    direction = numeric(),   # proportion of edge occurrences pointing from -> to
    stringsAsFactors = FALSE
  )
  
  # Iterate over all unordered node pairs (each pair processed exactly once)
  for (i in 1:(length(nodes) - 1)) {
    for (j in (i + 1):length(nodes)) {
      node1 <- nodes[i]
      node2 <- nodes[j]
      
      # Weighted occurrence counts in each direction across bootstrap replicates
      n12 <- 0   # node1 -> node2
      n21 <- 0   # node2 -> node1
      
      for (b in 1:length(boot_nets)) {
        
        i1 <- i2 <- 0
        arcs <- boot_nets[[b]]$arcs
        
        if (nrow(arcs) > 0) {
          if (any(arcs[, "from"] == node1 & arcs[, "to"] == node2)) i1 <- 1
          if (any(arcs[, "from"] == node2 & arcs[, "to"] == node1)) i2 <- 1
          
          if (i1 + i2 == 1) {
            # Unambiguous direction: full count assigned to the observed direction
            n12 <- n12 + i1
            n21 <- n21 + i2
          } else if (i1 + i2 == 2) {
            # Both directions present in this replicate (conflicting orientations):
            # split contribution equally
            n12 <- n12 + 0.5
            n21 <- n21 + 0.5
          }
          # i1 + i2 == 0: edge absent in this replicate; no contribution to counts
        }
      }
      
      total <- n12 + n21
      if (total > 0) {
        strength <- total / length(boot_nets)
        
        # Append row for node1 -> node2
        result <- rbind(result, data.frame(
          from      = node1,
          to        = node2,
          strength  = strength,
          direction = n12 / total,   # P(node1 -> node2 | edge present)
          stringsAsFactors = FALSE
        ))
        
        # Append row for node2 -> node1 (direction is the complement)
        result <- rbind(result, data.frame(
          from      = node2,
          to        = node1,
          strength  = strength,
          direction = n21 / total,   # P(node2 -> node1 | edge present)
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(result)
}
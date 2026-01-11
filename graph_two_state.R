# graph_two_state_basic.R


library(gfpop)
library(igraph)

gf_graph <- gfpop::graph(
  gfpop::Edge("normal", "normal", type = "null"),
  gfpop::Edge("normal", "normal", type = "std", penalty = 5)
)

nodes <- unique(c(gf_graph$state1, gf_graph$state2))

edges <- data.frame(
  from  = gf_graph$state1,
  to    = gf_graph$state2,
  label = gf_graph$type,
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(
  edges,
  vertices = data.frame(name = nodes),
  directed = TRUE
)

plot(
  g,
  vertex.size = 40,
  vertex.label = V(g)$name,
  edge.label = E(g)$label,
  main = "Basic two-state GFPOP graph"
)

out_path <- file.path(Sys.getenv("HOME"), "Desktop", "graph_two_state_basic.png")

png(out_path, width = 800, height = 600)
plot(
  g,
  vertex.size = 40,
  vertex.label = V(g)$name,
  edge.label = E(g)$label,
  main = "Basic two-state GFPOP graph"
)
dev.off()

cat("Graph saved to:", out_path, "\n")

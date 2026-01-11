# graph2.R
# GFPOP test graph 2 using default choices

library(gfpop)
library(igraph)

# Define GFPOP graph
gf_graph <- gfpop::graph(
  gfpop::Edge("start", "normal", type = "start"),
  gfpop::Edge("normal", "up", type = "up", penalty = 3),
  gfpop::Edge("up", "normal", type = "down", penalty = 2),
  gfpop::Edge("normal", "end", type = "end")
)

# Convert to igraph
nodes <- unique(c(gf_graph$state1, gf_graph$state2))
edges <- data.frame(
  from  = gf_graph$state1,
  to    = gf_graph$state2,
  label = gf_graph$type,
  stringsAsFactors = FALSE
)
g <- graph_from_data_frame(edges, vertices = data.frame(name=nodes), directed=TRUE)

# Plot in RStudio first
plot(
  g,
  vertex.size=40,
  vertex.label=V(g)$name,
  edge.label=E(g)$label,
  main="GFPOP Graph 2"
)

desktop_results <- file.path(Sys.getenv("HOME"), "Desktop", "GFPOP_results")
if(!dir.exists(desktop_results)) dir.create(desktop_results)

png(file.path(desktop_results, "graph2.png"), width=800, height=600)
plot(
  g,
  vertex.size=40,
  vertex.label=V(g)$name,
  edge.label=E(g)$label,
  main="GFPOP Graph 2"
)
dev.off()

cat("Graph 2 saved to Desktop/GFPOP_results/graph2.png\n")
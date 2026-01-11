
# GFPOP Enhanced Test Script


# Load libraries
if(!require(gfpop)) install.packages("gfpop", repos=NULL, type="source")
library(gfpop)

if(!require(igraph)) install.packages("igraph")
library(igraph)

# Create a GFPOP graph with more states and rules
graph_complex <- gfpop::graph(
  gfpop::Edge("normal", "normal", type="null"),
  gfpop::Edge("normal", "normal", type="std", penalty=5),
  gfpop::Edge("normal", "up", type="up", penalty=10),
  gfpop::Edge("normal", "down", type="down", penalty=10),
  gfpop::Edge("up", "up", type="std", penalty=3),
  gfpop::Edge("up", "normal", type="std"),
  gfpop::Edge("down", "down", type="std", penalty=3),
  gfpop::Edge("down", "normal", type="std"),
  gfpop::Edge("normal", "noChange", type="std"),
  gfpop::Edge("noChange", "noChange", type="null")
)

# Convert to igraph object
edges <- data.frame(
  from = graph_complex$state1,
  to   = graph_complex$state2,
  label = graph_complex$type,
  stringsAsFactors = FALSE
)

g <- igraph::graph_from_data_frame(edges, directed = TRUE)

# Plot the graph in R
plot(
  g,
  vertex.color="lightgreen",
  vertex.size=40,
  vertex.label=V(g)$name,
  vertex.label.color="black",
  edge.arrow.size=0.5,
  edge.label=E(g)$label,
  main="Enhanced GFPOP Graph"
)

# Save the plot as PNG on Desktop
png(filename="~/Desktop/gfpop_graph_enhanced.png", width=1000, height=800)
plot(
  g,
  vertex.color="lightgreen",
  vertex.size=40,
  vertex.label=V(g)$name,
  vertex.label.color="black",
  edge.arrow.size=0.5,
  edge.label=E(g)$label,
  main="Enhanced GFPOP Graph"
)
dev.off()

cat("Enhanced plot saved to Desktop as gfpop_graph_enhanced.png\n")
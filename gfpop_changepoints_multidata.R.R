# Load necessary packages
library(gfpop)
library(ggplot2)

# Define edges
e1 <- Edge("normal", "normal", "null", penalty = 0, decay = 1)
e2 <- Edge("normal", "noChange", "std", penalty = 1, decay = 1)
e3 <- Edge("noChange", "normal", "null", penalty = 0, decay = 1)

# Create graph
mygraph <- graph(e1, e2, e3)

# Example datasets
datasets <- list(
  data1 = c(1, 2, 2, 3, 5, 6, 5, 7),
  data2 = c(1, 1, 2, 2, 3, 7, 8, 8, 5, 4),
  data3 = c(5, 5, 5, 6, 7, 7, 6, 5, 5, 6)
)

# Desktop path (adjust if needed)
desktop_path <- "~/Desktop/"

# Loop over datasets
for(name in names(datasets)) {
  data <- datasets[[name]]
  
  # Run gfpop
  res <- gfpop(data, mygraph, type = "mean")
  
  # Create data frame for plotting
  plot_df <- data.frame(
    x = 1:length(data),
    y = data
  )
  
  seg_df <- data.frame(
    x = res$changepoints,
    y = res$parameters
  )
  
  # Plot
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_vline(xintercept = res$changepoints, color = "red", linetype = "dashed") +
    geom_point(data = seg_df, aes(x = x, y = y), color = "green", size = 3) +
    ggtitle(paste("GFPOP changepoints -", name)) +
    theme_minimal()
  
  # Save plot to Desktop
  ggsave(filename = paste0(desktop_path, "gfpop_plot_", name, ".png"), plot = p, width = 6, height = 4)
}
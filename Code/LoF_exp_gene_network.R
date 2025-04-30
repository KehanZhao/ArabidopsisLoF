# Libraries ---------------------------------------------------------------
library(data.table)
library(igraph)
library(pbapply)
library(ggrepel)

# Load data
candidates <- read.csv("/path/to/S5_LoF_burden_expression_significant_associations.csv")
candidates <- as.data.table(candidates)
candidates[, NLoF := .N, by = LoF.gene]
highlight_genes <- c(candidates$Exp.gene[which(candidates$LoF.gene=="AT4G00650")],"AT4G00650") # highlight genes associated with FRI

# Function to plot ---------------------------------------------------------------
# This function builds and plots a directed gene network from LoF-Exp candidate associations, optionally highlighting selected genes.
make_lof_exp_network <- function(candidates, highlight_genes = NULL) {
  
  # Create the links data frame
  links <- data.frame(
    source = candidates$LoF.gene,
    target = candidates$Exp.gene,
    importance = -log10(candidates$p.raw),
    dir = ifelse(candidates$Beta.coefficient > 0, "A", "R")
  )
  
  # Turn it into an igraph object
  network <- graph_from_data_frame(d = links, directed = TRUE)
  
  # Define vertex types: "source", "target", or "both"
  vertex_types <- sapply(V(network)$name, function(vertex) {
    is_source <- vertex %in% links$source
    is_target <- vertex %in% links$target
    if (is_source & is_target) {
      "both"
    } else if (is_source) {
      "source"
    } else {
      "target"
    }
  })
  
  # Set vertex shape based on vertex type (LoF vs Exp gene)
  V(network)$shape <- ifelse(vertex_types == "source", "circle",
                             ifelse(vertex_types == "target", "square", "circle"))
  
  # Highlight specific genes: set a larger size and a different color
  if (!is.null(highlight_genes)) {
    V(network)$color <- ifelse(V(network)$name %in% highlight_genes, "orange", "gray90")
    V(network)$size <- ifelse(V(network)$name %in% highlight_genes, 4, 0.5) # Larger size for highlighted genes
    V(network)$label <- ifelse(V(network)$name %in% highlight_genes, V(network)$name, NA) # Label highlighted genes
  } else {
    V(network)$color <- "gray90"
    V(network)$size <- 0.5
    V(network)$label <- NA
  }
  
  # Calculate and capture the default layout used by igraph
  layout <- layout_nicely(network)  # This captures the default layout that igraph would use
  
  # Open a PDF device to save the plots
  pdf("/path/to/LoF_exp_network.pdf", width=5, height=5)
  
  # Plot the network with labels
  plot(network,
       layout = layout,               # Use the captured layout
       vertex.label = V(network)$label,  # Display labels for highlighted genes
       vertex.size = V(network)$size,    # Adjust size
       vertex.color = V(network)$color,  # Adjust color
       edge.color = ifelse(links$dir == "A", adjustcolor("#e9a3c9", alpha.f = 1), adjustcolor("#04b7ec", alpha.f = 1)),  # Edge color
       edge.width = .25,                 # Edge width
       edge.arrow.size = 0,              # Arrow size
       edge.arrow.width = 0,             # Arrow width
       edge.curved = 0.3,                # Curved edges
       edge.alpha = 0.5                  # Edge transparency
  )
  
  # Plot the network without labels
  plot(network,
       layout = layout,               # Use the same captured layout
       vertex.label = NA,                # No labels
       vertex.size = V(network)$size,    # Adjust size
       vertex.color = V(network)$color,  # Adjust color
       edge.color = ifelse(links$dir == "A", adjustcolor("#e9a3c9", alpha.f = 1), adjustcolor("#04b7ec", alpha.f = 1)),  # Edge color
       edge.width = .25,                 # Edge width
       edge.arrow.size = 0,              # Arrow size
       edge.arrow.width = 0,             # Arrow width
       edge.curved = 0.3,                # Curved edges
       edge.alpha = 0.5                  # Edge transparency
  )
  
  # Close the PDF device
  dev.off()
}

# Make network graph ------------------------------------------------------
make_lof_exp_network(candidates[NLoF>=10], highlight_genes = highlight_genes)

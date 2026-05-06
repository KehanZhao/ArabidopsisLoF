# LoF_exp_gene_network.R
# Builds and plots a directed gene network from significant LoF-expression
# associations. Edges go from LoF gene (source) to expression gene (target),
# colored by direction of effect (activation vs. repression).
# Outputs two PDF plots: one with gene labels for highlighted genes, one without.

library(data.table)
library(igraph)
library(pbapply)
library(ggrepel)

# Load and annotate significant associations ------------------------------------

candidates <- as.data.table(read.csv("/path/to/Supplementary Table S2.csv"))

# Count how many expression genes are associated with each LoF gene
candidates[, NLoF := .N, by = LoF.gene]

# Genes to highlight: all expression targets of FRI (AT4G00650) + FRI itself
highlight_genes <- c(candidates$Exp.gene[candidates$LoF.gene == "AT4G00650"], "AT4G00650")

# Network plotting function -----------------------------------------------------
# Builds an igraph directed network and saves two versions to PDF:
#   1. With gene labels on highlighted nodes
#   2. Without any labels
#
# Arguments:
#   candidates      - data.table with columns: LoF.gene, Exp.gene, p.raw, Beta.coefficient
#   highlight_genes - character vector of gene IDs to highlight (optional)

make_lof_exp_network <- function(candidates, highlight_genes = NULL) {

  # Build edge table
  links <- data.frame(
    source     = candidates$LoF.gene,
    target     = candidates$Exp.gene,
    importance = -log10(candidates$p.raw),
    dir        = ifelse(candidates$Beta.coefficient > 0, "A", "R")  # A=activation, R=repression
  )

  network <- graph_from_data_frame(d = links, directed = TRUE)

  # Classify vertices as source-only, target-only, or both
  vertex_types <- sapply(V(network)$name, function(vertex) {
    is_source <- vertex %in% links$source
    is_target <- vertex %in% links$target
    if (is_source & is_target) "both"
    else if (is_source)        "source"
    else                       "target"
  })

  # Vertex shape: circle = LoF (source), square = expression (target)
  V(network)$shape <- ifelse(vertex_types == "source", "circle",
                      ifelse(vertex_types == "target", "square", "circle"))

  # Vertex appearance: highlight selected genes in orange
  if (!is.null(highlight_genes)) {
    V(network)$color <- ifelse(V(network)$name %in% highlight_genes, "orange", "gray90")
    V(network)$size  <- ifelse(V(network)$name %in% highlight_genes, 4, 0.5)
    V(network)$label <- ifelse(V(network)$name %in% highlight_genes, V(network)$name, NA)
  } else {
    V(network)$color <- "gray90"
    V(network)$size  <- 0.5
    V(network)$label <- NA
  }

  # Compute layout once and reuse for both plots (ensures identical positioning)
  layout <- layout_nicely(network)

  # Edge colors: pink = activation, blue = repression
  edge_colors <- ifelse(
    links$dir == "A",
    adjustcolor("#e9a3c9", alpha.f = 1),
    adjustcolor("#04b7ec", alpha.f = 1)
  )

  # Shared plot parameters
  plot_args <- list(
    x            = network,
    layout       = layout,
    vertex.size  = V(network)$size,
    vertex.color = V(network)$color,
    edge.color   = edge_colors,
    edge.width   = 0.25,
    edge.arrow.size  = 0,
    edge.arrow.width = 0,
    edge.curved  = 0.3,
    edge.alpha   = 0.5
  )

  pdf("/path/to/LoF_exp_network.pdf", width = 5, height = 5)

  # Plot 1: with labels on highlighted genes
  do.call(plot, c(plot_args, list(vertex.label = V(network)$label)))

  # Plot 2: no labels
  do.call(plot, c(plot_args, list(vertex.label = NA)))

  dev.off()
}

# Generate network plot ---------------------------------------------------------
# Filter to LoF genes with >= 10 significant expression associations

make_lof_exp_network(candidates[NLoF >= 10], highlight_genes = highlight_genes)

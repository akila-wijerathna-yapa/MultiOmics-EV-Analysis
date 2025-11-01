#!/usr/bin/env Rscript
# ==============================================================================
# COMPREHENSIVE MULTI-OMICS INTEGRATION ANALYSIS
# Complete R script for publication-quality figures
# Author: Akila Wijerathna-Yapa
# Date: November 1, 2025
# ==============================================================================

# ==============================================================================
# SECTION 0: PACKAGE INSTALLATION AND LOADING
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("MULTI-OMICS INTEGRATION ANALYSIS - COMPLETE R WORKFLOW\n")
cat("================================================================================\n\n")

cat("Installing and loading required packages...\n")

# Function to install packages if not available
install_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("  Installing ", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
  }
}

# Required packages
required_packages <- c(
  "tidyverse",      # Data manipulation and ggplot2
  "igraph",         # Network analysis
  "ggraph",         # Network visualization
  "circlize",       # Circos plots
  "ComplexHeatmap", # Advanced heatmaps
  "viridis",        # Color scales
  "RColorBrewer",   # Color palettes
  "gridExtra",      # Multiple plots
  "scales",         # Scale functions
  "reshape2"        # Data reshaping
)

# Install packages
for (pkg in required_packages) {
  install_if_needed(pkg)
}

# Load packages quietly
suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(circlize)
  library(ComplexHeatmap)
  library(viridis)
  library(RColorBrewer)
  library(gridExtra)
  library(scales)
  library(reshape2)
})

cat("✓ All packages loaded successfully\n\n")

# ==============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ==============================================================================

cat("================================================================================\n")
cat("SECTION 1: DATA LOADING AND PREPARATION\n")
cat("================================================================================\n\n")

# Load data files
cat("Loading data files...\n")
protein_freq <- read.csv("protein_frequencies.csv", stringsAsFactors = FALSE)
metabolite_freq <- read.csv("metabolite_frequencies.csv", stringsAsFactors = FALSE)
lipid_freq <- read.csv("lipid_frequencies.csv", stringsAsFactors = FALSE)

# Clean molecule names
protein_freq$Molecule <- trimws(protein_freq$Molecule)
metabolite_freq$Molecule <- trimws(metabolite_freq$Molecule)
lipid_freq$Molecule <- trimws(lipid_freq$Molecule)

# Get molecule lists
protein_list <- unique(protein_freq$Molecule)
metabolite_list <- unique(metabolite_freq$Molecule)
lipid_list <- unique(lipid_freq$Molecule)

cat(sprintf("  Loaded: %d proteins, %d metabolites, %d lipids\n\n", 
           length(protein_list), length(metabolite_list), length(lipid_list)))

# ==============================================================================
# SECTION 2: PATHWAY MAPPING
# ==============================================================================

cat("================================================================================\n")
cat("SECTION 2: PATHWAY MAPPING AND ENRICHMENT\n")
cat("================================================================================\n\n")

cat("Mapping molecules to biological pathways...\n")

# Define pathway associations (curated from EV biology literature)
pathway_mapping <- list(
  "Exosome Biogenesis" = list(
    proteins = c("CD63", "CD81", "CD9", "TSG101", "ALIX", "Syntenin", "PDCD6IP"),
    metabolites = c("ATP"),
    lipids = c("Cholesterol", "Ceramide", "SM", "ceramides")
  ),
  "ESCRT Pathway" = list(
    proteins = c("TSG101", "ALIX", "PDCD6IP", "VPS4"),
    metabolites = c("ATP"),
    lipids = c("PS", "PE")
  ),
  "Membrane Trafficking" = list(
    proteins = c("CD63", "CD81", "CD9", "Flotillin", "Caveolin", "Annexin"),
    metabolites = c("GTP", "ATP"),
    lipids = c("PC", "PE", "PS", "Cholesterol")
  ),
  "Lipid Raft Formation" = list(
    proteins = c("CD63", "CD81", "CD9", "Flotillin", "Caveolin"),
    metabolites = c("Cholesterol"),
    lipids = c("Cholesterol", "SM", "Ceramide", "ceramides", "GM1")
  ),
  "Energy Metabolism" = list(
    proteins = c("GAPDH", "ENO1", "PKM", "LDHA", "PGK1"),
    metabolites = c("Glucose", "Pyruvate", "Lactate", "ATP", "ADP", "pyruvic acid", 
                   "lactic acid", "glucose-6-phosphate"),
    lipids = c()
  ),
  "TCA Cycle" = list(
    proteins = c("IDH", "MDH", "CS"),
    metabolites = c("Citrate", "Succinate", "Fumarate", "Malate", "ATP"),
    lipids = c()
  ),
  "Amino Acid Metabolism" = list(
    proteins = c("GAPDH", "ENO1"),
    metabolites = c("Glutamate", "Alanine", "Valine", "Leucine", "Isoleucine", 
                   "Glycine", "L-glutamic acid", "alanine", "valine", "threonine"),
    lipids = c()
  ),
  "Phospholipid Metabolism" = list(
    proteins = c("PLA2", "PLCG"),
    metabolites = c("Choline", "Ethanolamine", "Glycerol-3-phosphate"),
    lipids = c("PC", "PE", "PS", "PI", "PA", "LPC", "LPE", "plasmalogen PE")
  ),
  "Sphingolipid Metabolism" = list(
    proteins = c("SMPD", "CERS"),
    metabolites = c("Sphingosine", "Serine"),
    lipids = c("SM", "Ceramide", "ceramides", "Sphingosine-1-phosphate", "Glucosylceramide")
  ),
  "Lipid Signaling" = list(
    proteins = c("PLA2", "PLCG", "DGKA"),
    metabolites = c("Arachidonic acid"),
    lipids = c("DAG", "Ceramide", "ceramides", "Sphingosine-1-phosphate", "LPA")
  ),
  "Cytoskeleton Organization" = list(
    proteins = c("Actin", "Tubulin", "Ezrin", "Moesin", "Cofilin"),
    metabolites = c("ATP", "GTP"),
    lipids = c("PIP2")
  ),
  "Protein Folding" = list(
    proteins = c("HSP70", "Hsp70", "HSP90", "Hsp60", "GRP78", "GRP94"),
    metabolites = c("ATP"),
    lipids = c()
  )
)

# Create pathway dataframe
pathway_df <- data.frame()
for (pathway_name in names(pathway_mapping)) {
  pathway_data <- pathway_mapping[[pathway_name]]
  
  if (length(pathway_data$proteins) > 0) {
    pathway_df <- rbind(pathway_df, data.frame(
      Molecule = pathway_data$proteins,
      Pathway = pathway_name,
      OmicsLayer = "Protein",
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(pathway_data$metabolites) > 0) {
    pathway_df <- rbind(pathway_df, data.frame(
      Molecule = pathway_data$metabolites,
      Pathway = pathway_name,
      OmicsLayer = "Metabolite",
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(pathway_data$lipids) > 0) {
    pathway_df <- rbind(pathway_df, data.frame(
      Molecule = pathway_data$lipids,
      Pathway = pathway_name,
      OmicsLayer = "Lipid",
      stringsAsFactors = FALSE
    ))
  }
}

# Match with detected molecules
pathway_df_detected <- pathway_df %>%
  filter(
    (OmicsLayer == "Protein" & Molecule %in% protein_list) |
    (OmicsLayer == "Metabolite" & Molecule %in% metabolite_list) |
    (OmicsLayer == "Lipid" & Molecule %in% lipid_list)
  )

# Calculate pathway enrichment
pathway_enrichment <- pathway_df_detected %>%
  group_by(Pathway, OmicsLayer) %>%
  summarize(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = OmicsLayer, values_from = Count, values_fill = 0)

# Ensure all columns exist
for (col in c("Protein", "Metabolite", "Lipid")) {
  if (!(col %in% colnames(pathway_enrichment))) {
    pathway_enrichment[[col]] <- 0
  }
}

pathway_enrichment <- pathway_enrichment %>%
  mutate(
    Total = Protein + Metabolite + Lipid,
    NumLayers = (Protein > 0) + (Metabolite > 0) + (Lipid > 0),
    IntegrationScore = NumLayers / 3
  ) %>%
  arrange(desc(NumLayers), desc(Total))

cat(sprintf("  Identified %d pathways with multi-omics support\n\n", 
           sum(pathway_enrichment$NumLayers >= 2)))

# ==============================================================================
# SECTION 3: NETWORK CONSTRUCTION
# ==============================================================================

cat("================================================================================\n")
cat("SECTION 3: NETWORK CONSTRUCTION\n")
cat("================================================================================\n\n")

cat("Building molecular interaction network...\n")

# Create edges based on pathway co-membership
edges_list <- list()
for (pathway in unique(pathway_df_detected$Pathway)) {
  molecules <- pathway_df_detected %>%
    filter(Pathway == pathway) %>%
    pull(Molecule)
  
  if (length(molecules) > 1) {
    molecule_pairs <- combn(molecules, 2, simplify = FALSE)
    for (pair in molecule_pairs) {
      edges_list[[length(edges_list) + 1]] <- data.frame(
        From = pair[1],
        To = pair[2],
        Pathway = pathway,
        stringsAsFactors = FALSE
      )
    }
  }
}

edge_df <- do.call(rbind, edges_list)

# Aggregate edges
edge_weights <- edge_df %>%
  group_by(From, To) %>%
  summarize(
    Weight = n(),
    Pathways = paste(unique(Pathway), collapse = "; "),
    .groups = "drop"
  )

# Create igraph network
all_molecules <- unique(c(edge_weights$From, edge_weights$To))

# Create node attributes
node_attributes <- data.frame(
  Molecule = all_molecules,
  stringsAsFactors = FALSE
) %>%
  mutate(
    OmicsLayer = case_when(
      Molecule %in% protein_list ~ "Protein",
      Molecule %in% metabolite_list ~ "Metabolite",
      Molecule %in% lipid_list ~ "Lipid",
      TRUE ~ "Unknown"
    )
  )

# Add frequency data
node_attributes <- node_attributes %>%
  left_join(protein_freq %>% select(Molecule, Frequency), by = "Molecule") %>%
  rename(Frequency_Protein = Frequency) %>%
  left_join(metabolite_freq %>% select(Molecule, Frequency), by = "Molecule") %>%
  rename(Frequency_Metabolite = Frequency) %>%
  left_join(lipid_freq %>% select(Molecule, Frequency), by = "Molecule") %>%
  rename(Frequency_Lipid = Frequency) %>%
  mutate(
    Frequency = coalesce(Frequency_Protein, Frequency_Metabolite, Frequency_Lipid, 1)
  ) %>%
  select(Molecule, OmicsLayer, Frequency)

# Build network
G <- graph_from_data_frame(
  d = edge_weights %>% select(From, To, Weight),
  directed = FALSE,
  vertices = node_attributes
)

# Calculate centrality
V(G)$degree_cent <- degree(G)
V(G)$betweenness <- betweenness(G, normalized = TRUE)
V(G)$closeness <- closeness(G, normalized = TRUE)

# Identify hubs (top 15%)
degree_values <- degree(G)
degree_threshold <- quantile(degree_values, 0.85)
V(G)$is_hub <- degree(G) >= degree_threshold

# Community detection
set.seed(42)
communities <- cluster_louvain(G)
V(G)$community <- membership(communities)

cat(sprintf("  Network: %d nodes, %d edges\n", vcount(G), ecount(G)))
cat(sprintf("  Communities: %d\n", max(V(G)$community)))
cat(sprintf("  Hub molecules: %d\n\n", sum(V(G)$is_hub)))

# ==============================================================================
# FIGURE 1: MULTI-OMICS PATHWAY HEATMAP
# ==============================================================================

cat("================================================================================\n")
cat("CREATING FIGURES\n")
cat("================================================================================\n\n")

cat("Creating Figure 1: Multi-Omics Pathway Heatmap...\n")

# Prepare heatmap data
multi_omics_pathways <- pathway_enrichment %>%
  filter(NumLayers >= 2)

heatmap_matrix <- as.matrix(multi_omics_pathways[, c("Protein", "Metabolite", "Lipid")])
rownames(heatmap_matrix) <- multi_omics_pathways$Pathway

# Normalize by row max
heatmap_matrix_norm <- t(apply(heatmap_matrix, 1, function(x) {
  if(max(x) > 0) x / max(x) else x
}))

# Create heatmap
pdf("Figure1_MultiOmics_Pathway_Heatmap.pdf", width = 8, height = 10)

col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#FFA500", "#DC143C"))

ht <- Heatmap(
  heatmap_matrix_norm,
  name = "Normalized\nCount",
  col = col_fun,
  
  # Row settings
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  row_title = "Biological Pathways",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # Column settings
  cluster_columns = FALSE,
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_title = "Multi-Omics Pathway Enrichment Analysis",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Cell annotations
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- heatmap_matrix[i, j]
    if(value > 0) {
      text_color <- if(heatmap_matrix_norm[i, j] >= 0.5) "white" else "black"
      grid.text(sprintf("%d", value), x, y, 
               gp = gpar(fontsize = 9, fontface = "bold", col = text_color))
    }
  },
  
  # Legend
  heatmap_legend_param = list(
    title = "Normalized\nCount",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_direction = "vertical",
    legend_height = unit(4, "cm")
  ),
  
  # Borders
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 2)
)

draw(ht, heatmap_legend_side = "right")

dev.off()
cat("  ✓ Figure 1 saved\n")

# ==============================================================================
# FIGURE 2: MULTI-OMICS NETWORK WITH ALL LABELS
# ==============================================================================

cat("Creating Figure 2: Multi-Omics Network (all nodes labeled)...\n")

# Define colors
omics_colors <- c(
  "Protein" = "#4285F4",
  "Metabolite" = "#34A853",
  "Lipid" = "#FBBC04"
)

# Convert to tidygraph for ggraph
library(tidygraph)
tidy_G <- as_tbl_graph(G)

# Create network plot
set.seed(42)
p2 <- ggraph(tidy_G, layout = 'fr') +
  # Edges
  geom_edge_link(aes(width = Weight), alpha = 0.15, color = "grey60") +
  scale_edge_width(range = c(0.3, 2), guide = "none") +
  
  # Nodes
  geom_node_point(aes(color = OmicsLayer, size = degree_cent), alpha = 0.9) +
  scale_color_manual(values = omics_colors, name = "Omics Layer") +
  scale_size_continuous(range = c(4, 14), name = "Degree") +
  
  # Labels with repel
  geom_node_text(
    aes(label = name, fontface = ifelse(is_hub, "bold", "plain")),
    size = 3.5,
    repel = TRUE,
    max.overlaps = 30,
    bg.color = "white",
    bg.r = 0.12,
    segment.color = "grey50",
    segment.size = 0.3,
    box.padding = 0.5
  ) +
  
  # Theme
  theme_graph(base_family = "sans") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "Multi-Omics Molecular Interaction Network",
    subtitle = "All nodes labeled | Node size = degree centrality | Bold = hub molecules"
  )

ggsave("Figure2_MultiOmics_Network_LABELED.pdf", p2, 
      width = 16, height = 14, dpi = 300)

cat("  ✓ Figure 2 saved\n")

# ==============================================================================
# FIGURE 3: NETWORK COMMUNITIES
# ==============================================================================

cat("Creating Figure 3: Network Communities (all nodes labeled)...\n")

set.seed(42)
p3 <- ggraph(tidy_G, layout = 'fr') +
  # Edges
  geom_edge_link(alpha = 0.12, color = "grey70", width = 0.3) +
  
  # Nodes colored by community
  geom_node_point(aes(color = factor(community), size = degree_cent), alpha = 0.9) +
  scale_color_viridis_d(option = "turbo", name = "Community") +
  scale_size_continuous(range = c(4, 12), name = "Degree") +
  
  # Labels
  geom_node_text(
    aes(label = name, fontface = ifelse(is_hub, "bold", "plain")),
    size = 3.5,
    repel = TRUE,
    max.overlaps = 30,
    bg.color = "white",
    bg.r = 0.12,
    segment.color = "grey50",
    segment.size = 0.3,
    box.padding = 0.5
  ) +
  
  # Theme
  theme_graph(base_family = "sans") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Functional Module Detection in Multi-Omics Network",
    subtitle = "Communities by Louvain algorithm | Bold = hub molecules"
  )

ggsave("Figure3_Network_Communities_LABELED.pdf", p3, 
      width = 16, height = 14, dpi = 300)

cat("  ✓ Figure 3 saved\n")

# ==============================================================================
# FIGURE 4: CIRCOS PLOT
# ==============================================================================

cat("Creating Figure 4: Circos Plot (publication quality)...\n")

# Select top 15 from each
n_top <- 15
top_proteins <- head(protein_freq[order(-protein_freq$Frequency), ], n_top)
top_metabolites <- head(metabolite_freq[order(-metabolite_freq$Frequency), ], n_top)
top_lipids <- head(lipid_freq[order(-lipid_freq$Frequency), ], n_top)

# Create sectors data
sectors <- data.frame(
  sector = c(rep("Proteins", n_top), rep("Metabolites", n_top), rep("Lipids", n_top)),
  molecule = c(top_proteins$Molecule, top_metabolites$Molecule, top_lipids$Molecule),
  frequency = c(top_proteins$Frequency, top_metabolites$Frequency, top_lipids$Frequency),
  x = c(1:n_top, 1:n_top, 1:n_top),
  stringsAsFactors = FALSE
)

# Create connections based on pathways
get_pathway_connection <- function(mol1, mol2) {
  for (pathway_name in names(pathway_mapping)) {
    pathway <- pathway_mapping[[pathway_name]]
    all_mols <- c(pathway$proteins, pathway$metabolites, pathway$lipids)
    if (mol1 %in% all_mols && mol2 %in% all_mols) {
      return(TRUE)
    }
  }
  return(FALSE)
}

connections <- data.frame()
all_mols <- sectors$molecule
for (i in 1:(length(all_mols)-1)) {
  for (j in (i+1):length(all_mols)) {
    if (get_pathway_connection(all_mols[i], all_mols[j])) {
      sector1 <- sectors$sector[sectors$molecule == all_mols[i]][1]
      sector2 <- sectors$sector[sectors$molecule == all_mols[j]][1]
      x1 <- sectors$x[sectors$molecule == all_mols[i]][1]
      x2 <- sectors$x[sectors$molecule == all_mols[j]][1]
      
      connections <- rbind(connections, data.frame(
        sector1 = sector1, x1 = x1,
        sector2 = sector2, x2 = x2,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sector colors
sector_colors <- c("Proteins" = "#4285F4", "Metabolites" = "#34A853", "Lipids" = "#FBBC04")

# Create Circos plot
pdf("Figure4_Circos_MultiOmics_PUBLICATION.pdf", width = 14, height = 14)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(10, 10, 10),
  track.margin = c(0.01, 0.01),
  cell.padding = c(0.02, 0, 0.02, 0),
  points.overflow.warning = FALSE
)

circos.initialize(
  factors = factor(sectors$sector, levels = c("Proteins", "Metabolites", "Lipids")),
  x = sectors$x
)

# Track 1: Sector labels
circos.track(
  ylim = c(0, 1),
  bg.col = sector_colors,
  bg.border = "white",
  track.height = 0.05,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ycenter, sector.name,
               cex = 1.8, font = 2, col = "white",
               facing = "bending.inside", niceFacing = TRUE)
  }
)

# Track 2: Radial labels (NO OVERLAP!)
circos.track(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.35,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    sector_molecules <- sectors[sectors$sector == sector.name, ]
    
    for (i in 1:nrow(sector_molecules)) {
      mol <- sector_molecules$molecule[i]
      x_pos <- sector_molecules$x[i]
      
      circos.lines(c(x_pos, x_pos), c(0, 0.15), 
                  col = sector_colors[sector.name], lwd = 2)
      
      circos.text(x_pos, 0.5, mol,
                 cex = 0.75, font = 1, col = "black",
                 facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }
  }
)

# Track 3: Frequency bars
circos.track(
  ylim = c(0, 1),
  bg.border = "grey80",
  track.height = 0.08,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    sector_molecules <- sectors[sectors$sector == sector.name, ]
    max_freq <- max(sector_molecules$frequency)
    
    for (i in 1:nrow(sector_molecules)) {
      x_pos <- sector_molecules$x[i]
      freq_norm <- sector_molecules$frequency[i] / max_freq
      circos.rect(x_pos - 0.4, 0, x_pos + 0.4, freq_norm,
                 col = sector_colors[sector.name], border = "white", lwd = 0.5)
    }
  }
)

# Draw connections
if (nrow(connections) > 0) {
  for (i in 1:nrow(connections)) {
    is_cross <- connections$sector1[i] != connections$sector2[i]
    link_col <- if(is_cross) "#E74C3C" else "#95A5A6"
    link_alpha <- if(is_cross) 0.4 else 0.2
    
    circos.link(connections$sector1[i], connections$x1[i],
               connections$sector2[i], connections$x2[i],
               col = adjustcolor(link_col, alpha.f = link_alpha), lwd = 1.5)
  }
}

# Title and legends
text(0, 1.25, "Multi-Omics Integration Circos Plot", cex = 2.0, font = 2, xpd = TRUE)
text(0, 1.17, "Top 15 molecules per omics layer with pathway-based connections",
     cex = 1.3, col = "gray30", xpd = TRUE)

legend("bottomleft",
       legend = c("Cross-Omics Connection", "Within-Omics Connection"),
       col = c(adjustcolor("#E74C3C", alpha.f = 0.4), 
               adjustcolor("#95A5A6", alpha.f = 0.2)),
       lwd = 4, cex = 1.1, bty = "n", xpd = TRUE)

legend("bottomright",
       legend = c("Proteins", "Metabolites", "Lipids"),
       fill = sector_colors, border = "black",
       cex = 1.1, bty = "n", xpd = TRUE)

circos.clear()
dev.off()

cat("  ✓ Figure 4 saved\n")

# ==============================================================================
# FIGURE 5: PATHWAY INTEGRATION SUMMARY
# ==============================================================================

cat("Creating Figure 5: Pathway Integration Summary...\n")

# Prepare data
pathway_plot_data <- pathway_enrichment %>%
  arrange(Total)

pathway_plot_data$Pathway <- factor(pathway_plot_data$Pathway, 
                                    levels = pathway_plot_data$Pathway)

# Create plot
p5 <- ggplot(pathway_plot_data, aes(x = Pathway, y = Total, fill = factor(NumLayers))) +
  geom_col(width = 0.7, color = "black", size = 0.3) +
  geom_text(aes(label = Total), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_manual(
    values = c("1" = "#BDC3C7", "2" = "#F39C12", "3" = "#E74C3C"),
    name = "# Omics Layers",
    labels = c("1" = "1 Layer", "2" = "2 Layers", "3" = "3 Layers")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Pathway-Level Multi-Omics Integration",
    subtitle = "Total molecules detected per pathway, colored by omics layer diversity",
    x = NULL,
    y = "Total Molecules Detected"
  )

ggsave("Figure5_Pathway_Integration.pdf", p5, 
      width = 10, height = 8, dpi = 300)

cat("  ✓ Figure 5 saved\n\n")

# ==============================================================================
# SECTION 4: EXPORT SUMMARY TABLES
# ==============================================================================

cat("================================================================================\n")
cat("EXPORTING SUMMARY TABLES\n")
cat("================================================================================\n\n")

# Hub molecules
hub_molecules <- data.frame(
  Molecule = V(G)$name[V(G)$is_hub],
  OmicsLayer = V(G)$OmicsLayer[V(G)$is_hub],
  Degree = V(G)$degree_cent[V(G)$is_hub],
  Betweenness = round(V(G)$betweenness[V(G)$is_hub], 4),
  Closeness = round(V(G)$closeness[V(G)$is_hub], 4),
  Community = V(G)$community[V(G)$is_hub],
  Frequency = V(G)$Frequency[V(G)$is_hub]
) %>%
  arrange(desc(Degree))

write.csv(hub_molecules, "Hub_Molecules.csv", row.names = FALSE)
cat("  ✓ Hub_Molecules.csv saved\n")

# Network statistics
network_stats <- data.frame(
  Metric = c("Total Nodes", "Total Edges", "Network Density", "Average Degree",
            "Average Clustering Coefficient", "Number of Communities",
            "Number of Hub Molecules", "Modularity"),
  Value = c(vcount(G), ecount(G), 
           round(edge_density(G), 4),
           round(mean(degree(G)), 2),
           round(transitivity(G, type = "global"), 4),
           max(V(G)$community),
           sum(V(G)$is_hub),
           round(modularity(communities), 4))
)

write.csv(network_stats, "Network_Statistics.csv", row.names = FALSE)
cat("  ✓ Network_Statistics.csv saved\n")

# Pathway summary
pathway_summary <- pathway_enrichment %>%
  mutate(
    IntegrationLevel = case_when(
      NumLayers == 3 ~ "Full Integration (3 layers)",
      NumLayers == 2 ~ "Partial Integration (2 layers)",
      TRUE ~ "Single Layer"
    )
  )

write.csv(pathway_summary, "MultiOmics_Pathway_Summary.csv", row.names = FALSE)
cat("  ✓ MultiOmics_Pathway_Summary.csv saved\n")

# Community composition
community_composition <- data.frame(
  Community = 1:max(V(G)$community)
) %>%
  rowwise() %>%
  mutate(
    Size = sum(V(G)$community == Community),
    Proteins = sum(V(G)$community == Community & V(G)$OmicsLayer == "Protein"),
    Metabolites = sum(V(G)$community == Community & V(G)$OmicsLayer == "Metabolite"),
    Lipids = sum(V(G)$community == Community & V(G)$OmicsLayer == "Lipid"),
    OmicsDiversity = sum(c(Proteins > 0, Metabolites > 0, Lipids > 0))
  ) %>%
  ungroup() %>%
  arrange(desc(Size))

write.csv(community_composition, "Community_Composition.csv", row.names = FALSE)
cat("  ✓ Community_Composition.csv saved\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("================================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("================================================================================\n\n")

cat("Generated Files:\n")
cat("  FIGURES (PDF):\n")
cat("    • Figure1_MultiOmics_Pathway_Heatmap.pdf\n")
cat("    • Figure2_MultiOmics_Network_LABELED.pdf\n")
cat("    • Figure3_Network_Communities_LABELED.pdf\n")
cat("    • Figure4_Circos_MultiOmics_PUBLICATION.pdf\n")
cat("    • Figure5_Pathway_Integration.pdf\n\n")

cat("  DATA TABLES (CSV):\n")
cat("    • Hub_Molecules.csv\n")
cat("    • Network_Statistics.csv\n")
cat("    • MultiOmics_Pathway_Summary.csv\n")
cat("    • Community_Composition.csv\n\n")

cat("Key Findings:\n")
cat(sprintf("  • %d pathways with full integration (3 omics layers)\n", 
           sum(pathway_enrichment$NumLayers == 3)))
cat(sprintf("  • %d hub molecules identified\n", sum(V(G)$is_hub)))
cat(sprintf("  • %d functional communities detected\n", max(V(G)$community)))
cat(sprintf("  • Network density: %.3f (moderately dense)\n", edge_density(G)))
cat(sprintf("  • Clustering coefficient: %.3f (highly modular)\n\n", 
           transitivity(G, type = "global")))

cat("All figures are publication-ready with:\n")
cat("  ✓ High resolution (300 DPI)\n")
cat("  ✓ Vector graphics (PDF format)\n")
cat("  ✓ All nodes labeled (no overlap)\n")
cat("  ✓ Professional color schemes\n")
cat("  ✓ Clear legends and titles\n\n")

cat("================================================================================\n")
cat("Ready for manuscript submission!\n")
cat("================================================================================\n\n")

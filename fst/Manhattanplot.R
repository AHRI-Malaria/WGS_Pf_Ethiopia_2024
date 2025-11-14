setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dplyr)
library(sf)
library(ggforce)
library(ggpattern) 
library(rnaturalearth)
library(geodata)
library(readr)
library(RColorBrewer)
library(ggrepel)

fst <- read.csv("all_fst_annotate_sites.csv", stringsAsFactors = FALSE)

# Prepare data - create continuous x-axis across all chromosomes
fst_plot <- fst %>%
  filter(chr != "Pf3D7_API_v3") %>%
  mutate(
    Population_1 = ifelse(Population_1 == "Gonder", "Gondar", Population_1),
    Population_2 = ifelse(Population_2 == "Gonder", "Gondar", Population_2),
    Population_1 = ifelse(Population_1 == "Metehera_Fentale", "Metehera Fentale", Population_1),
    Population_2 = ifelse(Population_2 == "Metehera_Fentale", "Metehera Fentale", Population_2),
    Population_1 = ifelse(Population_1 == "Arbaminch", "Arba Minch", Population_1),
    Population_2 = ifelse(Population_2 == "Arbaminch", "Arba Minch", Population_2),
  ) %>%
  filter(Population_1 %in% districts$site & Population_2 %in% districts$site) %>%
  mutate(
    chr_num = as.numeric(gsub("Pf3D7_", "", gsub("_v3", "", chr))),
    Comparison = paste(Population_1, "vs", Population_2)
  ) %>%
  arrange(chr_num, pos) %>%
  group_by(chr_num) %>%
  mutate(pos_index = row_number()) %>%
  ungroup()

# Calculate cumulative position for each chromosome
chr_lengths <- fst_plot %>%
  group_by(chr_num) %>%
  summarise(max_pos = max(pos_index), .groups = 'drop') %>%
  mutate(
    cum_pos = cumsum(as.numeric(max_pos)) - max_pos,
    chr_center = cum_pos + max_pos/2
  )

# Add cumulative position to data
fst_plot <- fst_plot %>%
  left_join(chr_lengths %>% select(chr_num, cum_pos), by = "chr_num") %>%
  mutate(x_pos = pos_index + cum_pos)

# First, identify genes that appear more than 4 times with Fst > 0.7
high_fst_genes <- fst_plot %>%
  filter(Fst > 0.7) %>%
  count(Gene) %>%
  filter(n > 5) %>%
  pull(Gene)

# Add a column to identify which points should be colored
fst_plot <- fst_plot %>%
  mutate(highlight = ifelse(Fst > 0.7 & Gene %in% high_fst_genes, Comparison, "Other"))

# Filter for plotting (Fst > 0.5)
fst_plot_filtered <- fst_plot %>% filter(Fst > 0.5)

# Then filter the data for labeling
labels_data <- fst_plot %>%
  filter(Fst > 0.7, Gene %in% high_fst_genes)


# Define gene regions from scratch
gene_regions <- data.frame(
  gene_name = c("CRT", "UT", "ACS8", "hyp12", "PF3D7_0705500", 
                "PF3D7_1321100", "PF3D7_1408700", "PF3D7_1448500", "PF3D7_1321300"),
  chr = c("Pf3D7_07_v3", "Pf3D7_07_v3", "Pf3D7_02_v3", "Pf3D7_07_v3", 
          "Pf3D7_07_v3", "Pf3D7_13_v3", "Pf3D7_14_v3", "Pf3D7_14_v3", "Pf3D7_13_v3"),
  start = c(402385, 214949, 628091, 83838, 269635, 
            861001, 324334, 1983780, 876120),
  end = c(406341, 230713, 632681, 86436, 280332, 
          872925, 347119, 1996262, 886756),
  stringsAsFactors = FALSE
)

# Extract chromosome number
gene_regions$chr_num <- as.numeric(gsub("Pf3D7_", "", gsub("_v3", "", gene_regions$chr)))

# Initialize x position columns
gene_regions$x_start <- NA_real_
gene_regions$x_end <- NA_real_


# Loop through each gene
for(i in 1:nrow(gene_regions)) {
  current_chr <- gene_regions$chr_num[i]
  current_start <- gene_regions$start[i]
  current_end <- gene_regions$end[i]
  
  cat(sprintf("\nProcessing %s (chr %d, %d-%d)...\n", 
              gene_regions$gene_name[i], current_chr, current_start, current_end))
  
  # Filter SNPs for this chromosome
  chr_snps <- fst_plot %>% filter(chr_num == current_chr)
  
  cat(sprintf("  Found %d SNPs on chromosome %d\n", nrow(chr_snps), current_chr))
  
  if(nrow(chr_snps) > 0) {
    # Find closest start position
    start_snp <- chr_snps %>%
      mutate(dist = abs(pos - current_start)) %>%
      slice_min(dist, n = 1)
    
    gene_regions$x_start[i] <- start_snp$x_pos[1]
    
    # Find closest end position
    end_snp <- chr_snps %>%
      mutate(dist = abs(pos - current_end)) %>%
      slice_min(dist, n = 1)
    
    gene_regions$x_end[i] <- end_snp$x_pos[1]
    
    cat(sprintf("  x_pos: %.2f to %.2f\n", gene_regions$x_start[i], gene_regions$x_end[i]))
  } else {
    cat("  WARNING: No SNPs found on this chromosome!\n")
  }
}

gene_regions_with_pos <- gene_regions

library(ggforce)  # For drawing half circles

# Define site colors - these match your data exactly
site_colors <- c(
  "Arba Minch" = "#f8766d",       # red
  "Gambella" = "#a3a500",          # orange
  "Gondar" = "#00bf7d",            # cyan
  "Metehera Fentale" = "#00b0f6",  # blue
  "Mizan" = "#e76bf3"              # magenta
)

# Create left and right datasets
fst_left <- fst_plot_filtered %>%
  filter(highlight != "Other") %>%
  mutate(
    x_adj = x_pos - 0.15,
    site_color = Population_1
  )

fst_right <- fst_plot_filtered %>%
  filter(highlight != "Other") %>%
  mutate(
    x_adj = x_pos + 0.15,
    site_color = Population_2
  )

# Filter for semicircles
fst_semicircles <- fst_plot_filtered %>%
  filter(highlight != "Other")

p <- ggplot(fst_plot_filtered, aes(x = x_pos, y = Fst)) +
  
  # --- Gene regions (transparent, colored by gene_name) ---
  geom_rect(
    data = gene_regions_with_pos,
    aes(xmin = x_start - 100,
        xmax = x_end + 100,
        ymin = 0.5,
        ymax = 1,
        fill = gene_name),
    inherit.aes = FALSE,
    color = NA,
    alpha = 0.25        # make much lighter
  ) +
  geom_rect(
    data = gene_regions_with_pos,
    aes(xmin = x_start - 100,
        xmax = x_end + 100,
        ymin = 0.5,
        ymax = 1,
        color = gene_name),
    fill = NA,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_d(
    option = "viridis",
    alpha = 0.25  # Add alpha here in the scale
  )+
  scale_fill_viridis_d(
    guide = "none", 
    option = "viridis",
    alpha = 0.25  # Add alpha here in the scale
  )+
  
  # --- Gene labels (repelled text) ---
  #geom_text_repel(
  #data = gene_regions_with_pos,
  #aes(x = (x_start + x_end) / 2, y = 0.98, label = gene_name),
  #inherit.aes = FALSE,
  #size = 3,
  #fontface = "bold",
  #color = "black",
  #bg.color = "white",
  #bg.r = 0.15,
  #max.overlaps = 100
  #) +
  
  # --- Reset fill/color scales for Fst points ---
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() +
  
  # --- Background points ---
  geom_point(
    data = fst_plot_filtered %>% filter(highlight == "Other"),
    color = "grey70", alpha = 0.3, size = 2
  ) +
  
  # --- Outer ring (Population_2) ---
  geom_point(
    data = fst_semicircles,
    aes(x = x_pos, y = Fst, color = Population_2),
    shape = 21, fill = NA, size = 4, stroke = 1.3
  ) +
  
  # --- Inner fill (Population_1) ---
  geom_point(
    data = fst_semicircles,
    aes(x = x_pos, y = Fst, fill = Population_1),
    shape = 21, color = "transparent", size = 3.2, stroke = 0
  ) +
  
  scale_fill_manual(values = site_colors, name = "Filling: Population 1") +
  scale_color_manual(values = site_colors, name = "Ring: Population 2") +
  
  # --- Axes & layout ---
  scale_x_continuous(
    breaks = chr_lengths$chr_center,
    labels = chr_lengths$chr_num,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(0.5, 1)) +
  coord_cartesian(
    xlim = c(min(fst_plot_filtered$x_pos), max(fst_plot_filtered$x_pos)),
    ylim = c(0.5, 1)
  ) +
  labs(x = "Chromosome", y = "Fst Value") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = chr_lengths$cum_pos[-1],
             linetype = "dotted", color = "gray80", alpha = 0.5)

print(p)
ggsave("manhattan_split_colors.pdf", p, width = 16, height = 10)

# Print color summary
cat("\n📊 Color distribution:\n")
cat("\nLeft side (Population 1):\n")
print(table(fst_left$site_color))
cat("\nRight side (Population 2):\n")
print(table(fst_right$site_color))

# Show example comparisons
cat("\n📋 Example comparisons and their colors:\n")
examples <- fst_left %>%
  select(Population_1, Population_2, site_color) %>%
  distinct() %>%
  arrange(Population_1, Population_2) %>%
  head(10)
print(examples)
################################################################################
# SNP-LEVEL SELECTION SIGNAL VISUALISATION - OLD VS NEW ETHIOPIA
# Groups: old (pre-2024) vs new (2024)
# Two plots per group:
#   1. Two-panel Manhattan (iHS top, XP-EHH bottom), significant SNPs highlighted
#   2. iHS logp vs XP-EHH logp scatter for co-detected SNPs, labelled by gene
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(cowplot)

BASE_DIR <- "/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection"
setwd(BASE_DIR)

IN_FILE        <- "old_vs_new/old_vs_new_integrated_selection_analysis_SNPs.tsv"
OUT_PREFIX     <- "old_vs_new/snp_plots_old_vs_new"
ANALYSIS_LABEL <- "Old vs New Ethiopia"
IHS_THRESH     <- -log10(2 * pnorm(-abs(2.0)))   # |iHS| >= 2.0 on logp scale
XPEHH_THRESH   <- -log10(2 * pnorm(-abs(2.5)))   # |XP-EHH| >= 2.5 on logp scale
TOP_N_GENES    <- 20

GROUP_COLORS <- c(
  "old" = "#4E79A7",
  "new" = "#F28E2B"
)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading:", IN_FILE, "\n")
snps <- read.delim(IN_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    Chr        = as.integer(Chr),
    Position   = as.numeric(Position),
    iHS_logp   = as.numeric(iHS_logp),
    XPEHH_logp = as.numeric(XPEHH_logp),
    n_datasets = as.integer(n_datasets)
  ) %>%
  filter(!is.na(Chr), !is.na(Position))

cat("  Total SNP rows:", nrow(snps), "\n")
cat("  Detected by both tools:", sum(snps$n_datasets == 2, na.rm = TRUE), "\n\n")

# ==============================================================================
# 2. CHROMOSOME OFFSETS
# ==============================================================================

chr_lengths <- data.frame(
  Chr    = 1:14,
  length = c(643292, 947102, 1067971, 1200490, 1343557, 1418242, 1445207,
             1472805, 1541735, 1687656, 2038340, 2271494, 2925236, 3291936)
)

chr_offsets <- chr_lengths %>%
  arrange(Chr) %>%
  mutate(offset = cumsum(lag(length, default = 0)) + (Chr - 1) * 2e5)

chr_midpoints <- chr_offsets %>% mutate(mid = offset + length / 2)

snps <- snps %>%
  left_join(chr_offsets %>% select(Chr, offset), by = "Chr") %>%
  mutate(
    pos_cum = Position + offset,
    detection = case_when(
      n_datasets == 2                         ~ "Both tools",
      !is.na(iHS_logp)  & is.na(XPEHH_logp)  ~ "iHS only",
      is.na(iHS_logp)   & !is.na(XPEHH_logp) ~ "XP-EHH only",
      TRUE                                    ~ "Both tools"
    )
  )

det_colors <- c("Both tools" = "#D62728", "iHS only" = "#AEC7E8", "XP-EHH only" = "#FFBB78")
det_sizes  <- c("Both tools" = 1.2,       "iHS only" = 0.6,       "XP-EHH only" = 0.6)
det_alpha  <- c("Both tools" = 1.0,       "iHS only" = 0.5,       "XP-EHH only" = 0.5)

x_breaks <- chr_midpoints$mid
x_labs   <- as.character(chr_midpoints$Chr)
x_lim    <- range(snps$pos_cum, na.rm = TRUE)

# ==============================================================================
# 3 & 4. PLOTS PER GROUP
# ==============================================================================

GROUP_COL <- "iHS_Region"

groups_in_data <- intersect(names(GROUP_COLORS), unique(snps[[GROUP_COL]]))
cat("Groups found in data:", paste(groups_in_data, collapse = ", "), "\n\n")

for (grp in groups_in_data) {

  cat("==============================\n")
  cat("Processing:", grp, "\n")
  cat("==============================\n")

  snps_c       <- filter(snps, .data[[GROUP_COL]] == grp)
  out_prefix_c <- paste0(OUT_PREFIX, "_", grp)

  if (nrow(snps_c) == 0) {
    cat("  No SNPs for", grp, "— skipping.\n\n")
    next
  }

  # ----------------------------------------------------------------------------
  # PLOT 1 — TWO-PANEL MANHATTAN
  # ----------------------------------------------------------------------------

  cat("  Generating Manhattan plot...\n")

  ihs_pts <- snps_c %>%
    filter(!is.na(iHS_logp)) %>%
    group_by(Chr, Position, pos_cum, iHS_Region, iHS_Gene, detection) %>%
    summarise(iHS_logp = max(iHS_logp, na.rm = TRUE), .groups = "drop")

  xp_pts <- snps_c %>%
    filter(!is.na(XPEHH_logp)) %>%
    group_by(Chr, Position, pos_cum, XPEHH_Comparison, XPEHH_Gene, detection) %>%
    summarise(XPEHH_logp = max(XPEHH_logp, na.rm = TRUE), .groups = "drop")

  ihs_fg <- ihs_pts %>% filter(iHS_logp   >= IHS_THRESH,   detection == "Both tools")
  xp_fg  <- xp_pts  %>% filter(XPEHH_logp >= XPEHH_THRESH, detection == "Both tools")

  ihs_labels <- ihs_fg %>%
    filter(!is.na(iHS_Gene) & iHS_Gene != "" & iHS_Gene != "NA") %>%
    group_by(iHS_Gene) %>%
    filter(n() > 1) %>%
    slice_max(iHS_logp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    slice_max(iHS_logp, n = TOP_N_GENES, with_ties = FALSE)

  xp_labels <- xp_fg %>%
    filter(!is.na(XPEHH_Gene) & XPEHH_Gene != "" & XPEHH_Gene != "NA") %>%
    group_by(XPEHH_Gene) %>%
    filter(n() > 1) %>%
    slice_max(XPEHH_logp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    slice_max(XPEHH_logp, n = TOP_N_GENES, with_ties = FALSE)

  p_ihs <- ggplot(mapping = aes(x = pos_cum, y = iHS_logp)) +
    geom_point(data = ihs_fg, colour = "#D62728", size = 2.0, stroke = 0) +
    geom_text_repel(data = ihs_labels, aes(label = iHS_Gene), size = 2.5, colour = "black",
                    max.overlaps = Inf, box.padding = 0.3, segment.colour = "grey50",
                    segment.size = 0.3, min.segment.length = 0, seed = 42) +
    geom_hline(yintercept = IHS_THRESH, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(breaks = x_breaks, labels = x_labs, limits = x_lim,
                       expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = expression(-log[10](p)~"(iHS)")) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

  p_xp <- ggplot(mapping = aes(x = pos_cum, y = XPEHH_logp)) +
    geom_point(data = xp_fg, colour = "#D62728", size = 2.0, stroke = 0) +
    geom_text_repel(data = xp_labels, aes(label = XPEHH_Gene), size = 2.5, colour = "black",
                    max.overlaps = Inf, box.padding = 0.3, segment.colour = "grey50",
                    segment.size = 0.3, min.segment.length = 0, seed = 42,
                    direction = "y") +
    geom_hline(yintercept = XPEHH_THRESH, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(breaks = x_breaks, labels = x_labs, limits = x_lim,
                       expand = expansion(mult = 0.01)) +
    scale_y_reverse() +
    labs(x = "Chromosome", y = expression(-log[10](p)~"(XP-EHH)")) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

  p_manhattan <- plot_grid(p_ihs, p_xp, ncol = 1, rel_heights = c(1, 1), align = "v")
  ggsave(paste0(out_prefix_c, "_manhattan.pdf"), plot = p_manhattan, width = 10, height = 6, device = "pdf")
  message("  Saved: ", out_prefix_c, "_manhattan.pdf")

  # ----------------------------------------------------------------------------
  # PLOT 2 — iHS vs XP-EHH SCATTER (co-detected SNPs)
  # ----------------------------------------------------------------------------

  cat("  Generating scatter plot...\n")

  scatter_df <- snps_c %>%
    filter(n_datasets == 2, !is.na(iHS_logp), !is.na(XPEHH_logp)) %>%
    group_by(Chr, Position, iHS_Region) %>%
    slice_max(XPEHH_logp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(gene_label = case_when(
      !is.na(iHS_Gene)   & iHS_Gene   != "" & iHS_Gene   != "NA" ~ iHS_Gene,
      !is.na(XPEHH_Gene) & XPEHH_Gene != "" & XPEHH_Gene != "NA" ~ XPEHH_Gene,
      TRUE ~ NA_character_
    ))

  cat("  Co-detected SNPs:", nrow(scatter_df), "\n")

  if (nrow(scatter_df) == 0) {
    cat("  No co-detected SNPs — skipping scatter.\n\n")
    next
  }

  sig_labels <- scatter_df %>%
    filter(!is.na(gene_label),
           iHS_logp >= IHS_THRESH, XPEHH_logp >= XPEHH_THRESH) %>%
    group_by(gene_label) %>%
    filter(n() > 1) %>%
    slice_max(iHS_logp + XPEHH_logp, n = 1, with_ties = FALSE) %>%
    ungroup()

  p_scatter <- ggplot(scatter_df, aes(x = iHS_logp, y = XPEHH_logp)) +
    geom_point(alpha = 0.6, size = 1.5, colour = GROUP_COLORS[grp]) +
    geom_vline(xintercept = IHS_THRESH,   linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = XPEHH_THRESH, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_text_repel(data = sig_labels, aes(label = gene_label), size = 2.5, colour = "black",
                    max.overlaps = Inf, box.padding = 0.3, segment.colour = "grey60",
                    segment.size = 0.3, min.segment.length = 0, seed = 42) +
    labs(x = expression(-log[10](p)~"(iHS)"), y = expression(-log[10](p)~"(XP-EHH)"),
         title    = paste0("Co-detected SNPs — ", grp),
         subtitle = paste0(nrow(scatter_df), " SNPs detected by both iHS and XP-EHH")) +
    theme_classic(base_size = 11) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

  ggsave(paste0(out_prefix_c, "_scatter.pdf"), plot = p_scatter, width = 8, height = 7, device = "pdf")
  message("  Saved: ", out_prefix_c, "_scatter.pdf")

  cat("\n")
}

cat("Done.\n")

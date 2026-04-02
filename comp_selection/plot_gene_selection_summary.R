################################################################################
# GENE SELECTION SIGNAL SUMMARY FIGURE
# For genes of interest, shows which populations have selection signals
# across all 6 comparison analyses (SNP-based and window-based)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

setwd("/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection")

# ==============================================================================
# DEFINE GENES OF INTEREST
# ==============================================================================

# Gene names as they appear in the data files
# SNP files use short names; window files use PF3D7_ID(name) format
genes_of_interest <- c("CRT", "DHPS", "PPPK-DHPS", "DHFR", "MDR1", "AAT1",
                        "Kelch13", "PX1", "MCA2",
                        "KIC4", "KIC5", "KIC7", "KIC8", "KIC9", "KIC10")

# Mapping to display names
gene_display <- c(
  "CRT"       = "crt",
  "PPPK-DHPS" = "dhps",
  "DHPS"      = "dhps",
  "DHFR"      = "dhfr",
  "MDR1"      = "mdr1",
  "AAT1"      = "aat1",
  "Kelch13"   = "kelch13",
  "PX1"       = "px1",
  "MCA2"      = "mca2",
  "KIC4"      = "kic4",
  "KIC5"      = "kic5",
  "KIC7"      = "kic7",
  "KIC8"      = "kic8",
  "KIC9"      = "kic9",
  "KIC10"     = "kic10"
)

# PF3D7 IDs for window file searching
gene_ids <- c(
  "PF3D7_0709000" = "crt",
  "PF3D7_0810800" = "dhps",
  "PF3D7_0417200" = "dhfr",
  "PF3D7_0523000" = "mdr1",
  "PF3D7_0629500" = "aat1",
  "PF3D7_1343700" = "kelch13",
  "PF3D7_0720700" = "px1",
  "PF3D7_1438400" = "mca2",
  "PF3D7_1246300" = "kic4",
  "PF3D7_1138700" = "kic5",
  "PF3D7_0813000" = "kic7",
  "PF3D7_1014900" = "kic8",
  "PF3D7_1442400" = "kic9",
  "PF3D7_1122500" = "kic10"
)

# ==============================================================================
# DEFINE COMPARISONS
# ==============================================================================

comparisons <- c("countries", "regions", "hrp", "k13", "north_vs_other", "old_vs_new")

# ==============================================================================
# FUNCTION: Extract gene signals from SNP file
# ==============================================================================

extract_snp_signals <- function(comp_name) {
  f <- paste0(comp_name, "/", comp_name, "_integrated_selection_analysis_SNPs.tsv")
  if (!file.exists(f)) return(NULL)

  dat <- read.delim(f, stringsAsFactors = FALSE)
  results <- list()

  for (i in seq_len(nrow(dat))) {
    row <- dat[i, ]

    # Check iHS gene column
    if (!is.na(row$iHS_Gene) && row$iHS_Gene != "") {
      gene_short <- row$iHS_Gene
      if (gene_short %in% names(gene_display)) {
        results[[length(results) + 1]] <- data.frame(
          gene       = gene_display[gene_short],
          population = row$iHS_Region,
          tool       = "iHS",
          analysis   = "SNP",
          comparison = comp_name,
          stringsAsFactors = FALSE
        )
      }
    }

    # Check XP-EHH gene column
    if (!is.na(row$XPEHH_Gene) && row$XPEHH_Gene != "") {
      gene_short <- row$XPEHH_Gene
      if (gene_short %in% names(gene_display)) {
        # For XP-EHH, signal is in pop1 (the population showing selection)
        results[[length(results) + 1]] <- data.frame(
          gene       = gene_display[gene_short],
          population = row$XPEHH_pop1,
          tool       = "XP-EHH",
          analysis   = "SNP",
          comparison = comp_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# ==============================================================================
# FUNCTION: Check if a gene is present in a semicolon-separated gene string
# ==============================================================================

find_gene_in_string <- function(gene_str, gene_ids, gene_display) {
  if (is.na(gene_str) || gene_str == "") return(NULL)

  found <- c()

  # Check for PF3D7 IDs
  for (gid in names(gene_ids)) {
    if (grepl(gid, gene_str, fixed = TRUE)) {
      found <- c(found, gene_ids[gid])
    }
  }

  # Check for short names (window files sometimes use short names like "CRT", "MDR1")
  for (gname in names(gene_display)) {
    # Use word boundary matching to avoid partial matches
    pattern <- paste0("\\b", gname, "\\b")
    if (grepl(pattern, gene_str)) {
      found <- c(found, gene_display[gname])
    }
  }

  unique(found)
}

# ==============================================================================
# FUNCTION: Extract gene signals from window file
# ==============================================================================

extract_window_signals <- function(comp_name) {
  f <- paste0(comp_name, "/", comp_name, "_integrated_selection_analysis_windows.tsv")
  if (!file.exists(f)) return(NULL)

  dat <- read.delim(f, stringsAsFactors = FALSE)
  results <- list()

  # Gene columns and their associated tool and population columns
  # Column names are consistent but order may vary
  gene_tool_map <- list(
    list(gene_col = "XPEHH_genes", tool = "XP-EHH", pop_col = "XPEHH_pop1"),
    list(gene_col = "IBD_genes",   tool = "IBD",     pop_col = "IBD_Region"),
    list(gene_col = "iHS_genes",   tool = "iHS",     pop_col = "iHS_Region"),
    list(gene_col = "TD_genes",    tool = "Tajima's D", pop_col = "TD_Region")
  )

  for (i in seq_len(nrow(dat))) {
    row <- dat[i, ]

    for (mapping in gene_tool_map) {
      gcol <- mapping$gene_col
      if (!(gcol %in% names(row))) next

      gene_str <- as.character(row[[gcol]])
      pop <- as.character(row[[mapping$pop_col]])

      if (is.na(gene_str) || gene_str == "" || is.na(pop) || pop == "") next

      matched_genes <- find_gene_in_string(gene_str, gene_ids, gene_display)

      for (g in matched_genes) {
        results[[length(results) + 1]] <- data.frame(
          gene       = g,
          population = pop,
          tool       = mapping$tool,
          analysis   = "window",
          comparison = comp_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# ==============================================================================
# COLLECT ALL SIGNALS
# ==============================================================================

cat("Extracting selection signals for genes of interest...\n")

all_signals <- list()

for (comp in comparisons) {
  cat("  Processing:", comp, "\n")

  snp_signals <- extract_snp_signals(comp)
  if (!is.null(snp_signals)) {
    cat("    SNP signals found:", nrow(snp_signals), "\n")
    all_signals[[paste0(comp, "_snp")]] <- snp_signals
  }

  win_signals <- extract_window_signals(comp)
  if (!is.null(win_signals)) {
    cat("    Window signals found:", nrow(win_signals), "\n")
    all_signals[[paste0(comp, "_win")]] <- win_signals
  }
}

signals_df <- bind_rows(all_signals)

# Deduplicate: unique combination of gene, population, tool, analysis, comparison
signals_df <- signals_df %>%
  distinct(gene, population, tool, analysis, comparison) %>%
  filter(!is.na(population) & population != "")

cat("\nTotal unique signals:", nrow(signals_df), "\n")
cat("Genes found:", paste(sort(unique(signals_df$gene)), collapse = ", "), "\n")
cat("Populations found:", paste(sort(unique(signals_df$population)), collapse = ", "), "\n")

# ==============================================================================
# ORDER GENES AND POPULATIONS
# ==============================================================================

gene_order <- c("crt", "dhps", "dhfr", "mdr1", "aat1", "kelch13", "px1", "mca2",
                 "kic4", "kic5", "kic7", "kic8", "kic9", "kic10")
signals_df$gene <- factor(signals_df$gene, levels = gene_order)

# Rename populations for display
pop_rename <- c(
  "new"   = "Newly sequenced Ethiopian isolates",
  "old"   = "Older Ethiopian isolates",
  "north" = "Northern Ethiopia",
  "rest"  = "Other parts of Ethiopia"
)
signals_df$population <- ifelse(
  signals_df$population %in% names(pop_rename),
  pop_rename[signals_df$population],
  signals_df$population
)

# Order populations: countries first, then regions, then other groups
pop_order <- c(
  # Countries
  "Ethiopia", "Eritrea", "Sudan", "Kenya", "Uganda", "Tanzania",
  "Mali", "Nigeria", "Ghana",
  # Regions
  "Horn_of_Africa", "East_Africa", "West_Africa",
  # HRP groups
  "hrp_wild", "hrp3_only", "hrp_double_deleted",
  # K13 groups
  "r622i_wt", "r622i_mut",
  # North vs other
  "Northern Ethiopia", "Other parts of Ethiopia",
  # Old vs new
  "Newly sequenced Ethiopian isolates", "Older Ethiopian isolates"
)

# Only include populations that are actually in the data
pop_present <- intersect(pop_order, unique(signals_df$population))
signals_df$population <- factor(signals_df$population, levels = rev(pop_present))

# ==============================================================================
# PLOT
# ==============================================================================

# Tool colours
tool_colors <- c(
  "iHS"        = "#E41A1C",
  "XP-EHH"    = "#377EB8",
  "IBD"        = "#4DAF4A",
  "Tajima's D" = "#984EA3"
)

# Drop genes with no signals in this plot
signals_df_plot1 <- signals_df %>% filter(!is.na(gene)) %>% mutate(gene = droplevels(gene))

p <- ggplot(signals_df_plot1, aes(x = gene, y = population)) +
  geom_point(aes(shape = analysis, colour = tool),
             size = 3.5, stroke = 0.8,
             position = position_jitter(width = 0.15, height = 0.15, seed = 42)) +
  scale_shape_manual(values = c("SNP" = 17, "window" = 16),
                     labels = c("SNP" = "SNP-based", "window" = "Window-based")) +
  scale_colour_manual(values = tool_colors) +
  labs(x = "Gene", y = "Population",
       shape = "Analysis type", colour = "Tool",
       title = "Selection signals in genes of interest across populations") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold")
  )

# Save
ggsave("gene_selection_summary.pdf", p,
       width = 10, height = 8, units = "in")

cat("\nFigure saved: gene_selection_summary.pdf\n")

# ==============================================================================
# PLOT 2: Colour = number of tools, shape = analysis type
# Split into 3 figures: geographical, R622I, HRP
# ==============================================================================

# Count number of distinct tools per gene x population x analysis
signals_count_all <- signals_df %>%
  group_by(gene, population, analysis, comparison) %>%
  summarise(n_tools = n_distinct(tool), .groups = "drop") %>%
  mutate(n_tools = factor(n_tools))

# Filtered version (2+ tools) for geographical plot
signals_count <- signals_count_all %>% filter(as.numeric(as.character(n_tools)) >= 2)

# Offset SNP left, window right (instead of jitter)
signals_count <- signals_count %>%
  mutate(
    x_num    = as.numeric(gene),
    x_offset = ifelse(analysis == "SNP", x_num - 0.15, x_num + 0.15)
  )

# Combine analysis + n_tools for colour
signals_count <- signals_count %>%
  mutate(analysis_ntools = paste0(analysis, ": ", n_tools, " tools"))

# Build colour palette: blues for SNP, green/red/purple for window
# Map colours by n_tools value (not position) so colours are consistent across plots
snp_color_map <- c("1" = "#4292C6", "2" = "#9ECAE1", "3" = "#08519C", "4" = "#08306B")
win_color_map <- c("1" = "#7F2704", "2" = "#E41A1C", "3" = "#984EA3", "4" = "#4DAF4A")

snp_levels <- sort(unique(signals_count_all$n_tools[signals_count_all$analysis == "SNP"]))
win_levels <- sort(unique(signals_count_all$n_tools[signals_count_all$analysis == "window"]))

snp_palette <- if (length(snp_levels) > 0) {
  setNames(snp_color_map[as.character(snp_levels)],
           paste0("SNP: ", snp_levels, " tools"))
} else { c() }

win_palette <- if (length(win_levels) > 0) {
  setNames(win_color_map[as.character(win_levels)],
           paste0("window: ", win_levels, " tools"))
} else { c() }
combined_palette <- c(snp_palette, win_palette)

# Shared theme
theme_selection <- theme_bw(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    plot.title         = element_text(size = 13, face = "bold")
  )

# ------------------------------------------------------------------------------
# Helper: build a plot for a subset, with horizontal separators between comparisons
# ------------------------------------------------------------------------------

make_selection_plot <- function(dat, title_text, pop_order_vec, comparison_groups,
                                subtitle_text = NULL) {
  # comparison_groups: list of list(label, pops) in display order (top to bottom)

  # Build ordered population vector (top to bottom on y-axis = last to first in factor)
  pop_vec <- unlist(lapply(comparison_groups, function(g) g$pops))
  dat <- dat %>% filter(as.character(population) %in% pop_vec)
  if (nrow(dat) == 0) return(NULL)

  dat$population <- factor(dat$population, levels = rev(pop_vec))

  # Drop genes with no signals in this subset
  dat <- dat %>% mutate(gene = droplevels(gene))
  local_gene_order <- levels(dat$gene)

  # Aggregate: take max n_tools per gene x population x analysis (across comparisons)
  dat <- dat %>%
    group_by(gene, population, analysis) %>%
    summarise(n_tools = max(as.numeric(as.character(n_tools))), .groups = "drop") %>%
    mutate(
      n_tools        = factor(n_tools),
      analysis_ntools = paste0(analysis, ": ", n_tools, " tools"),
      x_num          = as.numeric(gene),
      x_offset       = ifelse(analysis == "SNP", x_num - 0.18, x_num + 0.18)
    )

  # Horizontal separator lines between comparison groups
  hline_positions <- c()
  cumulative <- 0
  for (i in seq_along(comparison_groups)) {
    cumulative <- cumulative + length(comparison_groups[[i]]$pops)
    if (i < length(comparison_groups)) {
      hline_positions <- c(hline_positions,
                           length(pop_vec) - cumulative + 0.5)
    }
  }

  # Comparison group labels for left-hand side annotation
  label_df <- data.frame(
    label = character(), mid_y = numeric(), stringsAsFactors = FALSE
  )
  cumul <- 0
  for (g in comparison_groups) {
    n <- length(g$pops)
    mid_y <- length(pop_vec) - cumul - (n - 1) / 2
    label_df <- rbind(label_df, data.frame(label = g$label, mid_y = mid_y))
    cumul <- cumul + n
  }

  p <- ggplot(dat, aes(x = x_offset, y = population)) +
    geom_point(aes(shape = analysis, colour = analysis_ntools),
               size = 3.5, stroke = 0.8) +
    scale_shape_manual(values = c("SNP" = 17, "window" = 16),
                       labels = c("SNP" = "SNP-based", "window" = "Window-based")) +
    scale_colour_manual(values = combined_palette, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(local_gene_order),
      labels = local_gene_order,
      limits = c(0.5, length(local_gene_order) + 0.5)
    ) +
    labs(x = "Gene", y = NULL,
         shape = "Analysis type", colour = "Number of tools",
         title = title_text,
         subtitle = subtitle_text) +
    theme_selection

  # Add thick black separator lines between comparison groups
  for (ypos in hline_positions) {
    p <- p + geom_hline(yintercept = ypos, linewidth = 1, colour = "black")
  }

  # Add comparison labels on the left
  for (i in seq_len(nrow(label_df))) {
    p <- p + annotate("text", x = 0.3, y = label_df$mid_y[i],
                       label = label_df$label[i], angle = 90,
                       size = 3, fontface = "bold")
  }

  p <- p + coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5, 5, 5, 40))

  p
}

# ==============================================================================
# FIGURE 2A: Geographical comparisons (countries, regions, north_vs_other, old_vs_new)
# ==============================================================================

geo_groups <- list(
  list(label = "Countries",
       pops = c("Ethiopia", "Eritrea", "Sudan", "Kenya", "Uganda", "Tanzania",
                "Mali", "Nigeria", "Ghana")),
  list(label = "Regions",
       pops = c("Horn_of_Africa", "East_Africa", "West_Africa")),
  list(label = "North vs Other",
       pops = c("Northern Ethiopia", "Other parts of Ethiopia")),
  list(label = "Old vs New",
       pops = c("Newly sequenced Ethiopian isolates", "Older Ethiopian isolates"))
)

geo_data <- signals_count %>%
  filter(comparison %in% c("countries", "regions", "north_vs_other", "old_vs_new"))

p2a <- make_selection_plot(
  geo_data,
  "Selection signals — Geographical comparisons",
  NULL, geo_groups,
  subtitle_text = "Only signals detected by 2+ tools shown"
)

if (!is.null(p2a)) {
  ggsave("gene_selection_summary_geographical.pdf", p2a,
         width = 10, height = 8, units = "in")
  cat("Figure saved: gene_selection_summary_geographical.pdf\n")
}

# ==============================================================================
# FIGURE 2B: R622I (K13) comparisons
# ==============================================================================

k13_groups <- list(
  list(label = "R622I", pops = c("r622i_wt", "r622i_mut"))
)

k13_data <- signals_count_all %>%
  filter(comparison == "k13") %>%
  filter(!(analysis == "window" & as.numeric(as.character(n_tools)) == 1))

p2b <- make_selection_plot(
  k13_data,
  "Selection signals — R622I comparisons",
  NULL, k13_groups
)

if (!is.null(p2b)) {
  ggsave("gene_selection_summary_k13.pdf", p2b,
         width = 5, height = 4, units = "in")
  cat("Figure saved: gene_selection_summary_k13.pdf\n")
}

# ==============================================================================
# FIGURE 2C: HRP comparisons
# ==============================================================================

hrp_groups <- list(
  list(label = "HRP", pops = c("hrp_wild", "hrp3_only", "hrp_double_deleted"))
)

hrp_data <- signals_count_all %>%
  filter(comparison == "hrp") %>%
  filter(!(analysis == "window" & as.numeric(as.character(n_tools)) == 1))

p2c <- make_selection_plot(
  hrp_data,
  "Selection signals — HRP deletions",
  NULL, hrp_groups
)

if (!is.null(p2c)) {
  ggsave("gene_selection_summary_hrp.pdf", p2c,
         width = 6, height = 4, units = "in")
  cat("Figure saved: gene_selection_summary_hrp.pdf\n")
}

# ==============================================================================
# FIGURE 2D: KIC genes across all comparisons and populations
# ==============================================================================

kic_data <- signals_count_all %>%
  filter(grepl("^kic", gene))

# Build comparison groups from all populations present in KIC data
kic_pop_list <- unique(as.character(kic_data$population))

kic_groups <- list(
  list(label = "Countries",
       pops = intersect(c("Ethiopia", "Eritrea", "Sudan", "Kenya", "Uganda", "Tanzania",
                           "Mali", "Nigeria", "Ghana"), kic_pop_list)),
  list(label = "Regions",
       pops = intersect(c("Horn_of_Africa", "East_Africa", "West_Africa"), kic_pop_list)),
  list(label = "HRP",
       pops = intersect(c("hrp_wild", "hrp3_only", "hrp_double_deleted"), kic_pop_list)),
  list(label = "R622I",
       pops = intersect(c("r622i_wt", "r622i_mut"), kic_pop_list)),
  list(label = "North vs Other",
       pops = intersect(c("Northern Ethiopia", "Other parts of Ethiopia"), kic_pop_list)),
  list(label = "Old vs New",
       pops = intersect(c("Newly sequenced Ethiopian isolates", "Older Ethiopian isolates"), kic_pop_list))
)

# Remove empty groups
kic_groups <- kic_groups[sapply(kic_groups, function(g) length(g$pops) > 0)]

p2d <- make_selection_plot(
  kic_data,
  "Selection signals — KIC genes across all comparisons",
  NULL, kic_groups
)

if (!is.null(p2d)) {
  ggsave("gene_selection_summary_kic.pdf", p2d,
         width = 8, height = 8, units = "in")
  cat("Figure saved: gene_selection_summary_kic.pdf\n")
}

# ==============================================================================
# FIGURE 3: KIC genes — IBD fraction vs genome-wide distribution
# Violin of genome-wide IBD per population, with KIC gene values overlaid
# ==============================================================================

cat("\nBuilding KIC IBD violin plot...\n")

ibd_all <- read.delim("../hmmIBD/hmmibd_combined_win10kb_annotated.tsv",
                       stringsAsFactors = FALSE)
ibd_all <- ibd_all %>% filter(!is.na(fraction))

# Recode category names to match the display names used elsewhere
ibd_recode <- c(
  "Eritrea" = "Eritrea", "Ethiopia" = "Ethiopia", "Sudan" = "Sudan",
  "Kenya" = "Kenya", "Uganda" = "Uganda", "Tanzania" = "Tanzania",
  "Mali" = "Mali", "Nigeria" = "Nigeria", "Ghana" = "Ghana",
  "Horn_of_Africa" = "Horn_of_Africa", "East_Africa" = "East_Africa",
  "West_Africa" = "West_Africa",
  "South_Asia" = "South_Asia", "Southeast_Asia" = "Southeast_Asia",
  "wild" = "hrp_wild", "hrp3_only" = "hrp3_only",
  "mutant" = "hrp_double_deleted", "hrp_double_deleted" = "hrp_double_deleted",
  "wt" = "r622i_wt", "mut" = "r622i_mut",
  "north" = "Northern Ethiopia", "rest" = "Other parts of Ethiopia",
  "new" = "Newly sequenced Ethiopian isolates",
  "old" = "Older Ethiopian isolates"
)

ibd_all$category <- ifelse(
  ibd_all$category %in% names(ibd_recode),
  ibd_recode[ibd_all$category],
  ibd_all$category
)

# KIC gene IDs
kic_ids <- gene_ids[grepl("^kic", gene_ids)]

# Extract windows containing KIC genes
kic_ibd <- ibd_all %>%
  filter(!is.na(gene_id)) %>%
  rowwise() %>%
  mutate(kic_gene = {
    matched <- names(kic_ids)[sapply(names(kic_ids), function(gid) grepl(gid, gene_id, fixed = TRUE))]
    if (length(matched) > 0) kic_ids[matched[1]] else NA_character_
  }) %>%
  ungroup() %>%
  filter(!is.na(kic_gene))

# Population order matching the rest of the script
kic_ibd_pops <- unique(kic_ibd$category)
all_pop_order <- c(
  "Ethiopia", "Eritrea", "Sudan", "Kenya", "Uganda", "Tanzania",
  "Mali", "Nigeria", "Ghana",
  "South_Asia", "Southeast_Asia",
  "Horn_of_Africa", "East_Africa", "West_Africa",
  "hrp_wild", "hrp3_only", "hrp_double_deleted",
  "r622i_wt", "r622i_mut",
  "Northern Ethiopia", "Other parts of Ethiopia",
  "Newly sequenced Ethiopian isolates", "Older Ethiopian isolates"
)
pops_present <- intersect(all_pop_order, unique(ibd_all$category))

# Filter genome-wide data to populations that have KIC hits
ibd_bg <- ibd_all %>% filter(category %in% pops_present)
ibd_bg$category <- factor(ibd_bg$category, levels = pops_present)
kic_ibd$category <- factor(kic_ibd$category, levels = pops_present)

# Comparison group labels for faceting
comp_lookup <- c(
  "Ethiopia" = "Countries", "Eritrea" = "Countries", "Sudan" = "Countries",
  "Kenya" = "Countries", "Uganda" = "Countries", "Tanzania" = "Countries",
  "Mali" = "Countries", "Nigeria" = "Countries", "Ghana" = "Countries",
  "South_Asia" = "Countries", "Southeast_Asia" = "Countries",
  "Horn_of_Africa" = "Regions", "East_Africa" = "Regions", "West_Africa" = "Regions",
  "hrp_wild" = "HRP", "hrp3_only" = "HRP", "hrp_double_deleted" = "HRP",
  "r622i_wt" = "R622I", "r622i_mut" = "R622I",
  "Northern Ethiopia" = "North vs Other", "Other parts of Ethiopia" = "North vs Other",
  "Newly sequenced Ethiopian isolates" = "Old vs New",
  "Older Ethiopian isolates" = "Old vs New"
)

ibd_bg$comp_group <- comp_lookup[as.character(ibd_bg$category)]
kic_ibd$comp_group <- comp_lookup[as.character(kic_ibd$category)]

group_order <- c("Countries", "Regions", "HRP", "R622I", "North vs Other", "Old vs New")
ibd_bg$comp_group <- factor(ibd_bg$comp_group, levels = group_order)
kic_ibd$comp_group <- factor(kic_ibd$comp_group, levels = group_order)

# KIC gene colours
kic_gene_order <- c("kic4", "kic5", "kic7", "kic8", "kic9", "kic10")
kic_colors <- c(
  "kic4"  = "#E41A1C", "kic5"  = "#377EB8", "kic7"  = "#4DAF4A",
  "kic8"  = "#984EA3", "kic9"  = "#FF7F00", "kic10" = "#A65628"
)
kic_ibd$kic_gene <- factor(kic_ibd$kic_gene, levels = kic_gene_order)

p3 <- ggplot() +
  geom_violin(data = ibd_bg,
              aes(x = category, y = fraction),
              fill = "grey85", colour = "grey60", scale = "width",
              width = 0.7, linewidth = 0.3) +
  geom_point(data = kic_ibd,
             aes(x = category, y = fraction, colour = kic_gene),
             size = 2.5, alpha = 0.85,
             position = position_jitter(width = 0.15, seed = 42)) +
  scale_colour_manual(values = kic_colors, name = "KIC gene") +
  facet_grid(. ~ comp_group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "IBD fraction",
       title = "KIC gene IBD signals vs genome-wide distribution") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 13, face = "bold")
  )

ggsave("gene_selection_summary_kic_ibd_violin.pdf", p3,
       width = 14, height = 5, units = "in")
cat("Figure saved: gene_selection_summary_kic_ibd_violin.pdf\n")

cat("\nAll figures saved.\n")

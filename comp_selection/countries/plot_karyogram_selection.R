################################################################################
# KARYOGRAM: INTEGRATED SELECTION SIGNALS
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggforce)   # geom_shape for rounded pill chromosomes
library(ggrepel)   # geom_text_repel for non-overlapping labels

setwd("/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection")
SUBDIR <- "countries"

# XP-EHH pop1 is always the selected population after calculate_rehh_metrics.R
# filters to XPEHH > 0 before detecting candidate regions.
CHR_SPACING  <- 1.6     # x-distance between chromosomes (>1 = wider apart)
CHR_W        <- 0.28    # half-width of chromosome arms
CEN_W        <- 0.06    # half-width at centromere constriction
CEN_ZONE     <- 40000   # bp, width of constriction zone each side of centromere
SHAPE_RADIUS <- unit(3, "mm")   # corner rounding for pill shape

# ==============================================================================
# 1. CHROMOSOME REFERENCE DATA
# ==============================================================================

chr_lengths <- data.frame(
  Chr = 1:14,
  length = c(643292, 947102, 1067971, 1200490, 1343557, 1418242, 1445207,
             1472805, 1541735, 1687656, 2038340, 2271494, 2925236, 3291936)
) %>%
  mutate(x_pos = Chr * CHR_SPACING)   # actual x-axis position

# Centromere positions from PfCENH3 ChIP-seq (doi: 10.1186/gb-2010-11-10-r108)
centromeres <- data.frame(
  Chr     = 1:14,
  cen_mid = c(
    (458021 + 462661) / 2,
    (446771 + 450941) / 2,
    (593691 + 597951) / 2,
    (648001 + 652321) / 2,
    (454541 + 458791) / 2,
    (477751 + 482011) / 2,
    (863741 + 868251) / 2,
    (299011 + 303631) / 2,
    (1241071 + 1245441) / 2,
    (935681 + 937821) / 2,
    (830781 + 835431) / 2,
    (1281511 + 1285931) / 2,
    (1167271 + 1171921) / 2,
    (1070851 + 1075311) / 2
  )
) %>%
  left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")

# ==============================================================================
# 2. CHROMOSOME SHAPE (10-vertex pill with centromere constriction)
# ==============================================================================

make_chr_vertices <- function(chr_x, chr_len, cen_y,
                               full_w = CHR_W, narrow_w = CEN_W,
                               cen_zone = CEN_ZONE) {
  data.frame(
    group = as.character(chr_x),
    x = c(chr_x - full_w,   chr_x + full_w,
          chr_x + full_w,   chr_x + narrow_w,
          chr_x + full_w,   chr_x + full_w,
          chr_x - full_w,   chr_x - full_w,
          chr_x - narrow_w, chr_x - full_w),
    y = c(0,                0,
          cen_y - cen_zone, cen_y,
          cen_y + cen_zone, chr_len,
          chr_len,          cen_y + cen_zone,
          cen_y,            cen_y - cen_zone)
  )
}

chr_polys <- bind_rows(
  mapply(make_chr_vertices,
         chr_x   = chr_lengths$x_pos,
         chr_len = chr_lengths$length,
         cen_y   = centromeres$cen_mid,
         SIMPLIFY = FALSE)
)

# Approximate bp depth of the rounded end caps at default plot size (20x10 in).
# SHAPE_RADIUS = 3mm; y-scale ~437,500 bp/inch → ~52,000 bp per 3mm.
CAP_ZONE_BP <- 55000

# Returns the chromosome half-width (in x-axis units) at a given genomic position.
# Mirrors the chromosome polygon profile: full width in the body, narrowing at
# the centromere constriction and at the rounded end caps.
chr_half_width_at <- function(chr_val, y_pos) {
  cen_y   <- centromeres$cen_mid[centromeres$Chr == chr_val]
  chr_len <- chr_lengths$length[chr_lengths$Chr == chr_val]
  y_bp <- c(0,             CAP_ZONE_BP,
            cen_y - CEN_ZONE, cen_y,  cen_y + CEN_ZONE,
            chr_len - CAP_ZONE_BP, chr_len)
  w_bp <- c(0,             CHR_W,
            CHR_W,         CEN_W,  CHR_W,
            CHR_W,         0)
  pmax(0, approx(y_bp, w_bp, xout = y_pos, rule = 2)$y)
}

# ==============================================================================
# 3. LOAD & PARSE RAW DATA
# ==============================================================================

raw <- read.table(
  file.path(SUBDIR, "countries_integrated_selection_analysis_windows.tsv"),
  header = TRUE, sep = "\t", fill = TRUE, quote = "", comment.char = ""
)

raw_parsed <- raw %>%
  filter(Chr != "Chr") %>%
  mutate(
    Chr            = as.integer(Chr),
    combined_start = as.numeric(combined_start),
    combined_end   = as.numeric(combined_end),
    n_datasets     = as.integer(n_datasets)
  ) %>%
  filter(!is.na(Chr), !is.na(combined_start), !is.na(n_datasets))

# ==============================================================================
# 4. ALL-COUNTRIES WINDOWS  (n_datasets >= 2)
# ==============================================================================

windows_all <- raw_parsed %>%
  filter(n_datasets >= 2) %>%
  group_by(Chr, combined_start, combined_end) %>%
  summarise(n_datasets = max(n_datasets), .groups = "drop") %>%
  left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")

# ==============================================================================
# 5. PER-COUNTRY WINDOWS
# ==============================================================================

extract_country_long <- function(df) {
  ihs <- df %>%
    filter(!is.na(iHS_Region) & iHS_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              country = iHS_Region, tool = "iHS")
  ibd <- df %>%
    filter(!is.na(IBD_Region) & IBD_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              country = IBD_Region, tool = "IBD")
  td <- df %>%
    filter(!is.na(TD_Region) & TD_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              country = TD_Region, tool = "TD")
  # XPEHH_pop1 is always the selected population (positive XPEHH, focal pop)
  xp <- df %>%
    filter(!is.na(XPEHH_pop1) & XPEHH_pop1 != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              country = XPEHH_pop1, tool = "XP-EHH")
  bind_rows(ihs, ibd, td, xp) %>%
    filter(!is.na(country), country != "")
}

country_long <- extract_country_long(raw_parsed)

windows_by_country <- country_long %>%
  group_by(Chr, win_start, win_end, country) %>%
  summarise(n_datasets = n_distinct(tool), .groups = "drop") %>%
  filter(n_datasets >= 2) %>%
  left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")

# ==============================================================================
# 6. GENE LABEL EXTRACTION
# ==============================================================================

parse_gene_names <- function(df) {
  df %>%
    pivot_longer(cols = ends_with("_genes"), values_to = "genes") %>%
    filter(!is.na(genes), genes != "", genes != "NA(NA)") %>%
    # Extract everything inside parentheses, then drop placeholder values
    mutate(gene_name = str_extract_all(genes, "(?<=\\()[^)]+(?=\\))")) %>%
    unnest(gene_name) %>%
    filter(!gene_name %in% c("N/A", "NA", "N", "")) %>%
    mutate(gene_name = recode(gene_name, "PPPK-DHPS" = "DHPS"))
}

extract_gene_labels <- function(windows_df, raw_src,
                                 start_col, end_col, min_tools = 3) {
  sig <- windows_df %>% filter(n_datasets >= min_tools)
  if (nrow(sig) == 0) return(data.frame())
  raw_src %>%
    inner_join(sig, by = setNames(c("Chr", start_col, end_col),
                                   c("Chr", "combined_start", "combined_end"))) %>%
    select(Chr, combined_start, combined_end,
           XPEHH_genes, IBD_genes, iHS_genes, TD_genes) %>%
    parse_gene_names() %>%
    group_by(Chr, combined_start, combined_end) %>%
    summarise(label = paste(unique(gene_name), collapse = "/"), .groups = "drop") %>%
    mutate(pos_mid = (combined_start + combined_end) / 2) %>%
    group_by(Chr, pos_mid) %>% slice(1) %>% ungroup() %>%
    left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")
}

get_country_gene_labels <- function(country_name, windows_c, raw_parsed,
                                     min_tools = 2) {
  sig <- windows_c %>% filter(country == country_name, n_datasets >= min_tools)
  if (nrow(sig) == 0) return(data.frame())
  raw_parsed %>%
    inner_join(sig, by = c("Chr",
                           "combined_start" = "win_start",
                           "combined_end"   = "win_end")) %>%
    select(Chr, combined_start, combined_end,
           XPEHH_genes, IBD_genes, iHS_genes, TD_genes) %>%
    parse_gene_names() %>%
    group_by(Chr, combined_start, combined_end) %>%
    summarise(label = paste(unique(gene_name), collapse = "/"), .groups = "drop") %>%
    mutate(pos_mid = (combined_start + combined_end) / 2) %>%
    group_by(Chr, pos_mid) %>% slice(1) %>% ungroup() %>%
    left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")
}

gene_labels_all <- extract_gene_labels(
  windows_all, raw_parsed,
  start_col = "combined_start", end_col = "combined_end", min_tools = 3
)

# ==============================================================================
# 7. COLOR PALETTE & CENTROMERE POINTS
# ==============================================================================

# Distinct colours (not shades of red)
tool_colors <- c(
  "1" = "#AAAAAA",   # grey  (unused after n>=2 filter, kept for completeness)
  "2" = "#3288BD",   # blue
  "3" = "#F46D43",   # orange
  "4" = "#66C2A5"    # teal
)
tool_labels <- c("1" = "1 tool",  "2" = "2 tools",
                 "3" = "3 tools", "4" = "4 tools (all)")

cen_pts <- centromeres %>% rename(x = x_pos, y = cen_mid)

# ==============================================================================
# 8. PLOT FUNCTION
# ==============================================================================

make_karyogram <- function(windows_df, start_col, end_col,
                            gene_labs, title_str, subtitle_str) {

  win_segs <- windows_df %>%
    mutate(
      win_mid = (.data[[start_col]] + .data[[end_col]]) / 2,
      half_w  = mapply(chr_half_width_at, Chr, win_mid),
      n_fac   = factor(as.character(n_datasets), levels = as.character(1:4)),
      xmin    = x_pos - half_w,
      xmax    = x_pos + half_w,
      ymin    = .data[[start_col]],
      ymax    = .data[[end_col]]
    )

  # x-axis: spaced chromosome positions and labels
  x_breaks <- chr_lengths$x_pos
  x_labels <- as.character(chr_lengths$Chr)

  p <- ggplot() +

    # Chromosome white fill
    geom_shape(data = chr_polys, aes(x = x, y = y, group = group),
               fill = "white", color = NA, radius = SHAPE_RADIUS) +

    # Selection windows
    geom_rect(data = win_segs,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = n_fac),
              color = NA, alpha = 0.9) +

    # Chromosome outline (on top for clean border)
    geom_shape(data = chr_polys, aes(x = x, y = y, group = group),
               fill = NA, color = "black", linewidth = 0.5,
               radius = SHAPE_RADIUS) +

    # Centromere solid black disc
    geom_point(data = cen_pts, aes(x = x, y = y),
               shape = 16, color = "black", size = 4.5)

  if (!is.null(gene_labs) && nrow(gene_labs) > 0) {
    # Repelled labels: anchor point is just to the right of each chromosome,
    # direction = "y" keeps labels aligned to their chromosome's x column
    # Anchor exactly at chromosome right edge; repel only vertically
    lab_anchor <- gene_labs %>%
      mutate(x_anchor = x_pos + CHR_W, y_anchor = pos_mid)

    p <- p +
      geom_text_repel(
        data               = lab_anchor,
        aes(x = x_anchor, y = y_anchor, label = label),
        direction          = "y",       # only shift vertically to avoid overlap
        hjust              = 0,
        color              = "black",
        size               = 2.5,
        nudge_x            = 0.08,      # small fixed push right of chromosome edge
        segment.color      = "grey40",
        segment.size       = 0.3,
        segment.linetype   = "solid",
        min.segment.length = 0,
        box.padding        = 0.2,
        point.padding      = 0,
        max.overlaps       = Inf,
        seed               = 42
      )
  }

  p +
    scale_fill_manual(values = tool_colors, labels = tool_labels,
                      name = "Tools detecting window", drop = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels,
                       expand = expansion(add = c(0.8, 2.2))) +
    scale_y_reverse(labels = function(x) paste0(round(x / 1e6, 1), " Mb"),
                    expand = expansion(mult = 0.03)) +
    labs(x = "Chromosome", y = "Position",
         title = title_str, subtitle = subtitle_str) +
    theme_classic(base_size = 11) +
    theme(
      panel.grid      = element_blank(),
      axis.line       = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 9, color = "grey40")
    )
}

# ==============================================================================
# 9. ALL-COUNTRIES FIGURE
# ==============================================================================

p_all <- make_karyogram(
  windows_all,
  start_col    = "combined_start",
  end_col      = "combined_end",
  gene_labs    = gene_labels_all,
  title_str    = "Integrated Selection Signals - All Populations",
  subtitle_str = "Windows detected by >=2 tools | iHS, XP-EHH, Tajima's D, IBD"
)

ggsave(file.path(SUBDIR, "karyogram_all_countries.pdf"), plot = p_all,
       width = 20, height = 10, device = "pdf")
message("Saved: ", SUBDIR, "/karyogram_all_countries.pdf/.png")

# ==============================================================================
# 10. PER-COUNTRY FIGURES
# ==============================================================================

countries <- sort(unique(windows_by_country$country))
message("Generating per-country figures for: ", paste(countries, collapse = ", "))

for (ctry in countries) {
  win_c <- windows_by_country %>% filter(country == ctry)
  glabs <- get_country_gene_labels(ctry, windows_by_country, raw_parsed,
                                    min_tools = 2)

  p_c <- make_karyogram(
    win_c,
    start_col    = "win_start",
    end_col      = "win_end",
    gene_labs    = glabs,
    title_str    = paste0("Selection Signals - ", ctry),
    subtitle_str = "Windows detected by >=2 tools | iHS, XP-EHH, Tajima's D, IBD"
  )

  fname <- paste0("karyogram_", gsub(" ", "_", ctry))
  ggsave(file.path(SUBDIR, paste0(fname, ".pdf")), plot = p_c, width = 20, height = 10, device = "pdf")
  message("  Saved: ", SUBDIR, "/", fname)
}


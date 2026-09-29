################################################################################
# KARYOGRAM: INTEGRATED SELECTION SIGNALS - SITES GROUPING
# Groups: Co_endemic_sites (Gondar + Mizan + Metehera Fentale),
#         Gambella, Arba_Minch
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggforce)
library(ggrepel)

BASE_DIR <- "/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection"
setwd(BASE_DIR)

IN_FILE      <- "sites_grouping/sites_grouping_integrated_selection_analysis_windows.tsv"
CHR_SPACING  <- 1.6
CHR_W        <- 0.28
CEN_W        <- 0.06
CEN_ZONE     <- 40000
SHAPE_RADIUS <- unit(3, "mm")

# ==============================================================================
# 1. CHROMOSOME REFERENCE DATA
# ==============================================================================

chr_lengths <- data.frame(
  Chr    = 1:14,
  length = c(643292, 947102, 1067971, 1200490, 1343557, 1418242, 1445207,
             1472805, 1541735, 1687656, 2038340, 2271494, 2925236, 3291936)
) %>% mutate(x_pos = Chr * CHR_SPACING)

centromeres <- data.frame(
  Chr     = 1:14,
  cen_mid = c(
    (458021 + 462661) / 2, (446771 + 450941) / 2,
    (593691 + 597951) / 2, (648001 + 652321) / 2,
    (454541 + 458791) / 2, (477751 + 482011) / 2,
    (863741 + 868251) / 2, (299011 + 303631) / 2,
    (1241071 + 1245441) / 2, (935681 + 937821) / 2,
    (830781 + 835431) / 2, (1281511 + 1285931) / 2,
    (1167271 + 1171921) / 2, (1070851 + 1075311) / 2
  )
) %>% left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")

# ==============================================================================
# 2. CHROMOSOME SHAPE
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

CAP_ZONE_BP <- 55000

chr_half_width_at <- function(chr_val, y_pos) {
  cen_y   <- centromeres$cen_mid[centromeres$Chr == chr_val]
  chr_len <- chr_lengths$length[chr_lengths$Chr == chr_val]
  y_bp <- c(0, CAP_ZONE_BP, cen_y - CEN_ZONE, cen_y, cen_y + CEN_ZONE,
            chr_len - CAP_ZONE_BP, chr_len)
  w_bp <- c(0, CHR_W, CHR_W, CEN_W, CHR_W, CHR_W, 0)
  pmax(0, approx(y_bp, w_bp, xout = y_pos, rule = 2)$y)
}

# ==============================================================================
# 3. LOAD & PARSE DATA
# ==============================================================================

raw <- read.table(IN_FILE, header = TRUE, sep = "\t",
                  fill = TRUE, quote = "", comment.char = "")

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
# 4. ALL-GROUPS WINDOWS  (n_datasets >= 2)
# ==============================================================================

windows_all <- raw_parsed %>%
  filter(n_datasets >= 2) %>%
  group_by(Chr, combined_start, combined_end) %>%
  summarise(n_datasets = max(n_datasets), .groups = "drop") %>%
  left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")

# ==============================================================================
# 5. PER-GROUP WINDOWS
# ==============================================================================

extract_group_long <- function(df) {
  ihs <- df %>%
    filter(!is.na(iHS_Region) & iHS_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              group = iHS_Region, tool = "iHS")
  ibd <- df %>%
    filter(!is.na(IBD_Region) & IBD_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              group = IBD_Region, tool = "IBD")
  td <- df %>%
    filter(!is.na(TD_Region) & TD_Region != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              group = TD_Region, tool = "TD")
  # XPEHH_pop1 is always the selected population (positive XPEHH, focal pop)
  xp <- df %>%
    filter(!is.na(XPEHH_pop1) & XPEHH_pop1 != "") %>%
    transmute(Chr, win_start = combined_start, win_end = combined_end,
              group = XPEHH_pop1, tool = "XP-EHH")
  bind_rows(ihs, ibd, td, xp) %>%
    filter(!is.na(group), group != "")
}

group_long <- extract_group_long(raw_parsed)

windows_by_group <- group_long %>%
  group_by(Chr, win_start, win_end, group) %>%
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
           any_of(c("XPEHH_genes", "IBD_genes", "iHS_genes", "TD_genes"))) %>%
    parse_gene_names() %>%
    group_by(Chr, combined_start, combined_end) %>%
    summarise(label = paste(unique(gene_name), collapse = "/"), .groups = "drop") %>%
    mutate(pos_mid = (combined_start + combined_end) / 2) %>%
    group_by(Chr, pos_mid) %>% slice(1) %>% ungroup() %>%
    left_join(chr_lengths %>% select(Chr, x_pos), by = "Chr")
}

get_group_gene_labels <- function(grp, windows_g, raw_parsed, min_tools = 2) {
  sig <- windows_g %>% filter(group == grp, n_datasets >= min_tools)
  if (nrow(sig) == 0) return(data.frame())
  raw_parsed %>%
    inner_join(sig, by = c("Chr",
                           "combined_start" = "win_start",
                           "combined_end"   = "win_end")) %>%
    select(Chr, combined_start, combined_end,
           any_of(c("XPEHH_genes", "IBD_genes", "iHS_genes", "TD_genes"))) %>%
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

tool_colors <- c(
  "1" = "#AAAAAA",
  "2" = "#3288BD",
  "3" = "#F46D43",
  "4" = "#66C2A5"
)
tool_labels <- c("1" = "1 tool", "2" = "2 tools",
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

  x_breaks <- chr_lengths$x_pos
  x_labels <- as.character(chr_lengths$Chr)

  p <- ggplot() +
    geom_shape(data = chr_polys, aes(x = x, y = y, group = group),
               fill = "white", color = NA, radius = SHAPE_RADIUS) +
    geom_rect(data = win_segs,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = n_fac),
              color = NA, alpha = 0.9) +
    geom_shape(data = chr_polys, aes(x = x, y = y, group = group),
               fill = NA, color = "black", linewidth = 0.5,
               radius = SHAPE_RADIUS) +
    geom_point(data = cen_pts, aes(x = x, y = y),
               shape = 16, color = "black", size = 4.5)

  if (!is.null(gene_labs) && nrow(gene_labs) > 0) {
    lab_anchor <- gene_labs %>%
      mutate(x_anchor = x_pos + CHR_W, y_anchor = pos_mid)
    p <- p +
      geom_text_repel(
        data               = lab_anchor,
        aes(x = x_anchor, y = y_anchor, label = label),
        direction          = "y",
        hjust              = 0,
        color              = "black",
        size               = 2.5,
        nudge_x            = 0.08,
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
# 9. ALL-GROUPS FIGURE
# ==============================================================================

p_all <- make_karyogram(
  windows_all,
  start_col    = "combined_start",
  end_col      = "combined_end",
  gene_labs    = gene_labels_all,
  title_str    = "Integrated Selection Signals - Site Groups",
  subtitle_str = "Windows detected by >=2 tools | iHS, XP-EHH, Tajima's D, IBD"
)

ggsave("sites_grouping/karyogram_sites_grouping_all.pdf", plot = p_all,
       width = 20, height = 10, device = "pdf")
message("Saved: karyogram_sites_grouping_all.pdf")

# ==============================================================================
# 10. PER-GROUP FIGURES
# ==============================================================================

groups <- sort(unique(windows_by_group$group))
message("Generating per-group figures for: ", paste(groups, collapse = ", "))

for (grp in groups) {
  win_g <- windows_by_group %>% filter(group == grp)
  glabs <- get_group_gene_labels(grp, windows_by_group, raw_parsed, min_tools = 2)

  p_g <- make_karyogram(
    win_g,
    start_col    = "win_start",
    end_col      = "win_end",
    gene_labs    = glabs,
    title_str    = paste0("Selection Signals - ", gsub("_", " ", grp)),
    subtitle_str = "Windows detected by >=2 tools | iHS, XP-EHH, Tajima's D, IBD"
  )

  fname <- paste0("sites_grouping/karyogram_", grp)
  ggsave(paste0(fname, ".pdf"), plot = p_g, width = 20, height = 10, device = "pdf")
  message("  Saved: ", fname)
}

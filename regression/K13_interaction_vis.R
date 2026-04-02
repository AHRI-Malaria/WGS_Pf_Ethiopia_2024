
library(ggplot2)
library(ggrepel)
library(dplyr)
setwd("/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/k13_interaction")

populations <- c("EthEri", "Ethiopia", "EthNew", "HOA")

for (pop in populations) {

  data <- read.csv(paste0(pop, "_kelch13_interaction_results_ann.csv"), sep = ",")

  plot_data <- data %>%
    mutate(
      logp = -log10(p_adj),
      chr = factor(chr, levels = paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3")),
      pos_mb = pos / 1e6
    )

  plot_data <- plot_data %>%
    mutate(significant = ifelse(p_adj < 0.05, "Significant", "Not significant"))

  plot_data <- plot_data %>%
    mutate(
      Gene_name_clean = sub("\\s*\\(.*\\)", "", Gene_name)
    ) %>%
    filter(chr != "NA")

  # Compute cumulative positions for single-line Manhattan layout
  chr_lengths <- plot_data %>%
    group_by(chr) %>%
    summarise(max_pos = max(pos_mb)) %>%
    arrange(chr) %>%
    mutate(cum_offset = cumsum(lag(max_pos, default = 0)))

  plot_data <- plot_data %>%
    left_join(chr_lengths, by = "chr") %>%
    mutate(cum_pos = pos_mb + cum_offset)

  # Chromosome number for alternating colours and axis labels
  plot_data <- plot_data %>%
    mutate(chr_num = as.integer(gsub("Pf3D7_(\\d+)_v3", "\\1", chr)))

  # Assign colour: significant = red, others alternate blue shades by chromosome
  plot_data <- plot_data %>%
    mutate(point_col = case_when(
      significant == "Significant" ~ "sig",
      chr_num %% 2 == 1 ~ "odd",
      TRUE ~ "even"
    ))

  # Axis tick positions at chromosome midpoints
  axis_ticks <- plot_data %>%
    group_by(chr_num) %>%
    summarise(center = (min(cum_pos) + max(cum_pos)) / 2)

  # K13 position in cumulative coords
  k13_pos_mb <- 1725133 / 1e6
  k13_offset <- chr_lengths %>% filter(chr == "Pf3D7_13_v3") %>% pull(cum_offset)
  k13_cum <- k13_pos_mb + k13_offset

  k13_label <- data.frame(
    cum_pos = k13_cum,
    logp = max(plot_data$logp, na.rm = TRUE) + 0.3
  )

  # Label data for significant genes
  label_data <- plot_data %>%
    filter(
      p_adj < 0.05 &
        !is.na(Gene_name) &
        Gene_name_clean != "NA" &
        Gene_name_clean != "NA (NA)" &
        Gene_name_clean != "." &
        !grepl("^@", Gene_name) &
        !grepl("^PF3D7", Gene_name) &
        !grepl("^SURF", Gene_name) &
        !grepl("^RIF", Gene_name) &
        !grepl("^VAR", Gene_name, ignore.case = TRUE) &
        !grepl("^\\s*$", Gene_name)
    ) %>%
    mutate(Gene_name = sub("\\s*\\(.*\\)", "", Gene_name))

  plot <- ggplot(plot_data, aes(x = cum_pos, y = logp, color = point_col)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
    geom_text_repel(
      data = label_data,
      aes(label = Gene_name),
      size = 3,
      color = "red",
      max.overlaps = 30,
      min.segment.length = 0,
      segment.color = "grey60"
    ) +
    geom_point(data = k13_label, aes(x = cum_pos, y = logp),
               color = "blue", shape = 17, size = 3, inherit.aes = FALSE) +
    geom_text_repel(data = k13_label, aes(x = cum_pos, y = logp, label = "K13 R622I"),
                    color = "blue", size = 3, fontface = "bold",
                    min.segment.length = 0, segment.color = "blue",
                    inherit.aes = FALSE) +
    scale_color_manual(
      values = c("sig" = "red", "odd" = "steelblue4", "even" = "steelblue2")
    ) +
    scale_x_continuous(breaks = axis_ticks$center, labels = axis_ticks$chr_num) +
    theme_minimal() +
    labs(
      title = pop,
      x = "Chromosome",
      y = expression(-log[10](p[adj])),
      color = ""
    ) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

  ggsave(paste0(pop, "_R622I_interaction_manhatten.pdf"), plot,
         width = 30, height = 15, units = "cm", dpi = 600)

}

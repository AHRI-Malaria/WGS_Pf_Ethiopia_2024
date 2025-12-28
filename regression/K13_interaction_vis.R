

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# With PC components
data<-read.csv("kelch13_interaction_results_ann.csv",sep=",")

plot_data <- data %>%
  mutate(
    logp = -log10(p_adj),                # avoid log(0)
    chr = factor(chr, levels = paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3")),
    pos_mb = pos / 1e6                            # convert to Mb
  )

plot_data <- plot_data %>%
  mutate(significant = ifelse(p_adj < 0.05, "Significant", "Not significant"))

plot_data <- plot_data %>%
  mutate(
    Gene_name_clean = sub("\\s*\\(.*\\)", "", Gene_name)  # keep only name before '('
  ) %>%
  filter(chr!="NA")


plot<-ggplot(plot_data, aes(x = pos_mb, y = logp, color = significant)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_text_repel(
    data = subset(
      plot_data,
      p_adj < 0.05 &
        !is.na(Gene_name) &
        Gene_name_clean != "NA" &
        Gene_name_clean != "NA (NA)" &
        !grepl("^PF3D7", Gene_name)  &         # exclude PF3D7 genes
        !grepl("^\\s*$", Gene_name)            # exclude empty names
    ) %>%
      transform(Gene_name = sub("\\s*\\(.*\\)", "", Gene_name)),  # keep only name before '('
    aes(label = Gene_name),
    size = 3,
    max.overlaps = 20,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  facet_wrap(~ chr, nrow = 2, scales = "free_x") +
  scale_color_manual(
    values = c("Significant" = "red", "Not significant" = "grey70")
  ) +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](p[adj])),
    color = ""
  ) +
  theme(legend.position = "none")

ggsave("R622I_interaction_manhatten.pdf", plot,
       width = 30, height = 18, units = "cm", dpi = 600)


# Without PC components
data<-read.csv("kelch13_interaction_results_ann_noPC.csv",sep=",")

plot_data <- data %>%
  mutate(
    logp = -log10(p_adj),                # avoid log(0)
    chr = factor(chr, levels = paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3")),
    pos_mb = pos / 1e6                            # convert to Mb
  )

plot_data <- plot_data %>%
  mutate(significant = ifelse(p_adj < 0.05, "Significant", "Not significant"))

plot_data <- plot_data %>%
  mutate(
    Gene_name_clean = sub("\\s*\\(.*\\)", "", Gene_name)  # keep only name before '('
  ) %>%
  filter(chr!="NA")


plot<-ggplot(plot_data, aes(x = pos_mb, y = logp, color = significant)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_text_repel(
    data = subset(
      plot_data,
      p_adj < 0.05 &
        !is.na(Gene_name) &
        Gene_name_clean != "NA" &
        Gene_name_clean != "NA (NA)" &
        !grepl("^PF3D7", Gene_name)  &         # exclude PF3D7 genes
        !grepl("^\\s*$", Gene_name)            # exclude empty names
    ) %>%
      transform(Gene_name = sub("\\s*\\(.*\\)", "", Gene_name)),  # keep only name before '('
    aes(label = Gene_name),
    size = 3,
    max.overlaps = 20,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  facet_wrap(~ chr, nrow = 2, scales = "free_x") +
  scale_color_manual(
    values = c("Significant" = "red", "Not significant" = "grey70")
  ) +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](p[adj])),
    color = ""
  ) +
  theme(legend.position = "none")

ggsave("R622I_interaction_manhatten_noPC.pdf", plot,
       width = 30, height = 18, units = "cm", dpi = 600)

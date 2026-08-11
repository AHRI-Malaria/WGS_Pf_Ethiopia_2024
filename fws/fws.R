
library(ggplot2)
library(patchwork)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the data
meta <- read.csv("../input_files/Metadata_Popgen_Ethiopia_withfws.csv", stringsAsFactors = FALSE)

# Convert Site to factor (controls order)
meta$site <- factor(meta$site, levels = unique(meta$site))

# Count samples per site
site_counts <- meta %>%
  group_by(site) %>%
  summarise(n = n()) %>%
  mutate(label = paste0("n=", n))

# Calculate median fws per site
site_medians <- meta %>%
  group_by(site) %>%
  summarise(median_fws = median(fws, na.rm = TRUE))

# Calculate proportion of polyclonal samples (Fws < 0.95) per site
site_polyclonal <- meta %>%
  group_by(site) %>%
  summarise(
    n = n(),
    n_polyclonal = sum(fws < 0.95, na.rm = TRUE),
    prop_polyclonal = n_polyclonal / n
  )
  
site_colors <- c(
  "Arba Minch" = "#F8766D",       # Red/pink
  "Gambella" = "#B79F00",         # Yellow/gold
  "Gondar" = "#00BA38",           # Green
  "Metehera Fentale" = "#00BFC4", # Cyan
  "Mizan" = "#F564E3"             # Pink/magenta
)

meta$site<-factor(meta$site, levels=c("Gambella","Mizan","Gondar","Arba Minch","Metehera Fentale"))
# Plot
ggplot(meta, aes(x = site, y = fws, color = site)) +
  # Boxplot
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  
  # Individual datapoints
  geom_jitter(width = 0.25, alpha = 0.7, size = 2) +
  
  # Dashed gray line for overall median
  geom_hline(yintercept = median(meta$fws, na.rm = TRUE), 
             linetype = "dotted", color = "gray40", linewidth = 1) +
  
  # n= sample count labels
  geom_text(data = site_counts, aes(x = site, y = 1.02, label = label),
            inherit.aes = FALSE, size = 4, color = "black") +
  
  scale_color_manual(values = site_colors) +
  
  theme_test() +
  labs(y = expression(F[WS]), x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("figs/fws.pdf", width=3, height=5)


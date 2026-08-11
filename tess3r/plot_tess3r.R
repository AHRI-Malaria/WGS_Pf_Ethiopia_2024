library(tess3r)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(showtext)
library(unikn)
library(tidyr)
library(rworldmap)
showtext_auto()

## Load funcitons with themes for plots
# For PCA
theme_PCA <- function(Fsize = 10, Fcolor = "black", Fangle = 0, Fhjust = 0.5) {
    theme_classic() +
    theme(axis.text.x = element_text(size = Fsize, color = Fcolor, angle = Fangle,
          hjust= Fhjust),
          axis.text.y = element_text(size = Fsize, color = Fcolor),
          axis.title.x = element_text(size = Fsize, color = Fcolor),
          axis.title.y = element_text(size = Fsize, color = Fcolor),
          legend.text = element_text(color = Fcolor, size = Fsize - 1,
                                     margin = margin(l = 0, r = 10)),
          legend.title = element_text(color = Fcolor, size = Fsize,
                                      hjust = Fhjust),
          legend.key.size = unit(0.2, "mm"),
          legend.position = "bottom")
}

# For ADMIXTURE
theme_ADMIX <- function(Fsize = 10, Fcolor = "black", Fangle = 0, Lposition = "none", Fhjust = 1) {
    theme_minimal() +
        theme(strip.text.x = element_text(angle = Fangle, hjust = Fhjust,
                                         size = Fsize - 3, color = "black",
                                         margin = margin(t = 0.25, r = 0.5, b = 0.25, l = 0.5, unit = "lines")),
            panel.spacing.x = unit(0.2, "lines"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = Fsize, color = "black"),
            axis.text.y = element_text(size = Fsize, color = "black",
                                       angle = 90),
            legend.text = element_text(color = Fcolor, size = Fsize - 1,
                                     margin = margin(l = 0, r = 20)),
            panel.grid = element_blank(),
            legend.position = Lposition)
}
# For MAP
theme_tess3r_map <- function(Fsize = 10, Fcolor = "black", Lposition="none") {
    theme_classic() +
    theme(axis.line.x = element_line(color = Fcolor),
          axis.line.y = element_line(color = Fcolor),
          axis.text.x = element_text(size = Fsize, color = Fcolor),
          axis.text.y = element_text(size = Fsize - 2, color = Fcolor),
          axis.title.x = element_text(size = Fsize, color = Fcolor),
          axis.title.y = element_text(size = Fsize, color = Fcolor),
          plot.title = element_text(size = Fsize, color = Fcolor, hjust = 0.5),          
          legend.position = Lposition)
}

#################################################################################
workdir <- "/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/tess3r/run2"

metadata <- readr::read_csv(file.path(workdir, "metadata_db_tess3r_022026.csv")) %>%
    mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude)) 
matbinT <- readRDS(file.path(workdir, "Pf_fulldb_022026_matbinT.mis.rds"))
coordsM <- readRDS(file.path(workdir, "Pf_fulldb_022026_coordsM.rds"))
samples <- readRDS(file.path(workdir, "Pf_fulldb_022026.rds"))

# Load tess3 object
tess3_all <- list()
Ks <- seq(1,10)
tess3_all <- setNames(
  lapply(Ks, function(K) readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", K)))),
  paste0("K", Ks)
)
res1 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 1)))
res2 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 2)))
res3 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 3)))
res4 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 4)))
res5 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 5)))
res6 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 6)))
res7 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 7)))
res8 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 8)))
res9 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 9)))
res10 <- readRDS(file.path(workdir, sprintf("Pf_fulldb_022026_rep50_K%d.rds", 10)))

# Plot error
Ks <- 1:10
scores <- tibble::tibble(
  K = Ks,
  crossentropy = sapply(Ks, function(k) Gettess3res(tess3_all[[paste0("K", k)]], K = k)$crossentropy),
  rmse         = sapply(Ks, function(k) Gettess3res(tess3_all[[paste0("K", k)]], K = k)$rmse)
)

print(scores)

ggplot(scores, aes(K, crossentropy)) + geom_line() + geom_point() +
  labs(y = "Cross-entropy (CV error)")

#Diminishing-returns (Δ% rule)
#Choose the smallest K where the relative improvement from K to K+1 drops below a small threshold (e.g., 1%).
ce_mean <- sapply(Ks, \(k) Gettess3res(tess3_all[[paste0("K",k)]], K=k)$crossentropy)
scores <- data.frame(K = Ks, ce = ce_mean)
scores <- transform(scores,
  rel_impr = c(NA, diff(ce)/head(ce, -1) * -1)  # positive = improvement
)
K_delta <- with(scores, min(K[which(rel_impr < 0.01)], na.rm = TRUE))  # 1% cutoff
K_delta

print(scores)
# Percent improvement from K-1 -> K
scores$rel_impr_pct <- round(100 * scores$rel_impr, 2)
scores

ggplot(scores, aes(K, ce)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = K_delta, linetype = 2) +
  labs(y = "Cross-entropy (CV error)", subtitle = paste("Δ% cutoff:", 100*0.01, "%; chosen K =", K_delta))


# Colors
my_cols <- unikn::usecol(c(Bordeaux, "white", Petrol), n = 6)
print(my_cols)
cols <- c("#8E2043", "#BB798E", "#CDE2E6", "#6AA9B6", "#077187",
          "#F4A261", "#E76F51", "#2A9D8F", "#E9C46A", "#264653")
names(cols) <- paste0("K", c(1:10))

# Cumulative plots
kk <- 10
res <- res10
q.matrix <- qmatrix(res, K = kk)
p.values <- pvalue(res, K = kk)

df_matrix <- as.data.frame(q.matrix)
colnames(df_matrix) <- paste0("K", c(1:ncol(df_matrix)))
df_matrix$sample_id <- samples

df_matrix_G <- df_matrix %>% mutate(ind = row_number()) %>%
    tidyr::pivot_longer(cols = matches("K"), names_to = "K", values_to = "Q") %>%
    group_by(ind) %>%
    mutate(max.K = paste0("K", which.max(Q))) %>%
    left_join(metadata %>% select(sample_id, country, Site)) %>% ungroup() %>%
    mutate(country = factor(country, levels = c("Ethiopia", "Eritrea", "Sudan", "Uganda", "Kenya", "Tanzania","India","Mali","Ghana","Nigeria")))

# Cumulative barplot max.K
# Country
ggplot(df_matrix_G) +
    geom_bar(aes(x = country, y =..count.., fill = max.K), position="fill") +
    scale_fill_manual(breaks = cols,
                      values = cols) +
    labs(fill="") +
    theme_PCA(Fsize = 8, Fangle = 90)


cumR <- df_matrix_G %>% group_by(country, K) %>%
     summarise(avg = sum(Q)/n())


# Cumulative barplot all K
ggplot(cumR, aes(x = country, y = avg, fill = K, group=K)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(breaks = names(cols),
                      values = cols,
                      labels = names(cols)) +
    geom_text(aes(x = country, y = avg, label = round(avg, 3)*100),
                  position = position_dodge(width = 1),
                  hjust = -0.2, vjust = 0.4, size = 3.5, angle = 90) +
    labs(x = "country", y = "K% per country", fill="") +
    ylim(c(0, 1)) +
    theme_PCA(Fsize = 8, Fangle = 90)


### Admixture data
admix_gg <- df_matrix_G %>%
  mutate(K = as.numeric(gsub("K", "", K))) %>%
  select(sample_id, country, Site, max.K, K, Q) %>%
  mutate(Q = as.numeric(Q)) %>%
  mutate(K = as.factor(K))

admix_gg_order <- admix_gg %>%
     mutate(country = factor(country, levels = c("Ethiopia", "Eritrea", "Sudan", "Uganda", "Kenya", "Tanzania","India","Mali","Ghana","Nigeria"))) %>%
     arrange(country) %>%
     filter(Q != 0)  %>%
     arrange(country, max.K, desc(Q), K) %>%
    mutate(sample = factor(sample_id, levels = unique(sample_id))) %>%
    mutate(K = paste0("K", K))

## Admixture plot
site_colors <- c(
  "Arba Minch" = "#F4D03F",      # Yellow
  "Gambella" = "#52BE80",         # Green
  "Gondar" = "#5DADE2",           # Blue
  "Metehera Fentale" = "#AF7AC5", # Purple
  "Mizan" = "#EC7063"             # Red/Orange
)

site_ribbon <- admix_gg_order %>%
  filter(country == "Ethiopia") %>%
  select(sample, Site) %>%
  distinct()

## First, define your site colors (using the colors from your earlier screenshot)
site_colors <- c(
  "Arba Minch" = "#f8766e",      # Yellow
  "Gambella" = "#a3a500",         # Green
  "Gondar" = "#36bf7d",           # Blue
  "Metehera Fentale" = "#2fb0f6", # Purple
  "Mizan" = "#e76bf3"             # Red/Orange
)

## Create a separate data frame for the site ribbon
site_ribbon <- admix_gg_order %>%
  filter(country == "Ethiopia") %>%
  select(sample, Site,country) %>%
  distinct()

## Modified admixture plot with site ribbon
ggplot() +
  # Site ribbon layer (only for Ethiopia samples)
  geom_tile(data = site_ribbon, 
            aes(x = sample, y = 1.05, fill = Site), 
            height = 0.05, width = 1) +
  # Main ancestry plot
  geom_col(data = admix_gg_order, 
           aes(x = sample, y = Q, fill = K), 
           width = 1, alpha = 0.9) +
  facet_grid(~country, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     limits = c(0, 1.1),  # Extended to accommodate ribbon
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c(cols, site_colors),  # Combine both color schemes
    breaks = c(names(cols), names(site_colors)),  # Control legend order
    name = ""
  ) +
  scale_x_discrete(expand = expansion(add = 3)) +
  theme_ADMIX(Fangle = 90, Lposition = "bottom", Fhjust = 0.5, Fsize = 8) +
  guides(fill = guide_legend(nrow = 2))  # Optional: organize legend in 2 rows

options(repr.plot.width = 16, repr.plot.height = 4)
ggplot(data = admix_gg_order, aes(x = sample, y = Q, fill = K)) +
  geom_col(width = 1, alpha = 0.9) +
  facet_grid(~country, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry", fill = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(add = 3)) +
  theme_ADMIX(Fangle = 90, Lposition = "bottom", Fhjust = 0.5, Fsize = 8)

ggplot(data = admix_gg_order, aes(x = sample, y = Q, fill = K)) +
  geom_col(width = 1, alpha = 0.9) +
  facet_grid(~Site, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry", fill = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(add = 3)) +
  theme_ADMIX(Fangle = 90, Lposition = "bottom", Fhjust = 0.5, Fsize = 8)

### Map
my_palette <- CreatePalette(cols[1:kk])
coordsMM <- as.data.frame(coordsM) %>% left_join(metadata %>%
 dplyr::select(country, longitude, latitude) %>% distinct(), c("V1"="longitude", "V2"="latitude")) %>%
 distinct() %>% group_by(country) %>% top_n(1, V1)
map.polygon <- rworldmap::getMap(resolution = "coarse")

pl <- ggtess3Q(q.matrix, coordsM, map.polygon = map.polygon,
  col.palette = my_palette)

pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), color = "grey") +
    xlim(-30, 100) +
    ylim(-15, 40) +
    coord_equal() +
    geom_point(data = as.data.frame(coordsM), aes(x = V1, y = V2), color = "black", size = 0.5) + 
    geom_label_repel(data = as.data.frame(coordsMM), aes(x = V1, y = V2, label = country),
        size = 2, max.overlaps = 20) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_tess3r_map()

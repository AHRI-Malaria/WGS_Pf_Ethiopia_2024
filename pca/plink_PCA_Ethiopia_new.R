library(showtext)
library(argparse)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()

# Run plink
# plink --vcf myvcf.vcf.gz --distance square --double-id --allow-extra-chr --out prefix

# Run this script
# Rscript plink_PCA_Ethiopia_new.R --workdir . \
                    #--prefix prefix \
                    #--metadata Popgen_db_samples_072025_filt_FINAL.csv \
                    #--r622i ../drugres/R622I_longformat.csv \
                    #--a675v ../drugres/A675V_longformat.csv \
                    #--k76t ../drugres/K76T_longformat.csv \
                    #--hrpdel ../hrp/sample_meta_with_hrp_status_updated.csv


# Parse command-line arguments
parser <- ArgumentParser(description = 'Run PCA and tree construction from PLINK distance matrix')
parser$add_argument('--workdir', required=TRUE, help='Working directory containing PLINK distance files')
parser$add_argument('--prefix', required=TRUE, help='Prefix for PLINK .dist and .dist.id files')
parser$add_argument('--metadata', required=TRUE, help='Metadata file path (tab-separated)')
parser$add_argument('--r622i', required=FALSE, default=NULL, help='R622I genotype file path (optional)')
parser$add_argument('--a675v', required=FALSE, default=NULL, help='A675V genotype file path (optional)')
parser$add_argument('--k76t', required=FALSE, default=NULL, help='K76T genotype file path (optional)')
parser$add_argument('--hrpdel', required=FALSE, default=NULL, help='hrp deletions file path (optional)')

args <- parser$parse_args()

workdir <- args$workdir
prefix <- args$prefix
metadata <- args$metadata
r622i_file <- args$r622i
a675v_file <- args$a675v
k76t_file <- args$k76t
hrpdel_file <- args$hrpdel

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "wgs_id"))


dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA ##
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Export dist_m to .newick to make neighbour joining tree
tree <- nj(dist_m)
write.tree(phy = tree, file = file.path(workdir, paste0(prefix, ".newick")))

# Ethiopia-only: color by District

# normalise country strings
desc$Country_norm <- tolower(trimws(desc$Country))
eth_idx <- which(desc$Country_norm == "ethiopia")

if (length(eth_idx) == 0) {
  warning("No samples with Country == 'Ethiopia' found in metadata.")
} else {
  # PCA coordinates for Ethiopian samples (from the global MDS)
  df_eth <- as.data.frame(cmd$points[eth_idx, , drop = FALSE], stringsAsFactors = FALSE)
  df_eth$District <- gsub("_", " ", desc$District[eth_idx])
  df_eth$Year <- desc$Year[eth_idx]
  df_eth$wgs_id <- desc$V1[eth_idx]
  colnames(df_eth) <- sub("^V", "PC", colnames(df_eth))

  # Clean Year
  df_eth$Year <- ifelse(df_eth$Year == "2017/2018", 2017.5, as.numeric(df_eth$Year))

  # check missing District within Ethiopian subset
  na_admin <- is.na(df_eth$District) | df_eth$District == ""
  if (any(na_admin)) {
    missing_ids <- df_eth$wgs_id[na_admin]
    message("Ethiopian samples missing District: ", paste(missing_ids, collapse = ", "))
    write.table(
      missing_ids,
      file = file.path(workdir, paste0(prefix, ".ethiopia_missing_District_wgs_ids.txt")),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }

  # Remove samples with missing District
  df_eth <- df_eth[!na_admin, ]
  
  # PC1 vs PC2
  p_eth <- ggplot(df_eth, aes(x = PC1, y = PC2, color = District)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)")
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth, filename = paste0("figs/PC_12_", prefix, "_District.pdf"), width = 6, height = 5)

  # PC1 vs PC3
  q_eth <- ggplot(df_eth, aes(x = PC1, y = PC3, color = District)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)")
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth, filename = paste0("figs/PC_13_",prefix,"_District.pdf"), width = 6, height = 5)

  # PC1 vs PC2 coloured by Year
  p_eth_year <- ggplot(df_eth, aes(x = PC1, y = PC2, color = Year)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "Collection year"
    ) +
    scale_color_gradient(
      low = "#c6dbef",
      high = "#08306b",
      breaks = c(2013, 2016, 2019, 2022, 2024)
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_year, filename = file.path(workdir, "figs", paste0("PC_12_", prefix, "_Year.pdf")), width = 6, height = 5)

  # PC1 vs PC3 coloured by Year
  q_eth_year <- ggplot(df_eth, aes(x = PC1, y = PC3, color = Year)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "Collection year"
    ) +
    scale_color_gradient(
      low = "#c6dbef",
      high = "#08306b",
      breaks = c(2013, 2016, 2019, 2022, 2024)
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_year, filename = file.path(workdir, "figs", paste0("PC_13_", prefix, "_Year.pdf")), width = 6, height = 5)
}

# annotate by R622I (only if file provided)
if (!is.null(r622i_file)) {
  R622I <- read.csv(r622i_file)
  
  # Join R622I genotype data with Ethiopian PCA data
  df_eth <- df_eth %>%
  left_join(R622I %>% select(Sample, Genotype) %>% rename(R622I_Genotype = Genotype),
            by = c("wgs_id" = "Sample"))
  
  # Create PCA plot colored by R622I genotype
  # PC1 vs PC2 with R622I genotype
  p_eth_R622I <- ggplot(df_eth, aes(x = PC1, y = PC2, color = R622I_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "R622I Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_R622I, filename = file.path(workdir, paste0("figs/PC_12_", prefix, "_R622I.pdf")), width = 6, height = 5)
  
  # PC1 vs PC3 with R622I genotype
  q_eth_R622I <- ggplot(df_eth, aes(x = PC1, y = PC3, color = R622I_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "R622I Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_R622I, filename = file.path(workdir, paste0("figs/PC_13_", prefix, "_R622I.pdf")), width = 6, height = 5)
}

# annotate by A675V (only if file provided)
if (!is.null(a675v_file)) {
  A675V <- read.csv(a675v_file)
  
  # Join A675V genotype data with Ethiopian PCA data
  df_eth <- df_eth %>%
  left_join(A675V %>% select(Sample, Genotype) %>% rename(A675V_Genotype = Genotype),
            by = c("wgs_id" = "Sample"))
  
  # Create PCA plot colored by A675V genotype
  # PC1 vs PC2 with A675V genotype
  p_eth_A675V <- ggplot(df_eth, aes(x = PC1, y = PC2, color = A675V_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "A675V Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_A675V, filename = file.path(workdir, paste0("figs/PC_12_", prefix, "_A675V.pdf")), width = 6, height = 5)
  
  # PC1 vs PC3 with A675V genotype
  q_eth_A675V <- ggplot(df_eth, aes(x = PC1, y = PC3, color = A675V_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "A675V Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_A675V, filename = file.path(workdir, paste0("figs/PC_13_", prefix, "_A675V.pdf")), width = 6, height = 5)
}

# Combined R622I/A675V plot (only if both files provided)
if (!is.null(r622i_file) && !is.null(a675v_file)) {
  
  # Create combined classification
  df_eth <- df_eth %>%
    mutate(Combined_Genotype = case_when(
      R622I_Genotype == "1/1" & A675V_Genotype == "1/1" ~ "R622I+A675V",
      R622I_Genotype == "1/1"                            ~ "R622I",
      R622I_Genotype == "0/1"                            ~ "Mixed R622I + WT",
      A675V_Genotype == "1/1"                            ~ "A675V",
      TRUE                                               ~ "WT"
    ))

  colors_combined <- c(
    "WT"               = "#6cd669",
    "R622I"            = "red",
    "Mixed R622I + WT" = "#ff7f0e",
    "A675V"            = "yellow"
  )

  # PC1 vs PC2
  p_eth_combined <- ggplot(df_eth, aes(x = PC1, y = PC2, color = Combined_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "PfK13 Genotype"
    ) +
    scale_color_manual(
      values = colors_combined,
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_combined, filename = file.path(workdir, "figs", paste0("PC_12_", prefix, "_k13_combined.pdf")), width = 6, height = 5)

  # PC1 vs PC3
  q_eth_combined <- ggplot(df_eth, aes(x = PC1, y = PC3, color = Combined_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "PfK13 Genotype"
    ) +
    scale_color_manual(
      values = colors_combined,
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_combined, filename = file.path(workdir, "figs", paste0("PC_13_", prefix, "_k13_combined.pdf")), width = 6, height = 5)
}

# annotate by K76T (only if file provided)
if (!is.null(k76t_file)) {
  K76T <- read.csv(k76t_file)
  
  # Join k76t genotype data with Ethiopian PCA data
  df_eth <- df_eth %>%
  left_join(K76T %>% select(Sample, Genotype) %>% rename(K76T_Genotype = Genotype),
            by = c("wgs_id" = "Sample"))
  
  # Create PCA plot colored by k76t genotype
  # PC1 vs PC2 with k76t genotype
  p_eth_K76T <- ggplot(df_eth, aes(x = PC1, y = PC2, color = K76T_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "PfCRT K76T Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#6cd669", "0|1" = "#ff7f0e", "1|1" = "red"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_K76T, filename = file.path(workdir, paste0("figs/PC_12_", prefix, "_K76T.pdf")), width = 6, height = 5)
  
  # PC1 vs PC3 with k76t genotype
  q_eth_K76T <- ggplot(df_eth, aes(x = PC1, y = PC3, color = K76T_Genotype)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "PfCRT K76T Genotype"
    ) +
    scale_color_manual(
      values = c("0/0" = "#6cd669", "0|1" = "#ff7f0e", "1|1" = "red"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_K76T, filename = file.path(workdir, paste0("figs/PC_13_", prefix, "_K76T.pdf")), width = 6, height = 5)
}


# annotate by hrp del (only if file provided)
if (!is.null(hrpdel_file)) {
  hrpdel <- read.csv(hrpdel_file)
  
  # Join hrpdel genotype data with Ethiopian PCA data
  df_eth <- df_eth %>%
  left_join(hrpdel %>% select(Sample, hrp_deletion_call) %>% rename(hrp_deletion_call = hrp_deletion_call),
            by = c("wgs_id" = "Sample"))
  
  # Create PCA plot colored by k76t genotype
  # PC1 vs PC2 with k76t genotype
  p_eth_hrp_del <- ggplot(df_eth, aes(x = PC1, y = PC2, color = hrp_deletion_call)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      color = "Pfhrp"
    ) +
    scale_color_manual(
      values = c("hrp2+/hrp3+" = "#6cd669", "hrp2-/hrp3=" = "#ff7f0e", "hrp2+/hrp3-" = "yellow", "hrp2-/hrp3-" = "red"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth_hrp_del, filename = file.path(workdir, paste0("figs/PC_12_", prefix, "_hrpdel.pdf")), width = 6, height = 5)
  
  # PC1 vs PC3 with k76t genotype
  q_eth_hrp_del <- ggplot(df_eth, aes(x = PC1, y = PC3, color = hrp_deletion_call)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)"),
      color = "Pfhrp"
    ) +
    scale_color_manual(
      values = c("hrp2+/hrp3+" = "#6cd669", "hrp2-/hrp3=" = "#ff7f0e", "hrp2+/hrp3-" = "yellow", "hrp2-/hrp3-" = "red"),
      na.value = "grey50"
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth_hrp_del, filename = file.path(workdir, paste0("figs/PC_13_", prefix, "_hrpdel.pdf")), width = 6, height = 5)
}
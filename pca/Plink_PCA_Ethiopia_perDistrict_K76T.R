library(dplyr)
library(ggplot2)

# Per-district PCA for all Ethiopian samples, coloured by PfCRT K76T (codon 76).
#
# For EACH district we subset the Ethiopian distance matrix to that district's
# samples and recompute a fresh PCA (classical MDS, as in Plink_PCA_Ethiopia_new.R),
# so each plot shows the within-district structure. Points are coloured by the
# K76T genotype.
#
# Districts that carry only a single allele at codon 76 are skipped: with nothing
# to contrast, the colouring carries no information (see --min_alleles).
#
# Run:
#   Rscript Plink_PCA_Ethiopia_perDistrict_K76T.R \
#       --workdir . \
#       --prefix Aug_Eth_2025 \
#       --metadata Popgen_db_samples_072025_filt_FINAL.csv \
#       --k76t ../drugres/K76T_longformat.csv

library(argparse)
parser <- ArgumentParser(description = 'Per-district PCA (Ethiopia) coloured by PfCRT K76T')
parser$add_argument('--workdir',  required = TRUE)
parser$add_argument('--prefix',   required = TRUE, help = 'Prefix for PLINK .dist and .dist.id files')
parser$add_argument('--metadata', required = TRUE)
parser$add_argument('--k76t',     required = TRUE, help = 'K76T long-format genotype file')
parser$add_argument('--min_n',    required = FALSE, type = 'integer', default = 3,
                    help = 'Minimum samples per district to attempt a PCA (default 3)')
parser$add_argument('--min_alleles', required = FALSE, type = 'integer', default = 2,
                    help = 'Minimum distinct codon-76 alleles per district (default 2; set 1 to keep monomorphic districts)')
args <- parser$parse_args()

workdir  <- args$workdir
prefix   <- args$prefix
min_n    <- args$min_n

outdir <- file.path(workdir, "figs", "per_district_K76T")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

calc_variance_explained <- function(pc_points) {
  eig <- pc_points$eig
  eig[eig < 0] <- 0
  round(eig / sum(eig) * 100, 1)
}

## ---- Normalise K76T genotype into interpretable classes -------------------
# Alt at Pf3D7_07_v3:403625 is "C,G": allele 1 = C = 76T (main resistant),
# allele 2 = G = 76R (c.227A>G, p.Lys76Arg).
classify_k76t <- function(gt) {
  gt <- trimws(gt)
  alleles <- strsplit(gt, "[/|]")[[1]]
  if (length(alleles) < 2 || any(alleles == ".") || gt == "" || is.na(gt)) return(NA_character_)
  if (any(alleles == "2"))                     return("K76R")
  if (all(alleles == "0"))                     return("K76 (WT)")
  if (all(alleles == "1"))                     return("76T")
  if ("0" %in% alleles && "1" %in% alleles)    return("Het (K76/76T)")
  NA_character_
}

# Distinct alleles called at codon 76 across a set of samples, ignoring missing.
codon76_alleles <- function(gts) {
  gts <- gts[!is.na(gts)]
  a <- unlist(strsplit(trimws(gts), "[/|]"))
  sort(unique(a[a != "." & a != ""]))
}

k76t_levels <- c("K76 (WT)", "Het (K76/76T)", "76T", "K76R")
k76t_colors <- c("K76 (WT)"       = "#6cd669",
                 "Het (K76/76T)"  = "#ff7f0e",
                 "76T"            = "red",
                 "K76R"           = "#7570b3")

## ---- Load inputs -----------------------------------------------------------
met <- read.table(args$metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id   <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

dist_m <- as.matrix(dist)
rownames(dist_m) <- id$V1
colnames(dist_m) <- id$V1

desc <- id %>% left_join(met, by = c("V1" = "wgs_id"))

K76T <- read.csv(args$k76t, stringsAsFactors = FALSE)
K76T$K76T_class <- vapply(K76T$Genotype, classify_k76t, character(1))
gt_map  <- setNames(K76T$K76T_class, K76T$Sample)
raw_map <- setNames(trimws(K76T$Genotype), K76T$Sample)

## ---- Restrict to Ethiopian samples with a District ------------------------
desc$Country_norm <- tolower(trimws(desc$Country))
desc$District_disp <- gsub("_", " ", desc$District)
eth <- desc[desc$Country_norm == "ethiopia" &
            !is.na(desc$District_disp) & desc$District_disp != "", ]

districts <- sort(unique(eth$District_disp))
message("Ethiopian districts found: ", paste(districts, collapse = ", "))

summary_rows <- list()
combined_list <- list()

for (d in districts) {
  ids <- eth$V1[eth$District_disp == d]
  ids <- ids[ids %in% rownames(dist_m)]
  n <- length(ids)
  if (n < min_n) {
    message(sprintf("Skipping %s (n=%d < min_n=%d)", d, n, min_n))
    next
  }

  # Skip districts that are monomorphic at codon 76 - the colouring would
  # carry no contrast, so the plot says nothing about K76T.
  alleles_here <- codon76_alleles(raw_map[ids])
  if (length(alleles_here) < args$min_alleles) {
    message(sprintf("Skipping %s (n=%d, only %d allele(s) at codon 76: %s)",
                    d, n, length(alleles_here),
                    paste(alleles_here, collapse = ",")))
    next
  }

  sub_m <- dist_m[ids, ids, drop = FALSE]
  k <- min(10, n - 1)
  cmd  <- cmdscale(sub_m, k = k, eig = TRUE, x.ret = TRUE)
  vars <- calc_variance_explained(cmd)

  df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
  colnames(df) <- paste0("PC", seq_len(ncol(df)))
  df$wgs_id <- ids
  df$K76T   <- factor(gt_map[ids], levels = k76t_levels)

  combined_list[[d]] <- transform(df[, c("PC1", "PC2", "wgs_id", "K76T")], District = d)

  safe <- gsub("[^A-Za-z0-9]+", "_", d)
  subtitle <- sprintf("n = %d", n)

  p12 <- ggplot(df, aes(PC1, PC2, color = K76T)) +
    geom_point(size = 2.4) +
    scale_color_manual(values = k76t_colors, drop = FALSE, na.value = "grey70",
                       name = "PfCRT K76T") +
    labs(title = d, subtitle = subtitle,
         x = sprintf("PC1 (%.1f%%)", vars[1]),
         y = sprintf("PC2 (%.1f%%)", vars[2])) +
    theme_test() + theme(legend.position = "right")
  ggsave(file.path(outdir, sprintf("PC_12_%s_%s_K76T.pdf", prefix, safe)),
         p12, width = 6, height = 5)

  if (k >= 3) {
    p13 <- ggplot(df, aes(PC1, PC3, color = K76T)) +
      geom_point(size = 2.4) +
      scale_color_manual(values = k76t_colors, drop = FALSE, na.value = "grey70",
                         name = "PfCRT K76T") +
      labs(title = d, subtitle = subtitle,
           x = sprintf("PC1 (%.1f%%)", vars[1]),
           y = sprintf("PC3 (%.1f%%)", vars[3])) +
      theme_test() + theme(legend.position = "right")
    ggsave(file.path(outdir, sprintf("PC_13_%s_%s_K76T.pdf", prefix, safe)),
           p13, width = 6, height = 5)
  }

  cnts <- table(factor(df$K76T, levels = k76t_levels), useNA = "no")
  summary_rows[[d]] <- data.frame(
    District = d, n = n,
    PC1_pct = vars[1], PC2_pct = vars[2],
    K76_WT = cnts["K76 (WT)"], Het = cnts["Het (K76/76T)"],
    T76 = cnts["76T"], K76R = cnts["K76R"],
    Missing = sum(is.na(df$K76T)),
    row.names = NULL)
  message(sprintf("  %-18s n=%2d  PC1=%.1f%% PC2=%.1f%%", d, n, vars[1], vars[2]))
}

## ---- Combined multi-panel (PC1 vs PC2, one facet per district) -------------
if (length(combined_list) > 0) {
  comb <- do.call(rbind, combined_list)
  comb$K76T <- factor(comb$K76T, levels = k76t_levels)
  pfac <- ggplot(comb, aes(PC1, PC2, color = K76T)) +
    geom_point(size = 1.8) +
    facet_wrap(~ District, scales = "free") +
    scale_color_manual(values = k76t_colors, drop = FALSE, na.value = "grey70",
                       name = "PfCRT K76T") +
    labs(title = "Per-district PCA (Ethiopia), coloured by PfCRT K76T") +
    theme_test() + theme(legend.position = "bottom")
  ggsave(file.path(outdir, sprintf("PC_12_%s_perDistrict_K76T_facet.pdf", prefix)),
         pfac, width = 12, height = 9)
}

## ---- Write a small summary table ------------------------------------------
if (length(summary_rows) > 0) {
  summ <- do.call(rbind, summary_rows)
  write.csv(summ, file.path(outdir, sprintf("%s_perDistrict_K76T_summary.csv", prefix)),
            row.names = FALSE)
  message("\nSummary written. Figures in: ", outdir)
}

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
                    #--metadata Popgen_db_samples_072025_filt_FINAL.csv


# Parse command-line arguments
parser <- ArgumentParser(description = 'Run PCA and tree construction from PLINK distance matrix')
parser$add_argument('--workdir', required=TRUE, help='Working directory containing PLINK distance files')
parser$add_argument('--prefix', required=TRUE, help='Prefix for PLINK .dist and .dist.id files')
parser$add_argument('--metadata', required=TRUE, help='Metadata file path (tab-separated)')

args <- parser$parse_args()

workdir <- args$workdir
prefix <- args$prefix
metadata <- args$metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "wgs_id"))

# Columns you rely on for colouring
cols_check <- c("Country", "Year", "Region")

# Any missing in any of those columns?
na_mask <- apply(desc[cols_check], 1, function(r) any(is.na(r)))
na_samples <- desc$V1[na_mask]

if (length(na_samples)) {
  message("Samples with missing metadata (any of Country/Year/Region): ",
          paste(na_samples, collapse = ", "))
  write.table(
    na_samples,
    file = file.path(workdir, paste0(prefix, ".missing_metadata_wgs_ids.txt")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA ##
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Export dist_m to .newick to make neighbour joining tree
tree <- nj(dist_m)
write.tree(phy = tree, file = file.path(workdir, paste0(prefix, ".newick")))

# Ethiopia-only: color by Admin.level

# normalise country strings
desc$Country_norm <- tolower(trimws(desc$Country))
eth_idx <- which(desc$Country_norm == "ethiopia")

if (length(eth_idx) == 0) {
  warning("No samples with Country == 'Ethiopia' found in metadata.")
} else {
  # PCA coordinates for Ethiopian samples (from the global MDS)
  df_eth <- as.data.frame(cmd$points[eth_idx, , drop = FALSE], stringsAsFactors = FALSE)
  df_eth$Admin.level <- gsub("_", " ", desc$Admin.level[eth_idx])
  df_eth$wgs_id <- desc$V1[eth_idx]
  colnames(df_eth) <- sub("^V", "PC", colnames(df_eth))

  # check missing Admin.level within Ethiopian subset
  na_admin <- is.na(df_eth$Admin.level) | df_eth$Admin.level == ""
  if (any(na_admin)) {
    missing_ids <- df_eth$wgs_id[na_admin]
    message("Ethiopian samples missing Admin.level: ", paste(missing_ids, collapse = ", "))
    write.table(
      missing_ids,
      file = file.path(workdir, paste0(prefix, ".ethiopia_missing_Admin.level_wgs_ids.txt")),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }

  # PC1 vs PC2
  p_eth <- ggplot(df_eth, aes(x = PC1, y = PC2, color = Admin.level)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)")
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = p_eth, filename = "figs/PC_12_Ethiopia_new_Admin.level.pdf", width = 6, height = 5)

  # PC1 vs PC3
  q_eth <- ggplot(df_eth, aes(x = PC1, y = PC3, color = Admin.level)) +
    geom_point(size = 2) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC3 (", vars["PC3"], "%)")
    ) +
    theme_test() +
    theme(legend.position = "right")
  ggsave(plot = q_eth, filename = "figs/PC_13_Ethiopia_new_Admin.level.pdf", width = 6, height = 5)
}

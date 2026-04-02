################################################################################
# EPISTASIS ANALYSIS - KELCH13 INTERACTION SCREEN
# Please note this code has been cleaned using Claude
################################################################################

library(data.table)
library(dplyr)
library(parallel)
library(stringr)
library(tidyr)

setwd("/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/k13")

# ==============================================================================
# USER-DEFINED INPUTS
# ==============================================================================

n_cores   <- 16
kelch_row <- NULL   # Set to row number or NULL to auto-detect
kelch_id  <- 1725133   # SNP position used to auto-detect kelch13 row

GENO_PREFIX <- "012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps"

# Each scenario: label, genotype matrix, PCA file, output prefix
scenarios <- list(
  list(
    label       = "EthEri",
    geno_file   = paste0(GENO_PREFIX, ".EthEri.filt.mono.maf0.01.mat.bin"),
    pcs_file    = "PC_032026_pca_coordinates.csv",
    output_file = "EthEri_kelch13_interaction_results.csv"
  ),
  list(
    label       = "HOA",
    geno_file   = paste0(GENO_PREFIX, ".HOA.filt.mono.maf0.01.mat.bin"),
    pcs_file    = "PC_HOA_pca_coordinates.csv",
    output_file = "HOA_kelch13_interaction_results.csv"
  ),
  list(
    label       = "Ethiopia",
    geno_file   = paste0(GENO_PREFIX, ".Ethiopia.filt.mono.maf0.01.mat.bin"),
    pcs_file    = "PC_Ethiopia_pca_coordinates.csv",
    output_file = "Ethiopia_kelch13_interaction_results.csv"
  ),
  list(
    label       = "EthNew",
    geno_file   = paste0(GENO_PREFIX, ".EthNew.filt.mono.maf0.01.mat.bin"),
    pcs_file    = "PC_EthNew_pca_coordinates.csv",
    output_file = "EthNew_kelch13_interaction_results.csv"
  )
)

# ==============================================================================
# ANALYSIS FUNCTION
# ==============================================================================

run_analysis <- function(sc) {

  cat("\n")
  cat("##############################################################\n")
  cat("# Scenario:", sc$label, "\n")
  cat("##############################################################\n\n")

  # ----------------------------------------------------------------------------
  # Load data
  # ----------------------------------------------------------------------------

  cat("Loading genotype matrix...\n")
  geno <- fread(sc$geno_file, sep = "\t", header = TRUE, data.table = FALSE)
  cat("  Dimensions:", nrow(geno), "SNPs x", ncol(geno) - 3, "samples\n\n")

  cat("Loading PC covariates...\n")
  pcs <- fread(sc$pcs_file)
  colnames(pcs) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5")

  # ----------------------------------------------------------------------------
  # Align samples
  # ----------------------------------------------------------------------------

  cat("Aligning samples...\n")
  common_samples <- intersect(colnames(geno)[4:ncol(geno)], pcs$sample)
  cat("  Common samples:", length(common_samples), "\n\n")

  geno <- geno[, c("chr", "pos", "ref", common_samples)]
  pcs  <- pcs[match(common_samples, pcs$sample), ]
  stopifnot(all(pcs$sample == common_samples))

  # ----------------------------------------------------------------------------
  # Extract kelch13 row
  # ----------------------------------------------------------------------------

  kr <- kelch_row
  if (!is.null(kr)) {
    cat("Using user-specified kelch13 row:", kr, "\n")
  } else if (!is.null(kelch_id)) {
    cat("Auto-detecting kelch13 row from identifier:", kelch_id, "\n")
    kr <- which(apply(geno[, c("chr","pos","ref")], 1, function(r) any(grepl(kelch_id, r))))
    if (length(kr) == 0) stop("Could not find kelch_id in genotype matrix. Check identifier.")
    if (length(kr) > 1) {
      warning("Multiple rows match kelch_id — using first match.")
      kr <- kr[1]
    }
    cat("  Found kelch13 at row:", kr, "\n")
  } else {
    stop("Must specify either kelch_row (integer) or kelch_id to identify the kelch13 SNP.")
  }

  kelch_geno <- as.numeric(geno[kr, 4:ncol(geno)])
  cat("  Kelch13 genotype — 0s:", sum(kelch_geno == 0, na.rm = TRUE),
      "| 1s:", sum(kelch_geno == 1, na.rm = TRUE),
      "| NAs:", sum(is.na(kelch_geno)), "\n")

  kelch_mut_samples <- common_samples[which(kelch_geno == 1)]
  cat("  Samples with kelch13 mutation (n =", length(kelch_mut_samples), "):\n")
  cat("   ", paste(kelch_mut_samples, collapse = ", "), "\n\n")

  # ----------------------------------------------------------------------------
  # Filter genotype matrix
  # ----------------------------------------------------------------------------

  geno_test <- geno[-kr, ]
  cat("Testing", nrow(geno_test), "SNPs\n\n")

  # ----------------------------------------------------------------------------
  # Parallel association test
  # ----------------------------------------------------------------------------

  cat("Running parallel association tests on", n_cores, "cores...\n")

  results <- mclapply(seq_len(nrow(geno_test)), function(i) {

    row  <- geno_test[i, ]
    snp2 <- as.numeric(row[1, 4:ncol(row)])

    interaction_term <- kelch_geno * snp2
    if (sum(interaction_term, na.rm = TRUE) <= 1) return(NULL)

    df <- data.frame(
      y    = kelch_geno,
      snp2 = snp2,
      PC1  = pcs$PC1,
      PC2  = pcs$PC2,
      PC3  = pcs$PC3,
      PC4  = pcs$PC4,
      PC5  = pcs$PC5
    )
    df <- df[complete.cases(df), ]

    if (nrow(df) <= 10) return(NULL)

    fit <- tryCatch(
      glm(y ~ snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = df, family = "binomial"),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NULL)

    coef_table <- coef(summary(fit))
    if (!"snp2" %in% rownames(coef_table)) return(NULL)

    return(data.frame(
      chr       = row[["chr"]],
      pos       = row[["pos"]],
      ref       = row[["ref"]],
      beta_snp2 = coef_table["snp2", "Estimate"],
      p_snp2    = coef_table["snp2", "Pr(>|z|)"],
      n         = nrow(df)
    ))

  }, mc.cores = n_cores)

  # ----------------------------------------------------------------------------
  # Combine and save
  # ----------------------------------------------------------------------------

  cat("Combining results...\n")
  results_df <- do.call(rbind, Filter(Negate(is.null), results))

  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No SNPs passed the co-occurrence filter for scenario: ", sc$label)
    return(invisible(NULL))
  }

  cat("  SNPs passing filters:", nrow(results_df), "\n")

  results_df$p_adj <- p.adjust(results_df$p_snp2, method = "BH")
  results_df       <- results_df[order(results_df$p_adj), ]

  write.csv(results_df, sc$output_file, row.names = FALSE)
  cat("Results saved to:", sc$output_file, "\n")

  # ----------------------------------------------------------------------------
  # Annotate results
  # ----------------------------------------------------------------------------

  cat("Adding annotations...\n")
  annotations_final <- read.csv("annotations_final.tsv", sep = "\t")
  annotated_output  <- sub("\\.csv$", "_ann.csv", sc$output_file)
  results_df_ann    <- left_join(results_df, annotations_final,
                                 by = c("chr" = "chr", "pos" = "pos"))
  write.csv(results_df_ann, annotated_output, row.names = FALSE)
  cat("Annotated results saved to:", annotated_output, "\n")

  return(invisible(results_df_ann))
}

# ==============================================================================
# RUN ALL SCENARIOS
# ==============================================================================

for (sc in scenarios) {
  run_analysis(sc)
}

cat("\nAll scenarios complete.\n")

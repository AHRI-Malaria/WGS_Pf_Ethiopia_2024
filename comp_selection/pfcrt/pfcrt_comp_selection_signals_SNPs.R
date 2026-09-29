################################################################################
# GENOMIC SELECTION ANALYSIS INTEGRATION - SNP-BASED
# Combines iHS and XP-EHH analyses at the SNP level
# Groups: pfcrt K76 vs 76T, within the co-endemic sites
#         (Gondar, Mizan, Metehera Fentale; Ethiopia 2024)
# (cf. pfcrt_comp_selection_signals_windows.R for the window-based version)
#
# GROUP DEFINITION: the split is made on pfcrt K76T (Pf3D7_07_v3:403,625),
#       from input_files/K76T_longformat.csv — K76 = 0/0, 76T = 1/1. Named for
#       the codon-76 genotype rather than the CVIET/CVMNK haplotype because
#       M74I and N75E are not variant in this callset; see
#       ehh/pfcrt_genotype_full_pipeline.sh.
################################################################################

library(dplyr)
library(stringr)
library(data.table)

BASE_DIR <- "/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection"
setwd(BASE_DIR)

EHH_DIR  <- "../ehh/pfcrt_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem"
OUT_FILE <- "pfcrt/pfcrt_integrated_selection_analysis_SNPs.tsv"

# ==============================================================================
# PART 1: LOAD iHS DATA
# ==============================================================================

cat("Loading iHS data...\n")
ihs <- read.delim(file.path(EHH_DIR, "high_ihs_all_categories.tsv"))
cat("  iHS SNPs:", nrow(ihs), "\n\n")

# ==============================================================================
# PART 2: LOAD XP-EHH DATA
# ==============================================================================

cat("Loading XP-EHH data...\n")
xpehh <- read.delim(file.path(EHH_DIR, "high_xpehh_all_categories.tsv"))
cat("  XP-EHH SNPs:", nrow(xpehh), "\n\n")

# ==============================================================================
# PART 3: STANDARDIZE FORMATS
# ==============================================================================

cat("Standardizing formats...\n")

IHS <- ihs %>% dplyr::select(CHR, POSITION, IHS, LOGPVALUE, Gene_name, category_name)
colnames(IHS) <- c("Chr", "Position", "IHS", "iHS_logp", "iHS_Gene", "iHS_Region")
setDT(IHS)
IHS <- unique(IHS)
cat("  iHS unique SNPs after dedup:", nrow(IHS), "\n")

XP <- xpehh %>% dplyr::select(CHR, POSITION, XPEHH, LOGPVALUE, Gene_name, category_name)
colnames(XP) <- c("Chr", "Position", "XPEHH", "XPEHH_logp", "XPEHH_Gene", "XPEHH_Comparison")
setDT(XP)
XP[, c("XPEHH_pop1", "XPEHH_pop2") := tstrsplit(XPEHH_Comparison, "|", fixed = TRUE)]
XP[, XPEHH_pop1 := trimws(XPEHH_pop1)]
XP[, XPEHH_pop2 := trimws(XPEHH_pop2)]
XP <- unique(XP)
cat("  XP-EHH unique SNPs after dedup:", nrow(XP), "\n")

cat("  Done.\n\n")

# ==============================================================================
# PART 4: MERGE PER CHROMOSOME
# ==============================================================================

cat("Merging datasets...\n\n")

all_chrs <- sort(unique(c(IHS$Chr, XP$Chr)))
all_chrs <- all_chrs[!is.na(all_chrs)]

results_list <- list()

for (chr in all_chrs) {
  cat("  Chr", chr, "...\n")

  ihs_chr <- IHS[Chr == chr]
  xp_chr  <- XP[Chr == chr]

  ihs_chr[, .ihs_id := .I]
  xp_chr[,  .xp_id  := .I]

  cat("    iHS:", nrow(ihs_chr), "| XP-EHH:", nrow(xp_chr), "\n")

  result <- data.table()

  if (nrow(ihs_chr) > 0 && nrow(xp_chr) > 0) {
    ov_xp <- merge(ihs_chr, xp_chr, by = c("Chr", "Position"),
                   all = FALSE, allow.cartesian = TRUE)

    ov_xp <- ov_xp[
      !is.na(iHS_Region) & iHS_Region == XPEHH_pop1
    ]

    ihs_matched <- unique(ov_xp$.ihs_id)
    xp_matched  <- unique(ov_xp$.xp_id)

    ihs_only <- ihs_chr[!.ihs_id %in% ihs_matched]
    xp_only  <- xp_chr[!.xp_id  %in% xp_matched]

    result <- rbind(ov_xp, ihs_only, fill = TRUE)
    if (nrow(xp_only) > 0) result <- rbind(result, xp_only, fill = TRUE)

  } else if (nrow(ihs_chr) > 0) {
    result <- ihs_chr
  } else if (nrow(xp_chr) > 0) {
    result <- xp_chr
  }

  cat("    After iHS+XP-EHH:", nrow(result), "rows\n")

  temp_cols <- grep("^\\.", names(result), value = TRUE)
  if (length(temp_cols) > 0) result[, (temp_cols) := NULL]

  results_list[[as.character(chr)]] <- result

  rm(ihs_chr, xp_chr, result)
  gc(verbose = FALSE)
}

# ==============================================================================
# PART 5: COMBINE AND FINALIZE
# ==============================================================================

cat("\nCombining all chromosomes...\n")
result <- rbindlist(results_list, fill = TRUE, use.names = TRUE)
setorder(result, Chr, Position, na.last = TRUE)

# ==============================================================================
# PART 6: SUMMARY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("INTEGRATION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Final dataset: 'result'\n")
cat("  Total rows:", nrow(result), "\n")
cat("  Chromosomes:", paste(sort(unique(result$Chr)), collapse = ", "), "\n\n")

result[, n_datasets := (!is.na(iHS_Region)) + (!is.na(XPEHH_Comparison))]

cat("Dataset overlap distribution:\n")
print(table(result$n_datasets))
cat("\n")

result[, combo := paste(
  ifelse(!is.na(iHS_Region),       "iHS", ""),
  ifelse(!is.na(XPEHH_Comparison), "XP",  ""),
  sep = "+"
)]
result[, combo := gsub("\\+{2,}", "+", combo)]
result[, combo := gsub("^\\+|\\+$", "", combo)]

cat("Dataset combinations (top 10):\n")
print(head(sort(table(result$combo), decreasing = TRUE), 10))

cat("\nFirst 20 rows:\n")
print(head(result, 20))

cat("\n", rep("=", 80), "\n", sep = "")
cat("Result stored in 'result' data.table\n")
cat(rep("=", 80), "\n", sep = "")

# ==============================================================================
# ANNOTATION CHECK
# ==============================================================================

gene_cols <- grep("_Gene$", names(result), value = TRUE)
gene_strings <- unlist(result[, ..gene_cols])
gene_strings <- gene_strings[!is.na(gene_strings) & gene_strings != ""]

px1_hits <- grep("PF3D7_0720700", gene_strings, value = TRUE)
if (length(px1_hits) > 0) {
  for (col in gene_cols) {
    result[get(col) == "PF3D7_0720700", (col) := "PX1"]
  }
  cat("CHECK PASSED: PF3D7_0720700 found and relabelled as PX1 in", length(px1_hits), "row(s)\n")
} else {
  cat("CHECK NOTE: PF3D7_0720700 not present in this result (gene may not be in a significant SNP)\n")
}

setwd(BASE_DIR)
fwrite(result, OUT_FILE, sep = "\t")
cat("Saved:", OUT_FILE, "\n")


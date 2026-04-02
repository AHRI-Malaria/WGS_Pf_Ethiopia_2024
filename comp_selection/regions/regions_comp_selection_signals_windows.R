################################################################################
# GENOMIC SELECTION ANALYSIS INTEGRATION - WINDOW-BASED
# Combines Tajima's D, iHS, XP-EHH, and IBD analyses
# Groups: Horn_of_Africa, East_Africa, West_Africa
################################################################################

library(dplyr)
library(stringr)
library(data.table)

BASE_DIR <- "/Users/lsh2004630/Library/CloudStorage/GoogleDrive-leenvanheer@gmail.com/My Drive/Research_Fellow/Popgen/comp_selection"
setwd(BASE_DIR)

GROUPS <- c("Horn_of_Africa", "East_Africa", "West_Africa")
EHH_DIR <- "../ehh/regions_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem"
OUT_FILE <- "regions/regions_integrated_selection_analysis_windows.tsv"

# ==============================================================================
# PART 1: LOAD AND PREPARE TAJIMA'S D DATA
# ==============================================================================

cat("Loading Tajima's D data...\n")

tajima_d <- read.csv("../tajimasD/tajimasD_all_annotated_merged.tsv", sep = "\t")
tajima_d <- tajima_d %>% filter(!chr %in% c("Pf_M76611", "Pf3D7_API_v3"))

percentiles <- tajima_d %>%
  group_by(group) %>%
  summarise(p_low  = quantile(TajimaD, 0.15, na.rm = TRUE),
            p_high = quantile(TajimaD, 0.85, na.rm = TRUE),
            .groups = "drop")

tajima_d2    <- tajima_d %>% left_join(percentiles, by = "group")
extremes_df  <- tajima_d2 %>% filter(!is.na(TajimaD) & TajimaD <= p_low)

extremes_df_groups <- extremes_df[extremes_df$group %in% GROUPS, ]

cat("  Extreme negative Tajima's D windows:", nrow(extremes_df_groups), "\n\n")

# ==============================================================================
# PART 2: LOAD SELECTION DATA
# ==============================================================================

cat("Loading EHH data...\n")
cr_ihs   <- read.delim(file.path(EHH_DIR, "cr_ihs_all_categories_annot.tsv"))
cr_xpehh <- read.delim(file.path(EHH_DIR, "cr_xpehh_all_categories_annot.tsv"))
cat("  iHS windows:    ", nrow(cr_ihs), "\n")
cat("  XP-EHH windows:", nrow(cr_xpehh), "\n\n")

# ==============================================================================
# PART 3: LOAD IBD DATA
# ==============================================================================

cat("Loading IBD data...\n")
ibd <- read.delim("../hmmIBD/hmmibd_combined_win10kb_annotated.tsv")
ibd <- ibd[ibd$category %in% GROUPS, ]

percentiles_ibd <- ibd %>%
  group_by(category) %>%
  summarise(p_low  = quantile(fraction, 0.15, na.rm = TRUE),
            p_high = quantile(fraction, 0.85, na.rm = TRUE),
            .groups = "drop")

ibd <- ibd %>%
  left_join(percentiles_ibd, by = "category") %>%
  filter(!is.na(fraction), fraction >= p_high)

cat("  High IBD windows:", nrow(ibd), "\n\n")

# ==============================================================================
# PART 4: STANDARDIZE FORMATS
# ==============================================================================

cat("Standardizing formats...\n")

# Tajima's D
TD <- extremes_df_groups %>%
  dplyr::select(chr, start, end, N_SNPS, TajimaD, genes, products, group)
colnames(TD) <- c("Chr", "TD_start", "TD_end", "TD_SNP_n", "TD",
                  "TD_genes", "TD_products", "TD_Region")
TD$Chr <- gsub("Pf3D7_0|Pf3D7_|_v3", "", TD$Chr)
TD$Chr <- as.numeric(as.character(TD$Chr))

# iHS
cr_ihs <- cr_ihs %>% dplyr::select(chr, start, end, category_name, genes, products)
colnames(cr_ihs) <- c("Chr", "iHS_start", "iHS_end", "iHS_Region", "iHS_genes", "iHS_products")

# IBD
ibd <- ibd %>% dplyr::select(chr, start, end, fraction, category, gene_name, gene_product)
colnames(ibd) <- c("Chr", "IBD_start", "IBD_end", "IBD_fraction", "IBD_Region",
                   "IBD_genes", "IBD_products")

# XP-EHH
cr_xpehh <- cr_xpehh %>%
  dplyr::select(chr, start, end, category_name, PERC_EXTR_MRK, genes, products)
colnames(cr_xpehh) <- c("Chr", "XPEHH_start", "XPEHH_end", "XPEHH_Comparison",
                        "XPEHH_extra_marks_p", "XPEHH_genes", "XPEHH_products")

setDT(TD); setDT(cr_ihs); setDT(ibd); setDT(cr_xpehh)
TD <- unique(TD); cr_ihs <- unique(cr_ihs)
ibd <- unique(ibd); cr_xpehh <- unique(cr_xpehh)

cat("  Done.\n\n")

# ==============================================================================
# PART 5: SPLIT XP-EHH INTO TWO POPULATIONS
# ==============================================================================

cat("Splitting XP-EHH comparisons into population pairs...\n")
cr_xpehh[, c("XPEHH_pop1", "XPEHH_pop2") := tstrsplit(XPEHH_Comparison, "|", fixed = TRUE)]
cr_xpehh[, XPEHH_pop1 := trimws(XPEHH_pop1)]
cr_xpehh[, XPEHH_pop2 := trimws(XPEHH_pop2)]
cat("  Split", nrow(cr_xpehh), "XP-EHH regions into pop1 and pop2\n\n")

# ==============================================================================
# PART 6: MERGE STEP-BY-STEP PER CHROMOSOME
# ==============================================================================

cat("Merging datasets...\n\n")

all_chrs <- sort(unique(c(TD$Chr, cr_ihs$Chr, ibd$Chr, cr_xpehh$Chr)))
all_chrs <- all_chrs[!is.na(all_chrs)]

results_list <- list()

for (chr in all_chrs) {
  cat("  Chr", chr, "...\n")

  TD_chr  <- TD[Chr == chr]
  ihs_chr <- cr_ihs[Chr == chr]
  ibd_chr <- ibd[Chr == chr]
  xp_chr  <- cr_xpehh[Chr == chr]

  TD_chr[,  .td_id  := .I]
  ihs_chr[, .ihs_id := .I]
  ibd_chr[, .ibd_id := .I]
  xp_chr[,  .xp_id  := .I]

  cat("    TD:", nrow(TD_chr), "| iHS:", nrow(ihs_chr),
      "| IBD:", nrow(ibd_chr), "| XP-EHH:", nrow(xp_chr), "\n")

  # === STEP 1: TD + iHS (same region, overlapping intervals) ===
  result      <- data.table()
  td_matched  <- integer(0)
  ihs_matched <- integer(0)

  if (nrow(TD_chr) > 0 && nrow(ihs_chr) > 0) {
    setkeyv(TD_chr,  c("Chr", "TD_start",  "TD_end"))
    setkeyv(ihs_chr, c("Chr", "iHS_start", "iHS_end"))

    ov <- foverlaps(TD_chr, ihs_chr, type = "any", nomatch = NULL)
    ov <- ov[TD_start < iHS_end & iHS_start < TD_end]
    ov <- ov[TD_Region == iHS_Region]

    if (nrow(ov) > 0) {
      result      <- rbind(result, ov, fill = TRUE)
      td_matched  <- ov$.td_id
      ihs_matched <- ov$.ihs_id
    }
  }

  if (nrow(TD_chr)  > 0) {
    td_only <- TD_chr[!.td_id %in% td_matched]
    if (nrow(td_only) > 0) result <- rbind(result, td_only, fill = TRUE)
  }
  if (nrow(ihs_chr) > 0) {
    ihs_only <- ihs_chr[!.ihs_id %in% ihs_matched]
    if (nrow(ihs_only) > 0) result <- rbind(result, ihs_only, fill = TRUE)
  }

  cat("    After TD+iHS:", nrow(result), "rows\n")

  # === STEP 2: + IBD ===
  if (nrow(result) > 0 && nrow(ibd_chr) > 0) {
    start_cols2 <- intersect(c("TD_start", "iHS_start"), names(result))
    end_cols2   <- intersect(c("TD_end",   "iHS_end"),   names(result))
    result[, unified_start := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = start_cols2]
    result[, unified_end   := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = end_cols2]
    result[, .res_id := .I]

    res_valid <- result[!is.na(unified_start) & !is.na(unified_end)]
    res_na    <- result[is.na(unified_start)  | is.na(unified_end)]

    if (nrow(res_valid) > 0) {
      setkeyv(res_valid, c("Chr", "unified_start", "unified_end"))
      setkeyv(ibd_chr,  c("Chr", "IBD_start",     "IBD_end"))

      ov_ibd <- foverlaps(res_valid, ibd_chr, type = "any", nomatch = NULL)
      ov_ibd <- ov_ibd[unified_start < IBD_end & IBD_start < unified_end]
      for (col in c("TD_Region", "iHS_Region")) {
        if (!col %in% names(ov_ibd)) ov_ibd[, (col) := NA_character_]
      }
      ov_ibd <- ov_ibd[
        (!is.na(TD_Region)  & TD_Region  == IBD_Region) |
        (!is.na(iHS_Region) & iHS_Region == IBD_Region)
      ]

      res_matched <- ov_ibd[!is.na(.res_id), unique(.res_id)]
      ibd_matched <- ov_ibd[!is.na(.ibd_id), unique(.ibd_id)]

      if (nrow(ov_ibd) > 0) {
        ov_ibd[, .match_pop := fcoalesce(TD_Region, iHS_Region)]
        first_nonNA <- function(x) { nna <- x[!is.na(x)]; if (length(nna)) nna[1] else x[NA_integer_] }
        ov_ibd <- ov_ibd[, lapply(.SD, first_nonNA), by = .(.ibd_id, .match_pop)]
        ov_ibd[, .match_pop := NULL]
      }

      res_unmatched <- res_valid[!.res_id %in% res_matched]
      result <- rbind(ov_ibd, res_unmatched, res_na, fill = TRUE)

      ibd_only <- ibd_chr[!.ibd_id %in% ibd_matched]
      if (nrow(ibd_only) > 0) result <- rbind(result, ibd_only, fill = TRUE)
    } else {
      result <- rbind(result, ibd_chr, fill = TRUE)
    }
    result[, c("unified_start", "unified_end", ".res_id") := NULL]
  } else if (nrow(ibd_chr) > 0) {
    result <- rbind(result, ibd_chr, fill = TRUE)
  }

  cat("    After +IBD:", nrow(result), "rows\n")

  # === STEP 3: + XP-EHH ===
  if (nrow(result) > 0 && nrow(xp_chr) > 0) {
    start_cols <- intersect(c("TD_start", "iHS_start", "IBD_start"), names(result))
    end_cols   <- intersect(c("TD_end",   "iHS_end",   "IBD_end"),   names(result))

    result[, unified_start := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = start_cols]
    result[, unified_end   := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = end_cols]
    result[, .res_id := .I]

    res_valid <- result[!is.na(unified_start) & !is.na(unified_end)]
    res_na    <- result[is.na(unified_start)  | is.na(unified_end)]

    if (nrow(res_valid) > 0) {
      setkeyv(res_valid, c("Chr", "unified_start", "unified_end"))
      setkeyv(xp_chr,   c("Chr", "XPEHH_start",   "XPEHH_end"))

      ov_xp <- foverlaps(res_valid, xp_chr, type = "any", nomatch = NULL)
      ov_xp <- ov_xp[unified_start < XPEHH_end & XPEHH_start < unified_end]

      for (col in c("TD_Region", "iHS_Region", "IBD_Region")) {
        if (!col %in% names(ov_xp)) ov_xp[, (col) := NA_character_]
      }

      ov_xp <- ov_xp[
        (!is.na(TD_Region)  & TD_Region  == XPEHH_pop1) |
          (!is.na(iHS_Region) & iHS_Region == XPEHH_pop1) |
          (!is.na(IBD_Region) & IBD_Region == XPEHH_pop1)
      ]

      res_matched <- ov_xp[!is.na(.res_id), unique(.res_id)]
      xp_matched  <- ov_xp[!is.na(.xp_id),  unique(.xp_id)]

      if (nrow(ov_xp) > 0) {
        ov_xp[, .match_pop := fcoalesce(TD_Region, iHS_Region, IBD_Region)]
        first_nonNA <- function(x) { nna <- x[!is.na(x)]; if (length(nna)) nna[1] else x[NA_integer_] }
        ov_xp <- ov_xp[, lapply(.SD, first_nonNA), by = .(.xp_id, .match_pop)]
        ov_xp[, .match_pop := NULL]
      }

      res_unmatched <- res_valid[!.res_id %in% res_matched]
      result <- rbind(ov_xp, res_unmatched, res_na, fill = TRUE)

      xp_only <- xp_chr[!.xp_id %in% xp_matched]
      if (nrow(xp_only) > 0) result <- rbind(result, xp_only, fill = TRUE)
    } else {
      result <- rbind(result, xp_chr, fill = TRUE)
    }
    result[, c("unified_start", "unified_end", ".res_id") := NULL]
  } else if (nrow(xp_chr) > 0) {
    result <- rbind(result, xp_chr, fill = TRUE)
  }

  cat("    After +XP-EHH:", nrow(result), "rows\n")

  temp_cols <- grep("^\\.|^i\\.", names(result), value = TRUE)
  if (length(temp_cols) > 0) result[, (temp_cols) := NULL]

  result <- unique(result)
  cat("    After dedup:", nrow(result), "rows\n")

  results_list[[as.character(chr)]] <- result

  rm(TD_chr, ihs_chr, ibd_chr, xp_chr, result)
  gc(verbose = FALSE)
}

# ==============================================================================
# PART 7: COMBINE AND FINALIZE
# ==============================================================================

cat("\nCombining all chromosomes...\n")
result <- rbindlist(results_list, fill = TRUE, use.names = TRUE)
result <- unique(result)
cat("  Rows after global dedup:", nrow(result), "\n")

result[, `:=`(
  combined_start = pmin(TD_start, iHS_start, IBD_start, XPEHH_start, na.rm = TRUE),
  combined_end   = pmax(TD_end,   iHS_end,   IBD_end,   XPEHH_end,   na.rm = TRUE)
)]

setorder(result, Chr, combined_start, na.last = TRUE)

# ==============================================================================
# PART 8: SUMMARY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("INTEGRATION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Final dataset: 'result'\n")
cat("  Total rows:", nrow(result), "\n")
cat("  Chromosomes:", paste(sort(unique(result$Chr)), collapse = ", "), "\n\n")

result[, n_datasets := (!is.na(TD_Region)) + (!is.na(iHS_Region)) +
         (!is.na(IBD_Region)) + (!is.na(XPEHH_Comparison))]

cat("Dataset overlap distribution:\n")
print(table(result$n_datasets))
cat("\n")

result[, combo := paste(
  ifelse(!is.na(TD_Region),         "TD",  ""),
  ifelse(!is.na(iHS_Region),        "iHS", ""),
  ifelse(!is.na(IBD_Region),        "IBD", ""),
  ifelse(!is.na(XPEHH_Comparison),  "XP",  ""),
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

result[!is.na(TD_Region), .(N = .N,
                             min_TD  = min(TD,  na.rm = TRUE),
                             max_TD  = max(TD,  na.rm = TRUE),
                             mean_TD = mean(TD, na.rm = TRUE)),
       by = TD_Region][order(TD_Region)]

# ==============================================================================
# ANNOTATION CHECK
# ==============================================================================

gene_cols <- grep("_genes$", names(result), value = TRUE)
gene_strings <- unlist(result[, ..gene_cols])
gene_strings <- gene_strings[!is.na(gene_strings) & gene_strings != ""]

px1_hits <- grep("PF3D7_0720700", gene_strings, value = TRUE)
if (length(px1_hits) > 0) {
  for (col in gene_cols) {
    result[grepl("PF3D7_0720700", get(col)),
           (col) := gsub("PF3D7_0720700\\([^)]*\\)", "PF3D7_0720700(PX1)", get(col))]
  }
  cat("CHECK PASSED: PF3D7_0720700 found and relabelled as PX1 in", length(px1_hits), "instance(s)\n")
} else {
  # Manually annotate windows overlapping PX1 gene (Pf3D7_07_v3:889883..899213)
  px1_chr <- 7
  px1_start <- 889883
  px1_end   <- 899213
  px1_gene  <- "PF3D7_0720700(PX1)"
  px1_prod  <- "phosphoinositide-binding protein"

  px1_rows <- result[Chr == px1_chr & combined_start < px1_end & combined_end > px1_start]

  if (nrow(px1_rows) > 0) {
    prod_cols <- grep("_products$", names(result), value = TRUE)
    idx <- result[, which(Chr == px1_chr & combined_start < px1_end & combined_end > px1_start)]

    for (i in seq_along(gene_cols)) {
      gcol <- gene_cols[i]
      pcol <- prod_cols[i]
      for (r in idx) {
        cur_gene <- result[[gcol]][r]
        if (is.na(cur_gene) || cur_gene == "") {
          set(result, r, gcol, px1_gene)
          set(result, r, pcol, px1_prod)
        } else if (!grepl("PF3D7_0720700", cur_gene)) {
          set(result, r, gcol, paste0(cur_gene, ",", px1_gene))
          set(result, r, pcol, paste0(result[[pcol]][r], ",", px1_prod))
        }
      }
    }
    cat("CHECK PASSED: PX1 manually annotated in", length(idx), "overlapping window(s)\n")
  } else {
    cat("CHECK NOTE: No windows on chr 7 overlap PX1 coordinates (889883-899213)\n")
  }
}

setwd(BASE_DIR)
fwrite(result, OUT_FILE, sep = "\t")
cat("Saved:", OUT_FILE, "\n")

library(dplyr)
library(stringr)
library(readr)

# ==============================================================================
# SETTINGS
# Usage: cd <data_dir> && Rscript /path/to/repo/tajimasD/annotate.R
#
# data_dir: directory containing tajimasD_by_* folders (default: working directory)
# ==============================================================================
data_dir      <- getwd()
gene_product_file <- "/mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv"
th_quantile   <- 0.1                      # annotate bottom 10%
rm_chr        <- "Pf_M76611,Pf3D7_API_v3"

# Directories to process: input folder -> output folder (relative to data_dir)
analyses <- list(
  tajimasD_by_country    = "tajimasD_annotated_by_country",
  tajimasD_by_region     = "tajimasD_annotated_by_region",
  tajimasD_old_vs_new    = "tajimasD_annotated_old_vs_new",
  tajimasD_north_vs_rest = "tajimasD_annotated_north_vs_rest",
  tajimasD_r622i         = "tajimasD_annotated_r622i",
  tajimasD_hrp_status    = "tajimasD_annotated_hrp_status",
  tajimasD_by_site       = "tajimasD_annotated_by_site",
  tajimasD_sites_grouping = "tajimasD_annotated_sites_grouping",
  tajimasD_pfcrt          = "tajimasD_annotated_pfcrt"
)

# ==============================================================================
# LOAD ANNOTATION
# ==============================================================================
annotation <- read_tsv(gene_product_file, col_types = cols())

rm_chr_ <- strsplit(rm_chr, ",")[[1]]
annotation <- annotation %>% filter(!chr %in% rm_chr_)
annotation$chr <- as.numeric(str_match(annotation$chr, "Pf3D7_(\\d+)_v3")[, 2])

# ==============================================================================
# ANNOTATION FUNCTION
# ==============================================================================
annotate_regions <- function(regions, annotation) {
  result <- regions %>%
    left_join(annotation, by = "chr", relationship = "many-to-many") %>%
    filter(pos_start <= end & pos_end >= start) %>%
    group_by(chr, start, end) %>%
    summarise(
      genes    = paste0(unique(na.omit(paste0(gene_id, "(", gene_name, ")"))), collapse = "; "),
      products = paste0(unique(na.omit(gene_product)), collapse = "; "),
      .groups  = "drop"
    )

  regions %>%
    dplyr::select(chr, start, end) %>%
    left_join(result, by = c("chr", "start", "end"))
}

# ==============================================================================
# PROCESS EACH ANALYSIS
# ==============================================================================
for (tajima_dir in names(analyses)) {
  outdir <- file.path(data_dir, analyses[[tajima_dir]])
  tajima_dir_full <- file.path(data_dir, tajima_dir)

  if (!dir.exists(tajima_dir_full)) {
    message("Skipping ", tajima_dir, " (directory not found)")
    next
  }

  dir.create(outdir, showWarnings = FALSE)
  files <- list.files(tajima_dir_full, pattern = "\\.Tajima\\.D$", full.names = TRUE)

  if (length(files) == 0) {
    message("No .Tajima.D files found in ", tajima_dir, ", skipping.")
    next
  }

  message("\n=== ", tajima_dir, " ===")

  for (f in files) {
    group <- gsub("tajimasd_|_10kb\\.Tajima\\.D", "", basename(f))
    message("Processing: ", group)

    tajima <- read_tsv(f, col_types = cols())
    tajima <- tajima %>% filter(!CHROM %in% rm_chr_)
    tajima$chr   <- as.numeric(str_match(tajima$CHROM, "Pf3D7_(\\d+)_v3")[, 2])
    tajima$start <- tajima$BIN_START
    tajima$end   <- tajima$BIN_START + 9999

    tajima <- tajima %>%
      mutate(q_low = quantile(TajimaD, th_quantile, na.rm = TRUE)) %>%
      filter(!is.na(TajimaD))

    extreme <- tajima %>% filter(TajimaD <= q_low)
    message("  Extreme windows: ", nrow(extreme))

    # Annotate all windows
    all_regions   <- tajima %>% dplyr::select(chr, start, end) %>% distinct()
    all_annotated <- annotate_regions(all_regions, annotation)

    result_all <- tajima %>%
      left_join(all_annotated, by = c("chr", "start", "end")) %>%
      dplyr::select(chr, start, end, N_SNPS, TajimaD, genes, products) %>%
      mutate(group = group)

    write.table(result_all,
      file.path(outdir, paste0("tajimasd_", group, "_10kb_annotated.tsv")),
      sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    message("  Written: tajimasd_", group, "_10kb_annotated.tsv")

    # Annotate filtered (bottom quantile) windows
    extreme_regions   <- extreme %>% dplyr::select(chr, start, end) %>% distinct()
    extreme_annotated <- annotate_regions(extreme_regions, annotation)

    result_extreme <- extreme %>%
      left_join(extreme_annotated, by = c("chr", "start", "end")) %>%
      dplyr::select(chr, start, end, N_SNPS, TajimaD, genes, products) %>%
      mutate(group = group)

    write.table(result_extreme,
      file.path(outdir, paste0("tajimasd_", group, "_10kb_annotated_q", th_quantile, ".tsv")),
      sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    message("  Written: tajimasd_", group, "_10kb_annotated_q", th_quantile, ".tsv")
  }
}

message("\nAll done.")

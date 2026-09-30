# Chromosome 7 (PF3D7_0720700 / PX1) EHH / furcation plots, new (2024) Ethiopian samples only.
#
# WHY THIS SCRIPT EXISTS IN THIS FORM
# -----------------------------------
# Three of the four PX1 codons of interest are NOT in the WGS-derived rehh input:
#   7_893360  (D384A)   -> ABSENT  (dropped by the biallelic filter, see below)
#   7_895874  (L1222P)  -> present in snp.info.inp.new
#   7_897312  (M1701I)  -> ABSENT  (never called in the filtered VCF)
#   7_897322  (D1705N)  -> ABSENT  (never called in the filtered VCF)
# 897312/897322 are lost upstream of the hypervariable/Tajima's-D filters (the
# region is indel-rich and short-read calling drops it), so the previous version
# of this script silently skipped 2 of its 3 SNPs via tryCatch.
#
# SeekDeep genotyped 1222/1701/1705 successfully from the amplicon
# Pf3D7_07_v3-897170-897388. This script therefore splices the SeekDeep calls
# into the chr7 haplotype matrix as extra markers before running rehh.
#
# D384A is a different problem. It is NOT in the SeekDeep amplicon, and WGS
# genotypes it fine - but Pf3D7_07_v3:893360 is multiallelic (A>C,G; C = D384A,
# G = D384G), so the ".bi." step of the pipeline that builds the binary matrix
# throws the whole site away before rehh ever sees it. Those calls are therefore
# pulled back out of the pre-biallelic VCF with bcftools (multiallelic_vcf
# below) and spliced in the same way as the SeekDeep ones. In the 2024 ("new")
# Ethiopian set 9/120 samples are C/C; the fws >= 0.95 monoclonal filter leaves
# 6 of those in the haplotype matrix, which is what the EHH curve is built from.
#
# Usage:
#   Rscript ehh_plot_PF3D7_0720700.R            # defaults to category "new"
#   Rscript ehh_plot_PF3D7_0720700.R old        # or any other category prefix

library(rehh)

# ---------------------------------------------------------------------------
# Paths / configuration
# ---------------------------------------------------------------------------
proj_dir <- "/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen"
base_dir <- file.path(proj_dir, "ehh/old_vs_new_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem")
sd_dir   <- file.path(proj_dir, "px1_genotype_result_metmiku")

bin_mat_file  <- file.path(base_dir, "012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.mat.bin")
metadata_file <- file.path(base_dir, "metadata_temporal_group.csv")   # tab-separated
sd_geno_file  <- file.path(sd_dir, "px1_genotype_result.tsv")

# The VCF one step BEFORE the ".bi." biallelic split, i.e. the same sample set
# and the same miss0.2/VQSLOD/SNP filtering as the binary matrix, but with the
# multiallelic sites still intact. This is the only place D384A survives.
multiallelic_vcf <- file.path(proj_dir, "012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.GT.miss0.2.vqslod.filt.snps.vcf.gz")
bcftools_bin <- Sys.which("bcftools")
if (!nzchar(bcftools_bin)) bcftools_bin <- "/mnt/storage13/leen/miniforge3/bin/bcftools"

# these must match what full_pipeline.sh passed to prepare_input_rehh_per_category.R
label_id       <- "wgs_id"
label_category <- "temporal_group"
label_fws      <- "fws"
fws_th         <- 0.95

chr <- 7
vcf_chrom <- sprintf("Pf3D7_%02d_v3", chr)   # contig name as spelled in the VCF

# EHH x-axis. NULL = auto: widen until both the ancestral and the derived curve
# have decayed to zero, so the point where each hits 0 is actually on the plot.
# Set to e.g. c(0.80E+6, 1.00E+6) to go back to a fixed window.
ehh_xlim <- NULL
# Furcation trees get unreadable over a wide span, so that window stays fixed.
furcation_xlim <- c(0.85E+6, 0.95E+6)

args <- commandArgs(trailingOnly = TRUE)
category <- if (length(args) >= 1) args[1] else "new"
country  <- ifelse(category == "new", "Ethiopia_new", category)

hap_file <- file.path(base_dir, sprintf("hap_chr%d_%s", chr, category))
map_file <- file.path(base_dir, sprintf("snp.info.inp.%s", category))

# augmented copies are written next to this script, the pipeline outputs are untouched
out_hap <- sprintf("hap_chr%d_%s.px1", chr, category)
out_map <- sprintf("snp.info.inp.%s.px1", category)

# ---------------------------------------------------------------------------
# The four PX1 markers.
# ref/alt verified against Pfalciparum.genome.fasta and against the bcftools csq
# annotation of the filtered VCF; PF3D7_0720700 is on the plus strand and these
# codons are intron-free, so:
#   893359-893361 gat (D384)  -> gCt  A   A>C   (893359G>T = D384Y, 893360A>G = D384G)
#   895873-895875 ctt (L1222) -> cCt  P   T>C
#   897310-897312 atg (M1701) -> atA  I   G>A
#   897322-897324 gat (D1705) -> Aat  N   G>A
# source  = where the genotypes come from, "seekdeep" or "vcf".
# sd_col  = column name in px1_genotype_result.tsv holding that codon call
#           (NA for markers that the SeekDeep amplicon does not cover).
# ---------------------------------------------------------------------------
px1 <- data.frame(
  pos      = c(893360, 895874, 897312, 897322),
  ref      = c("A", "T", "G", "G"),
  alt      = c("C", "C", "A", "A"),
  source   = c("vcf", "seekdeep", "seekdeep", "seekdeep"),
  sd_col   = c(NA, "1222", "1701", "1705"),
  mutation = c("D384A", "L1222P", "M1701I", "D1705N"),
  stringsAsFactors = FALSE
)
px1$name <- sprintf("%d_%d", chr, px1$pos)

# ---------------------------------------------------------------------------
# 1. Rebuild the sample order of the rows of the hap file.
#
# prepare_input_rehh_per_category.R built the matrix as
#   fread(bin_mat_file, select = c("chr","pos","ref", <metadata ids>))   # file order
#   ... %>% select(<ids with fws >= fws_th>)  %>% t()
# so row i of hap_chr<chr>_<category> is the i-th sample of the binary matrix
# header that is (a) in the metadata for this category and (b) fws >= threshold.
# ---------------------------------------------------------------------------
mat_samples <- scan(bin_mat_file, what = "", nlines = 1, sep = "\t", quiet = TRUE)
mat_samples <- mat_samples[-(1:3)]                       # drop chr, pos, ref

metadata <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c(label_id, label_category, label_fws) %in% colnames(metadata)))

keep <- metadata[[label_category]] == category &
        !is.na(metadata[[label_fws]]) &
        metadata[[label_fws]] >= fws_th
eligible <- metadata[[label_id]][keep]

hap_samples <- mat_samples[mat_samples %in% eligible]    # binary-matrix column order

# ---------------------------------------------------------------------------
# 2. Read the existing rehh inputs
# ---------------------------------------------------------------------------
hap <- read.delim(hap_file, header = FALSE, stringsAsFactors = FALSE)
hap_alleles <- as.matrix(hap[, -1])                      # col 1 is the haplotype id
storage.mode(hap_alleles) <- "integer"

map <- read.table(map_file, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("name", "chr", "pos", "anc", "der"))
map_chr <- map[map$chr == chr, ]

if (nrow(hap_alleles) != length(hap_samples)) {
  stop(sprintf("Reconstructed %d sample ids but hap file has %d rows - the sample order could not be recovered.",
               length(hap_samples), nrow(hap_alleles)))
}
if (ncol(hap_alleles) != nrow(map_chr)) {
  stop(sprintf("hap file has %d markers but the map has %d for chr %d.",
               ncol(hap_alleles), nrow(map_chr), chr))
}

# ---------------------------------------------------------------------------
# 3. External calls -> rehh allele coding, for exactly these samples in row order
#    rehh "12" coding: 1 = ancestral/reference, 2 = derived/mutant, 0 = missing
#
#    This matches prepare_input_rehh_per_category.R, which (with neither
#    --forced_recode nor --forced_mixed) maps ref -> 1, alt -> 2 and both mixed
#    and uncalled -> 0.
# ---------------------------------------------------------------------------
sd <- read.delim(sd_geno_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(sd) <- sd$Sample

sd_code <- function(sd_col, mutation) {
  raw <- sd[hap_samples, sd_col]                         # NA if sample absent
  out <- rep(0L, length(raw))                            # default: missing
  out[!is.na(raw) & raw == "WT"]      <- 1L
  out[!is.na(raw) & raw == mutation]  <- 2L
  out
}

# %TGT prints the genotype as bases ("A/A", "C/C", "./."), so there is no need
# to track which ALT index this site's mutation happens to be at.
vcf_genotypes <- function(pos) {
  cache <- sprintf("vcf_gt_%s_%d.tsv", vcf_chrom, pos)
  if (!file.exists(cache)) {
    message("Extracting ", vcf_chrom, ":", pos, " from ", basename(multiallelic_vcf), " ...")
    status <- system2(bcftools_bin,
                      c("query",
                        "-f", shQuote("[%SAMPLE\\t%TGT\\n]"),
                        "-r", shQuote(sprintf("%s:%d-%d", vcf_chrom, pos, pos)),
                        shQuote(multiallelic_vcf)),
                      stdout = cache)
    if (status != 0 || file.size(cache) == 0) {
      unlink(cache)
      stop("bcftools failed to return genotypes at ", vcf_chrom, ":", pos)
    }
  }
  gt <- read.delim(cache, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("sample", "tgt"))
  setNames(gt$tgt, gt$sample)
}

# Anything that is not a clean homozygous ref or homozygous alt call - missing,
# mixed, or the third allele at a multiallelic site (D384G at 893360) - is
# missing as far as a polarized EHH on this mutation is concerned.
vcf_code <- function(pos, ref, alt) {
  tgt <- vcf_genotypes(pos)[hap_samples]                 # NA if sample absent
  out <- rep(0L, length(tgt))                            # default: missing
  out[!is.na(tgt) & tgt == paste(ref, ref, sep = "/")] <- 1L
  out[!is.na(tgt) & tgt == paste(alt, alt, sep = "/")] <- 2L
  out
}

px1_calls <- lapply(seq_len(nrow(px1)), function(i) {
  if (px1$source[i] == "seekdeep") {
    sd_code(px1$sd_col[i], px1$mutation[i])
  } else {
    vcf_code(px1$pos[i], px1$ref[i], px1$alt[i])
  }
})
names(px1_calls) <- px1$name

n_missing <- vapply(px1_calls, function(x) sum(x == 0L), integer(1))
message("Calls recovered for ", length(hap_samples), " haplotypes:")
for (i in seq_len(nrow(px1))) {
  cnt <- table(factor(px1_calls[[i]], levels = 0:2, labels = c("missing", "WT", px1$mutation[i])))
  message("  ", px1$name[i], " (", px1$mutation[i], ", ", px1$source[i], "): ",
          paste(sprintf("%s=%d", names(cnt), as.integer(cnt)), collapse = "  "))
}

# ---------------------------------------------------------------------------
# 4. QC: for any PX1 marker already in the WGS map, the external calls must
#    agree with it. This is what proves the row order above is right - do not
#    remove it. (Only the SeekDeep markers are keyed by position rather than by
#    sample name, so this is the check that matters for them; the VCF ones are
#    joined on sample id and cannot be mis-ordered.)
# ---------------------------------------------------------------------------
for (i in seq_len(nrow(px1))) {
  j <- match(px1$name[i], map_chr$name)
  if (is.na(j)) next
  wgs <- hap_alleles[, j]
  ext <- px1_calls[[px1$name[i]]]
  both <- wgs != 0L & ext != 0L
  agree <- sum(wgs[both] == ext[both])
  message(sprintf("QC %s: %d/%d comparable calls agree between the matrix and %s",
                  px1$name[i], agree, sum(both), px1$source[i]))
  if (sum(both) > 0 && agree / sum(both) < 0.95) {
    stop("The binary matrix and the ", px1$source[i], " calls disagree at ", px1$name[i],
         " - the haplotype/sample mapping is probably wrong; stopping.")
  }
}

# ---------------------------------------------------------------------------
# 5. Splice the missing PX1 markers into the map and the haplotype matrix.
#    rehh needs markers in ascending position order within a chromosome, so the
#    new columns go in at their sorted index, in both objects at once.
# ---------------------------------------------------------------------------
for (i in seq_len(nrow(px1))) {
  nm <- px1$name[i]
  if (nm %in% map_chr$name) {
    message("Marker ", nm, " already present in the WGS map - keeping the WGS calls.")
    next
  }
  at <- sum(map_chr$pos < px1$pos[i])                    # insert after this many markers
  new_row <- data.frame(name = nm, chr = chr, pos = px1$pos[i],
                        anc = px1$ref[i], der = px1$alt[i],
                        stringsAsFactors = FALSE)
  map_chr <- rbind(map_chr[seq_len(at), ], new_row,
                   map_chr[seq_len(nrow(map_chr) - at) + at, ])
  hap_alleles <- cbind(hap_alleles[, seq_len(at), drop = FALSE],
                       px1_calls[[nm]],
                       hap_alleles[, seq_len(ncol(hap_alleles) - at) + at, drop = FALSE])
  message("Inserted ", nm, " (", px1$mutation[i], ") at chr", chr, " marker index ", at + 1)
}
colnames(hap_alleles) <- NULL
stopifnot(!is.unsorted(map_chr$pos))

# rebuild the whole-genome map with the chr7 block replaced
map_out <- rbind(map[map$chr < chr, ], map_chr, map[map$chr > chr, ])

# ---------------------------------------------------------------------------
# 6. Write the augmented rehh inputs.
#    Using the real sample ids as haplotype ids so furcation plots are readable.
# ---------------------------------------------------------------------------
write.table(cbind(hap_samples, hap_alleles), out_hap,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(map_out, out_map,
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Wrote ", out_hap, " (", nrow(hap_alleles), " haplotypes x ",
        ncol(hap_alleles), " markers) and ", out_map)

# ---------------------------------------------------------------------------
# 7. EHH / furcation plots
# ---------------------------------------------------------------------------
hap_chr_pop <- data2haplohh(hap_file = out_hap, map_file = out_map,
                            recode.allele = FALSE,
                            chr.name = chr, min_perc_geno.hap = 0.8, min_perc_geno.mrk = 0.8)

for (i in seq_len(nrow(px1))) {
  snp <- px1$name[i]
  pos <- px1$pos[i]
  mut <- px1$mutation[i]

  tryCatch({
    res <- calc_ehh(hap_chr_pop,
                    mrk = snp,
                    include_nhaplo = TRUE,
                    include_zero_values = TRUE,
                    discard_integration_at_border = FALSE,
                    polarized = TRUE)

    # Work out how far out the curves have to be drawn for the decay to zero to
    # be visible on both sides of the focal SNP, then pad a little.
    if (is.null(ehh_xlim)) {
      ed <- res$ehh
      ecols <- grep("^EHH", names(ed), value = TRUE)
      nz <- ed$POSITION[rowSums(ed[, ecols, drop = FALSE] > 0, na.rm = TRUE) > 0]
      if (length(nz) > 0) {
        # first marker left of, and first right of, the last non-zero one
        left  <- suppressWarnings(max(ed$POSITION[ed$POSITION < min(nz)]))
        right <- suppressWarnings(min(ed$POSITION[ed$POSITION > max(nz)]))
        if (!is.finite(left))  left  <- min(ed$POSITION)
        if (!is.finite(right)) right <- max(ed$POSITION)
        pad <- 0.03 * (right - left)
        this_xlim <- c(left - pad, right + pad)
      } else {
        this_xlim <- range(ed$POSITION)
      }
    } else {
      this_xlim <- ehh_xlim
    }
    message("  EHH x-range for ", snp, ": ",
            sprintf("%.0f - %.0f", this_xlim[1], this_xlim[2]))

    pdf(paste0("pf3d7_0720700_ehh_plot_", country, "_", pos, ".pdf"),
        width = 9, height = 6.5)
    plot(res, xlim = this_xlim,
         main = paste0("EHH - ", country, "\n", mut, " (pos ", pos, ")"))
    dev.off()

    f <- calc_furcation(hap_chr_pop, mrk = snp)

    pdf(paste0("pf3d7_0720700_furcation_plot_", country, "_", pos, ".pdf"),
        width = 7, height = 7)
    plot(f, lwd = 0.2,
         hap.names = hap.names(hap_chr_pop),
         legend.xy.coords = "none",
         xlim = furcation_xlim,
         main = paste0("Furcation - ", country, "\n", mut, " (pos ", pos, ")"))
    dev.off()
    message("Plotted ", snp, " (", mut, ")")
  }, error = function(e) {
    message("Skipping SNP ", snp, " (", country, "): ", conditionMessage(e))
  })
}

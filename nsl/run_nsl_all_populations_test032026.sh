#!/bin/bash
# TEST VERSION (March 2026) — do not overwrite production outputs.
# Differs from run_nsl_all_populations.sh in the filtering applied before selscan:
#   - Hypervariable regions are taken from hypervarregions_032026.csv (new CSV)
#     instead of hyper_variable_regions.txt
#   - Top 10 Tajima's D windows are excluded as before
# All outputs are written to <POP>_only_test032026/ directories.
#
# Usage:
#   ./run_nsl_all_populations_test032026.sh [populations_file]

source /mnt/storage13/leen/miniforge3/etc/profile.d/conda.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${DATA_DIR:-$PWD}"

POPULATIONS_FILE="${1:-$DATA_DIR/populations.txt}"
METADATA="$DATA_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
INPUT_VCF="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"

# New hypervariable regions CSV (tab-separated, has header row: Chromosome Start End Gene_ID)
NEW_HYPER_VAR_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/nsl/hypervarregions_032026.csv"

# Additional metadata for special population groups (same sources as hmmibd)
R622I_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/R622I_longformat.csv"
HRP_META="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/sample_meta_with_hrp_status_updated.csv"

if [[ ! -f "$POPULATIONS_FILE" ]]; then
    echo "ERROR: Populations file '$POPULATIONS_FILE' not found."
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Metadata file '$METADATA' not found."
    exit 1
fi

if [[ ! -f "$NEW_HYPER_VAR_CSV" ]]; then
    echo "ERROR: Hypervariable regions CSV '$NEW_HYPER_VAR_CSV' not found."
    exit 1
fi

if [[ ! -f "$R622I_CSV" ]]; then
    echo "ERROR: R622I CSV '$R622I_CSV' not found."
    exit 1
fi

if [[ ! -f "$HRP_META" ]]; then
    echo "ERROR: HRP metadata '$HRP_META' not found."
    exit 1
fi

while IFS= read -r POP || [[ -n "$POP" ]]; do
    [[ -z "$POP" || "$POP" == \#* ]] && continue

    echo ""
    echo "=========================================="
    echo "Processing population: $POP  [TEST 032026]"
    echo "=========================================="

    # Use a separate output directory so production outputs are untouched
    POP_DIR="$DATA_DIR/${POP}_only_test032026"
    PREFIX="${POP}_only_test032026"

    mkdir -p "$POP_DIR/vcfs_by_chromosome/imputed"
    mkdir -p "$POP_DIR/selscan_output"

    # ------------------------------------------------------------------
    # Generate sample list from metadata
    # ------------------------------------------------------------------
    echo "[0/8] Generating sample list for $POP from metadata..."
    conda activate fastq2matrix
    python3 "$SCRIPT_DIR/extract_samples.py" "$POP" "$METADATA" \
        --r622i-csv "$R622I_CSV" \
        --hrp-meta "$HRP_META" > "$POP_DIR/${POP}.txt"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to extract samples for $POP. Skipping."
        continue
    fi
    N_SAMPLES=$(wc -l < "$POP_DIR/${POP}.txt")
    echo "  $N_SAMPLES samples found for $POP"

    cd "$POP_DIR" || { echo "ERROR: Cannot cd to $POP_DIR. Skipping $POP."; continue; }

    # ------------------------------------------------------------------
    # 1. Extract population samples from the annotated bi-allelic VCF
    # ------------------------------------------------------------------
    echo "[1/8] Extracting samples for $POP..."
    bcftools view -S "${POP}.txt" "$INPUT_VCF" \
        -O z -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"

    # ------------------------------------------------------------------
    # 2. Remove non-variant sites
    # ------------------------------------------------------------------
    echo "[2/8] Removing non-variant sites..."
    bcftools view -c 1 "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz" \
        -Oz -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.filtered.vcf.gz"

    # ------------------------------------------------------------------
    # 3. Calculate Tajima's D on pre-exclusion VCF (after non-variant removal),
    #    then build a combined exclusion BED from:
    #      (a) hypervarregions_032026.csv regions
    #      (b) top 10% of Tajima's D windows
    #    Apply both exclusions together in a single bcftools step.
    # ------------------------------------------------------------------
    echo "[3/8] Computing Tajima's D (pre-exclusion) and building combined exclusion BED..."
    vcftools --gzvcf "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.filtered.vcf.gz" \
        --out tajimasd --TajimaD 10000

    # Top 10% of windows by Tajima's D (ceiling division)
    tail -n +2 tajimasd.Tajima.D | sort -k4 -rn > tajimasd_sorted.txt
    n_windows=$(wc -l < tajimasd_sorted.txt)
    top10pct=$(( (n_windows + 9) / 10 ))
    echo "  Tajima's D: $n_windows windows total, excluding top $top10pct (10%)"
    head -n "$top10pct" tajimasd_sorted.txt \
        | awk '{print $1"\t"$2"\t"$2+10000}' > top10pct_tajima_windows.bed

    # Hypervariable regions from new CSV (strip header, cols 1-3)
    tail -n +2 "$NEW_HYPER_VAR_CSV" | cut -f1-3 > hypervar_regions_032026.bed

    # Merge into one exclusion BED
    cat hypervar_regions_032026.bed top10pct_tajima_windows.bed \
        | sort -k1,1 -k2,2n > combined_exclusion_032026.bed
    echo "  Combined exclusion BED: $(wc -l < combined_exclusion_032026.bed) regions"

    # Single exclusion step (replaces separate steps 3 & 4 from the original)
    echo "[4/8] Excluding combined regions (hypervar + top 10% Tajima's D)..."
    n_snps=$(bcftools view -H -v snps "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.filtered.vcf.gz" | wc -l)
    echo "  SNPs before exclusion: $n_snps"

    bcftools view --targets-file "^combined_exclusion_032026.bed" \
        "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.filtered.vcf.gz" \
        -Oz -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz"
    bcftools index "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz"

    n_snps=$(bcftools view -H -v snps "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz" | wc -l)
    echo "  SNPs after exclusion: $n_snps"

    # ------------------------------------------------------------------
    # 5. Exclude sites with > 10% missingness
    # ------------------------------------------------------------------
    echo "[5/8] Filtering by missingness..."
    bcftools view -e 'F_MISSING > 0.1' \
        "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz" \
        -Oz -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.filtered.vcf.gz"
    bcftools index "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.filtered.vcf.gz"

    n_snps=$(bcftools view -H -v snps "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.filtered.vcf.gz" | wc -l)
    echo "  SNPs after missingness filtering: $n_snps"

    # ------------------------------------------------------------------
    # 6. Split VCF by chromosome
    # ------------------------------------------------------------------
    echo "[6/8] Splitting by chromosome..."
    for chr in Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 \
               Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 \
               Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3; do
        bcftools view -r "$chr" \
            "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.filtered.vcf.gz" \
            -O z -o "vcfs_by_chromosome/${chr}.vcf.gz" \
            && bcftools index "vcfs_by_chromosome/${chr}.vcf.gz"
    done

    # ------------------------------------------------------------------
    # 7. Impute with Beagle, then run selscan nSL and normalise
    # ------------------------------------------------------------------
    echo "[7/8] Imputing and running nSL..."
    conda activate imputation
    for chr in Pf3D7_0{1..9}_v3 Pf3D7_{10..14}_v3; do
        java -Xmx8g -jar /mnt/storage13/leen/miniforge3/envs/imputation/share/beagle-5.5_27Feb25.75f-0/beagle.jar \
            gt="vcfs_by_chromosome/${chr}.vcf.gz" \
            out="vcfs_by_chromosome/imputed/${chr}.imputed" \
            chrom="${chr}" \
            nthreads=4 \
            iterations=20 \
            ne=1000000
        bcftools index "vcfs_by_chromosome/imputed/${chr}.imputed.vcf.gz"
    done

    conda activate selscan

    for VCF in vcfs_by_chromosome/imputed/*.imputed.vcf.gz; do
        CHR=$(basename "$VCF" .vcf.gz)
        echo "  selscan: $CHR"
        selscan --nsl --vcf "$VCF" --unphased --maf 0.01 --out "selscan_output/NSL_${CHR}"
    done

    for FILE in selscan_output/NSL_Pf3D7*.imputed.nsl.out; do
        echo "  normalising: $(basename "$FILE")"
        norm --nsl --bins 100 --files "$FILE"
    done

    MERGED_FILE="NSL_allchromosomes_norm_merged_${PREFIX}.txt"
    echo -e "File\tChrom\tPos\tFreqBin\tIHH1\tIHH0\tRaw_nSL\tNorm_nSL\tCritFlag" > "$MERGED_FILE"
    for FILE in selscan_output/NSL_Pf3D7*imputed*.norm; do
        BASENAME=$(basename "$FILE")
        awk -v f="$BASENAME" '{print f "\t" $0}' "$FILE" >> "$MERGED_FILE"
    done
    echo "  Merged output: $MERGED_FILE"

    # ------------------------------------------------------------------
    # 8. Annotate hits
    # ------------------------------------------------------------------
    echo "[8/8] Annotating..."
    conda activate fastq2matrix

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\n' \
        "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz" > annotations.tsv

    python3 << 'PYEOF'
import pandas as pd

anno = pd.read_csv("annotations.tsv", sep="\t", header=None)
gene_name1 = anno[4].str.split('|').str[1]
gene_name2 = anno[4].str.split('|').str[2]
anno_final = anno.iloc[:, 0:4].copy()
anno_final['Gene_name'] = gene_name1 + " (" + gene_name2 + ")"
anno_final.columns = ['chr', 'pos', 'ref', 'alt', 'Gene_name']
anno_final.to_csv("annotations_final.tsv", index=False, sep="\t")
PYEOF

    Rscript -e "
library(dplyr); library(stringr)
setwd('$POP_DIR')
output <- read.csv('$MERGED_FILE', sep='\t')
output <- output %>% mutate(Chrom = str_extract(File, 'Pf3D7_[0-9]+_v3'))
anno   <- read.csv('annotations_final.tsv', sep='\t')
out_ann <- left_join(output, anno, by = c('Chrom' = 'chr', 'Pos' = 'pos'))
write.csv(out_ann, 'NSL_allchromosomes_norm_merged_${PREFIX}_annotated.csv', row.names = FALSE)
"

    echo "Done: $POP -> $POP_DIR/NSL_allchromosomes_norm_merged_${PREFIX}_annotated.csv"

done < "$POPULATIONS_FILE"

# ------------------------------------------------------------------
# Merge all population annotated outputs into one CSV
# ------------------------------------------------------------------
echo ""
echo "=========================================="
echo "Merging all population outputs..."
echo "=========================================="

cd "$DATA_DIR"
MERGED_ALL="$DATA_DIR/NSL_allchromosomes_norm_merged_ALL_populations_test032026_annotated.csv"
echo "Population,File,Chrom,Pos,FreqBin,IHH1,IHH0,Raw_nSL,Norm_nSL,CritFlag,ref,alt,Gene_name" > "$MERGED_ALL"

while IFS= read -r POP || [[ -n "$POP" ]]; do
    [[ -z "$POP" || "$POP" == \#* ]] && continue
    CSV="$DATA_DIR/${POP}_only_test032026/NSL_allchromosomes_norm_merged_${POP}_only_test032026_annotated.csv"
    if [[ ! -f "$CSV" ]]; then
        echo "  WARNING: output not found for $POP, skipping: $CSV"
        continue
    fi
    tail -n +2 "$CSV" | awk -v pop="$POP" 'BEGIN{FS=","; OFS=","} {print pop, $0}' >> "$MERGED_ALL"
    echo "  Added: $POP"
done < "$POPULATIONS_FILE"

N=$(tail -n +2 "$MERGED_ALL" | wc -l)
echo "Done. Total rows: $N"
echo "Output: $MERGED_ALL"

echo ""
echo "All populations complete (test032026 run)."

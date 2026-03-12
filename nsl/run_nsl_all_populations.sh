#!/bin/bash
# Master wrapper: runs the full nsl pipeline for each population listed in a file.
#
# Usage:
#   ./run_nsl_all_populations.sh [populations_file]
#
# populations_file (default: populations.txt) — one population name per line.
# Sample lists are auto-generated from the metadata file via extract_samples.py.
#
# Outputs per population <POP>:
#   <POP>_only/NSL_allchromosomes_norm_merged_<POP>_only_annotated.csv

source /mnt/storage13/leen/miniforge3/etc/profile.d/conda.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"  # directory containing this script (repo)
DATA_DIR="${DATA_DIR:-$PWD}"                 # working data directory (default: current dir)

POPULATIONS_FILE="${1:-$DATA_DIR/populations.txt}"
METADATA="$DATA_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
INPUT_VCF="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
HYPER_VAR_REGIONS="$DATA_DIR/hyper_variable_regions.txt"

if [[ ! -f "$POPULATIONS_FILE" ]]; then
    echo "ERROR: Populations file '$POPULATIONS_FILE' not found."
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Metadata file '$METADATA' not found."
    exit 1
fi

while IFS= read -r POP || [[ -n "$POP" ]]; do
    # Skip empty lines and comments
    [[ -z "$POP" || "$POP" == \#* ]] && continue

    echo ""
    echo "=========================================="
    echo "Processing population: $POP"
    echo "=========================================="

    POP_DIR="$DATA_DIR/${POP}_only"
    PREFIX="${POP}_only"

    mkdir -p "$POP_DIR/vcfs_by_chromosome/imputed"
    mkdir -p "$POP_DIR/selscan_output"

    # ------------------------------------------------------------------
    # Generate sample list from metadata
    # ------------------------------------------------------------------
    echo "[0/8] Generating sample list for $POP from metadata..."
    conda activate fastq2matrix
    python3 "$SCRIPT_DIR/extract_samples.py" "$POP" "$METADATA" > "$POP_DIR/${POP}.txt"
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
    # 3. Exclude hypervariable regions
    # ------------------------------------------------------------------
    echo "[3/8] Excluding hypervariable regions..."
    bcftools view --targets-file "^$HYPER_VAR_REGIONS" \
        "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.filtered.vcf.gz" \
        -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.vcf.gz"
    bcftools index "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.vcf.gz"

    n_snps=$(bcftools view -H -v snps "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.vcf.gz" | wc -l)
    echo "  SNPs before Tajima's D filtering: $n_snps"

    # ------------------------------------------------------------------
    # 4. Exclude top 10% Tajima's D windows
    # ------------------------------------------------------------------
    echo "[4/8] Filtering by Tajima's D..."
    vcftools --gzvcf "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.vcf.gz" \
        --out tajimasd --TajimaD 10000

    tail -n +2 tajimasd.Tajima.D | sort -k4 -rn | head -10 \
        | awk '{print $1"\t"$2"\t"$2+10000}' > top10_tajima_windows.bed

    bcftools view \
        --targets-file ^top10_tajima_windows.bed \
        "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.vcf.gz" \
        -o "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz" \
        -O z
    bcftools index "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz"

    n_snps=$(bcftools view -H -v snps "${PREFIX}.filt.csq.bi.GT.miss0.2.vqslod.filt.hypervarrem.snps.filtered.tajrem.vcf.gz" | wc -l)
    echo "  SNPs after Tajima's D filtering: $n_snps"

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

    # Run selscan nSL per chromosome
    for VCF in vcfs_by_chromosome/imputed/*.imputed.vcf.gz; do
        CHR=$(basename "$VCF" .vcf.gz)
        echo "  selscan: $CHR"
        selscan --nsl --vcf "$VCF" --unphased --maf 0.01 --out "selscan_output/NSL_${CHR}"
    done

    # Normalise
    for FILE in selscan_output/NSL_Pf3D7*.imputed.nsl.out; do
        echo "  normalising: $(basename "$FILE")"
        norm --nsl --bins 100 --files "$FILE"
    done

    # Merge normalised outputs
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

echo ""
echo "All populations complete."

#!/bin/bash
# Merge nSL annotated outputs for populations listed in populations.txt into a single CSV.
# Adds a "Population" column.
#
# Usage:
#   ./merge_all_populations_nsl.sh [populations_file] [output_file]
#
# Output (default): NSL_allchromosomes_norm_merged_ALL_populations_annotated.csv

DATA_DIR="${DATA_DIR:-$PWD}"                 # working data directory (default: current dir)
POPULATIONS_FILE="${1:-$DATA_DIR/populations.txt}"
OUT="${2:-$DATA_DIR/NSL_allchromosomes_norm_merged_ALL_populations_annotated.csv}"

if [[ ! -f "$POPULATIONS_FILE" ]]; then
    echo "ERROR: Populations file '$POPULATIONS_FILE' not found."
    exit 1
fi

echo "Merging populations from: $POPULATIONS_FILE"
echo "Output: $OUT"

# Write header with Population column prepended
echo "Population,File,Chrom,Pos,FreqBin,IHH1,IHH0,Raw_nSL,Norm_nSL,CritFlag,ref,alt,Gene_name" > "$OUT"

while IFS= read -r POP || [[ -n "$POP" ]]; do
    # Skip empty lines and comments
    [[ -z "$POP" || "$POP" == \#* ]] && continue

    CSV="$DATA_DIR/${POP}_only/NSL_allchromosomes_norm_merged_${POP}_only_annotated.csv"

    if [[ ! -f "$CSV" ]]; then
        echo "  WARNING: output not found for $POP, skipping: $CSV"
        continue
    fi

    # Skip header line, prepend population name
    tail -n +2 "$CSV" | awk -v pop="$POP" 'BEGIN{FS=","; OFS=","} {print pop, $0}' >> "$OUT"
    echo "  Added: $POP"

done < "$POPULATIONS_FILE"

N=$(tail -n +2 "$OUT" | wc -l)
echo ""
echo "Done. Total rows: $N"
echo "Output: $OUT"

#!/usr/bin/env python3
"""
Extract sample IDs from the metadata file for a given population name.

Usage:
    python3 extract_samples.py <population_name> <metadata_file>
                               [--r622i-csv PATH] [--hrp-meta PATH]

Prints one wgs_id per line to stdout.

Population name resolution order:
  1. Special Ethiopia subgroups (Ethiopia_new, Ethiopia_old,
     Ethiopia_north, Ethiopia_east_west_south)
  2. Kelch13 R622I groups (r622i_wt, r622i_mut) — requires --r622i-csv
  3. HRP status groups (hrp_wild, hrp3_only, hrp_double_deleted) — requires --hrp-meta
  4. Exact match on Country column
  5. Exact match on Region column
"""

import sys
import argparse
import pandas as pd

NORTH_DISTRICTS = {
    "Atseged Tsimbila",
    "West Armachiho",
    "K/Humera",
    "Gondar",
    "Shewa Robit Town",
}

# Mapping from populations.txt names to genotype values in R622I_longformat.csv
R622I_GROUPS = {
    "r622i_wt":  "0/0",
    "r622i_mut": "1/1",
}

# Mapping from populations.txt names to hrp_status values in HRP metadata
HRP_GROUPS = {
    "hrp_wild":           "wild",
    "hrp3_only":          "hrp3_only",
    "hrp_double_deleted": "mutant",
}


def get_samples(pop, meta_file, r622i_csv=None, hrp_meta=None):
    df = pd.read_csv(meta_file, sep="\t")

    # Restrict to monoclonal samples (FWS >= 0.95) — required for haplotype-based
    # statistics (nSL/iHS); polyclonal samples corrupt phasing and EHH signals.
    n_before = len(df)
    df = df[df["fws"] >= 0.95]
    n_after = len(df)
    print(f"  Monoclonal filter (FWS >= 0.95): {n_before} -> {n_after} samples retained", file=sys.stderr)

    if pop == "Ethiopia_new":
        mask = (df["Country"] == "Ethiopia") & (df["Year"].astype(str) == "2024")
    elif pop == "Ethiopia_old":
        mask = (df["Country"] == "Ethiopia") & (df["Year"].astype(str) != "2024")
    elif pop == "Ethiopia_north":
        mask = (df["Country"] == "Ethiopia") & (df["District"].isin(NORTH_DISTRICTS))
    elif pop == "Ethiopia_east_west_south":
        mask = (df["Country"] == "Ethiopia") & (~df["District"].isin(NORTH_DISTRICTS))
    elif pop in R622I_GROUPS:
        if not r622i_csv:
            print(f"ERROR: '{pop}' requires --r622i-csv.", file=sys.stderr)
            sys.exit(1)
        genotype = R622I_GROUPS[pop]
        r622i_df = pd.read_csv(r622i_csv)
        # Exclude heterozygous (0/1) and missing (./.) — same as hmmibd
        r622i_df = r622i_df[r622i_df["Genotype"].isin(["0/0", "1/1"])]
        samples_with_geno = r622i_df[r622i_df["Genotype"] == genotype]["Sample"].unique()
        mask = (df["Country"] == "Ethiopia") & df["wgs_id"].isin(samples_with_geno)
    elif pop in HRP_GROUPS:
        if not hrp_meta:
            print(f"ERROR: '{pop}' requires --hrp-meta.", file=sys.stderr)
            sys.exit(1)
        hrp_status_val = HRP_GROUPS[pop]
        hrp_df = pd.read_csv(hrp_meta)
        samples_with_status = hrp_df[hrp_df["hrp_status"] == hrp_status_val]["Sample"].unique()
        mask = (df["Country"] == "Ethiopia") & df["wgs_id"].isin(samples_with_status)
    elif pop in df["Country"].values:
        mask = df["Country"] == pop
    elif pop in df["Region"].values:
        mask = df["Region"] == pop
    else:
        print(f"ERROR: '{pop}' not found in Country or Region columns.", file=sys.stderr)
        sys.exit(1)

    samples = df.loc[mask, "wgs_id"]
    if samples.empty:
        print(f"ERROR: No samples found for population '{pop}'.", file=sys.stderr)
        sys.exit(1)

    print("\n".join(samples.tolist()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("population", help="Population name")
    parser.add_argument("metadata", help="Metadata TSV file (tab-separated)")
    parser.add_argument("--r622i-csv", default=None,
                        help="R622I long-format CSV (columns: Chrom,Position,Ref,Alt,Sample,Genotype)")
    parser.add_argument("--hrp-meta", default=None,
                        help="HRP status metadata CSV (columns include: Sample, hrp_status)")
    args = parser.parse_args()
    get_samples(args.population, args.metadata, args.r622i_csv, args.hrp_meta)

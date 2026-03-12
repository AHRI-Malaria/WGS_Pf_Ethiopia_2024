#!/usr/bin/env python3
"""
Extract sample IDs from the metadata file for a given population name.

Usage:
    python3 extract_samples.py <population_name> <metadata_file>

Prints one wgs_id per line to stdout.

Population name resolution order:
  1. Special Ethiopia subgroups (Ethiopia_new, Ethiopia_old,
     Ethiopia_north, Ethiopia_east_west_south)
  2. Exact match on Country column
  3. Exact match on Region column
"""

import sys
import pandas as pd

NORTH_DISTRICTS = {
    "Atseged Tsimbila",
    "West Armachiho",
    "K/Humera",
    "Gondar",
    "Shewa Robit Town",
}

def get_samples(pop, meta_file):
    df = pd.read_csv(meta_file, sep="\t")

    if pop == "Ethiopia_new":
        mask = (df["Country"] == "Ethiopia") & (df["Year"].astype(str) == "2024")
    elif pop == "Ethiopia_old":
        mask = (df["Country"] == "Ethiopia") & (df["Year"].astype(str) != "2024")
    elif pop == "Ethiopia_north":
        mask = (df["Country"] == "Ethiopia") & (df["District"].isin(NORTH_DISTRICTS))
    elif pop == "Ethiopia_east_west_south":
        mask = (df["Country"] == "Ethiopia") & (~df["District"].isin(NORTH_DISTRICTS))
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
    if len(sys.argv) != 3:
        print("Usage: extract_samples.py <population_name> <metadata_file>", file=sys.stderr)
        sys.exit(1)
    get_samples(sys.argv[1], sys.argv[2])

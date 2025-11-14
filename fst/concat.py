#! /usr/bin/env python
import argparse
import itertools
import pandas as pd
import numpy as np

def compile_all_fst(meta, pop, out_dir, selected_sites=None):
    """
    Compile all Fst values from existing paired Fst files
    """
    geo_r = meta[pop].unique()
    geo_r = [x for x in geo_r if str(x) != 'nan']
    
    # Filter for selected sites if provided
    if selected_sites:
        geo_r = [x for x in geo_r if x in selected_sites]
        print(f"📍 Filtered to {len(geo_r)} sites: {', '.join(geo_r)}\n")
    
    if len(geo_r) == 0:
        print("❌ No matching sites found in metadata!")
        return None
    
    all_fst_data = []
    
    for comb in itertools.combinations(geo_r, 2):
        pop1 = comb[0]
        pop2 = comb[1]
        filename1 = out_dir + "/" + pop1 + "_vs_" + pop2 + ".paired.snp_fst_df.csv"
        filename2 = out_dir + "/" + pop2 + "_vs_" + pop1 + ".paired.snp_fst_df.csv"
        
        try:
            site_fst = pd.read_csv(filename1, header=0)
        except FileNotFoundError:
            try:
                site_fst = pd.read_csv(filename2, header=0)
            except FileNotFoundError:
                print(f"⚠️ No file found for {pop1} vs {pop2}, skipping...")
                continue
        
        # Clean Fst values
        site_fst["Fst"] = pd.to_numeric(site_fst["Fst"], errors="coerce").fillna(0).clip(lower=0)
        
        # Add population comparison columns
        site_fst["Population_1"] = pop1
        site_fst["Population_2"] = pop2
        
        # Append to list
        all_fst_data.append(site_fst)
        
        print(f"✓ {pop1} vs {pop2}: {len(site_fst)} SNPs")
    
    if len(all_fst_data) == 0:
        print("❌ No Fst data found for selected sites!")
        return None
    
    # Combine all dataframes
    all_fst_df = pd.concat(all_fst_data, ignore_index=True)
    
    print(f"\n✅ Total: {len(all_fst_df)} Fst values across all comparisons")
    
    return all_fst_df

def main(args):
    # Define selected sites
    selected_sites = ['Arbaminch', 'Gambella', 'Gonder', 'Metehera_Fentale', 'Mizan']
    
    print(f"🎯 Filtering for sites: {', '.join(selected_sites)}\n")
    
    # Read metadata
    meta = pd.read_csv(args.metadata, sep=',')
    
    # Print available sites in metadata
    available_sites = meta[args.pop].unique()
    available_sites = [x for x in available_sites if str(x) != 'nan']
    print(f"📋 Available sites in metadata: {', '.join(available_sites)}\n")
    
    # Compile all Fst values
    print("📊 Compiling all Fst values...\n")
    all_fst_df = compile_all_fst(meta, args.pop, args.fst_dir, selected_sites)
    
    if all_fst_df is None or len(all_fst_df) == 0:
        print("\n❌ No data to process. Exiting.")
        return
    
    # Save all Fst values
    output_file = args.out_dir + "/all_fst_values.csv"
    all_fst_df.to_csv(output_file, index=False)
    print(f"\n✅ All Fst values saved to {output_file}")
    
    # Print summary statistics
    print(f"\n📈 Summary Statistics:")
    print(f"   Total SNPs: {len(all_fst_df)}")
    print(f"   Mean Fst: {all_fst_df['Fst'].mean():.4f}")
    print(f"   Median Fst: {all_fst_df['Fst'].median():.4f}")
    print(f"   Max Fst: {all_fst_df['Fst'].max():.4f}")
    print(f"   SNPs with Fst > 0.8: {(all_fst_df['Fst'] > 0.8).sum()}")
    
    # Summary by comparison
    print(f"\n📊 Summary by Comparison:")
    comparison_summary = all_fst_df.groupby(['Population_1', 'Population_2']).agg({
        'Fst': ['count', 'mean', 'median', 'max']
    }).round(4)
    print(comparison_summary)
    
    print("\n✨ Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Compile all Fst values from pairwise comparisons',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--metadata', 
                        help='CSV file containing metadata', 
                        required=True)
    parser.add_argument('--pop', 
                        help='Population column name in metadata', 
                        required=True)
    parser.add_argument('--fst_dir', 
                        help='Directory containing paired Fst files', 
                        required=True)
    parser.add_argument('--out_dir', 
                        help='Output directory for results', 
                        required=True)
    
    args = parser.parse_args()
    main(args)
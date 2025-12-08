#! /usr/bin/env python
import argparse
import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def create_custom_colormap():
    colors = ['#fde724', '#b0d743','#90d743', '#21918c', '#31688e', '#443983'] 
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)
    return cmap

def calculate_mean_fst(meta, pop, out_dir, selected_sites=None):
    """
    Calculate mean Fst values from existing paired Fst files
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
    
    mean_fst_data = []
    
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
        
        # Calculate mean Fst
        mean_val = np.mean(site_fst['Fst'])
        
        mean_fst_data.append({
            'Population_1': pop1,
            'Population_2': pop2,
            'Mean_Fst': mean_val,
            'N_SNPs': len(site_fst)
        })
        
        print(f"✓ {pop1} vs {pop2}: Mean Fst = {mean_val:.4f} ({len(site_fst)} SNPs)")
    
    if len(mean_fst_data) == 0:
        print("❌ No Fst data found for selected sites!")
        return None
    
    mean_fst_df = pd.DataFrame(mean_fst_data)
    return mean_fst_df

def create_fst_matrix(mean_fst_df):
    """
    Convert pairwise Fst values into a symmetric matrix
    """
    # Get all unique populations
    all_pops = sorted(list(set(mean_fst_df['Population_1'].tolist() + 
                                mean_fst_df['Population_2'].tolist())))
    
    # Create empty matrix
    fst_matrix = pd.DataFrame(0.0, index=all_pops, columns=all_pops)
    
    # Fill the matrix with mean Fst values (symmetric)
    for _, row in mean_fst_df.iterrows():
        pop1 = row['Population_1']
        pop2 = row['Population_2']
        mean_fst = row['Mean_Fst']
        
        fst_matrix.loc[pop1, pop2] = mean_fst
        fst_matrix.loc[pop2, pop1] = mean_fst
    
    # Diagonal is 0 (population vs itself)
    np.fill_diagonal(fst_matrix.values, 0)
    
    return fst_matrix

def plot_fst_heatmap(fst_matrix, out_dir, max_fst=None, figsize=(10, 8)):
    """
    Create and save heatmap of Fst values with custom colormap
    """
    # Create custom colormap
    custom_cmap = create_custom_colormap()
    
    # If max_fst not provided, use actual max from data
    if max_fst is None:
        max_fst = fst_matrix.values[~np.eye(len(fst_matrix), dtype=bool)].max()  # Exclude diagonal
        print(f"📊 Using max Fst value from data: {max_fst:.4f}")
    
    # Create a mask for the upper triangle
    mask = np.triu(np.ones_like(fst_matrix, dtype=bool))

    # Create figure with better proportions
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap with custom colormap and fixed scale
    sns.heatmap(fst_matrix, 
                mask=mask,
                annot=True, 
                fmt='.3f',  
                cmap=custom_cmap, 
                vmin=0, 
                vmax=max_fst,
                cbar_kws={
                    'label': 'Mean Fst',
                    'shrink': 1, 
                    'pad': 0.02
                },
                square=True, 
                linewidths=1, 
                linecolor='white',
                annot_kws={
                    'fontsize': 12,
                    'weight': 'bold'  
                },
                ax=ax)

    
    cbar = ax.collections[0].colorbar
    cbar.set_label('Mean Fst', fontsize=14, fontweight='bold')
    
    # Improve tick labels - make them larger and more readable
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(rotation=0, fontsize=14)
    
    # Remove axis spines for cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.tight_layout()
    
    # Save plots
    plt.savefig(out_dir + f"/fst_heatmap.png", dpi=300, bbox_inches='tight')
    plt.savefig(out_dir + f"/fst_heatmap.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"\n✅ Heatmap saved to {out_dir}/fst_heatmap.png")
    print(f"✅ Heatmap saved to {out_dir}/fst_heatmap.pdf")

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
    
    # Calculate mean Fst values
    print("📊 Calculating mean Fst values...\n")
    mean_fst_df = calculate_mean_fst(meta, args.pop, args.fst_dir, selected_sites)
    
    if mean_fst_df is None or len(mean_fst_df) == 0:
        print("\n❌ No data to process. Exiting.")
        return
    
    # Save mean Fst summary
    mean_fst_df.to_csv(args.out_dir + "/mean_fst_summary.csv", index=False)
    print(f"\n✅ Mean Fst summary saved to {args.out_dir}/mean_fst_summary.csv")
    
    # Create Fst matrix
    fst_matrix = create_fst_matrix(mean_fst_df)
    
    # Save matrix
    fst_matrix.to_csv(args.out_dir + "/fst_matrix.csv")
    print(f"✅ Fst matrix saved to {args.out_dir}/fst_matrix.csv")
    
    # Print matrix
    print(f"\n📊 Fst Matrix:\n")
    print(fst_matrix.round(4))
    
    # Create heatmap
    print("\n🎨 Creating heatmap...")
    plot_fst_heatmap(fst_matrix, args.out_dir, 
                     max_fst=args.max_fst, 
                     figsize=(args.figsize, args.figsize))
    
    print("\n✨ Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate mean Fst values and create heatmap for selected sites',
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
    parser.add_argument('--max_fst', 
                        type=float,
                        help='Maximum Fst value for color scale', 
                        default=0.1)
    parser.add_argument('--figsize', 
                        type=int,
                        help='Figure size (will be square)', 
                        default=10)
    
    args = parser.parse_args()
    main(args)
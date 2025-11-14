import itertools
import numpy as np
import pandas as pd
import allel
import os
import zarr
# --- Load data ---
callset = zarr.open('072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps_ethiopia.zarr', mode='r')
genotype_all = allel.GenotypeArray(callset['calldata/GT'])
samples = callset['samples'][:]  # sample IDs in Zarr
df_samples = pd.read_csv('ethiopia_samples.csv', index_col=0)  # index must match sample IDs

# --- Parameters ---
geo_r = df_samples['site'].unique()  # column name for population
min_call_rate = 2/3
out_dir = "./"
os.makedirs(out_dir, exist_ok=True)

# --- Iterate through all pairwise population comparisons ---
for pop1, pop2 in itertools.combinations(geo_r, 2):
    print(f"\nComparing {pop1} vs {pop2}")

    # --- Map population sample IDs to integer positions in Zarr ---
    subpops = {}
    for pop in [pop1, pop2]:
        sample_ids = df_samples.loc[df_samples['site'] == pop].index.to_list()
        positions = np.where(np.isin(samples, sample_ids))[0]
        subpops[pop] = positions

    n_samples_pop1 = len(subpops[pop1])
    n_samples_pop2 = len(subpops[pop2])
    print(f"{pop1}: {n_samples_pop1} samples, {pop2}: {n_samples_pop2} samples")

    # --- Compute allele counts ---
    acs = genotype_all.count_alleles_subpops(subpops)

    # --- Compute missingness per variant ---
    genotypes_pop1 = genotype_all[:, subpops[pop1]]
    genotypes_pop2 = genotype_all[:, subpops[pop2]]

    call_rate_pop1 = genotypes_pop1.is_called().mean(axis=1)
    call_rate_pop2 = genotypes_pop2.is_called().mean(axis=1)

    flt_call_both = (call_rate_pop1 >= min_call_rate) & (call_rate_pop2 >= min_call_rate)

    # --- Filter for biallelic, segregating SNPs ---
    acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
    flt_biallelic = acu.is_segregating() & (acu.max_allele() == 1)

    # --- Combine filters ---
    flt = flt_biallelic & flt_call_both
    n_retained = np.count_nonzero(flt)
    print(f"Retaining {n_retained} SNPs after filters")

    if n_retained == 0:
        print("No SNPs pass filters, skipping this population pair")
        continue

    # --- Apply filters to allele counts ---
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])

    # --- Compute Hudson FST ---
    num, den = allel.hudson_fst(ac1, ac2)
    snp_fst_hudson = num / den
    fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
    print(f"Hudson FST: {fst_hudson:.4f} ± {se_hudson:.4f}")

    # --- Create output DataFrames ---
    sites = f"{pop1}_vs_{pop2}"

    # Per-SNP FST
    snp_fst_df = pd.DataFrame({
        'CHR': callset['variants/CHROM'][:][flt],
        'POS': callset['variants/POS'][:][flt],
        'Fst': snp_fst_hudson
    })

    # Average Hudson FST
    fst_hudson_df = pd.DataFrame({
        'Site': [sites],
        'Hudson_Fst': [fst_hudson],
        'SE': [se_hudson],
        'n_retained': [n_retained]
    })

    # --- Save CSVs ---
    f_name1 = os.path.join(out_dir, f"{sites}.paired.snp_fst_df.csv")
    f_name2 = os.path.join(out_dir, f"{sites}.paired.hudson_fst_df.csv")
    snp_fst_df.to_csv(f_name1, index=False)
    fst_hudson_df.to_csv(f_name2, index=False)
    print(f"Saved: {f_name1}")
    print(f"Saved: {f_name2}")
# --- Loop over each focal population ---
for pop1 in geo_r:
    print(f"\nComparing {pop1} vs all other populations")

    # --- Define subpopulations ---
    # Focal group
    sample_ids_pop1 = df_samples.loc[df_samples['site'] == pop1].index.to_list()
    idx_pop1 = np.where(np.isin(samples, sample_ids_pop1))[0]

    # All other groups
    sample_ids_others = df_samples.loc[df_samples['site'] != pop1].index.to_list()
    idx_others = np.where(np.isin(samples, sample_ids_others))[0]

    n1, n2 = len(idx_pop1), len(idx_others)
    print(f"{pop1}: {n1} samples, All others: {n2} samples")

    subpops = {pop1: idx_pop1, "Others": idx_others}

    # --- Compute allele counts ---
    acs = genotype_all.count_alleles_subpops(subpops)

    # --- Compute missingness per variant ---
    genotypes_pop1 = genotype_all[:, idx_pop1]
    genotypes_others = genotype_all[:, idx_others]
    call_rate_pop1 = genotypes_pop1.is_called().mean(axis=1)
    call_rate_others = genotypes_others.is_called().mean(axis=1)
    flt_call_both = (call_rate_pop1 >= min_call_rate) & (call_rate_others >= min_call_rate)

    # --- Filter for biallelic segregating SNPs ---
    acu = allel.AlleleCountsArray(acs[pop1][:] + acs["Others"][:])
    flt_biallelic = acu.is_segregating() & (acu.max_allele() == 1)

    # --- Combine filters ---
    flt = flt_biallelic & flt_call_both
    n_retained = np.count_nonzero(flt)
    print(f"Retaining {n_retained} SNPs after filters")

    if n_retained == 0:
        print("No SNPs pass filters, skipping this population")
        continue

    # --- Apply filters ---
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs["Others"].compress(flt, axis=0)[:, :2])

    # --- Compute Hudson FST ---
    num, den = allel.hudson_fst(ac1, ac2)
    snp_fst_hudson = num / den
    fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
    print(f"{pop1} vs Others: FST = {fst_hudson:.4f} ± {se_hudson:.4f}")

    # --- Create per-SNP output ---
    snp_fst_df = pd.DataFrame({
        'CHR': callset['variants/CHROM'][:][flt],
        'POS': callset['variants/POS'][:][flt],
        'Fst': snp_fst_hudson
    })

    # --- Save per-pop results ---
    f_name_snp = os.path.join(out_dir, f"{pop1}_vs_all_others.snp_fst.csv")
    snp_fst_df.to_csv(f_name_snp, index=False)

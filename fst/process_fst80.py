#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
import itertools
import pandas as pd
import numpy as np


### Functions
def clean_data(metadata):
  meta = pd.read_csv(metadata, sep=',')
  return(meta)


def site_fst_pair_analysis(meta, pop, out_dir, thresh):
  geo_r = meta[pop].unique()
  geo_r  = [x for x in geo_r if str(x) != 'nan'] 
  high_fst_snps =[]
  fst_site =[]
  for comb in itertools.combinations(geo_r,2): #check combs vs permutations
    pop1 = comb[0]
    pop2 = comb[1]
    filename1 = out_dir + "/" + pop1 + "_vs_" + comb[1] + ".paired.snp_fst_df.csv"
    filename2 = out_dir + "/" + pop2 + "_vs_" + comb[0] + ".paired.snp_fst_df.csv"
    try:
        site_fst = pd.read_csv(filename1, header=0)
    except FileNotFoundError:
        try:
            site_fst = pd.read_csv(filename2, header=0)
        except FileNotFoundError:
            print(f"⚠️ No file found for {pop1} vs {pop2}, skipping...")
            continue
    #site_fst = pd.read_csv(filename,header=0)
    site_fst["Fst"] = pd.to_numeric(site_fst["Fst"], errors="coerce").fillna(0).clip(lower=0)

    site_fst.loc[site_fst['Fst'] < 0, 'Fst'] = 0
    mean_val = np.mean(site_fst['Fst'])

    fst_80 = site_fst.loc[site_fst['Fst'] > thresh].copy()  

    fst_80["Population 1"] = pop1
    fst_80["Population 2"] = pop2

    num_80 = len(fst_80)
    fst_site.append(pd.DataFrame({'Population 1': pop1, 'Population 2': pop2, 'Fst>threshold':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_site_df = pd.concat(fst_site, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_site_df, high_fst_snps_df)
  
def site_fst_one_analysis(meta, pop, out_dir, thresh):
  geo_r = meta[pop].unique()
  geo_r  = [x for x in geo_r if str(x) != 'nan']
  high_fst_snps =[]
  fst_site =[]  
  for site in geo_r:
    filename = out_dir + "/" + site + "_vs_all" + "_others.snp_fst.csv"
    site_fst = pd.read_csv(filename,header=0)
    site_fst["Fst"] = pd.to_numeric(site_fst["Fst"], errors="coerce").fillna(0).clip(lower=0)
    fst_80 = site_fst[site_fst['Fst']>thresh]
    fst_80["site"] = site
    num_80 = len(fst_80)
    fst_site.append(pd.DataFrame({'site': site,  'Fst>thresh':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_site_df = pd.concat(fst_site, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_site_df, high_fst_snps_df)


def main(args):
  #params = {"metadata": args.metadata}
  meta_clean = clean_data(args.metadata)
  if "paired" or "both" in args.comparison:
    fst_site_df1, high_fst_snps_df1 = site_fst_pair_analysis(meta_clean, args.pop, args.out_dir, args.thresh)
    fst_site_df1.to_csv(args.out_dir + "/" + "fst_site_pairs.csv")
    high_fst_snps_df1.to_csv(args.out_dir + "/" + "high_fst_snps_site_pairs.csv")
  if "one_against_all" or "both" in args.comparison:
    fst_site_df2, high_fst_snps_df2 = site_fst_one_analysis(meta_clean, args.pop, args.out_dir, args.thresh)
    fst_site_df2.to_csv(args.out_dir + "/" + "fst_site_single.csv")
    high_fst_snps_df2.to_csv(args.out_dir + "/" + "high_fst_snps_site_single.csv")

parser = argparse.ArgumentParser(description='Fst statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='csv file containing metadata',required=True)
parser.add_argument('--comparison',help='paired or one_against_all or both',required=True)
parser.add_argument('--thresh',type=float, default=80, help='paired or one_against_all or both',required=True)
parser.add_argument('--out_dir',help='out_dir',required=True)
parser.add_argument('--pop',help='population column',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

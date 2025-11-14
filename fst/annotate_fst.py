#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
import itertools
import pandas as pd
import fastq2matrix as fm


### Functions
def annotate_fst(high_fst_file, anno_file):
  fst = pd.read_csv(high_fst_file, header=0)
  anno = pd.read_table(anno_file, header=0)
  fst = fst.rename(columns={'CHR': 'chr', 'POS': 'pos'})
  fst_ann = fst.merge(anno, how='left')
  return(fst_ann)

def get_product(anno, product_file):
  prod = pd.read_table(product_file, header=0)
  prod = prod.rename(columns={'gene_id': 'gene'})
  fst_prod = anno.merge(prod, how='left')
  return(fst_prod)
  
def main(args):
  fst_ann = annotate_fst(args.high_fst_file, args.anno_file)
  fst_ann.to_csv(args.output_ann)
  fst_prod = get_product(fst_ann, args.product_file)
  fst_prod.to_csv(args.output_prod)

parser = argparse.ArgumentParser(description='Annotate high Fst',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--high_fst_file',help='csv file containing snps with high fst',required=True)
parser.add_argument('--anno_file',help='tsv file containing annotations',required=True)
parser.add_argument('--product_file',help='tsv file containing products',required=True)
parser.add_argument('--output_ann',help='name of output file',required=True)
parser.add_argument('--output_prod',help='name of output file',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from tqdm import tqdm
import time
import os.path
import pandas as pd
import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
print('scikit-allel', allel.__version__)

vcf = "../072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
out_zarr = "072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.zarr"
allel.vcf_to_zarr(vcf, out_zarr)

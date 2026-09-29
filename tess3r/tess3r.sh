#!/bin/bash

install.packages("devtools")
devtools::install_github("bcm-uga/TESS3_encho_sen")

Rscript prepare_tess3r_files.R

xargs -I {} Rscript run_tess3r_general.R -b Pf_mgen_ethiopia_matbinT.mis.rds -c Pf_mgen_ethiopia_coordsM.rds --prefix Pf_mgen_ethiopia --K {} --rep 50 --threads 10 < K_runs.txt

Rscript plot_tess3r.R


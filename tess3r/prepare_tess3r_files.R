library(tess3r)
library(readr)
library(dplyr)

wd <- "/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen"

# METADATA
# Make GPS coordinates uniform across site/location
# Keep columns numeric
metadata <- readr::read_csv(file.path(wd, "/tess3r/run2/metadata_db_tess3r_022026.csv"), col_names = T) %>%
    select(sample_id, latitude, longitude, Site, country) %>%
     mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# MATRIX
# Binary matrix
# Extract header and save sample order
matbin <- data.table::fread(file.path(wd, "/tess3r/run2/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.tess3r.miss0.2.vqslod.filt.snps.filt.MAF0.001.mat.bin"))
matbin_s <- matbin[ , -c(1:3)]
samples <- colnames(matbin_s)

# Check overlap
# Sample save
metadata <- metadata %>% filter(sample_id %in% samples)
matbin_s <- matbin_s %>% select(metadata$sample_id)
saveRDS(colnames(matbin_s), file.path(wd, "/tess3r/run2/Pf_fulldb_022026.rds"))

# Matrix save
# Transposed recoded mixed calls to alternative - 0.5 to 1
matbin_s[matbin_s == 0.5] <- 1
all(colnames(matbin_s) == metadata$sample_id)
matbinT <- t(matbin_s)
saveRDS(matbinT, file.path(wd, "/tess3r/run2/Pf_fulldb_022026_matbinT.mis.rds"))

# Coordinates
# Save as matrix
coords <- metadata %>% select(longitude, latitude) %>%
      mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude))
colnames(coords) <- c("V1", "V2")
coordsM <- data.matrix(coords, rownames.force=FALSE)
saveRDS(coordsM, file.path(wd, "/tess3r/run2/Pf_fulldb_022026_coordsM.rds"))

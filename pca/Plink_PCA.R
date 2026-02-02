library(showtext)
library(argparse)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()

# Filter out multiclonal samples
# awk -F',' 'NR==1 {for(i=1;i<=NF;i++) {if($i=="fws") fws_col=i; if($i=="wgs_id") id_col=i}} NR>1 && $fws_col > 0.95 {print $id_col}' metadata.csv > samples_0.95.csv

# bcftools view -S samples_0.95.csv -o allsamples_fws_0.95.vcf.gz -O z /mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz

# # Index the output
# bcftools index allsamples_fws_0.95.vcf.gz

# Run plink
# plink --vcf allsamples_fws_0.95.vcf.gz --distance square --double-id --allow-extra-chr --out allsamples_fws

# Run this script
# Rscript Plink_PCA.R --workdir . \
                    #--prefix Dec2025 \
                    #--metadata metadata.csv

# Parse command-line arguments
parser <- ArgumentParser(description = 'Run PCA and tree construction from PLINK distance matrix')
parser$add_argument('--workdir', required=TRUE, help='Working directory containing PLINK distance files')
parser$add_argument('--prefix', required=TRUE, help='Prefix for PLINK .dist and .dist.id files')
parser$add_argument('--metadata', required=TRUE, help='Metadata file path (tab-separated)')

args <- parser$parse_args()

workdir <- args$workdir
prefix <- args$prefix
metadata <- args$metadata
color_by <- args$color_by

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "wgs_id"))
# desc <- id %>% left_join(met, by = c("V1" = "sample_id"))

# Columns you rely on for colouring
cols_check <- c("Country", "Year", "Region")

# Any missing in any of those columns?
na_mask <- apply(desc[cols_check], 1, function(r) any(is.na(r)))
na_samples <- desc$V1[na_mask]

if (length(na_samples)) {
  message("Samples with missing metadata (any of Country/Year/Region): ",
          paste(na_samples, collapse = ", "))
  write.table(
    na_samples,
    file = file.path(workdir, paste0(prefix, ".missing_metadata_wgs_ids.txt")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA ##
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling - might take a while
# saveRDS(cmd, paste0(prefix, ".dist.rds") # save to RDS format
#cmd <- readRDS(file.path(workdir, paste0(prefix, ".dist.rds"))
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Export dist_m to .newick to make neighbour joining tree
tree <- nj(dist_m)
write.tree(phy = tree, file = file.path(workdir, paste0(prefix, ".newick")))



# Overlay COUNTRY info
df <- as.data.frame(cmd$points, stringsAsFactors = F)
df$Country <- gsub("_", " ", desc$Country)
#df$region <- gsub("_", " ", desc$region)
#df$sample_id <- rownames(matrix)
colnames(df) <- gsub("V", "PC", colnames(df))

color_by <- "Country" # specify if coloured by region or country
country_colors <- c(
  "Ethiopia" = "#E41A1C", 
  "Kenya" = "#377EB8", 
  "Tanzania" = "#4DAF4A",
  "Uganda" = "#984EA3",   
  "Sudan" = "#FF7F00",  
  "Ghana" = "#FFFF33", 
  "Nigeria" = "#A65628",
  "Mali" = "#eb88d4",
  "Eritrea" = "darkblue",
  "India" = "black",
  "Thailand" = "darkgreen",
  "Myanmar" = "lightsalmon"
)

# Graph with PC1 an PC2
p <- ggplot(data = df, aes(x = PC1, y = PC2,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC2", " (", vars["PC2"], "%)")) +
    theme_test() +
    scale_color_manual(values = country_colors, name = "Countries") +
    theme(legend.position = "bottom")
ggsave(plot = p, filename="figs/PC_12_country.pdf")


# Graph with PC1 and PC3
q <- ggplot(data = df, aes(x = PC1, y = PC3,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC3", " (", vars["PC3"], "%)")) +
    theme_test() +
    scale_color_manual(values = country_colors, name = "Countries") +
    theme(legend.position = "bottom")
ggsave(plot = q, filename = "figs/PC_13_country.pdf")


# Overlay YEAR info
df <- as.data.frame(cmd$points, stringsAsFactors = F)
df$Year <- gsub("_", " ", desc$Year)
#df$region <- gsub("_", " ", desc$region)
#df$sample_id <- rownames(matrix)
colnames(df) <- gsub("V", "PC", colnames(df))

color_by <- "Year" # specify if coloured by region or country

# Graph with PC1 an PC2
p <- ggplot(data = df, aes(x = PC1, y = PC2,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC2", " (", vars["PC2"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = p, filename="figs/PC_12_year.pdf")


# Graph with PC1 and PC3
q <- ggplot(data = df, aes(x = PC1, y = PC3,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC3", " (", vars["PC3"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = q, filename = "figs/PC_13_year.pdf")



# Overlay REGION info
df <- as.data.frame(cmd$points, stringsAsFactors = F)
df$Region <- gsub("_", " ", desc$Region)
#df$region <- gsub("_", " ", desc$region)
#df$sample_id <- rownames(matrix)
colnames(df) <- gsub("V", "PC", colnames(df))

color_by <- "Region" # specify if coloured by region or country

# Graph with PC1 an PC2
p <- ggplot(data = df, aes(x = PC1, y = PC2,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC2", " (", vars["PC2"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = p, filename="figs/PC_12_region.pdf")


# Graph with PC1 and PC3
q <- ggplot(data = df, aes(x = PC1, y = PC3,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC3", " (", vars["PC3"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = q, filename = "figs/PC_13_region.pdf")


# Overlay fws info
df <- as.data.frame(cmd$points, stringsAsFactors = F)
df$fws <- as.numeric(desc$fws)
colnames(df) <- gsub("V", "PC", colnames(df))

color_by <- "fws" 

# Graph with PC1 an PC2
p <- ggplot(data = df, aes(x = PC1, y = PC2,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    scale_color_viridis_c(name = "Fws", limits = c(0, 1)) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC2", " (", vars["PC2"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = p, filename="figs/PC_12_fws.pdf")


# Graph with PC1 and PC3
q <- ggplot(data = df, aes(x = PC1, y = PC3,
       color = !!sym(color_by))) +
    geom_point(size=2) +
    scale_color_viridis_c(name = "Fws", limits = c(0, 1)) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC3", " (", vars["PC3"], "%)")) +
    theme_test() +
    theme(legend.position = "bottom")
ggsave(plot = q, filename = "figs/PC_13_fws.pdf")


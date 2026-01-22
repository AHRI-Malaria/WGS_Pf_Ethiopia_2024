# Load required libraries
library(data.table)
library(dplyr)
library(parallel)
library(stringr)
library(tidyr)

# Set working directory and number of cores
setwd("/mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/interaction_analysis")
n_cores <- 16

# Load genotype matrix
geno <- fread("/mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.mat.bin", sep = "\t", header = TRUE, data.table = FALSE)

# Mask 0.5 with NA
geno[geno == 0.5] <- NA

# Load PC covariates
pcs <- fread("/mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/pca/Sep2025_pca_coordinates.csv")
colnames(pcs) <- c("sample", "PC1", "PC2", "PC3", "PC4","PC5")

# Keep only samples present in both genotype and PCs
common_samples <- intersect(colnames(geno)[4:ncol(geno)], pcs$sample)
geno <- geno[, c("chr", "pos", "ref", common_samples)]
pcs <- pcs %>% filter(sample %in% common_samples)
#write.csv(geno,"geno.csv",row.names=TRUE)

# Extract kelch13 mutation row (e.g., by annotation)
kelch_row <- which(geno$chr == "Pf3D7_13_v3" & geno$pos == 1725133)
kelch_geno <- as.numeric(geno[kelch_row, 4:ncol(geno)])

# Filter genotype matrix to exclude kelch13 itself
geno_test <- geno 

# Parallel association test
geno_rows <- split(geno_test, seq(nrow(geno_test)))


# With PC components
results <- lapply(geno_rows, function(row) {
  snp2 <- as.numeric(row[1, 4:ncol(row)])
  interaction_term <- kelch_geno * snp2
  if (sum(interaction_term, na.rm = TRUE) > 1) {
    print("passed interaction")
    chr <- row["chr"]
    pos <- row["pos"]
    ref <- row["ref"]
    df <- data.frame(
      y = kelch_geno,
      snp2 = snp2,
      PC1 = pcs$PC1,
      PC2 = pcs$PC2,
      PC3 = pcs$PC3,
      PC4 = pcs$PC4,
      PC5 = pcs$PC5
    )
    df <- df[complete.cases(df), ]
    if (nrow(df) > 10) {
      print("exceeds sample number")
      
      # Check if snp2 has variation
      if (length(unique(df$snp2[!is.na(df$snp2)])) < 2) {
      print("snp2 has no variation, skipping")
      return(NULL)
      }
      
      # Fit model with error handling
      fit <- tryCatch({
        glm(y ~ snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = df, family = "binomial") 
      }, error = function(e) {
        print(paste("Model fitting error:", e$message))
        return(NULL)
      })
      
      if (is.null(fit)) return(NULL)
      
      coef_table <- coef(summary(fit))
      
      # Check if snp2 coefficient exists
      if (!"snp2" %in% rownames(coef_table)) {
        print("snp2 coefficient not in model, skipping")
        return(NULL)
      }
      
      return(data.frame(
        chr = chr,
        pos = pos,
        ref = ref,
        beta_snp2 = coef_table["snp2", "Estimate"],
        p_snp2 = coef_table["snp2", "Pr(>|z|)"],
        n = nrow(df)
      ))
    }
  }
  return(NULL)
})


# Combine and save results
results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_snp2, method="BH")
write.csv(results_df, "kelch13_interaction_results.csv", row.names = FALSE)


# Run bcftools command from R
system("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%BCSQ\\n' /mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz > annotations.tsv")



make_anno <- function(annotations) {
  anno <- read.table(annotations, sep = "\t", header = FALSE)
  
  # Split the 5th column by '|' and extract gene names
  anno_split <- anno %>%
    separate(V5, into = c("field1", "gene_name1", "gene_name2"), 
             sep = "\\|", extra = "drop", fill = "right")
  
  # Create combined gene name
  anno_final <- anno_split %>%
    mutate(Gene_name = paste0(gene_name1, " (", gene_name2, ")")) %>%
    select(V1, V2, V3, V4, Gene_name)
  
  # Rename columns
  colnames(anno_final) <- c("chr", "pos", "ref", "alt", "Gene_name")
  
  # Write to file
  write.table(anno_final, "annotations_final.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Run the function
make_anno("annotations.tsv")

# R analysis
annotations_final <- read.csv("annotations_final.tsv", sep = "\t")
results_df_ann <- left_join(results_df, annotations_final, 
                                        by = c("chr" = "chr", "pos" = "pos"))
write.csv(results_df_ann, "kelch13_interaction_results_ann.csv", row.names = FALSE)





#Without PC components
# Load genotype matrix
geno <- fread("/mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.mat.bin", sep = "\t", header = TRUE, data.table = FALSE)

# Mask 0.5 with NA
geno[geno == 0.5] <- NA

# Extract kelch13 mutation row (e.g., by annotation)
kelch_row <- which(geno$chr == "Pf3D7_13_v3" & geno$pos == 1725133)
kelch_geno <- as.numeric(geno[kelch_row, 4:ncol(geno)])

# Filter genotype matrix to exclude kelch13 itself
geno_test <- geno 

# Parallel association test
geno_rows <- split(geno_test, seq(nrow(geno_test)))

results <- lapply(geno_rows, function(row) {
  snp2 <- as.numeric(row[1, 4:ncol(row)])
  interaction_term <- kelch_geno * snp2
  if (sum(interaction_term, na.rm = TRUE) > 1) {
    print("passed interaction")
    chr <- row["chr"]
    pos <- row["pos"]
    ref <- row["ref"]
    df <- data.frame(
      y = kelch_geno,
      snp2 = snp2#,
      #PC1 = pcs$PC1,
      #PC2 = pcs$PC2,
      #PC3 = pcs$PC3,
      #PC4 = pcs$PC4,
      #PC5 = pcs$PC5
    )
    df <- df[complete.cases(df), ] #removes rows with any NA values
    if (nrow(df) > 10) {
      print("exceeds sample number")
      
      # Check if snp2 has variation
      if (length(unique(df$snp2[!is.na(df$snp2)])) < 2) {
      print("snp2 has no variation, skipping")
      return(NULL)
      }
      
      # Fit model with error handling
      fit <- tryCatch({
        glm(y ~ snp2, data = df, family = "binomial") # + PC1 + PC2 + PC3 + PC4 + PC5
      }, error = function(e) {
        print(paste("Model fitting error:", e$message))
        return(NULL)
      })
      
      if (is.null(fit)) return(NULL)
      
      coef_table <- coef(summary(fit))
      
      # Check if snp2 coefficient exists
      if (!"snp2" %in% rownames(coef_table)) {
        print("snp2 coefficient not in model, skipping")
        return(NULL)
      }
      
      return(data.frame(
        chr = chr,
        pos = pos,
        ref = ref,
        beta_snp2 = coef_table["snp2", "Estimate"],
        p_snp2 = coef_table["snp2", "Pr(>|z|)"],
        n = nrow(df)
      ))
    }
  }
  return(NULL)
})

# Combine and save results
results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_snp2, method="BH")
write.csv(results_df, "kelch13_interaction_results_noPC.csv", row.names = FALSE)


# Run bcftools command from R
system("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%BCSQ\\n' /mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz > annotations.tsv")



make_anno <- function(annotations) {
  anno <- read.table(annotations, sep = "\t", header = FALSE)
  
  # Split the 5th column by '|' and extract gene names
  anno_split <- anno %>%
    separate(V5, into = c("field1", "gene_name1", "gene_name2"), 
             sep = "\\|", extra = "drop", fill = "right")
  
  # Create combined gene name
  anno_final <- anno_split %>%
    mutate(Gene_name = paste0(gene_name1, " (", gene_name2, ")")) %>%
    select(V1, V2, V3, V4, Gene_name)
  
  # Rename columns
  colnames(anno_final) <- c("chr", "pos", "ref", "alt", "Gene_name")
  
  # Write to file
  write.table(anno_final, "annotations_final.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Run the function
make_anno("annotations.tsv")

# R analysis
annotations_final <- read.csv("annotations_final.tsv", sep = "\t")
results_df_ann <- left_join(results_df, annotations_final, 
                                        by = c("chr" = "chr", "pos" = "pos"))
write.csv(results_df_ann, "kelch13_interaction_results_ann_noPC.csv", row.names = FALSE)
#options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("optparse") 
library(tidyverse)
library(isoRelate)
library(SeqArray)
library(moimix)
library(ggplot2)
library(data.table)
library(readr)
library(optparse)

#Rscript test.R ../072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps_ethiopia.vcf.gz  ../isorelate_metadata.csv ibd_proportion_metadata.csv test

#cat("Enter the path to the VCF file:", arg1, "\n")
#cat("isorelate_metadata.csv:", arg2, "\n")
#cat("Enter ibd_proportion_metadata.csv:", arg3, "\n")
#cat("prefix:", arg4, "\n")
#Enter the path to the VCF file
#Pf241_samples <- args[1]





option_list <- list(

  make_option(
    c( "--fws_metadata"),
    type = "character",
    help = "METADATA FILE WHICH CONTAIN THE FWS RESULT(CSV/TSV)",
    metavar = "FILE"
  ),

  make_option(
    c("-v", "--vcf"),
    type = "character",
    help = "VCF file containing variants",
    metavar = "FILE"
  ),

  make_option(
    c("-p", "--prefix"),
    type = "character",
    help = "prefix for the output files",
    metavar = "text"
  ),

  make_option(
    c( "--ir_metadata"),
    type = "character",
    help = "IR metadata which contain the iid , fid and group ,the fid needs to be similar with the previous fws_metadata - in .csv file",
    metavar = "FILE"
  ),

  make_option(
    c( "--annotation"),
    type = "character",
    help = "annotation .tsv file generated from the snpeff output which contain the variants with gene name",
    metavar = "FILE"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

Pf241_samples <- opt$vcf
fws_metadata <- opt$fws_metadata
prefix <- opt$prefix
metadata <- opt$ir_metadata
annotation <- opt$annotation 


#enter the path to the metadata file
#fws_metadata <- args[2] ~ isorelate_metadata.csv
fws_metadata <- read_delim(fws_metadata,delim=NULL,trim_ws=TRUE)%>%
	mutate(moi = ifelse(fws >= 0.95 , 1, 2))
  
vcf_header <- seqVCF_Header(Pf241_samples)
# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "EFF")
# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

Pf241_samples_gds<-seqVCF2GDS(Pf241_samples,
                              "gds_file.gds",genotype.var.name="GT",
                              header=vcf_header, info.import="ANN",
                              fmt.import=format.import)
# Use isolates file for downstream analysis 
isolates241samples_coregenome <- seqOpen(Pf241_samples_gds) 

ethiopia_ped_map_file_241samples<-extractPED(isolates241samples_coregenome, 
                                           moi.estimates = fws_metadata$moi, 
                                           outfile = NULL) 

#ethiopia_ped_map_file_241samples$map$genetic_distance <-ethiopia_ped_map_file_241samples$map$pos/13500

ethiopia_ped_map_file_241samples$ped <- ethiopia_ped_map_file_241samples$ped %>% 
  select(-famID) %>% 
  mutate(famID = case_when(
    grepl(pattern = "AC" , indID) ~ "Gambella",
    grepl(pattern = "GM", indID) ~ "Gonder",  #such kind of mis spelling will result in error for ir statistics so it would be great to change it to Gondar or similar with metadata
    grepl(pattern = "MM", indID) ~ "Mizan",
    grepl(pattern = "WG24AC", indID) ~ "Gambella",
    grepl(pattern = "WG24ALHC", indID) ~ "Arbaminch",
    grepl(pattern = "WG24GMHC", indID) ~ "Gonder",
    grepl(pattern = "WG24MFHC", indID) ~ "Metehera Fentale",
    grepl(pattern = "WG24MHC", indID) ~ "Mizan"


  ))
  

#to move the famID to the first column

ethiopia_ped_map_file_241samples$ped <- ethiopia_ped_map_file_241samples$ped %>% 
  select(famID , everything())

saveRDS(ethiopia_ped_map_file_241samples, file = paste0(prefix,"wgs_ped_map_file.rds"))

#isorelate_script2.R
#ethiopia_ped_map_file_241samples <- readRDS("wgs_ped_map_file.rds")
samples241_final_genotypes<-
  getGenotypes(ped.map = ethiopia_ped_map_file_241samples,
  reference.ped.map = NULL,
  maf = 0.02, # you can modify the threshold based on you data type
  isolate.max.missing = 0.2,# you can modify the threshold based on you data type
  snp.max.missing = 0.1,# you can modify the threshold based on you data type
  chromosomes = NULL,
  input.map.distance = "cM",
  reference.map.distance = "cM")

  samples241_final_parameters <- getIBDparameters(ped.genotypes = samples241_final_genotypes,
                                                number.cores = 8) 
  #saveRDS(samples241_final_genotypes, file = paste0(prefix,"wgs_final_genotypes.rds"))
  #write.csv(samples241_final_parameters, paste0(prefix,"wgs_final_parameters.csv"),row.names=F) 

  samples241_final_ibd_segments<-
  getIBDsegments(ped.genotypes= samples241_final_genotypes,                    
   parameters = samples241_final_parameters, 
   number.cores = 8, #The number of cores used for parallel execution
   minimum.snps = 20, #the minimum number of SNPs in an IBD segment for it to be reported
   minimum.length.bp = 50000,# The minimum length of a reported IBD segment
   error = 0.001) #The genotyping error rate
#write.csv(samples241_final_ibd_segments, file = paste0(prefix,"wgs_final_ibd_segments.csv"),row.names=F)
pdf(paste0(prefix,"wgs_final_ibd_segments.pdf"), width = 12, height = 10)
plotIBDsegments(ped.genotypes = samples241_final_genotypes , ibd.segments = samples241_final_ibd_segments ,    interval = NULL,
                annotation.genes = NULL,
                annotation.genes.color = NULL,
                highlight.genes = NULL, 
                highlight.genes.labels = TRUE,
                highlight.genes.color = NULL,
                highlight.genes.alpha = 0.1,
                segment.height = 0.6,
                number.per.page = NULL, 
                fid.label = TRUE, 
                iid.label = FALSE, 
                ylabel.size = 9, 
                add.rug = FALSE,
                plot.title = "Distribution of IBD segments in PNG", 
                add.legend = TRUE, 
                segment.color = NULL)
dev.off()

#isorelate_script3.R

#samples241_final_genotypes <- readRDS("wgs_final_genotypes.rds")
#samples241_final_ibd_segments <- read.csv("wgs_final_ibd_segments.csv")
samples241_final_genotypes$genotypes <- samples241_final_genotypes$genotypes %>% 
  filter(chr != "Pf3D7_API_v3")%>% 
  filter(chr != "Pf_M76611")
samples241_final_ibd_segments <- samples241_final_ibd_segments %>%
  filter(chr != "Pf3D7_API_v3") %>% 
  filter(chr != "Pf_M76611")

# This metadata is for the ibd proportion calculation
#Enter ibd_proportion_metadata.csv
#metadata <- read.csv(args[3])
#	filter(group == "Gonder")
    #  mutate(group = case_when(
     #   grepl("ERR", fid) ~ "old",
     #   grepl("SRR", fid) ~ "old",
      #  TRUE ~ "new"
   # ))

samples241_final_IBD_matrix<- getIBDmatrix(ped.genotypes = samples241_final_genotypes, 
                                           ibd.segments = samples241_final_ibd_segments)
#samples241_final_proportion <- getIBDproportion(ped.genotypes = samples241_final_genotypes,
#                                               ibd.matrix = samples241_final_IBD_matrix, 
#                                                groups = metadata)
write.csv(samples241_final_IBD_matrix , paste0(prefix,"wgs_final_ibd_matrix.csv"),row.names=F)
write.csv(samples241_final_parameters, paste0(prefix,"wgs_final_parameters.csv"),row.names=F) 
write.csv(samples241_final_ibd_segments, paste0(prefix,"wgs_final_ibd_segments.csv"),row.names=F)
saveRDS(samples241_final_genotypes, paste0(prefix,"wgs_final_genotypes.rds"))


metadata <- read.csv(metadata,header=T) # this file should contain iid,fid ,and group , 
# the other thing is that the fid column (names of sites) should be similar with the fws_metadata or the one created above

y  <- function (ped.genotypes, ibd.matrix, groups) 
{
    if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) 
        stop("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
    if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) 
        stop("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
    ped.genotypes <- samples241_final_genotypes
    
    pedigree <- ped.genotypes[["pedigree"]]
    genotypes <- ped.genotypes[["genotypes"]]
    groups <- groups[,c("fid","iid","group")]
    if (!is.data.frame(pedigree)) 
        stop("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")
    if (!is.data.frame(genotypes)) 
        stop("'ped.genotypes' has incorrect format - 'genotypes' is not a data.frame")
    if (ncol(pedigree) != 6) 
        stop("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
    if (nrow(pedigree) < 3) 
        stop("cannot calculate iR statistic with only 2 isolates")
    if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", 
        "moi", "aff"))) 
        stop("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")
    number.isolates <- nrow(pedigree)
    number.pairs <- (number.isolates * (number.isolates - 1))/2
    if (ncol(genotypes) <= 6) 
        stop("'ped.genotypes' has incorrect format - minimum of 2 isolates required for analysis")
    if (nrow(genotypes) <= 1) 
        stop("'ped.genotypes' has incorrect format - too few SNPs for analysis")
    if (any(colnames(genotypes)[1:5] != c("chr", "snp_id", "pos_M", 
        "pos_bp", "freq"))) 
        stop("'ped.genotypes' has incorrect format - 'genotypes' must have columns labelled: chr, snp_id, pos_M, pos_bp and freq")
    if (ncol(genotypes) != (number.isolates + 5)) 
        stop("'ped.genotypes' has incorrect format - some isolates without genotype data")
    if (ncol(ibd.matrix) < 5) 
        stop("'ibd.matrix' has incorrect format - must have minimum 5 columns")
    if (any(colnames(ibd.matrix)[1:4] != c("chr", "snp_id", "pos_M", 
        "pos_bp"))) 
        stop("'ibd.matrix' has incorrect format - must have first 4 columns labelled: chr, snp_id, pos_M, pos_bp")
    if (ncol(ibd.matrix) != (number.pairs + 4)) 
        stop("'ibd.matrix' has incorrect format - somes pairs are missing")
    if (!is.null(groups)) {
        if (!is.data.frame(groups)) 
            stop("'groups' has incorrect format - must be a data.frame")
        if (ncol(groups) != 3) 
            stop("'groups' has incorrect format - must have 3 columns: fid, iid, group")
        colnames(groups)[1:2] <- c("fid", "iid")
        group.names <- paste(groups[, "fid"], groups[, "iid"],  sep = "/")
        isolate.names <- paste(pedigree[, "fid"], pedigree[, "iid"], sep = "/")
        print(setdiff(isolate.names, group.names) )
        if (!all(isolate.names %in% group.names)) 
            stop("'groups' has incorrect format - some isolates 'ped.genotypes' are missing from 'groups'")
        pedigree.0 <- data.frame(num.id = 1:nrow(pedigree), pedigree)
        pedigree.group <- merge(pedigree.0, groups, by = c("fid", 
            "iid"))
        pedigree.group <- pedigree.group[order(pedigree.group[, 
            "num.id"]), ]
        group.pairs <- groupPairs(as.character(pedigree.group[, 
            8]))
        groups.unique <- as.character(unique(pedigree.group[, 
            8]))
        groups.unique.pairs <- groupPairs(groups.unique)
        if (nrow(groups.unique.pairs) > 0) {
            for (i in 1:nrow(groups.unique.pairs)) {
                change.pair <- which(group.pairs[, 1] == groups.unique.pairs[i, 
                  2] & group.pairs[, 2] == groups.unique.pairs[i, 
                  1])
                group.pairs[change.pair, 1] <- groups.unique.pairs[i, 
                  1]
                group.pairs[change.pair, 2] <- groups.unique.pairs[i, 
                  2]
            }
        }
        npairs.groups <- NULL
        group.pairs.1 <- paste(group.pairs[, 1], group.pairs[, 
            2], sep = "/")
        for (i in unique(group.pairs.1)) {
            id_1 <- unlist(strsplit(i, "/"))[1]
            id_2 <- unlist(strsplit(i, "/"))[2]
            npairs <- dim(group.pairs[group.pairs[, 1] == id_1 & 
                group.pairs[, 2] == id_2, ])[1]
            if (is.null(npairs)) 
                npairs <- 1
            npairs.groups <- c(npairs.groups, npairs)
        }
    }
    #library(isoRelate)
    locus.pairs <- ibd.matrix[, 5:ncol(ibd.matrix)]
    if (!is.null(groups)) {
        locus.prop <- NULL
        locus.pvalue <- NULL
        na.iR <- NULL
        locus.prop.colnames <- NULL
        for (g in 1:length(unique(group.pairs.1))) {
            if (npairs.groups[g] == 1) {
                locus.g <- rep(NA, nrow(ibd.matrix))
            }
            else {
                locus.pairs.g <- locus.pairs[, group.pairs.1 == 
                  unique(group.pairs.1)[g]]
                locus.g <- isoRelate:::iRfunction(locus.pairs.g, genotypes[, 
                  "freq"])
            }
            if (all(is.na(locus.g[, 2]))) 
                na.iR <- c(na.iR, g)
            locus.prop <- cbind(locus.prop, locus.g[, 1])
            locus.pvalue <- cbind(locus.pvalue, locus.g[, 2])
            locus.prop.colnames <- c(locus.prop.colnames, unique(group.pairs.1)[g])
        }
        if (length(na.iR) > 0) 
            warning("NA values generated for: ", paste(unique(group.pairs.1)[na.iR], 
                collapse = ", "))
        if (length(unique(group.pairs.1)) != 1) {
            colnames(locus.prop) <- unique(group.pairs.1)
            colnames(locus.pvalue) <- unique(group.pairs.1)
        }
    }
    else {
        locus.prop <- isoRelate:::iRfunction(locus.pairs, genotypes[, "freq"])
        if (all(is.na(locus.prop[, 2]))) 
            warning("NA values generated")
        locus.pvalue <- locus.prop[, 2]
        locus.prop <- locus.prop[, 1]
    }
    return.locus.prop <- cbind(ibd.matrix[, 1:4], pop = rep(1, 
        nrow(ibd.matrix)), locus.prop)
    return.locus.pvalue <- cbind(ibd.matrix[, 1:4], pop = rep(1, 
        nrow(ibd.matrix)), locus.pvalue)
    return.locus.prop <- as.data.table(return.locus.prop)
    return.locus.pvalue <- as.data.table(return.locus.pvalue)
    if (!is.null(groups) & ncol(return.locus.prop) > 6) {
        return.locus.prop <- data.table::melt(return.locus.prop, 
            id.vars = colnames(return.locus.prop)[1:5])
        return.locus.pvalue <- data.table::melt(return.locus.pvalue, 
            id.vars = colnames(return.locus.pvalue)[1:5])
        return.locus.prop <- data.frame(return.locus.prop, return.locus.pvalue[, 
            7])
    }
    else {
        return.locus.prop <- data.frame(ibd.matrix[, 1:4], pop = rep(1, 
            nrow(ibd.matrix)), subpop = 1, locus.prop, locus.pvalue)
    }
    colnames(return.locus.prop) <- c("chr", "snp_id", "pos_M", 
        "pos_bp", "pop", "subpop", "iR", "log10_pvalue")
    return(data.frame(return.locus.prop))
}

ibd <- y(ped.genotypes =samples241_final_genotypes,
                                    ibd.matrix =  samples241_final_IBD_matrix, 
                                    groups = metadata)
write.csv(ibd, file = paste0(prefix,"wgs_final_IR.csv"), row.names = FALSE)
print("iR calculation completed and saved to wgs_final_IR.csv")

jpeg(paste0(prefix,"wgs_final_IR_plot.jpeg"), width = 1200, height = 800)
plotIBDiR(ibd.iR = ibd, 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels = FALSE,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.1,
          point.size = 1,
          point.color = NULL,
          add.rug = FALSE, 
          plot.title = "Significance of IBD sharing", 
          add.legend = FALSE,
          facet.label = TRUE, 
          facet.scales = "fixed")
dev.off()
print("ir plot is completed and saved successfully")

ir_file <- ibd
#annotation <- read.csv("../../annotation/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.ann.wide.tsv", header = TRUE, stringsAsFactors = FALSE,sep="")
annotation <- read.csv(annotation, header = TRUE, stringsAsFactors = FALSE,sep="")
print(colnames(annotation))
#colnames(annotation) <- c("chrom","pos","ref","alt","gene")
ir_annotated <- ir_file %>%
  left_join(annotation, by = c("pos_bp" = "pos","chr"="chr"))
write.csv(ir_annotated, file = paste0(prefix,"wgs_final_IR_annotated_snpeff.csv"), row.names = FALSE)  
print("succesfuully annotated")
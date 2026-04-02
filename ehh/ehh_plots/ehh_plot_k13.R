# the command can be run 
# Rscript ehh_plot_k13.R /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/hap_chr13_Ethiopia /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/snp.info.inp.Ethiopia
# Rscript ehh_plot_k13.R /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/old_vs_new_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/hap_chr13_new /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/old_vs_new_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/snp.info.inp.new
library(R.utils)
#install.packages("vcfR" , "data.table")
library(vcfR)
library(data.table)
library(qqman)
library(rehh)
library(R.utils)
args = commandArgs(trailingOnly=TRUE)
pop_name <- sub("hap_chr[0-9]+_", "", basename(args[1]))
 hap_chrhap_chrhap_chr_pop <- data2haplohh(hap_file= args[1] ,map_file = args[2],
                                recode.allele = FALSE,
                                chr.name = 13, min_perc_geno.hap = 0.8 , min_perc_geno.mrk = 0.8)
  
  res_k13 <- calc_ehh(hap_chrhap_chrhap_chr_pop, 
                mrk = "13_1725133", 
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE , polarized = TRUE   )

  res_675 <- calc_ehh(hap_chrhap_chrhap_chr_pop,
                mrk = "13_1724974",
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE , polarized = TRUE   )
pdf(paste0("k13_ehh_plot_675_", pop_name, ".pdf"), width = 7, height = 7)
plot(res_675, xlim = c(1.6E+6, 1.85E+6) )
dev.off()

pdf(paste0("k13_ehh_plot_R622I_", pop_name, ".pdf"), width = 7, height = 7)
plot(res_k13, xlim = c(1.4E+6, 2.5E+6) )
dev.off()

pdf(paste0("k13_furcation_plot_R622I_", pop_name, ".pdf"),
    width = 7, height = 7)
f <- calc_furcation(hap_chrhap_chrhap_chr_pop, mrk = "13_1725133")
plot(f,lwd = 0.2,
hap.names = hap.names(hap_chrhap_chrhap_chr_pop),
# no place for a legend inside the plot
legend.xy.coords = "none", xlim = c(1.65E+6, 1.8E+6) )
dev.off()

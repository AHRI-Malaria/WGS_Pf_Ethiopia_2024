# the command can be run
# Rscript ehh_plot_MCA2.R /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/hap_chr14_Ethiopia /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/snp.info.inp.Ethiopia
# Rscript ehh_plot_MCA2.R /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/hap_chr14_Eritrea /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/snp.info.inp.Eritrea
#Rscript ehh_plot_MCA2.R /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/hap_chr14_Nigeria /mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/ehh/countries_MAF0.01_hypervarrem_ihsth_2_xpehhth_2.5_tajimaDrem/snp.info.inp.Nigeria
library(rehh)
library(R.utils)
library(vcfR)
library(data.table)
library(qqman)

args = commandArgs(trailingOnly=TRUE)
country <- sub(".*hap_chr[0-9]+_", "", basename(args[1]))
hap_chr_pop <- data2haplohh(hap_file = args[1], map_file = args[2],
                             recode.allele = FALSE,
                             chr.name = 14, min_perc_geno.hap = 0.8, min_perc_geno.mrk = 0.8)

snps <- c("14_1552256", "14_1551150", "14_1551105", "14_1551095", "14_1549843")

for (snp in snps) {
  pos <- sub("14_", "", snp)

  tryCatch({
    # EHH plot
    res <- calc_ehh(hap_chr_pop,
                    mrk = snp,
                    include_nhaplo = TRUE,
                    include_zero_values = TRUE,
                    discard_integration_at_border = FALSE,
                    polarized = TRUE)

    pdf(paste0("mca2_ehh_plot_", country, "_", pos, ".pdf"),
        width = 7, height = 7)
    plot(res, xlim = c(1.4E+6, 1.7E+6),
         main = paste0("EHH - ", country, " (pos ", pos, ")"))
    dev.off()

    # Furcation plot
    f <- calc_furcation(hap_chr_pop, mrk = snp)

    pdf(paste0("mca2_furcation_plot_", country, "_", pos, ".pdf"),
        width = 7, height = 7)
    plot(f, lwd = 0.2,
         hap.names = hap.names(hap_chr_pop),
         legend.xy.coords = "none",
         xlim = c(1.5E+6, 1.6E+6),
         main = paste0("Furcation - ", country, " (pos ", pos, ")"))
    dev.off()
  }, error = function(e) {
    message("Skipping SNP ", snp, " (", country, "): ", conditionMessage(e))
  })
}

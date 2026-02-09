## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

# setwd("/Users/warri/OneDrive/Documents/Denver/Research/Data/ashad_n110")
sessionInfo()

# RDAC R packages location '/home/rcorcoran/R/x86_64-pc-linux-gnu-library/4.1'



## ----Libraries and Data------------------------------------------------------------------------------------------
suppressWarnings({
   library(vcfR)
   library(data.table)
   library(tibble)
   library(stringr)
   library(tidyr)
   library(dplyr)
})


## ----Setting up the Other Inputs------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]

## Making a table of chromosome name and chromosome length
genome_fai <- read.csv("SCRIPTS_DIR/chrom_length.tsv", sep = "\t")
chrom_list <- genome_fai$chrom


## ----Finding Fixed Alternative SNPs---------------------------------------------------------------------------------------------

## All individuals in the same file
all_snps <- read.table(paste(file_name, ".frq", sep=""),
                       sep="\t", header=FALSE, fill=TRUE, row.names=NULL, skip=1, 
                       col.names=c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF_ref", "AF_alt", "AF_alt2", "AF_alt3", "AF_alt4"))
print(paste0("all_snps: nrow=", nrow(all_snps))
head(all_snps)

## Populations tested separately
## for ( lake in c("amos", "bride", "long", "pat", "quon")) {
##    print(paste("Loading in", lake, "file", sep=" "))
##    pop <- read.table(paste("./all_chrom_scaffolds_bylake-vc_allsites_snps_gatkrec_maxMissing75_rmA28Q17_hard_softfiltered.recode_", lake, ".frq", sep=""),
##                      sep="\t", header=FALSE, fill=TRUE, row.names=NULL, col.names=c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF_ref", "AF_alt"), skip=1)[,c(1,2,4,5)]
##    assign(value=pop, x=paste(lake, "pop", sep="_"))
## }

## all_pops1 <- merge(amos_pop, bride_pop, by=c("CHROM", "POS"), all=FALSE, suffixes=c("_amos", "_bride"))
## all_pops2 <- merge(long_pop, pat_pop, by=c("CHROM", "POS"), all=FALSE, suffixes=c("_long", "_pat"))
## all_pops3 <- merge(all_pops2, quon_pop, by=c("CHROM", "POS"), all=FALSE, suffixes=c("", "_quon"))
## all_snps <- merge(all_pops1, all_pops3, by=c("CHROM", "POS"), all=FALSE)

## colnames(all_snps)[c(11,12)] <- c("N_CHR_quon", "AF_ref_quon")


## Separate out the allele from the frequency
all_snps_af <- all_snps %>%
   separate(AF_ref, into = c("allele_ref", "freq_ref"), sep = c(":")) %>%
   separate(AF_alt, into = c("allele_alt_major", "freq_alt_major"), sep = c(":"))
print(paste0("all_snps_af: nrow=", nrow(all_snps_af))
head(all_snps_af)


## If population-separate input:
## all_snps_af <- all_snps[,c(1:2,4,6,8,10,12)] %>%
##    separate(AF_ref_amos, into = c("allele_amos", "freq_amos"), sep = c(":")) %>%
##    separate(AF_ref_bride, into = c("allele_bride", "freq_bride"), sep = c(":")) %>%
##    separate(AF_ref_long, into = c("allele_long", "freq_long"), sep = c(":")) %>%
##    separate(AF_ref_pat, into = c("allele_pat", "freq_pat"), sep = c(":")) %>%
##    separate(AF_ref_quon, into = c("allele_quon", "freq_quon"), sep = c(":"))

## Pull out just the frequency and pos columns
## all_snps_freq <- all_snps_af[,c(1:2,4,6,8,10,12)]


### Separating the final step depending on if I'm making a positions file for reference replacement (refhom_zero) or site exclusion (fixedAlt)
if ( file_name == "all_chrom_scaffolds_bylake-vc_allsites_snps_gatkrec_maxMissing75_rmA28Q17_FARM_hard_softfiltered" ) {
   ### Getting list of SNPs to replace in the reference - if there are no individuals with the reference allele, pull that SNP ###
   refhom_zero_snps <- subset(all_snps_af, freq_ref == 0)
   print(paste0("refhom_zero_snps: nrow=", nrow(refhom_zero_snps))
   head(refhom_zero_snps)

   print(paste("nrow(subset(refhom_zero_snps, N_ALLELES == 2)):", nrow(subset(refhom_zero_snps, N_ALLELES == 2)))
   print(paste("nrow(subset(refhom_zero_snps, N_ALLELES == 3)):", nrow(subset(refhom_zero_snps, N_ALLELES == 3)))
   print(paste("nrow(subset(refhom_zero_snps, N_ALLELES == 4)):", nrow(subset(refhom_zero_snps, N_ALLELES == 4)))
   
   write(file=paste(file_name, "_refhom_zero_snps_vcftools.txt", sep=""),
         x=paste(refhom_zero_snps$CHROM, refhom_zero_snps$POS, sep="\t"), sep=", ", ncolumns=1)
} else {
   ### Getting list of SNPs to exclude in vcftools for analyses ### 
   ## Remove any line where all populations have the same AF
   fixedAlt_snps <- subset(all_snps_af, freq_alt_major == 1)
   print(paste0("fixedAlt_snps: nrow=", nrow(fixedAlt_snps))
   head(fixedAlt_snps)

   ## intergenic_snps_freq <- intergenic_snps_freq[apply(intergenic_snps_freq[,c(3:7)],
   ##    MARGIN = 1, # rowwise
   ##    FUN = function(x) length(unique(x)) > 1), ]


   ## vcftools positions file for removing fixedAlt SNPs
   write(file=paste(file_name, "_fixedAlt_snps_vcftools.txt", sep=""),
         x=paste(fixedAlt_snps$CHROM, fixedAlt_snps$POS, sep="\t"), sep=", ", ncolumns=1)
}

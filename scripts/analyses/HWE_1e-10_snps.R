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
   library(readr)
})


## ----Setting up Other Inputs-------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]


## ----Subsetting Neutral Allele Frequencies------------------------------------------------------------------------------

HWE_file <- read.table(paste(file_name, ".hwe", sep=""), sep="\t", header=TRUE)

write.table(file=paste(file_name, "HWE1e-10_snps.txt", sep="_"), x=subset(HWE_file, P_HWE <= 1e-10)[,c(1:2)], 
   quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

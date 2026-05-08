## ----setup, include=FALSE----------------------------------------------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

# setwd("/Users/warri/OneDrive/Documents/Denver/Research/Data/ashad_n110")
#sessionInfo()

# RDAC R packages location '/home/rcorcoran/R/x86_64-pc-linux-gnu-library/4.1'



## ----Libraries and Data------------------------------------------------------------------------------------------
suppressWarnings({
   library(vcfR)
   library(data.table)
   library(ggplot2)
   library(cowplot)
   library(tibble)
   library(stringr)
   library(tidyr)
   library(dplyr)
})


## ----Setting up the Other Inputs---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
pop <- args[2]

#chroms <- data.frame(
#  chr=c("NC_014690.1", "NC_055980.1", "NC_055979.1", "NC_055978.1", "NC_055977.1", "NC_055976.1", "NC_055975.1", "NC_055974.1", "NC_055973.1", "NC_055972.1", "NC_055971.1", "NC_055970.1", "NC_055969.1", "NC_055968.1", "NC_055967.1", "NC_055966.1", "NC_055965.1", "NC_055964.1", "NC_055963.1", "NC_055962.1", "NC_055961.1", "NC_055960.1", "NC_055959.1", "NC_055958.1", "NC_055957.1"),
#  chr_num=seq(25,1))
#chrom_length <- read.table("./../../../../8_filtervcf/chrom_length.tsv", header=TRUE, sep="\t")

chroms <- data.frame(chr="A", chr_num=1)
chrom_length <- data.frame(chrom="A", length=30000000)


## ----Running the Sliding Window Code-----------------------------------------------------------------------------
all_ohana_windows <- c()
ohana_pos <- cbind(read.table(paste0(in_file, "_scan_lle-ratios.txt"), header=TRUE, sep="\t")[,-1], 
                   read.table(paste0(in_file, ".map"), header=FALSE, sep="\t")[,-c(2:3)]) %>%
  rename_with(~c("lle_ratio", "global_lle", "local_lle", "f_1", "f_2", "chr", "pos"))


## 50kb sliding window (by 10kb)
num_windows <- ceiling(sum(chrom_length[,2]) / 10000) + 5
ohana_windows <- data.frame("chr"=rep(0, num_windows), "window_pos_1"=rep(0, num_windows), "window_pos_2"=rep(0, num_windows), 
                            "mean_lle_ratio"=rep(0, num_windows), "cum_lle_ratio"=rep(0, num_windows), "num_snps"=rep(0, num_windows))
ohana_windows_row <- 1
slide <- 10000

for (chrom in 1) {

  print(paste("Running", pop, "chr", chrom))

  window_start <- 1
  window_end <- 50000

  while ( window_start + 40000 <= chrom_length[chrom,2]) {
    subset_tmp <- subset(ohana_pos, chr %in% chrom_length[chrom,1] & pos >= window_start & pos <= window_end)
    ohana_windows$chr[ohana_windows_row] <- chrom_length[chrom,1]
    ohana_windows$window_pos_1[ohana_windows_row] <- window_start
    ohana_windows$window_pos_2[ohana_windows_row] <- window_end
    ohana_windows$mean_lle_ratio[ohana_windows_row] <- mean(subset_tmp$lle_ratio, na.rm=TRUE)
    ohana_windows$cum_lle_ratio[ohana_windows_row] <- sum(subset_tmp$lle_ratio)
    ohana_windows$num_snps[ohana_windows_row] <- nrow(subset_tmp)
    ohana_windows$top_mean_lle_ratio[ohana_windows_row] <- mean(subset(subset_tmp, lle_ratio >= quantile(subset_tmp$lle_ratio, 0.8))$lle_ratio, na.rm=TRUE)
    ohana_windows$max_lle_ratio[ohana_windows_row] <- ifelse(ohana_windows$num_snps[ohana_windows_row] != 0, max(subset_tmp$lle_ratio, na.rm=TRUE), NA)

    window_start <- window_start + slide
    if ( window_end + 10000 <= chrom_length[chrom,2] ) { window_end <- window_end + slide }
    else if ( window_end + 10000 > chrom_length[chrom,2] ) { window_end <- chrom_length[chrom,2] }
    ohana_windows_row <- ohana_windows_row + 1
  }
}

print(summary(ohana_windows$num_snps))

ohana_windows <- subset(ohana_windows, !is.na(cum_lle_ratio)) %>%
  merge(., chroms, by="chr") %>%
  mutate(chr_num = as.factor(chr_num))
write.table(file=paste0(in_file, "_selscan_50kb_10kb.txt"), x=ohana_windows, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## ----setup, include=FALSE----------------------------------------------------------------------------------------
## ----Libraries and Data------------------------------------------------------------------------------------------
suppressWarnings({
   library(vcfR)
   library(data.table)
   library(ggplot2)
   library(cowplot)
   library(tibble)
   library(stringr)
   library(tidyr)
})


## ----Setting up the Other Inputs------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

run_table <- read.table(paste0(out_dir, "/ohana_pop_structure_lls.txt"), header=FALSE, sep="\t") %>%
   mutate(K = c(rep(2, 20), rep(3, 20), rep(4, 20), rep(5, 20), rep(6, 20)))
colnames(run_table)[3] <- "ll"

mins <- data.frame("K"=c(2, 3, 4, 5, 6), 
                   "log_likelihood"=c(min(subset(run_table, K %in% 2)$ll), min(subset(run_table, K %in% 3)$ll), 
                                      min(subset(run_table, K %in% 4)$ll), min(subset(run_table, K %in% 5)$ll), 
                                      min(subset(run_table, K %in% 6)$ll)))

ggsave(filename=paste0(out_dir, "best_K_2_6.pdf"), 
       plot=(ggplot(mins, aes(x=K, y=log_likelihood)) + 
                geom_point() + 
                geom_line() +
                ylab("Minimum Log-Likelihood for each K (20 repititions each)")), 
       width=6, height=6)


## Which run do the mins correspond to?
for ( K_val in c(2,3,4,5,6)) {
   print(paste(K_val, row.names(subset(run_table, ll %in% subset(mins, K %in% K_val)$log_likelihood)), sep=": "))
}

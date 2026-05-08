suppressWarnings({
   library(afvaper)
   library(vcfR)
   library(data.table)
   library(ggplot2)
   library(cowplot)
   library(tibble)
   library(stringr)
   library(tidyr)
})

args <- commandArgs(trailingOnly=TRUE)
out_full <- args[1]
k_val <- args[2]
out_subset < paste0(out_full, "_LD50kb") 

run_table_all <- c()
for ( run in 1:20 ) {
   run_table <- read.table(paste0("./pop_structure/ohana_fish_admixture_infer_K", k_val, "_", run, ".out"), header=TRUE, sep="\t", skip=3)
   run_table$run <- run
   if ( run == 1 ) {
      run_table_all <- run_table
   } else {
      run_table_all <- rbind(run_table_all, run_table)
   }
}

print(paste("K =", k_val, "highest loglikelihood run:", subset(run_table_all, -log_likelihood == max(-run_table_all$log_likelihood))$run)

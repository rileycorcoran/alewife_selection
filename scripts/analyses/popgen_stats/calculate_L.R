amos_ns <- 543893
bride_ns <- 543893
long_ns <- 543893
pat_ns <- 543893
quon_ns <- 543893
amos_i <- 20
bride_i <- 30
long_i <- 20
pat_i <- 19
quon_i <- 18
amos_het <- 0.14329768056384
bride_het <- 0.255855625571989
long_het <- 0.156552960008373
pat_het <- 0.187106414828695
quon_het <- 0.17847314076479


for ( lake in c("amos", "bride", "long", "pat", "quon") ){
  ns <- get(paste0(lake, "_ns"))
  ni <- get(paste0(lake, "_i"))
  het <- get(paste0(lake, "_het"))
  
  l <- log(0.05/(ns*ni)) / log(1-het)
  print(paste(lake, "# of SNPs:", l))

  N_out <- 4
  t <- round((N_out + 1)/round(l, digits=3), digits=3)  
  print(paste(lake, "threshold:", t))
}

# [1] "amos: 129.3465839998"
# [1] "bride: 175.330974570368"
# [1] "long: 211.517686942076"
# [1] "pat: 188.565016456386"
# [1] "quon: 182.078049571537"

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(plyr)
  library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
test <- args[1]
chr_num <- args[2]

chroms <- data.frame(
  chr=c("NC_014690.1", "NC_055980.1", "NC_055979.1", "NC_055978.1", "NC_055977.1", "NC_055976.1", "NC_055975.1", "NC_055974.1", "NC_055973.1", "NC_055972.1", "NC_055971.1", "NC_055970.1", "NC_055969.1", "NC_055968.1", "NC_055967.1", "NC_055966.1", "NC_055965.1", "NC_055964.1", "NC_055963.1", "NC_055962.1", "NC_055961.1", "NC_055960.1", "NC_055959.1", "NC_055958.1", "NC_055957.1"),
  chr_num=seq(25,1))

all_lakes <- c("bride", "amos", "long", "quon", "pat")
l_lakes <- c("amos", "long", "quon", "pat")

if ( test %in% "pi" ) {
  names_from <- c("pop")
  col_name <- "avg_pi"
  lake_list <- c("amos", "bride", "long", "quon", "pat")
  nuc_div_threshold <- 0.01
} else if ( test %in% "dxy" ) {
  names_from <- c("pop1", "pop2")
  lake_list <- c("amos", "long", "quon", "pat")
  col_name <- "avg_dxy"
  nuc_div_threshold <- 0.99
} else if ( test %in% "fst" ) {
  names_from <- c("pop1", "pop2")
  lake_list <- c("amos", "long", "quon", "pat")
  col_name <- "avg_wc_fst"
  nuc_div_threshold <- 0.99
}
  
## Load in site-level data
chrom_length <- read.table("./../../../../8_filtervcf/chrom_length.tsv", header=TRUE, sep="\t")
# file_list <- list.files(path=paste("./sitelevel/", test, "/subset", sep=""), pattern="\\.txt$", full.names=TRUE)
# data_list <- lapply(file_list, function(file) { read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE) })
# data_long <- do.call(rbind, data_list) %>%
#   #rename_with(~c("avg"), get(colname)) %>%
#   filter(!is.na(get(col_name)))

file_list <- list.files(path=paste("./sitelevel/", test, "/subset", sep=""), pattern=paste0("all_", test, "_rbind_chr", chr_num, "_\\d+\\.txt$"), full.names=TRUE)
data_list <- lapply(file_list, function(file) { read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE) })
data_long <- do.call(rbind, data_list) %>%
   filter(!is.na(get(col_name))) %>%
   distinct(.keep_all = TRUE)
    
data_wide <- pivot_wider(data_long, id_cols=c("chromosome", "window_pos_1", "window_pos_2"), names_from=all_of(names_from), values_from=all_of(col_name)) %>%
  select(-window_pos_2) %>%
  merge(., chroms, by.x="chromosome", by.y="chr") %>%
  mutate(chr_num = as.factor(chr_num),
         window_pos_1 = as.numeric(window_pos_1))
colnames(data_wide) <- tolower(colnames(data_wide)) 
    
## Subset to just the comparisons of interest if necessary (i.e. for fst and dxy)
if ( test %in% "fst" | test %in% "dxy" ) {
  data_wide <- data_wide %>% select(contains(c("chromosome", "window", "bride"))) %>%
    rename_with(~c("amos", "long", "quon", "pat"), c("amos_bride", "bride_long", "bride_pat", "bride_quon"))
}
  
for ( lake in lake_list ) {
#for ( lake in c("bride", "long", "pat", "quon") ) {
  lake_test <- cbind(data_wide[,1:2], data_wide[,grepl(lake, colnames(data_wide))]) %>%
    rename_with(~c("chr", "pos", "value")) %>%
    filter(!is.na(value), chr != "NC_014690.1") %>%
    mutate(value = as.numeric(value))
  assign(x=paste(lake, test, sep="_"), value=lake_test, pos=1)
      
  num_windows <- ceiling(sum(chrom_length[,2]) / 10000) + 5
  windows <- data.frame("chr"=rep(0, num_windows), "window_pos_1"=rep(0, num_windows), "window_pos_2"=rep(0, num_windows),
                            "mean_value"=rep(0, num_windows), "num_snps"=rep(0, num_windows))
  windows_row <- 1
      
  for (chrom in chr_num) {
 
    print(paste("Running", lake, "chr", chrom, sep=" "))
 
    chr_subset <- subset(lake_test, chr %in% chrom_length[chrom,1]) %>%
      arrange(pos) %>%
      mutate(snp_count = 1:nrow(.))
  
    window_start <- 1
    window_end <- 50000
    slide <- 10000
  
    while ( window_start + 40000 <= chrom_length[chrom,2] ) {
      subset_tmp <- subset(chr_subset, pos >= window_start & pos <= window_end)
    
      windows$chr[windows_row] <- chrom_length[chrom,1]
      windows$window_pos_1[windows_row] <- window_start
      windows$window_pos_2[windows_row] <- window_end
      windows$mean_value[windows_row] <- mean(subset_tmp$value, na.rm=TRUE)
      windows$num_snps[windows_row] <- nrow(subset_tmp)
    
      window_start <- window_start + slide
      if ( window_end + slide <= chrom_length[chrom,2] ) { window_end <- window_end + slide }
      else if ( window_end + slide > chrom_length[chrom,2] ) { window_end <- chrom_length[chrom,2] }
      windows_row <- windows_row + 1
    }
  } 
  
  print(summary(windows$num_snps))
  windows <- subset(windows, !is.na(mean_value)) %>%
    merge(., chroms, by="chr") %>%
    mutate(chr_num = as.factor(chr_num))
      
#  if ( test %in% "fst" | test %in% "dxy" ) {
#    windows <- windows %>%
#      mutate(outlier = ifelse(mean_value >= quantile(mean_value, probs = nuc_div_threshold, na.rm=TRUE), 1, 0)) %>%
#      rename_with(~paste("outlier", lake, sep="_"), "outlier")
#  } else {
#    windows <- windows %>%
#      mutate(outlier = ifelse(mean_value <= quantile(mean_value, probs = nuc_div_threshold, na.rm=TRUE), 1, 0)) %>%
#      rename_with(~paste("outlier", lake, sep="_"), "outlier")
#  }
  assign(x=paste(lake, "windows", sep="_"), value=windows, pos=1)
}
if ( test %in% "fst" | test %in% "dxy" ) {
  data_wide <- merge(amos_windows, long_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), suffixes=c("_amos", "_long"), all=TRUE) %>%
    merge(., quon_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), all=TRUE) %>%
    merge(., pat_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), suffixes=c("_quon", "_pat"), all=TRUE) %>%
    mutate(window_pos_1 = as.numeric(window_pos_1),	
           window_pos_2 = as.numeric(window_pos_2)) # ,
#           total_outliers = rowSums(.[,c("outlier_amos","outlier_long","outlier_quon","outlier_pat")]), 
#           outlier_lakes = NA) %>%
#    filter(!is.na(total_outliers))
} else {
  data_wide <- merge(amos_windows, bride_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), suffixes=c("_amos", "_bride"), all=TRUE) %>%
    merge(., long_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), all=TRUE) %>%
    merge(., quon_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), suffixes=c("_long", "_quon"), all=TRUE) %>%
    merge(., pat_windows, by=c("chr", "window_pos_1", "window_pos_2", "chr_num"), all=TRUE) %>%
    rename_with(~c("mean_value_pat", "num_snps_pat"), c("mean_value", "num_snps")) %>%
    mutate(window_pos_1 = as.numeric(window_pos_1),
           window_pos_2 = as.numeric(window_pos_2)) # ,
#           total_outliers = rowSums(.[,c("outlier_amos", "outlier_bride", "outlier_long","outlier_quon","outlier_pat")]), 
#           outlier_lakes = NA)) %>%
#    filter(!is.na(total_outliers))
}

write.table(file=paste0("./sitelevel/", test, "/subset/all_", test, "_chr", chr_num, "_wide_50kb_10kb_windows.txt"), data_wide, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
# write.table(file=paste("all", test, "wide_50kb_10kb_windows.txt", sep="_"), data_wide, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

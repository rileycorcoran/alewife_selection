## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

# setwd("/Users/warri/OneDrive/Documents/Denver/Research/Data/ashad_n110")
sessionInfo()

# RDAC R packages location '/home/rcorcoran/R/x86_64-pc-linux-gnu-library/4.1'



## ----Libraries and Data------------------------------------------------------------------------------------------
suppressWarnings({
   library(afvaper)
   library(vcfR)
   library(data.table)
   library(ggplot2)
   library(cowplot)
   library(tibble)
   library(stringr)
   library(tidyr)
   library(dplyr)
})


## ----General Inputs------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
in_folder <- args[2]

## Colors
amos_col <- "#1D3203"
bride_col <- "#1855F2"
long_col <- "#496E12"
pat_col <- "#A0E72F"
quon_col <- "#74AB20"
colors <- c("amos" = amos_col, "bride" = bride_col, "long" = long_col, "pat" = pat_col, "quon" = quon_col)
afvaper_col <- "#ef476f"

chrom1 <- "grey75"
chrom2 <- "grey25"

chroms <- data.frame(
  chr=c("NC_014690.1", "NC_055980.1", "NC_055979.1", "NC_055978.1", "NC_055977.1", "NC_055976.1", "NC_055975.1", "NC_055974.1", "NC_055973.1", "NC_055972.1", "NC_055971.1", "NC_055970.1", "NC_055969.1", "NC_055968.1", "NC_055967.1", "NC_055966.1", "NC_055965.1", "NC_055964.1", "NC_055963.1", "NC_055962.1", "NC_055961.1", "NC_055960.1", "NC_055959.1", "NC_055958.1", "NC_055957.1"),
  chr_num=c(25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1))


## ----F - Heterozygosity-------------------------------------------------------------------------------------

vcftools_het <- c()
for ( lake in c("amos", "bride", "long", "pat", "quon")) {
  example <- read.table(paste(in_folder, in_file, "_", lake, ".het", sep=""), header=TRUE, sep="\t")
  example_mean<-mean(subset(example, !is.na(example$F))$F)
  example$lake <- lake
  
  example_hist <- ggplot(example, aes(x=F, color=lake, fill=lake)) +
    geom_histogram(binwidth = 0.005) +
    geom_vline(xintercept = example_mean, color="red") +
    scale_color_manual(values=c(get(paste(lake,"_col", sep="")), "red"), aesthetics=c("color", "fill")) +
    xlab("F - Heterozygosity") +
    ylab(lake) +
    xlim(-1, 0) +
    theme(legend.position = "none")

  assign(value=example, x=paste(lake, "F", sep="_"))
  assign(value=example_hist, x=paste(lake, "F_hist", sep="_"))
  
  vcftools_het <- rbind(vcftools_het, example)
}
vcftools_het$lake <- factor(vcftools_het$lake)


## Boxplots
vcftools_f_plot <- ggplot(vcftools_het, aes(x=lake, y=F)) +
  geom_violin(mapping=aes(color=lake, fill=lake, alpha=0.5), width=0.8) +
  geom_boxplot(mapping=aes(color=lake), width=0.15, linewidth=0.8) +
  scale_color_manual(values=colors, aesthetics=c("color", "fill")) +
  xlab("Lake") +
  ylab("F - Heterozygosity") +
  ylim(-1,0) +
  theme_cowplot() +
  theme(legend.position = "none")
vcftools_f_plot
ggsave(paste(in_folder, "vcftools_fis_boxplot.pdf", sep=""), vcftools_f_plot, height=2, width=4)

vcftools_f_anova <- aov(F ~ lake, data = vcftools_het)
print(vcftools_f_anova)
vcftools_f_tukey <- TukeyHSD(vcftools_f_anova)
print(vcftools_f_tukey)


## vcftools Heterozygosity counts
perc_df <- c()
lake_column <- 5
for (lake in c("amos", "bride", "long", "pat", "quon")) {
  het_example <- read.table(paste(in_folder, in_file, "_", lake, ".hwe", sep=""), header=TRUE, sep="\t")
  het_example <- het_example %>%
    separate(OBS.HOM1.HET.HOM2., into = c("OHOM1", "OHET", "OHOM2"), sep = c("/")) %>%
    separate(E.HOM1.HET.HOM2., into = c("EHOM1", "EHET", "EHOM2"), sep = c("/"))
  het_example$OHOM1 <- as.numeric(het_example$OHOM1)
  het_example$OHOM2 <- as.numeric(het_example$OHOM2)
  het_example$OHET <- as.numeric(het_example$OHET)
  het_example$EHOM1 <- as.numeric(het_example$EHOM1)
  het_example$EHOM2 <- as.numeric(het_example$EHOM2)
  het_example$EHET <- as.numeric(het_example$EHET)

  ## There are some rows with NA's in the expected columns (any site where no individual in that pop had a retained genotype) so I'm going to remove any rows with them
  het_example <- subset(het_example, !is.na(EHOM1) & !is.na(EHOM2) & !is.na(EHET)) 
  
  ## Calculate fraction heterozygous
  example_ohet_fraction <- sum(het_example$OHET) / sum(sum(het_example$OHOM1), sum(het_example$OHOM2), sum(het_example$OHET))
  example_ehet_fraction <- sum(het_example$EHET) / sum(sum(het_example$EHOM1), sum(het_example$EHOM2), sum(het_example$EHET))
   
  print(paste(lake, "Fraction Expected Heterozygous Loci:", example_ehet_fraction, sep=" "))
  print(paste(lake, "Fraction Observed Heterozygous Loci:", example_ohet_fraction, sep=" "))


  ## Calculating and graphing site-level hom/het differences
  het_example$n_ind <- rowSums(het_example[,c(3:5)])
  het_row <- het_example %>% rowwise() %>%
    mutate(OHOM1_frac = (OHOM1 / n_ind)) %>%
    mutate(OHOM2_frac = (OHOM2 / n_ind)) %>%
    mutate(OHET_frac = (OHET / n_ind)) %>%
    mutate(EHOM1_frac = (EHOM1 / n_ind)) %>%
    mutate(EHOM2_frac = (EHOM2 / n_ind)) %>%
    mutate(EHET_frac = (EHET / n_ind))

  print(paste(lake, "Mean OHET_frac Across All Loci:", mean(het_row$OHET_frac), sep=" "))
  print(paste(lake, "SE:", mean_se(het_row$OHET_frac)$ymax - mean_se(het_row$OHET_frac)$y, sep=" "))
  print(paste(lake, "Mean EHET_frac Across All Loci:", mean(het_row$EHET_frac), sep=" "))
  print(paste(lake, "SE:", mean_se(het_row$EHET_frac)$ymax - mean_se(het_row$EHET_frac)$y, sep=" "))

  for ( i in 1:nrow(het_row) ) {
    if ( i %% 10000 == 0 ) { print(paste("Row", i, "out of", nrow(het_row))) }
    if ( het_row$OHOM1_frac[i] == "1" ) {
      het_row$obs_group[i] <- "all_hom1"
      het_row$obs_value[i] <- 1
    } else if ( het_row$OHOM2_frac[i] == "1" ) {
      het_row$obs_group[i] <- "all_hom2"
      het_row$obs_value[i] <- 1
    } else if ( het_row$OHET_frac[i] == "1" ) {
      het_row$obs_group[i] <- "all_het"
      het_row$obs_value[i] <- 1
    } else {
      het_row$obs_group[i] <- "mixed"
      het_row$obs_value[i] <- het_row$OHET_frac[i]
    }
    if ( het_row$EHOM1_frac[i] == "1" ) {
      het_row$exp_group[i] <- "all_hom1"
      het_row$exp_value[i] <- 1
    } else if ( het_row$EHOM2_frac[i] == "1" ) {
      het_row$exp_group[i] <- "all_hom2"
      het_row$exp_value[i] <- 1
    } else if ( het_row$EHET_frac[i] == "1" ) {
      het_row$exp_group[i] <- "all_het"
      het_row$exp_value[i] <- 1
    } else {
      het_row$exp_group[i] <- "mixed"
      het_row$exp_value[i] <- het_row$EHET_frac[i]
    }
  }
  write.table(x=het_row, file=paste(in_folder, lake, "_hom_het_calculations.txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

  het_row_longer <- pivot_longer(het_row, cols=c(obs_group, exp_group), 
                                 names_to=c("type", NA), names_sep="_", values_to="group") ## Pivot the "group" columns to put obs and exp into the same column for graphing
  het_row_longer <- pivot_longer(het_row_longer, cols=c(obs_value, exp_value), 
                                 names_to=c("type2", NA), names_sep="_", values_to="value") ## Pivot the "value" columns to put obs and exp into the same column for graphing
  het_row_longer <- subset(het_row_longer, type == type2) ## Keep only rows in which both the "group" and "value" refer to the same "type"
  het_row_longer <- het_row_longer[,-18]  ## Remove the duplicate "type" column from the double pivoting
  
  het_row_longer$group <- factor(het_row_longer$group, levels=c("all_hom1", "all_hom2", "all_het", "mixed"))
  
  het_row_box <- ggplot(het_row_longer, aes(x=group, group=type, color=type, fill=type)) + 
    geom_bar(position=position_dodge()) +
    ylab(paste("Total Sites in", lake, "out of", nrow(het_row), sep=" ")) +
    xlab("Population characterization") +
    scale_color_manual(values=c("grey75", "grey25")) +
    scale_fill_manual(values=c("grey75", "grey25")) +
    theme_cowplot() +
    theme(legend.position = "null")
  
  het_mixed_bar <- ggplot(subset(het_row_longer, group %in% "mixed" | group %in% "all_het"), aes(x=value, group=type, color=type, fill=type)) +
    geom_histogram(binwidth = 1/max(het_row$n_ind), position="identity", alpha=0.2) + 
    ylab(paste("Total Mixed or All Heterozygous Sites in", lake, sep=" ")) +
    xlab("Fraction of Individuals Heterozygous") +
    scale_color_manual(values=c("grey75", "grey25")) +
    scale_fill_manual(values=c("grey75", "grey25")) +
    theme_cowplot()
  
  if ( lake != "pat" ) {
    het_row_box <- het_row_box + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    het_mixed_bar <- het_mixed_bar + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  het_alleles_grid <- plot_grid(het_row_box, het_mixed_bar, nrow=1, rel_widths=c(0.75,1))
  ggsave(paste(in_folder, lake, "_hom_het_allele_distribution.pdf", sep=""), het_alleles_grid, width=12, height=6)
  assign(value=het_alleles_grid, x=paste(lake, "het_alleles_grid", sep="_"))

  
  ## Calculate fraction heterozygous in 50kb windows
  het_windows <- c()
  for (i in 25:2) {
    chr <- chroms$chr[i]
    het_frame <- read.table(paste("./../pixy/pixy_dxy_50kb_", chr, ".txt", sep=""), header=TRUE, sep="\t")
    het_frame <- pivot_wider(het_frame, id_cols=c("chromosome", "window_pos_1", "window_pos_2"), names_from=c("pop1", "pop2"), values_from="avg_dxy")
    het_frame <- het_frame[,c(1:3)]
    het_frame$EHET <- 0
    het_frame$OHET <- 0
    het_frame$total <- 0
    het_windows <- rbind(het_windows, het_frame)
  }
  het_windows$row_num <- row_number(het_windows)

  ## Adding the observed and expected heterozygosity up for each 50kb window to make a graph like van der Zee  
  het_example <- subset(het_example, CHR != "NC_014690.1")
  for (i in 1:nrow(het_example)) {
    het_window <- subset(het_windows, chromosome %in% het_example$CHR[i] & 
                           window_pos_1 <= het_example$POS[i] & 
                           window_pos_2 >= het_example$POS[i]) ## Pull the window that the row of het_example would fall into
    het_window$EHET <- sum(het_window$EHET, het_example$EHET[i]) ## Add the expected # of het loci in row [i] to the running window sum
    het_window$OHET <- sum(het_window$OHET, het_example$OHET[i]) ## Add the observed # of het loci in row [i] to the running window sum
    het_window$total <- sum(het_window$total, het_example$n_ind[i]) ## Add the total # of reads at loci [i] to the running window sum
    het_windows[het_windows$row_num == het_window$row_num, ] <- het_window ## Replace the window row with the edited window pulled before
  }
  het_windows <- het_windows %>% rowwise() %>% 
    mutate(OHET_frac = (OHET / total)) %>%
    mutate(EHET_frac = (EHET / total))
  write.table(x=het_windows, file=paste(in_folder, lake, "_het_50kb_windows.txt", sep=""), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  het_windows$lake <- lake
  assign(value=het_windows, x=paste(lake, "het_windows", sep="_"))

  ## Create df with % observed and % expected heterozygosity per loci
  ## I can add code to make it calculate this in 50kb windows (make POS fit within increments of 50000)
  het_example <- het_example %>% rowwise() %>% 
    mutate(obs = OHET/(OHET + OHOM1 + OHOM2)) %>%
    mutate(exp = EHET/(EHET + EHOM1 + EHOM2))
  if (lake == "amos") {
    perc_df$CHR <- het_example$CHR
    perc_df$POS <- het_example$POS
    perc_df$amos_obs <- het_example$obs
    perc_df$amos_exp <- het_example$exp
    perc_df <- data.frame(perc_df)
  } else {
    perc_df <- merge(perc_df, het_example[,c(1:2,9:10)], by=c("CHR", "POS"), all=TRUE)
    colnames(perc_df)[lake_column] <- paste(lake, "obs", sep="_")
    colnames(perc_df)[lake_column+1] <- paste(lake, "exp", sep="_")
    lake_column <- lake_column + 2
  }
}

all_het_alleles_grid <- plot_grid(bride_het_alleles_grid, amos_het_alleles_grid, long_het_alleles_grid, quon_het_alleles_grid, pat_het_alleles_grid, nrow=5)
ggsave(paste(in_folder, "all_lakes_hom_het_alleles_distributions.pdf", sep=""), all_het_alleles_grid, width=12, height=18)


## Boxplots and Histograms 
test <- "windows"
for (type in c("obs", "exp")) local({
  
  if (test == "sitelevel") { 
    if (type == "obs") { count_subset <- perc_df[,c(3,5,7,9,11)] } 
    else { count_subset <- perc_df[,c(4,6,8,10,12)] } 
  }
  
  if (test == "windows") { 
    het_windows <- rbind(get("bride_het_windows", pos=1), get("amos_het_windows", pos=1), get("long_het_windows", pos=1), get("pat_het_windows", pos=1), get("quon_het_windows", pos=1))
    count_subset <- pivot_wider(het_windows, id_cols = c("chromosome", "chr_num", "window_pos_1", "window_pos_2"), names_from = "lake", values_from = c("EHET_frac", "OHET_frac"))
    if (type == "obs") { count_subset <- count_subset[,c(10:14)] }
    else { count_subset <- count_subset[,c(5:9)] } 
  }

  colnames(count_subset) <- c("bride", "amos", "long", "pat", "quon")
  count_pivot <- pivot_longer(count_subset, cols=c("amos", "bride", "long", "pat", "quon"), names_to="lake", values_to="frac")
  count_pivot <- na.omit(count_pivot)
  
  ## Boxplots
  count_pivot$lake <- factor(count_pivot$lake, levels=c("bride", "amos", "long", "quon", "pat"))
  type <- type
  count_plot <<- ggplot(count_pivot, aes(x=lake, y=frac)) +
    geom_violin(mapping=aes(color=lake, fill=lake, alpha=0.5), width=0.8) +
    geom_boxplot(mapping=aes(color=lake), width=0.15, linewidth=0.8) +
    scale_color_manual(values=colors, aesthetics=c("color", "fill")) +
    xlab("Lake") +
    ylab(paste("Fraction", type, "Heterozygosity", sep=" ")) +
    ylim(0, 0.2) +
    theme_cowplot() +
    theme(legend.position = "none")
  #count_plot
  assign(x=paste(type, "het_boxplot", sep="_"), value=count_plot, pos=1)
  
  # One large boxplot
  count_pivot2 <- pivot_longer(het_windows, cols = c("EHET_frac", "OHET_frac"), names_to = "test", values_to = "frac")
  count_pivot2$lake <- factor(count_pivot2$lake, levels=c("bride", "amos", "long", "quon", "pat"))
  count_plot2 <<- ggplot(count_pivot2, aes(x=lake, y=frac, fill=test)) +
    geom_violin(mapping=aes(color=lake, fill=lake, alpha=0.5), width=0.8) +
    geom_boxplot(mapping=aes(color=lake, fill=test), width=0.4, linewidth=0.8) +
    scale_color_manual(values=c(colors, "EHET_frac"=chrom1, "OHET_frac"=chrom2), aesthetics=c("color", "fill")) +
    xlab("Lake") +
    ylab(paste("Fraction", type, "Heterozygosity", sep=" ")) +
    ylim(0, 0.2) +
    theme_cowplot() +
    theme(legend.position = "none")
  assign(x="het_boxplot2", value=count_plot2, pos=1)
  
  
  ## Histograms
  for ( lake in c("amos", "bride", "long", "pat", "quon") ) {
    het_lake <- count_subset[lake]
    het_lake <- na.omit(het_lake)
    colnames(het_lake)[1] <- "count"
  
    het_bin <- c()
    row_num <- 1
    increment <- 0.002
    for (i in seq(0, 1, by = increment)) {
      het_bin$bin[row_num] <- i
      het_bin$count[row_num] <- nrow(subset(het_lake, count > (i - increment) & count <= i))
      row_num <- row_num + 1
    }
    het_bin <- data.frame(het_bin)
    colnames(het_bin)[2] <- lake
    assign(x=paste(lake, "het_bin", sep="_"), value=het_bin, pos=1)
  }

  het_bins <- cbind(amos_het_bin, bride_het_bin$bride, long_het_bin$long, pat_het_bin$pat, quon_het_bin$quon)
  colnames(het_bins) <- c("bins", "amos", "bride", "long", "pat", "quon")
  het_bins <- pivot_longer(het_bins, cols=c("amos", "bride", "long", "pat", "quon"), names_to = "lake", values_to = "value")
  het_bins$lake <- as.factor(het_bins$lake)
  
  if (type == "exp") { het_bins$value <- if_else(het_bins$bins > 0.5, NA, het_bins$value)}
  
  het_bins_plot <<- ggplot(het_bins) +
    geom_line(aes(x=bins, y=value, col=lake), linewidth=0.6) +
    scale_color_manual(values=colors) +
    scale_x_continuous(name=paste("Fraction", type, "Heterozygous", sep=" "), limits=c(0,1)) +
    scale_y_continuous(name="Count", limits=c(0, 1000)) +
    theme_cowplot() +
    theme(legend.position = "none")
  het_bins_plot
  assign(x=paste(type, "het_hist", sep="_"), value=het_bins_plot, pos=1)
  
  ## Zoom
  het_bins_zoom <<- ggplot(het_bins) +
    geom_line(aes(x=bins, y=value, col=lake), linewidth=0.6) +
    scale_color_manual(values=colors) +
    scale_x_continuous(name=paste("Fraction", type, "Heterozygous", sep=" "), limits=c(0,0.4)) +
    scale_y_continuous(name="Count", limits=c(0, 1200)) +
    theme_cowplot() +
    theme(legend.position = "none")
  het_bins_zoom_end <<- ggplot(het_bins) +
    geom_line(aes(x=bins, y=value, col=lake), linewidth=0.6) +
    scale_color_manual(values=colors) +
    scale_x_continuous(limits=c(0.95,1)) +
    scale_y_continuous(limits=c(0, 1200)) +
    theme_cowplot() +
    theme(legend.position = "none", 
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  het_bins_zoom_all <- plot_grid(het_bins_zoom, het_bins_zoom_end, nrow=1, rel_widths = c(1,0.25))
  assign(x=paste(type, "het_hist_zoom", sep="_"), value=het_bins_zoom_all, pos=1)
  
  ## Are the lakes significantly different from each other?
  assign(value=aov(frac ~ lake, count_pivot), x=paste(type, "aov", sep="_"), pos=1)
  assign(value=TukeyHSD(get(paste(type, "aov", sep="_"))), x=paste(type, "tukey", sep="_"), pos=1)
})
ggsave(paste(in_folder, "obs_het_hist_zoom_extended.pdf", sep=""), obs_het_hist_zoom, width=2, height = 2)
ggsave(paste(in_folder, "exp_het_hist_zoom_extended.pdf", sep=""), exp_het_hist_zoom, width=2, height = 2)
het_zoom_grid <- plot_grid(obs_het_hist_zoom, exp_het_hist_zoom, nrow=2, align="hv")
ggsave(paste(in_folder,"heterozygosity_histograms_zoom.pdf", sep=""), het_zoom_grid, width=12, height=6)

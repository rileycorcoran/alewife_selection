library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tidyr)

in_dir <- "/Users/corcorri/Documents/Denver/Manuscript/popgen_stats"
out_dir <- "/Users/corcorri/Documents/Denver/Manuscript/r_figures/"

bootstrap_all <- data.table("lake"=c("Bride", "Amos", "Long", "Quon", "Pat"), 
                            "fst"=c(rep("NA",5)), 
                            "pi"=c(rep("NA",5)), 
                            "td"=c(rep("NA",5)), 
                            "ho"=c(rep("NA",5)), 
                            "he"=c(rep("NA",5)))

# Load data (in separate files for each chromosome)
for ( stat_type in c("pi", "fst", "td", "ho", "he")) {
  if ( stat_type %in% "pi" ) {
    folder_path <- paste0(in_dir, "pixy_files/pi")
    file_pattern <- "\\.txt$"
    col_name <- "avg_pi"
  } else if ( stat_type %in% "fst") {
    folder_path <- paste0(in_dir, "pixy_files/fst")
    file_pattern <- "\\.txt$"
    col_name <- "avg_wc_fst"
  } else if ( stat_type %in% "td" ) {
    folder_path <- in_dir
    file_pattern <- "\\.Tajima.D$"
    col_name <- "TajimaD"
  } else if ( stat_type %in% "ho" ) {
    folder_path <- in_dir
    file_pattern <- "het_50kb_windows\\.txt$"
    col_name <- "OHET_frac"
  } else if ( stat_type %in% "he" ) {
    folder_path <- in_dir
    file_pattern <- "het_50kb_windows\\.txt$"
    col_name <- "EHET_frac"
  }
  
  setwd(folder_path) #set directory
  file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)
  data_list <- lapply(file_list, function(file) {
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)})
  
  combined_data <- do.call(rbind, data_list)
  
  head(combined_data)
  data_df <- combined_data
  
  if ( stat_type %in% "fst" ) {
    data_df <- data_df %>% filter(pop1 == "Bride" | pop2 == "Bride")
    data_df$pop <- ifelse(data_df$pop1 == "Bride", data_df$pop2, data_df$pop1)
  } else if ( stat_type %in% "td" ) {
    i <- 1
    data_df$pop <- NA
    td_row_groups <- c(as.numeric(row.names(subset(data_df, CHROM %in% "NC_055957.1" & BIN_START == 0))), nrow(data_df))
    for ( lake in c("Amos", "Bride", "Long", "Pat", "Quon")) { 
      data_df$pop[td_row_groups[i]:td_row_groups[i+1]] <- lake
      i <- i + 1
    }
  } else if ( stat_type %like% "h" ) {
    data_df <- data_df %>%
      mutate(pop = c(rep("Amos", 17935), rep("Bride", 17935), rep("Long", 17935), rep("Pat", 17935), rep("Quon", 17935)))
  }
  
  ### boostrapping to produce 95% CIs #####
  
  n_windows <- 10000  # Number of windows per population (40k is half the number of observed windows)
  #bootstrap function
  bootstrap_data <- function(data, n_boot = 1000, sample_size = NULL) {
    results <- data.frame()

    for (pop_name in unique(data$pop)) {
      pop_data <- data %>% filter(pop == pop_name)

      if (nrow(pop_data) == 0) {
        warning(paste("No data found for population:", pop_name))
        next
      }

      n <- if (is.null(sample_size)) nrow(pop_data) else sample_size

      boot_means <- replicate(n_boot, {
        sample_data <- sample(pop_data[,col_name], size = n, replace = TRUE)
        mean(sample_data, na.rm = TRUE)
      })

      boot_df <- data.frame(
        pop = pop_name,
        mean_data = mean(boot_means, na.rm = TRUE),
        lower_ci = quantile(boot_means, 0.025, na.rm = TRUE),
        upper_ci = quantile(boot_means, 0.975, na.rm = TRUE)
      )

      results <- rbind(results, boot_df)
    }

    return(results)
  }
  
  # Run bootstrap code #
  bootstrap_results <- bootstrap_data(data_df, n_boot = 1000, sample_size = n_windows)
  print(bootstrap_results)
  
  # Add to bootstrap_all table:
  for ( for_lake in bootstrap_results$pop ) {
    bootstrap_all[which(bootstrap_all$lake == for_lake), stat_type] <- paste(signif(bootstrap_results$mean_data[which(bootstrap_results$pop == for_lake)], 3), 
                                                                         " (", signif(bootstrap_results$lower_ci[which(bootstrap_results$pop == for_lake)], 3), "-",
                                                                         signif(bootstrap_results$upper_ci[which(bootstrap_results$pop == for_lake)], 3), ")", sep="")
  }
  
  # write to table
  write.csv(bootstrap_results, paste0(out_dir, stat_type, "_boostrap_results.csv"), row.names=FALSE)
  
  # plot the means and confidence intervals
  bootstrap_mean_ci <- ggplot(bootstrap_results, aes(x = pop, y = mean_data)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    ylab(paste("Mean", stat_type, "(with 95% CI)")) +
    theme_minimal()
  print(bootstrap_mean_ci)
  
  ###################################################
  ###### testing differences between means, pairwise, with permutation tests ####
  ###############################################
  
  # ---- Define permutation test function ----
  permutation_test_data <- function(data, pop1, pop2, nperm = 10000) {
    g1 <- data[,col_name][data$pop == pop1]
    g2 <- data[,col_name][data$pop == pop2]
    g1 <- g1[is.finite(g1)]
    g2 <- g2[is.finite(g2)]
    
    if (length(g1) == 0 | length(g2) == 0) {
      warning(paste("Missing data for", pop1, pop2))
      return(NULL)
    }
    
    obs_diff <- mean(g1) - mean(g2)
    combined <- c(g1, g2)
    n1 <- length(g1)
    
    perm_diffs <- replicate(nperm, {
      shuffled <- sample(combined)
      mean(shuffled[1:n1]) - mean(shuffled[(n1 + 1):length(combined)])
    })
    
    p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
    
    tibble(
      pop1 = pop1,
      pop2 = pop2,
      observed_diff = obs_diff,
      mean_perm_diff = mean(perm_diffs),
      sd_perm_diff = sd(perm_diffs),
      lower_95_perm = quantile(perm_diffs, 0.025),
      upper_95_perm = quantile(perm_diffs, 0.975),
      p_value = p_val,
      perm_diff = list(perm_diffs)  # store full distribution for plotting
      )
  }
  
  # ----  Run all pairwise population comparisons ----
  pops <- unique(data_df$pop)
  pairs <- combn(pops, 2, simplify = FALSE)
  
  perm_results <- map_dfr(pairs, ~permutation_test_data(data_df, .x[1], .x[2], nperm = 5000))
  
  # ----  Add multiple testing corrections ----
  perm_results <- perm_results %>%
    mutate(
      p_adj_fdr = p.adjust(p_value, method = "fdr"),
      p_adj_bonf = p.adjust(p_value, method = "bonferroni"),
      direction = paste(pop1, ifelse(observed_diff > 0, " > ", " < "), pop2),
      comparison = paste(pop1, "vs", pop2))
  
  # ----  Expand permutation distributions for plotting ----
  perm_plot_df <- perm_results %>%
    select(pop1, pop2, observed_diff, perm_diff, comparison) %>%
    unnest(cols = c(perm_diff))
  
  # ---- Plot permutation distributions ----
  ggplot(perm_plot_df, aes(x = perm_diff)) +
    geom_histogram(bins = 50, fill = "gray85", color = "white") +
    geom_vline(aes(xintercept = observed_diff), color = "red", size = 1) +
    facet_wrap(~ comparison, scales = "free", ncol = 2) +
    theme_bw(base_size = 12) +
    labs(
      x = paste("Permuted difference in mean", stat_type, "(pop1 − pop2)"),
      y = "Frequency",
      title = paste("Permutation Distributions of Mean", stat_type, "Differences"),
      subtitle = "Red line = observed difference")
  
  # --- Save summary results ----
  perm_results_export <- perm_results %>% 
    select(-c(perm_diff, comparison))
  write.csv(perm_results_export, paste0(out_dir, "pairwise_", stat_type, "_permutation_results.csv"), row.names = FALSE)
}


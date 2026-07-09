setwd("/tempor/merlat/exogap")

## Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ape, dplyr, ggplot2, patchwork, cowplot, ggtext, Matrix,
               circlize, phylogram, tidyverse, ggtree, RColorBrewer, glue)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.22", ask = FALSE)
BiocManager::install("ComplexHeatmap", ask = FALSE)

library(ComplexHeatmap)

## Load custom functions (copied next to this script by the Nextflow module)
functions_file = "/enadisk/tempor/merlat/exogap/bin/exostats_functions.R"
if (file.exists(functions_file)) {
  source(functions_file)
}

# Files definition
newick_file <- "/enadisk/tempor/merlat/exogap/cache/preprocessing/download_newick/newick/tree.nwk"
stats_file <- "/tempor/merlat/exogap/cache/preprocessing/merge_stats/all/all.stats"
summary_all_file <- "/enadisk/tempor/merlat/exogap/cache/repeats_annotation/merge_summaries/repeats_summary_all_only_repetitive_elements/repeats_summary_all_only_repetitive_elements.tsv"
summary_order_file <- "/enadisk/tempor/merlat/exogap/cache/repeats_annotation/merge_summaries/repeats_summary_by_order/repeats_summary_by_order.tsv"
summary_superfamily_file <- "/enadisk/tempor/merlat/exogap/cache/repeats_annotation/merge_summaries/repeats_summary_by_superfamily/repeats_summary_by_superfamily.tsv"

# Load phylogeny
phylo <- read.tree(newick_file)
phylo$tip.label <- sapply(phylo$tip.label, format_genome_name)

# Correlation plots
## Load stats
stats <- read_summary(stats_file, phylo)

## Load overall summary
summary_all <- read_summary(summary_all_file, phylo)
summary_all <- summary_all %>%
  dplyr::left_join(stats, by = "genome_id", select(genome_id, sum_len, N50))  %>%
  dplyr::rename(size = sum_len)

### Add units
res <- select_unit(summary_all, list('size' = 'pb', 'load' = 'n', 'N50' = 'pb'))
summary_all <- res$df
units <- res$units

draw_correlation <- function(df, x_column, y_column, tree, units, out_dir) {
  subtrees <- subtrees(as.phylo(tree))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (subtree in subtrees) {
    parent <- subtree[["node.label"]][[1]]
    childs <- subtree[["tip.label"]]
    parent <- sub("\\|.*$", "", parent)
    childs <- sub("\\|.*$", "", childs)
    if (length(childs) > 3) {
      print(parent)
      sub_df <- subset(df, df$genome_id %in% childs)
      print(sub_df)
      # Skip if y column missing
      if (!y_column %in% colnames(sub_df)) next
      # Compute correlations
      spearman <- suppressWarnings(cor.test(sub_df[[x_column]], sub_df[[y_column]], method = 'spearman', exact = FALSE))
      pearson <- suppressWarnings(cor.test(sub_df[[x_column]], sub_df[[y_column]], method = "pearson"))
      sub_df <- sub_df %>% mutate(
        label = glue("*r* = {round(pearson$estimate, 3)}, *p* = {round(pearson$p.value, 4)}\n\n*&rho;* = {round(spearman$estimate, 3)}, *p* = {round(spearman$p.value, 4)}")
      )
      x_label <- ifelse(x_column == 'coverage', "Coverage (%)", glue("Load{units$load}"))
      y_label <- ifelse(y_column == 'size', glue("Assembly size{units$size}"), "N50")
      
      p <- ggplot(sub_df, aes(x = .data[[x_column]], y = .data[[y_column]])) +
        geom_point() +
        geom_smooth(method = 'lm', formula = y ~ x) +
        ggtitle(parent) +
        ylab(y_label) +
        xlab(x_label) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.5, color = "grey80"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 1, color = "black"),
          plot.margin = margin(10, 10, 10, 10)
        ) +
        geom_richtext(data = sub_df, 
                      aes(label = label, x = Inf, y = -Inf), 
                      hjust = 1, 
                      vjust = -0.2,
                      fill = scales::alpha("white", 0.5))
      
      plot_name <- glue("correlation_{parent}_{tolower(x_column)}_vs_{tolower(y_column)}")
      save_plot(plot_name, p, 4, 4, out_dir = out_dir)
      # file_name <- glue("{out_dir}/correlation_{parent}_{tolower(x_column)}_vs_{tolower(y_column)}.pdf")
    }
  }
}

run_correlations <- function(df, y_columns, phylo, units, out_prefix) {
  y_columns = c("size", "N50")
  for (y in y_columns) {
    draw_correlation(df, "coverage", y, phylo, units, paste0(out_prefix, "_", tolower(y)))
    draw_correlation(df, "load",     y, phylo, units, paste0(out_prefix, "_", tolower(y)))
  }
}

run_correlations(
  df = summary_all,
  phylo = phylo,
  units = units,
  out_prefix = "correlation"
)

## Overview plot
## Load summary by TE orders
summary_order <- read_summary(summary_order_file, phylo)
summary_order <- filter_out_non_transposable_elements(summary_order)

### Add normalization of coverage and load 
summary_order_total <- summary_order %>%
  group_by(genome_id) %>%
  summarize(
    total_coverage = sum(coverage, na.rm = TRUE),
    total_load = sum(load, na.rm = TRUE)
  )

summary_order <- summary_order %>%
  left_join(summary_order_total %>%
              select(genome_id, total_coverage, total_load),
            by = c("genome_id"))

summary_order$norm_load <- summary_order$load * 100 / summary_order$total_load
summary_order$Norm_coverage <- summary_order$coverage * 100 / summary_order$total_coverage

### Add units
res <- select_unit(summary_order, list('load' = 'n'))
summary_order <- res$df

### Factorize by orders
summary_order <- factorize_by_order(summary_order)

## Assembly plots
assembly <- ggplot(summary_all, aes(x = genome_id, y = size)) +
  geom_col(fill = "gray30") +  # Set the color of the bars to dark gray
  labs(y = glue("Assembly size{units$size}")) +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.6, color = "gray"),  # Keep the major y grid lines
    panel.grid.major.x = element_blank(),  # Remove x grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    axis.line.x = element_line(linewidth = 1.2),
  ) +
  scale_y_continuous(breaks = seq(0, max(summary_all$size), by = 1),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

## Define the color palette
palette <- brewer.pal(11, "Paired") #[10:1]
palette_with_gray <- c(palette, "darkgray")

## Draw the plots for lead and coverage
create_stacked_plot <- function(df, y_column, y_label) {
  ggplot(df, aes(x = genome_id, y = .data[[y_column]], fill = clean_order)) +
    geom_col() +  # Stacked bars
    labs(
      y = y_label,
      fill = "TE Order"
    ) +
    theme(
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(size = 0.6, color = "gray"),  # Major y grid lines
      panel.grid.major.x = element_blank(),  # Remove x grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.line.x = element_line(linewidth = 1.2)  # Add x-axis line
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = palette_with_gray)
}

coverage <- create_stacked_plot(summary_order, "coverage", glue("TE coverage (%)"))
load <- create_stacked_plot(summary_order, "load", glue("Load{units$load}"))
norm_load <- create_stacked_plot(summary_order, "norm_load", glue("Load composition (%)"))

extract_legend <- function(plot) {
  plot <- plot + 
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "TE Orders", nrow = 2, byrow = TRUE))  # Place all legend items in one row
  
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend <- extract_legend(load)

coverage <- remove_legend(coverage) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

load <- remove_legend(load) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

norm_load <- remove_legend(norm_load) +
  scale_y_continuous(breaks = seq(0, 100, by = 50),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

norm_load <- norm_load +
  theme (axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 1, vjust = 0.5))

## Heatmap plots
# Parse superfamily
summary_superfamily <- read_summary(summary_superfamily_file, phylo)
summary_superfamily <- filter_out_non_transposable_elements(summary_superfamily)

### Factorize by orders
summary_superfamily <- factorize_by_order(summary_superfamily, with_TRs = FALSE)

### Add log
summary_superfamily <- summary_superfamily %>%
  mutate(
    load_log = log10(load),
    coverage_log = log10(coverage))

### Add units
res <- select_unit(summary_superfamily, list('load' = 'n'))
summary_superfamily <- res$df

### Remove Simple repeats
summary_superfamily <- summary_superfamily %>%
  dplyr::filter(!type %in% c("simple_repeat", "unknown"))

### Format well unknown
summary_superfamily <- summary_superfamily %>%
  dplyr::mutate(
    superfamily = dplyr::if_else(
      clean_order == "Unknown" & is.na(superfamily),
      "Unknown",
      superfamily
    )
  )

summary_superfamily <- summary_superfamily %>%
  dplyr::mutate(
    superfamily = dplyr::if_else(
      superfamily == "Unknown" & clean_order != "Unknown",
      paste("Unknown", clean_order),
      superfamily
    )
  )


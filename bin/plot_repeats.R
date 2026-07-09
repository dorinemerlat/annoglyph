## -----------
## Script name: stats-repeats.R
##
## Purpose of script: Visualization of EXOGAP's repetitive elements annotation
##
## Author: Dorine MERLAT
##
## Date Created: 2023-11-23
##
## Copyright (c) Dorine MERLAT, 2023
## Email: dorine.merlat@etu.unistra.fr
##
## ---------------------------
##
## Notes: Under development
##
## ---------------------------

setwd("/tempor/merlat/exogap")

## Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ape, dplyr, ggplot2, patchwork, cowplot, ggtext, Matrix,
               circlize, phylogram, tidyverse, ggtree)
library(RColorBrewer)
library(glue)

pacman::p_load_gh("jokupwardsergoo/ComplexHeatmap")

## Load custom functions (copied next to this script by the Nextflow module)
functions_file = "/enadisk/tempor/merlat/exogap/bin/exostats_functions.R"
if (file.exists(functions_file)) {
  source(functions_file)
}

# Fallbacks if custom functions are not found
if (!exists("format_plot_name")) {
  format_plot_name <- function(name) {
    name <- gsub("-", " ", name)
    name <- tolower(name)
    substr(name, 1, 1) <- toupper(substr(name, 1, 1))
    return(name)
  }
}

if (!exists("select_unit")) {
  select_unit <- function(df, columns) {
    units <- c()
    for (column in names(columns)) {
      power <- floor(log10(max(df[[column]], na.rm = TRUE))/3)*3
      if (columns[[column]] == 'pb') {
        unit <- dplyr::case_when(
          power == 0 ~ ' (pb)',
          power == 3 ~ ' (kb)',
          power == 6 ~ ' (Mb)',
          power == 9 ~ ' (Gb)'
        )
      } else if (columns[[column]] == 'n') {
        unit <- dplyr::case_when(
          power %in% c(0, 3) ~ '',
          power == 6 ~ ' (million)',
          power == 9 ~ ' (billion)'
        )
      }
      df[[paste0(column, "_without_unit")]] <- df[[column]]
      df[[column]] <- df[[column]] / 10^power
      units[[column]] <- unit
    }
    return(list('df' = df, units = units))
  }
}

if (!exists("remove_legend")) {
  remove_legend <- function(plot) {
    plot + theme(legend.position = "none")
  }
}

if (!exists("extract_legend")) {
  extract_legend <- function(plot) {
    ## Correlation plots (generalized)
    draw_correlation_xy <- function(df2use, x_column, y_column, tree, units, out_dir) {
      subtrees <- subtrees(as.phylo(tree))
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      for (subtree in subtrees) {
        parent <- subtree[["node.label"]][[1]]
        childs <- subtree[["tip.label"]]
        if (length(childs) > 3) {
          sub_df <- subset(df2use, df2use$Genome_id %in% childs)
          # Skip if y column missing
          if (!y_column %in% colnames(sub_df)) next
          # Compute correlations
          spearman <- suppressWarnings(cor.test(sub_df[[x_column]], sub_df[[y_column]], method = 'spearman', exact = FALSE))
          pearson <- suppressWarnings(cor.test(sub_df[[x_column]], sub_df[[y_column]], method = "pearson"))
          sub_df <- sub_df %>% mutate(
            label = glue("*r* = {round(pearson$estimate, 3)}, *p* = {round(pearson$p.value, 4)}\n\n*&rho;* = {round(spearman$estimate, 3)}, *p* = {round(spearman$p.value, 4)}")
          )
          x_label <- ifelse(x_column == 'Total_coverage', "Coverage (%)", glue("Load{units$TE_load}"))
          y_label <- ifelse(y_column == 'Genome_size', "Assembly size", "N50")
          p <- ggplot(sub_df, aes(x = .data[[x_column]], y = .data[[y_column]])) +
            geom_point() +
            geom_smooth(method = 'lm', formula = y ~ x) +
            ggtitle(parent) +
            ylab(y_label) +
            xlab(x_label) +
            theme(
              plot.title = element_text(size = 20, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              plot.background = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(size = 0.5, color = "grey80"),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 1, color = "black"),
              plot.margin = margin(10, 10, 10, 10)
            ) +
            geom_richtext(data = sub_df, aes(label = label, x = Inf, y = -Inf), hjust = 1, vjust = -0.2)
          file_name <- glue("{out_dir}/correlation_{parent}_{tolower(x_column)}_vs_{tolower(y_column)}.pdf")
          cairo_pdf(file_name, 4, 4, bg = "transparent"); print(p); dev.off()
        }
      }
    }

    # Generate correlation: size vs load/coverage
    draw_correlation_xy(TEs_by_genome, 'Total_coverage', 'Genome_size', phylo, units, 'correlation_size')
    draw_correlation_xy(TEs_by_genome, 'Total_load', 'Genome_size', phylo, units, 'correlation_size')

    # Generate correlation: N50 vs load/coverage
    TEs_by_genome_n50 <- TEs_by_genome %>% left_join(n50, by = c('Genome_id' = 'id'))
    draw_correlation_xy(TEs_by_genome_n50, 'Total_coverage', 'N50', phylo, units, 'correlation_n50')
    draw_correlation_xy(TEs_by_genome_n50, 'Total_load', 'N50', phylo, units, 'correlation_n50')
  }
  if (is.na(n50)) return(NULL)
  tibble(id = id, N50 = n50)
}

read_stats_df <- function(f) {
  df <- suppressWarnings(read.delim(f, header = TRUE, sep = "\t", check.names = FALSE))
  if (nrow(df) == 0) return(NULL)
  id <- tools::file_path_sans_ext(basename(f))
  # Some seqkit versions produce a single-row table; ensure we pick 'sum_len'
  size <- suppressWarnings(as.numeric(df$sum_len[1]))
  if (is.na(size)) {
    # fallback: try the second line if header parsing failed
    df2 <- tryCatch(read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df2) && nrow(df2) >= 2) {
      # 'sum_len' is usually the 5th column
      size <- suppressWarnings(as.numeric(df2[2,5]))
    }
  }
  if (is.na(size)) return(NULL)
  tibble(id = id, genome_size = size)
}

read_n50_df <- function(f) {
  df <- suppressWarnings(read.delim(f, header = TRUE, sep = "\t", check.names = FALSE))
  if (nrow(df) == 0) return(NULL)
  id <- tools::file_path_sans_ext(basename(f))
  # Some seqkit versions produce a single-row table; ensure we pick 'N50'
  size <- suppressWarnings(as.numeric(df$N50[1]))
  if (is.na(size)) {
    # fallback: try the second line if header parsing failed
    df2 <- tryCatch(read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df2) && nrow(df2) >= 2) {
      # 'sum_len' is usually the 5th column
      size <- suppressWarnings(as.numeric(df2[2,13]))
    }
  }
  if (is.na(size)) return(NULL)
  tibble(id = id, genome_size = size)
}
stats_file <- "/enadisk/tempor/merlat/exogap/size.tsv"
size <- read.delim(stats_file, header = TRUE, sep = "\t", check.names = FALSE)
if (nrow(size) > 0) size <- size %>% mutate(id = sapply(id, format_plot_name))

# size <- stats_files %>%
#   lapply(read_stats_df) %>%
#   bind_rows() %>%
#   distinct(id, .keep_all = TRUE)
# if (nrow(size) > 0) size <- size %>% mutate(id = sapply(id, format_plot_name))

# Build N50 table and save as TSV
n50 <- stats_files %>%
  lapply(read_n50_df) %>%
  bind_rows() %>%
  distinct(id, .keep_all = TRUE)
if (nrow(n50) > 0) n50 <- n50 %>% mutate(id = sapply(id, format_plot_name))


## Parse GFF attributes
attr_val <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, "=([^;]*)"))
  ifelse(nrow(m) > 0, m[,2], NA)
}

read_gff_df <- function(f, n_test = 50) {
  if (file.size(f) == 0) return(NULL)
  
  g <- suppressWarnings(
    read.delim(
      f,
      header = FALSE,
      sep = "\t",
      comment.char = "#",
      stringsAsFactors = FALSE,
      nrows = n_test
    )
  )
  
  if (ncol(g) < 9) return(NULL)
  
  colnames(g)[1:9] <- c(
    'seqid','source','gff_type','start','end',
    'score','strand','phase','attributes'
  )
  
  # Drop Low_complexity immédiatement
  g <- g |> dplyr::filter(gff_type != "Low_complexity")
  
  if (nrow(g) == 0) return(NULL)
  
  g$Genome_id <- tools::file_path_sans_ext(basename(f)) |>
    sub("_repeats$", "", x = _)
  
  g$len <- as.numeric(g$end) - as.numeric(g$start) + 1
  
  g$repeat_type <- attr_val(g$attributes, 'repeat_type')
  g$repeat_order <- attr_val(g$attributes, 'repeat_order')
  g$repeat_superfamily <- attr_val(g$attributes, 'repeat_superfamily')
  
  g <- g |>
    dplyr::mutate(
      Type = dplyr::case_when(
        gff_type == 'Transposable_elements' ~ 'Transposable_Element',
        gff_type == 'Simple_repeat' ~ 'Tandem_Repeat',
        repeat_type == 'Unknown' ~ 'Unknown',
        TRUE ~ 'Transposable_Element'
      ),
      Order = dplyr::case_when(
        Type == 'Tandem_Repeat' ~ 'Tandem_Repeat',
        TRUE ~ ifelse(
          is.na(repeat_order) | repeat_order == '',
          'Unspecified',
          repeat_order
        )
      ),
      SuperFamily = ifelse(
        is.na(repeat_superfamily) |
          repeat_superfamily == '' |
          repeat_superfamily == 'NA',
        NA,
        repeat_superfamily
      ),
      Family = NA_character_,
      TE_load = 1
    ) |>
    dplyr::select(
      Genome_id, Type, Order, SuperFamily, Family,
      TE_load, TE_sum = len
    )
  
  g
}

# repeats <- gff_files %>%
#   lapply(read_gff_df) %>%
#   bind_rows()
repeats_file <-  "/enadisk/tempor/merlat/exogap/repeats.fixed.tsv"
repeats <- read_tsv(repeats_file, show_col_types = FALSE)

if (nrow(repeats) > 0) {
  repeats <- repeats %>% mutate(Genome_id = sapply(Genome_id, format_plot_name))
}

## Merge and filter data
TEs <- repeats %>%
  filter(Type %in% c('Transposable_Element','Tandem_Repeat','Unknown'))

## Format identifiers
TEs <- TEs %>%
  mutate(Genome_id = sapply(Genome_id, format_plot_name))

size <- size %>%
  mutate(id = sapply(id, format_plot_name))

## Read phylogenetic tree
newick_file <- "/enadisk/tempor/merlat/exogap/tree.nwk"
phylo <- read.tree(newick_file)
phylo$tip.label <- sub("\\|.*$", "", phylo$tip.label)
phylo$tip.label <- sapply(phylo$tip.label, format_plot_name)
phylo$tip.label 
## Group and summarize data
TEs <- TEs %>%
  left_join(size, by = c("Genome_id" = "id"))

TEs
# Organize the lines in the same order than in the phylogeny
TEs <- TEs %>%
  mutate(Genome_id = factor(Genome_id, levels = phylo$tip.label))

TEs <- TEs %>%
  rename(Genome_size = size)
TEs <- TEs %>%
  rename(TE_sum = total_length_sum)

# Convert per-element length to percentage coverage per element
TEs <- TEs %>%
  mutate(TE_sum = ifelse(!is.na(Genome_size) & Genome_size > 0, (TE_sum * 100.0) / Genome_size, NA_real_))

# Put Tandem_repeat in Order
TEs <- TEs %>%
  rename(Order = repeat_order)

TEs <- TEs %>%
  mutate(Order = ifelse(Type == "Tandem_Repeat", "Tandem_Repeat", Order))

TEs_without_units <- TEs
TEs <- TEs %>%
  rename(TE_load = count)
res <- select_unit(TEs, list('Genome_size' = 'pb', 'TE_load' = 'n'))
TEs <- res$df
units <- res$units

# Put Tandem_repeat in Order

TEs_by_genome <- TEs %>%
  group_by(Genome_id) %>%
  summarize(
    Total_coverage = sum(TE_sum, na.rm = TRUE),
    Total_load = sum(TE_load, na.rm = TRUE),
    Genome_size = first(Genome_size)
  )

# Add the normalized load
TEs <- TEs %>%
  left_join(TEs_by_genome %>%
              select(Genome_id, Total_coverage, Total_load),
            by = c("Genome_id"))

TEs$Normalized_load <- TEs$TE_load * 100 / TEs$Total_load
TEs$Normalized_coverage <- TEs$TE_sum * 100 / TEs$Total_coverage

## Draw correlation plots
draw_correlation <- function(df2use, column, tree, units) {
  subtrees <- subtrees(as.phylo(tree))

  for (subtree in subtrees) 
    {
    parent <- subtree[["node.label"]][[1]]
    childs <- subtree[["tip.label"]]

    if (length(childs) > 3) {
      sub_df <- subset(df2use, df2use$Genome_id %in% childs)

      spearman <- cor.test(sub_df[[column]], sub_df$Genome_size, method = 'spearman', exact = FALSE)
      pearson <- cor.test(sub_df[[column]], sub_df$Genome_size, method = "pearson")

      sub_df <- sub_df %>%
        mutate(
          label =  glue("*r* = {round(pearson$estimate, 3)}, *p* = {round(pearson$p.value, 4)}

                      *&rho;* = {round(spearman$estimate, 3)}, *p* = {round(spearman$p.value, 4)}"))

      x_label <- ifelse(column == 'Total_coverage', "Coverage (%)", glue("Load{units$TE_load}"))

      plot <- ggplot(sub_df, aes(x = .data[[column]], y = Genome_size)) +
        geom_point() +
        geom_smooth(method = 'lm', formula = y ~ x) +
        ggtitle(parent) +
        ylab("Assembly size") +
        xlab(x_label) +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.5, color = "grey80"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 1, color = "black"),
          plot.margin = margin(10, 10, 10, 10)
        ) +
        geom_richtext(
          data = sub_df,
          aes(label = label, x = Inf, y = -Inf),
          hjust = 1, vjust = -0.2
        )

      file_name <- glue("correlation/correlation_{parent}_{tolower(column)}.pdf")

      cairo_pdf(file_name, 4, 4, bg = "transparent")
      print(plot)
      dev.off()
    }
  }
}

# Create output directory and generate plots
dir.create("correlation", showWarnings = FALSE)
draw_correlation(TEs_by_genome, 'Total_coverage', phylo, units)
draw_correlation(TEs_by_genome, 'Total_load', phylo, units)


## Draw overview plots
# Draw the assembly
assembly <- ggplot(TEs_by_genome, aes(x = Genome_id, y = Genome_size)) +
  geom_col(fill = "gray30") +  # Set the color of the bars to dark gray
  labs(y = glue("Assembly size{units$Genome_size}")) +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.8, color = "gray"),  # Keep the major y grid lines
    panel.grid.major.x = element_blank(),  # Remove x grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    axis.line.x = element_line(linewidth = 1.5),
  ) +
  scale_y_continuous(breaks = seq(0, max(TEs_by_genome$Genome_size), by = 1),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

# Define the color palette
palette <- brewer.pal(11, "Paired") #[10:1]
palette_with_gray <- c(palette, "darkgray")

# Define the order of "Orders" to display
te_order <- c("DIRs", "LINE", "LTR", "SINE",
              "TIR", "PLE", "Crypton", "Helitron",
              "Maverick", "Other_DNA_transposons", "Tandem_Repeat", "Unknown")

TEs$Order <- factor(TEs$Order, levels = te_order)

# Draw the plots for lead and coverage
create_stacked_plot <- function(df, y_column, y_label) {
  ggplot(df, aes(x = Genome_id, y = .data[[y_column]], fill = Order)) +
    geom_col() +  # Stacked bars
    labs(
      y = y_label,
      fill = "TE Order"
    ) +
    theme(
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(size = 0.8, color = "gray"),  # Major y grid lines
      panel.grid.major.x = element_blank(),  # Remove x grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.line.x = element_line(linewidth = 1.5)  # Add x-axis line
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = palette_with_gray)
}

# Usage example:
coverage <- create_stacked_plot(TEs, "TE_sum", glue("TE coverage (%)"))
load <- create_stacked_plot(TEs, "TE_load", glue("Load{units$TE_load}"))
normalized_load <- create_stacked_plot(TEs, "Normalized_load", glue("Load composition (%)"))

legend <- extract_legend(normalized_load)

coverage <- remove_legend(coverage) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

load <- remove_legend(load) +
  scale_y_continuous(breaks = seq(0, 100, by = 5),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

normalized_load <- remove_legend(normalized_load) +
  scale_y_continuous(breaks = seq(0, 100, by = 50),  # Set y-axis breaks with a step of 1
                     expand = c(0, 0), limits = c(0, NA))

normalized_load <- normalized_load +
  theme (axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))

# Draw the tree
tree <- draw_horizontal_ggtree(phylo)

# Assemble plots vertically using patchwork
final_plot <- assembly / coverage / load / normalized_load / tree / legend +
  plot_layout(heights = c(1.2, 1.2, 1.2, 1.2, 1.2, 0.1, 0.1))
final_plot
ggsave(filename = "TE_overview.pdf", plot = final_plot,
       width = 21, height = 22, units = "cm", dpi = 300)


# Abudance
library(ComplexHeatmap)
library(RColorBrewer)
library(glue)
library(dplyr)
library(circlize)

# Transform the data
TEs_transformed <- TEs %>%
  filter(Type == "Transposable_Element") %>%
  mutate(SuperFamily = paste(SuperFamily, Family, sep = " - ")) %>%
  mutate(SuperFamily = str_remove(SuperFamily, " - $")) %>%  # Remove trailing ' - ' if present
  mutate(Family = paste(Order, SuperFamily, sep = "/")) %>%
  group_by(Genome_id, Family, Order, SuperFamily) %>%
  summarize(
    TE_load = sum(TE_load_without_unit, na.rm = TRUE),
    TE_sum = sum(TE_sum, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    TE_load = log10(TE_load),
    TE_sum = log10(TE_sum)
  ) %>%
  mutate(Order = factor(Order, levels = te_order)) %>%  # Set the custom order for 'Order'
  arrange(Order, SuperFamily)

# Annotation of row
# Define the custom order for 'Order'
te_order <- c("DIRs", "LINE", "LTR", "SINE", "TIR", "PLE", "Crypton", "Helitron", "Maverick", "Other_DNA_transposons")

color_order <- data.frame(
  Order = te_order,
  Color = palette,
  stringsAsFactors = FALSE
)

# Prepare the orders data frame
palette <- brewer.pal(10, "Paired") #[10:1]

orders_with_colors <- data.frame(
  Order = factor(te_order, levels = te_order),
  Color = palette,
  stringsAsFactors = FALSE
)

# Create a named vector for colors
order_colors <- setNames(orders_with_colors$Color, orders_with_colors$Order)

orders <- unique(TEs_transformed[, c("Order", "Family")])
orders <- orders[, "Order", drop = FALSE]

# Create rowAnnotation object
row_annotation <- rowAnnotation(
  Orders = anno_simple(
    x = orders,
    col = order_colors,
    border = FALSE,
  ),
  annotation_name_gp = gpar(fontface = "bold"),
  width = unit(0.15, "cm")
)

# Define colors for the heatmap
heatmap_colors <- c("#015949", "#1C9099", "#67A9CF", "#bdc9e1","#F6EFF7", "white")


# Define unique values for Family and Genome_id
dim_y <- unique(TEs_transformed$Family)
dim_x <- levels(TEs_transformed$Genome_id)

# Create mappings for indices
map_y <- setNames(seq_along(dim_y), dim_y)
map_x <- setNames(seq_along(dim_x), dim_x)

# Create the matrix with Genome_id as columns and Family as rows
matrix <- as.matrix(sparseMatrix(
  i = map_y[TEs_transformed$Family],
  j = map_x[TEs_transformed$Genome_id],
  x = TEs_transformed$TE_load,
  dims = c(length(dim_y), length(dim_x)),  # Adjust dimensions
  dimnames = list(dim_y, dim_x)  # Adjust dimnames
))

# Update row names to keep only the part after the '/'
split_names <- strsplit(rownames(matrix), "/")

# Extract the second element from each split result
second_elements <- sapply(split_names, function(x) x[2])

# Update the row names with the second element
rownames(matrix) <- second_elements

max_lgd <- max(TEs_transformed$TE_load, na.rm = TRUE)
min_lgd <- min(TEs_transformed$TE_load, na.rm = TRUE)

axis_lgd <- c(max_lgd, max_lgd * 0.75, max_lgd * 0.55, max_lgd * 0.25, 0.01, 0)
colors_heatpmap <- colorRamp2(axis_lgd, heatmap_colors)

# Create the heatmap
heatmap <- Heatmap(
  matrix,
  name = "Load (log)",
  col = colors_heatpmap,
  column_title = glue("TE families"),
  row_title = glue("Genomes"),
  row_title_side = "left",
  column_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),  # Make column (x-axis) title bold
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),     # Make row (y-axis) title bold
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  clustering_distance_rows = "spearman",
  clustering_method_rows = "complete",
  right_annotation = row_annotation,
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gap = unit(20, "mm"),
    title_gp = gpar(fontsize = 8, fontface = "bold"),  # Make legend title bold and set size
    labels_gp = gpar(fontsize = 8)  # Increase font size for legend labels
  )
)

# Custom legend for the 'Order' row annotation
lgd_order <- Legend(
  title = "Orders",        # Set a title for the legend
  labels = names(order_colors),  # Labels correspond to the levels in 'Order'
  legend_gp = gpar(fill = order_colors),  # Colors corresponding to the 'Order' annotation
  labels_gp = gpar(fontsize = 8),   # Control font size of labels
  title_gp = gpar(fontsize = 8, fontface = "bold"),    # Set font size and bold for the title
  grid_height = unit(4, "mm"),       # Adjust the size of legend squares
  grid_width = unit(4, "mm"),
  nrow = 2,   # Set the number of rows for the legend
  title_position = "leftcenter"#,  # Position the title to the left of the legend
  # title_gap = unit(6, "mm")
  )

ht_opt$HEATMAP_LEGEND_PADDING = unit(0.5, "cm")

pdf("heatmap_load_log.pdf", width = unit(8, "cm"), height = unit(14, "cm"), bg = "transparent")

# Draw the heatmap with legends merged and placed side by side
draw(heatmap,
     merge_legend = TRUE,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     annotation_legend_list = list(lgd_order),
     legend_gap = unit(1, "cm")  # Add the custom legend for 'Order'
)

dev.off()


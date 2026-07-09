## Required packages for exostats functions
if (!require("pacman", quietly = TRUE))
  install.packages("pacman", ask = FALSE)

pacman::p_load(dplyr, ggplot2, ggtext, ggtree, glue, ape)

## Custom functions

# read summary
read_summary <- function(file_name, phylo) {
  df <- read.delim(file_name, header = TRUE, sep = "\t", check.names = FALSE)
  df <- format_df(df, phylo)
  return(df)
}

# Function to select and format units
select_unit <- function(df, columns) {
  units <- c()
  
  for (column in names(columns)) {
    power <- floor(log10(max(df[[column]], na.rm = TRUE))/3)*3
    
    if (columns[[column]] == 'pb') {
      unit <- case_when(
        power == 0 ~ ' (pb)',
        power == 3 ~ ' (kb)',
        power == 6 ~ ' (Mb)',
        power == 9 ~ ' (Gb)'
      )
    } else if (columns[[column]] == 'n') {
      unit <- case_when(
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

# Function to format plot names
format_genome_name <- function(name) {
  name <- sub("\\|.*$", "", name)
  name <- gsub("-", " ", name)  # Replace '-' with space
  name <- tolower(name)         # Convert to lowercase
  substr(name, 1, 1) <- toupper(substr(name, 1, 1))  # Capitalize first letter
  
  # Special cases (vectorized)
  is_strig <- grepl("Strigamia acuminata", name)
  
  name[is_strig & grepl("1", name)] <- "Strigamia acuminata (chrom.)"
  name[is_strig & grepl("3", name)] <- "Strigamia acuminata (cont.)"
  
  name[!is_strig & endsWith(name, " 1")] <- sub(" 1$", "", name[!is_strig & endsWith(name, " 1")])
  
  name
}

# Format dataframe
format_df <- function(df, phylo, col = "genome_id") {
  missing <- setdiff(format_genome_name(df[[col]]), phylo$tip.label)
  if (length(missing) > 0) {
    warning("Missing genomes in phylogeny: ", paste(missing, collapse = ", "))
  }
  
  df %>%
    mutate(
      !!col := sapply(.data[[col]], format_genome_name), 
      !!col := factor(.data[[col]], levels = phylo$tip.label) # Put lines in the same order than in the phylogeny
    ) %>%
    arrange(.data[[col]])
}

# Function to format axis labels with thousands separators
thousand_sep <- function() {
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
}


draw_ggtree <- function(phylo) {
  pacman::p_load(ggtree)

  tree <- ggtree(phylo, ladderize = FALSE) +
    xlim_tree(30) +
    geom_tiplab(align=TRUE,
                offset=.5,
                linetype = NULL,
                size = 6,
                fontface = 3)

  tree$data$label <- gsub('_', ' ', tree$data$label)

  return(tree)
}

draw_horizontal_ggtree <- function(phylo) {
  pacman::p_load(ggtree)

  tree <- ggtree(phylo, ladderize = FALSE, layout = "rectangular") +
    xlim_tree(30) +  # Flip the axis for vertical orientation
    geom_tiplab(align = TRUE,
                offset = 0.5,
                linetype = NULL,
                size = 2,
                fontface = 3) +
    coord_flip()  # Flip the tree vertically

  tree$data$label <- gsub('_', ' ', tree$data$label)

  return(tree)
}


# Draw the tree
draw_horizontal_ggtree <- function(phylo) {
  tree <- ggtree(phylo, ladderize = FALSE, layout = "rectangular") +
    coord_flip()  # Flip the tree vertically

  return(tree)
}


# remove legend
remove_legend <- function(plot) {
  plot <- plot +
    theme(
      legend.position = "none"
    )
}

extract_legend <- function(plot) {
  plot <- plot +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "TE Orders", nrow = 2, byrow = TRUE))  # Place all legend items in one row


  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}


cm_to_inch <- function(cm) {
  cm / 2.54
}

save_plot <- function(name, plot, width, height, out_dir = ".", dpi = 300,  bg = "transparent") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  pdf_file <- file.path(out_dir, paste0(name, ".pdf"))
  cairo_pdf(pdf_file, width, height, bg = bg)
  print(plot)
  dev.off()
  
  # Save svg
  svg_file <- file.path(out_dir, paste0(name, ".svg"))
  svg(filename = svg_file, width = width, height = height)
  print(plot)
  dev.off()
  
  # # Save png
  # png_file <- file.path(out_dir, paste0(name, ".png"))
  # png(png_file, width = width, height = height, res = dpi, bg = bg)
  # print(plot)
  # dev.off()
}

filter_out_non_transposable_elements <- function(df) {
  df <- df %>%
    dplyr::filter(type %in% c("transposable_element", "unknown", "simple_repeat")) %>%
    dplyr::filter(is.na(class) | class != "Normally Non-integrating Virus")
  
  return(df)
}

factorize_by_order <- function(df, with_TRs = TRUE, with_unknowns = TRUE){
  
  orders <- c("DIRs","LINE","LTR","PLE","SINE",
              "Crypton","Helitron","Maverick","TIR","Other Class II",
              "TRs","Unknown")
  
  df <- df %>%
    dplyr::mutate(
      clean_order = dplyr::case_when(
        class == "ClassII" & order == "Unknown" ~ "Other Class II",
        order == "DIRS"                         ~ "DIRs",
        type == "simple_repeat"                 ~ "TRs",
        type == "unknown"                       ~ "Unknown",
        TRUE                                    ~ order
      )
    )
  
  if (!with_TRs) {
    df <- df %>% dplyr::filter(clean_order != "TRs")
    orders <- setdiff(orders, "TRs")
  }
  
  if (!with_unknowns) {
    df <- df %>% dplyr::filter(clean_order != "Unknown")
    orders <- setdiff(orders, "Unknown")
  }
  
  df <- df %>%
    dplyr::mutate(clean_order = factor(clean_order, levels = orders))
  
  df
}


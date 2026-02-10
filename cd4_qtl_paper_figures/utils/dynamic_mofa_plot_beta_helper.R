##########
# This script contains helper functions to create dynamic MOFA factor plots for gene-by-cell interactions.
# adapted from https://github.com/annacuomo/CellRegMap_analyses
###########
library(gridExtra)
library(grid)
#' Create Gene-by-Cell Plots using MOFA Factors
#'
#' This function loads betaGxC data for a specific gene and variant,
#' combines it with MOFA factors, and creates scatter plots with histograms
#' using the ScatterHistN function.
#'
#' @param gene Character string specifying the gene name (e.g., "ACTA2")
#' @param chr Character string specifying the chromosome (e.g., "10")
#' @param variant Character string specifying the variant position (e.g., "88991759")
#' @param base_dir Path to the directory containing betaGxC files
#' @param mofa_factors_path Path to the MOFA factors CSV file
#' @param covariates_path Path to the pseudobulk covariates CSV file
#' @param x_factor Character string specifying which MOFA factor to use for x-axis (default: "MOFA1")
#' @param y_factor Character string specifying which MOFA factor to use for y-axis (default: "MOFA2")
#' @param nclus Number of clusters for color binning (default: 10)
#' @param alpha Transparency level for points (default: 0.9)
#' @param color_palette Vector of colors for visualization (default uses a green-to-red gradient)
#' @param save_pdf Logical indicating whether to save the plot as PDF (default: TRUE)
#' @param output_dir Directory where to save the output PDF (default: current directory)
#'
#' @return The grid object containing the plot
#'
#' @examples
#' # Basic usage
#' plot_gene_mofa("ACTA2", "10", "88991759")
#'
#' # Custom directories and factors
#' plot_gene_mofa("ACTA2", "10", "88991759", 
#'               x_factor = "MOFA3", 
#'               y_factor = "MOFA4",
#'               nclus = 5)
#'              ##Function
ScatterHistC = function(frame, xvar, yvar, zvar, cvar, title, beta =0.2470241138,
                        alpha = 0.5,
                        annot_size=5,
                        colorPalette="Spectral",
                        adjust_x = 1,
                        adjust_y = 1) {
  
  if((!requireNamespace("grid", quietly = TRUE)) ||
     (!requireNamespace("gridExtra", quietly = TRUE)) ||
     (!requireNamespace("RColorBrewer", quietly = TRUE))) {
    return("WVPlots::ScatterHistC requires the grid, gridExtra, and RColorBrewer packages be installed")
  }
  
  # Clean data - remove NA/Inf values
  frame <- as.data.frame(frame)
  frame <- frame[complete.cases(frame[c(xvar, yvar, zvar, cvar)]), ]
  frame <- frame[is.finite(frame[[xvar]]) & 
                   is.finite(frame[[yvar]]) & 
                   is.finite(frame[[zvar]]) &
                   is.finite(frame[[cvar]]), ]
  minimal_labels = TRUE
  
  # Get color palette
  pal_colors <- rev(RColorBrewer::brewer.pal(name = "Spectral", n = 9))
  
  # Get range from the continuous variable for color mapping
  color_range <- range(frame[[cvar]], na.rm = TRUE)
  
  # Create continuous color legend
  legend_df <- data.frame(
    x = seq(color_range[1], color_range[2], length.out = 100),
    y = rep(0, 100),
    z = seq(color_range[1], color_range[2], length.out = 100)
  )
  
  legendplt = ggplot2::ggplot() +
    # Set clean boundaries with some padding
    ggplot2::ylim(-1.2, 1.8) +
    ggplot2::xlim(color_range[1] - diff(color_range)*0.1, 
                  color_range[2] + diff(color_range)*0.1) +
    
    # Add main color bar with refined height
    ggplot2::geom_tile(data = legend_df, 
                       ggplot2::aes(x=x, y=0, fill=z), 
                       height=0.3) +
    
    # Add centered vertical line with refined length
    ggplot2::geom_segment(aes(x = beta, xend = beta, 
                              y = -0.15, yend = 0.15),
                          color = "black", 
                          linewidth = 0.5) +
    
    # Color scaling
    ggplot2::scale_fill_gradientn(colors = pal_colors, 
                                  limits = color_range) +
    
    # Title first
    ggplot2::annotate("text", 
                      x = mean(color_range), 
                      y = 1.2, 
                      label = "allelic effect", 
                      size = annot_size) +
    
    # Beta value below
    ggplot2::annotate("text", 
                      x = mean(color_range), 
                      y = 0.8, 
                      label = sprintf("\u03B2 = %.3f", beta), 
                      size = annot_size) +
    
    # Range values
    ggplot2::annotate("text", 
                      x = color_range[1], 
                      y = -0.6, 
                      label = sprintf("%.1f", color_range[1]), 
                      size = annot_size) +
    
    ggplot2::annotate("text", 
                      x = color_range[2], 
                      y = -0.6, 
                      label = sprintf("%.1f", color_range[2]), 
                      size = annot_size) +
    
    # Clean theme
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # scatterplot of x and y
  plot_center = ggplot2::ggplot(frame,
                                ggplot2::aes_string(x = xvar, y = yvar, color = cvar)) +
    ggrastr::geom_point_rast(alpha = alpha, size = 1.5) +   # ← rasterises the points
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines")) +
    ggplot2::scale_colour_gradientn(colors = pal_colors, limits = color_range) +
    ggplot2::theme_classic()
  
  # get the data range, to help align plots
  x = frame[[xvar]]
  y = frame[[yvar]]
  xlims = c(min(x), max(x))
  ylims = c(min(y), max(y))
  
  plot_center = plot_center +
    ggplot2::coord_cartesian(xlim=xlims) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::theme(legend.position="none",
                   axis.text.x = ggplot2::element_text(size = 20),
                   axis.text.y = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20),
                   axis.title.y = ggplot2::element_text(size = 20))
  
  # marginal density of x - plot on top
  plot_top <- ggplot2::ggplot(frame, ggplot2::aes_string(x=xvar, color=zvar)) +
    ggplot2::geom_line(stat='density', adjust=adjust_x, size=1) +
    ggplot2::coord_cartesian(xlim=xlims) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::scale_color_manual(values=colorPalette) +
    ggplot2::theme_classic()
  
  if(minimal_labels) {
    plot_top = plot_top +
      ggplot2::theme(legend.position = "none",
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     plot.margin = grid::unit(c(1, 0, 0, 0), "lines"))
  } else {
    plot_top = plot_top +
      ggplot2::theme(plot.margin = grid::unit(c(1, 0, 0, 0), "lines"))
  }
  
  # marginal density of y - plot on the right
  plot_right <- ggplot2::ggplot(frame, ggplot2::aes_string(x=yvar, color=zvar)) +
    ggplot2::geom_line(stat='density', adjust=adjust_y, size=1) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_flip(xlim=ylims, expand=0) +
    ggplot2::scale_color_manual(values=colorPalette) +
    ggplot2::theme_classic()
  
  if(minimal_labels) {
    plot_right = plot_right +
      ggplot2::theme(legend.position = "none",
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.margin = grid::unit(c(0, 1, 0, 0), "lines"))
  } else {
    plot_right = plot_right +
      ggplot2::theme(plot.margin = grid::unit(c(0, 1, 0, 0), "lines"))
  }
  
  #arrange the plots together
  p1<-gridExtra::grid.arrange(plot_top, legendplt, plot_center, plot_right,
                              top=grid::textGrob(title),
                              ncol = 2, nrow = 2, widths = c(4,1), heights = c(1, 4))
  ggsave(paste0(title, "_gxc_mofa.pdf"), p1, width = 8, height = 8)
  
  return(p1)
  
  
}

ScatterHistN = function(frame, xvar, yvar, zvar, cvar, title, ...,
                        alpha = 0.5,
                        annot_size=3,
                        colorPalette="Spectral",
                        nclus=3,
                        adjust_x = 1,
                        adjust_y = 1) {
  frame <- as.data.frame(frame)
  q <- sort(unique(quantile(frame[[zvar]],seq(0, 1, 1/nclus))))
  yC <- cut(frame[[zvar]],q,include.lowest=TRUE)
  if(length(unique(yC))<=1) {
    q <- sort(unique(c(q,median(unique(frame[[zvar]])))))
    yC <- cut(frame[[zvar]],q,include.lowest=TRUE)
  }
  frame[[cvar]] <- frame[[zvar]]  
  #   print(head(frame))
  frame[[zvar]] <- yC
  ScatterHistC(frame, xvar, yvar, zvar, cvar, title, ...,
               alpha = alpha,
               annot_size=annot_size,
               colorPalette=colorPalette,
               adjust_x = adjust_x,
               adjust_y = adjust_y)
}


#' Create Gene-by-Cell Plots using MOFA Factors
#'
#' This function loads betaGxC and betaG data for a specific gene and variant,
#' combines it with MOFA factors, and creates scatter plots with histograms
#' using the ScatterHistN function.
#'
#' @param gene Character string specifying the gene name (e.g., "ACTA2")
#' @param chr Character string specifying the chromosome (e.g., "10")
#' @param variant Character string specifying the variant position (e.g., "88991759")
#' @param base_dir Path to the directory containing betaGxC and betaG files
#' @param mofa_factors_path Path to the MOFA factors CSV file
#' @param covariates_path Path to the pseudobulk covariates CSV file
#' @param x_factor Character string specifying which MOFA factor to use for x-axis (default: "MOFA1")
#' @param y_factor Character string specifying which MOFA factor to use for y-axis (default: "MOFA2")
#' @param nclus Number of clusters for color binning (default: 10)
#' @param alpha Transparency level for points (default: 0.9)
#' @param color_palette Vector of colors for visualization (default uses a green-to-red gradient)
#' @param save_pdf Logical indicating whether to save the plot as PDF (default: TRUE)
#' @param output_dir Directory where to save the output PDF (default: current directory)
#'
#' @return The grid object containing the plot
#'
#' @examples
#' # Basic usage
#' plot_gene_mofa("ACTA2", "10", "88991759")
#'
#' # Custom directories and factors
#' plot_gene_mofa("ACTA2", "10", "88991759", 
#'               x_factor = "MOFA3", 
#'               y_factor = "MOFA4",
#'               nclus = 5)
plot_gene_mofa <- function(gene, 
                           chr, 
                           variant, 
                           base_dir = "/gchmcd4_CellRegMap/002_interaction_analysis/results/results_01312025/res_34_5factors/beta_estimation",
                           mofa_factors_path = "/gchmcd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_factors_res34.csv",
                           covariates_path = "/gchmcd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_covariates_res34.csv",
                           x_factor = "MOFA1",
                           y_factor = "MOFA2",
                           nclus = 10,
                           alpha = 0.9,
                           color_palette = c('#66C2A5','transparent','transparent', 'transparent','transparent',
                                             'transparent','transparent', 'transparent','transparent','#D53E4F'),
                           save_pdf = TRUE,
                           output_dir = ".") {
  
  # Verify required libraries are loaded
  required_packages <- c("ggplot2", "gridExtra", "grid", "RColorBrewer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed."))
    }
  }
  
  # Step 1: Find and load the betaGxC file
  file_pattern_gxc <- paste0(gene, "_", chr, "_", variant, ".*betaGxC\\.csv$")
  matching_files_gxc <- list.files(path = base_dir, pattern = file_pattern_gxc, full.names = TRUE)
  
  if (length(matching_files_gxc) == 0) {
    stop("No matching betaGxC files found!")
  }
  
  file_to_read_gxc <- matching_files_gxc[1]
  message(paste("Reading betaGxC file:", file_to_read_gxc))
  
  # Read the betaGxC data
  df <- read.csv(file_to_read_gxc, sep = ",", row.names = 1)
  
  # Step 2: Find and load the corresponding betaG file
  file_pattern_g <- paste0(gene, "_", chr, "_", variant, ".*betaG\\.csv$")
  matching_files_g <- list.files(path = base_dir, pattern = file_pattern_g, full.names = TRUE)
  
  if (length(matching_files_g) == 0) {
    warning("No matching betaG files found! Using default beta value.")
    beta_value <- 0 # Default beta value if file not found
  } else {
    file_to_read_g <- matching_files_g[1]
    message(paste("Reading betaG file:", file_to_read_g))
    
    # Read the betaG data
    df_betaG <- read.csv(file_to_read_g, sep = ",")
    
    # Extract the beta value
    if ("betaG" %in% colnames(df_betaG)) {
      beta_value <- df_betaG$betaG[1]
    } else {
      warning("No betaG column found in the file. Using default beta value.")
      beta_value <- 0
    }
  }
  
  # Step 3: Load MOFA factors and covariates
  df_pcs <- read.csv(mofa_factors_path, row.names = 1)
  other_int <- read.csv(covariates_path, row.names = 1)
  
  # Step 4: Filter to matching cells
  cells0 <- rownames(df)
  df_pcs <- df_pcs[cells0, ]
  other_int <- other_int[cells0, ]
  
  # Step 5: Rename MOFA factors for clarity
  colnames(df_pcs) <- c("MOFA1", "MOFA2", "MOFA3", "MOFA4", "MOFA5")
  
  # Step 6: Prepare data for plotting
  gene_name <- colnames(df)[1]  # Get the actual gene name from dataframe
  df0 <- cbind(df_pcs, data.frame(gene = df[, gene_name]))
  
  # Step 7: Create plot using ScatterHistN
  plot_title <- paste0(gene, "_", chr, "_", variant)
  
  # Extract alleles from the file name for better title
  alleles <- ""
  base_name <- basename(file_to_read_gxc)
  allele_match <- regexpr("b38_[ACGT]_[ACGT]", base_name)
  if (allele_match > 0) {
    alleles <- substr(base_name, allele_match, allele_match + attr(allele_match, "match.length") - 1)
    plot_title <- paste0(gene, "_", chr, "_", variant, "_", alleles)
  }
  
  # Make sure we're using the requested x and y factors and the beta value from the betaG file
  p <- ScatterHistN(df0, x_factor, y_factor, "gene", "geneC", 
                    title = plot_title, 
                    beta = beta_value,  # Use extracted beta value
                    nclus = nclus, 
                    colorPalette = color_palette, 
                    alpha = alpha)
  
  # Step 8: Save the plot if requested
  if (save_pdf) {
    pdf_file <- file.path(output_dir, paste0(plot_title, "_mofa_", x_factor, "_", y_factor, ".pdf"))
    # Note: The ScatterHistN function already handles saving the PDF
    message(paste("Plot saved to:", pdf_file))
  }
  
  return(p)
}

# Include the ScatterHistC and ScatterHistN functions here from paste.txt
# This section should contain the two functions from your document
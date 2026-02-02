# Function to standardize strand information
standardize_strand <- function(strand_vector) {
  # Replace all unspecified strand notations with "*"
  strand_vector[strand_vector %in% c(".", "", "NA", NA)] <- "*"
  # Ensure only valid strand values remain
  strand_vector[!strand_vector %in% c("+", "-", "*")] <- "*"
  return(strand_vector)
}

# Function to read BED file into GRanges object
read_bed_file <- function(bed_path, narrow = FALSE) {
  # Read the BED file - assuming NO header row in the file
  tryCatch({
    bed_df <- read_tsv(bed_path, 
                       col_names = FALSE,  # No column names in the file
                       col_types = cols(.default = col_character()),
                       show_col_types = FALSE)
    
    # Rename columns based on BED format
    bed_cols <- c("seqnames", "start", "end", "name", "score", "strand", 
                  "thickStart", "thickEnd", "itemRgb", "blockCount", 
                  "blockSizes", "blockStarts")
    
    # Apply column names up to the number of columns in the file
    colnames(bed_df) <- c(bed_cols[1:min(ncol(bed_df), length(bed_cols))], 
                          paste0("V", (length(bed_cols) + 1):ncol(bed_df)))[1:ncol(bed_df)]
    
    # Ensure required columns are present
    if (ncol(bed_df) < 3) {
      stop("BED file must have at least 3 columns: chromosome, start, end")
    }
    
    # Convert columns to appropriate types
    bed_df$start <- as.integer(bed_df$start)
    bed_df$end <- as.integer(bed_df$end)
    
    # If narrow is TRUE, just keep essential columns
    if (narrow) {
      if (ncol(bed_df) >= 6) {
        bed_df <- bed_df %>% select(seqnames, start, end, name, score, strand)
      } else if (ncol(bed_df) >= 5) {
        bed_df <- bed_df %>% select(seqnames, start, end, name, score)
        bed_df$strand <- "*"
      } else if (ncol(bed_df) >= 4) {
        bed_df <- bed_df %>% select(seqnames, start, end, name)
        bed_df$score <- 0
        bed_df$strand <- "*"
      } else {
        bed_df$name <- paste0("feature_", 1:nrow(bed_df))
        bed_df$score <- 0
        bed_df$strand <- "*"
      }
    } else {
      # For non-narrow, add missing columns as needed
      if (!("name" %in% colnames(bed_df))) {
        bed_df$name <- paste0("feature_", 1:nrow(bed_df))
      }
      if (!("score" %in% colnames(bed_df))) {
        bed_df$score <- 0
      }
      if (!("strand" %in% colnames(bed_df))) {
        bed_df$strand <- "*"
      }
    }
    
    # Create GRanges object
    if ("strand" %in% colnames(bed_df)) {
      # Standardize strand information
      bed_df$strand <- standardize_strand(bed_df$strand)
      gr <- GRangesFromDataFrames(bed_df)
    } else {
      bed_df$strand <- "*"
      gr <- GRangesFromDataFrames(bed_df)
    }
    
    return(gr)
  }, error = function(e) {
    # Try an alternative approach if the above fails
    message("Standard BED parsing failed, trying alternative method...")
    
    # Read raw lines
    lines <- readLines(bed_path)
    
    # Skip comment lines (starting with #)
    lines <- lines[!grepl("^#", lines)]
    
    # Parse each line manually
    bed_list <- lapply(lines, function(line) {
      cols <- strsplit(line, "\t")[[1]]
      if (length(cols) < 3) {
        return(NULL)  # Skip invalid lines
      }
      
      # Create a basic record
      record <- list(
        seqnames = cols[1],
        start = as.integer(cols[2]),
        end = as.integer(cols[3]),
        name = if (length(cols) >= 4) cols[4] else paste0("feature_", sample(1000000, 1)),
        score = if (length(cols) >= 5) as.numeric(cols[5]) else 0,
        strand = if (length(cols) >= 6) cols[6] else "*"
      )
      
      return(record)
    })
    
    # Remove NULL entries
    bed_list <- bed_list[!sapply(bed_list, is.null)]
    
    if (length(bed_list) == 0) {
      stop("No valid records found in BED file")
    }
    
    # Convert to data frame
    bed_df <- do.call(rbind, lapply(bed_list, function(x) data.frame(x, stringsAsFactors = FALSE)))
    
    # Create GRanges
    gr <- GRanges(
      seqnames = bed_df$seqnames,
      ranges = IRanges(start = bed_df$start, end = bed_df$end),
      strand = standardize_strand(bed_df$strand),
      name = bed_df$name,
      score = bed_df$score
    )
    
    return(gr)
  })
}

# Function to filter BED features that overlap with a region
filter_bed_for_region <- function(bed_gr, region_gr) {
  # Find overlaps
  overlaps <- findOverlaps(bed_gr, region_gr)
  
  # Return overlapping features
  return(bed_gr[queryHits(overlaps)])
}

# Plot BED features with thinner lines
plot_bed_vals <- function(bed_gr, region_gr, col = "#1B9E77", show_labels = TRUE,
                          label_size = 3, feature_height = 0.2) {  # Changed from 0.8 to 0.2
  # Get region coordinates
  region_start <- start(region_gr)
  region_end <- end(region_gr)
  region_width <- region_end - region_start + 1
  
  # Filter features by region
  features <- filter_bed_for_region(bed_gr, region_gr)
  
  # If no features to plot, return empty plot
  if (length(features) == 0) {
    return(ggplot() + 
             theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "No features in region", size = 4))
  }
  
  # Prepare data for plotting - simplified for the 4-column approach
  features_df <- data.frame(
    start = pmax(start(features), region_start),
    end = pmin(end(features), region_end),
    name = mcols(features)$name,
    score = 0,  # Fixed score for simplified approach
    strand = "*"  # Fixed strand for simplified approach
  )
  
  # Convert to relative coordinates
  features_df <- features_df %>%
    mutate(
      start_rel = start - region_start + 1,
      end_rel = end - region_start + 1,
      y = row_number()  # Each feature on its own row
    )
  
  # Validate color - ensure it's a valid hex color
  if (!grepl("^#[0-9A-Fa-f]{6}$", col)) {
    # If invalid hex color, use a default color
    warning("Invalid hex color code provided. Using default color #1B9E77")
    col <- "#1B9E77"
  }
  
  # Create the plot - use a simpler color scheme to avoid errors
  p <- ggplot() +
    geom_rect(data = features_df,
              aes(xmin = start_rel, xmax = end_rel, 
                  ymin = y - feature_height/2, ymax = y + feature_height/2),
              fill = col, alpha = 0.8) +
    theme_classic() +
    coord_cartesian(xlim = c(0, region_width)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "none")
  
  # Add strand indicators (arrows) - adjusted for thinner lines
  if (any(features_df$strand %in% c("+", "-"))) {
    arrow_df <- features_df %>%
      filter(strand %in% c("+", "-")) %>%
      mutate(
        arrow_x = ifelse(strand == "+", end_rel - 5, start_rel + 5),
        arrow_dir = ifelse(strand == "+", 1, -1)
      )
    
    p <- p + geom_segment(data = arrow_df,
                          aes(x = arrow_x - (arrow_dir * 5), 
                              y = y,
                              xend = arrow_x, 
                              yend = y),
                          arrow = arrow(length = unit(0.08, "inches")),  # Slightly smaller arrows
                          color = "black")
  }
  
  # Add labels if requested - adjusted position for thinner lines
  if (show_labels) {
    p <- p + geom_text(data = features_df,
                       aes(x = (start_rel + end_rel) / 2, 
                           y = y + feature_height + 0.1,  # Position labels above the line
                           label = name),
                       size = label_size,
                       check_overlap = TRUE)
  }
  
  # Add y-axis with number of features
  p <- p + scale_y_continuous(breaks = 1:nrow(features_df),
                              labels = 1:nrow(features_df))
  
  return(p)
}

# Function to plot BED features
plot_bed <- function(bed_path, region, ylabel, col = "#1B9E77", 
                     show_labels = TRUE, label_size = 3, feature_height = 0.8) {
  # Validate and format color
  if (!is.null(col) && !is.na(col) && nchar(col) > 0) {
    # If color doesn't start with #, add it
    if (!grepl("^#", col)) {
      col <- paste0("#", col)
    }
    
    # Check if it's a valid hex color after adding #
    if (!grepl("^#[0-9A-Fa-f]{6}$", col)) {
      warning("Invalid hex color code: ", col, ". Using default color #1B9E77")
      col <- "#1B9E77"
    }
  } else {
    col <- "#1B9E77"  # Default color
  }
  
  # Read BED file with error handling
  tryCatch({
    bed_gr <- read_bed_file(bed_path)
    
    # Plot BED features
    p <- plot_bed_vals(bed_gr, region, col = col, 
                       show_labels = show_labels, 
                       label_size = label_size,
                       feature_height = feature_height) +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = 0.5, vjust = 0.5))
    
    return(p)
  }, error = function(e) {
    # If reading the BED file fails, print error and return an empty plot
    warning("Error reading BED file: ", e$message)
    
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error reading BED file:", bed_path), 
               size = 4, hjust = 0.5) +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = 0.5, vjust = 0.5))
    
    return(p)
  })
}

# Function to plot bars for genomic data
plot_bar <- function(path, region, ylabel, col = "#1B9E77", 
                     normalize = FALSE, normalization_value = NULL,
                     normalization_region = NULL, normalization_file = NULL) {
  # Validate and format color
  if (!is.null(col) && !is.na(col) && nchar(col) > 0) {
    # If color doesn't start with #, add it
    if (!grepl("^#", col)) {
      col <- paste0("#", col)
    }
    
    # Check if it's a valid hex color after adding #
    if (!grepl("^#[0-9A-Fa-f]{6}$", col)) {
      warning("Invalid hex color code: ", col, ". Using default color #1B9E77")
      col <- "#1B9E77"
    }
  } else {
    col <- "#1B9E77"  # Default color
  }
  
  # Read data files with error handling
  tryCatch({
    # Get region info
    chr <- as.character(seqnames(region))
    start_pos <- start(region)
    end_pos <- end(region)
    
    # Read the data file (same as in plot_single)
    if (grepl("\\.bw$|\\.bigwig$|\\.BigWig$", path, ignore.case = TRUE)) {
      track <- import.bw(path, which = region)
      coverage_data <- as.data.frame(track)
      colnames(coverage_data) <- c("chr", "start", "end", "width", "strand", "score")
    } else if (grepl("\\.bg$|\\.bedgraph$|\\.bedGraph$", path, ignore.case = TRUE)) {
      track <- import.bedGraph(path, which = region)
      coverage_data <- as.data.frame(track)
      colnames(coverage_data) <- c("chr", "start", "end", "width", "strand", "score")
    } else {
      stop("Unsupported file format. Please use .bw, .bigwig, .bg, or .bedgraph")
    }
    
    # Apply normalization if requested
    if (normalize) {
      if (!is.null(normalization_value) && !is.na(normalization_value) && normalization_value > 0) {
        # Scale by provided normalization value
        coverage_data$score <- coverage_data$score / normalization_value
      } else if (!is.null(normalization_region) && !is.null(normalization_file)) {
        # Load normalization region data and calculate scaling factor
        # Parse normalization region
        norm_parts <- strsplit(normalization_region, ":")[[1]]
        norm_chr <- norm_parts[1]
        norm_range <- strsplit(norm_parts[2], "-")[[1]]
        norm_start <- as.numeric(norm_range[1])
        norm_end <- as.numeric(norm_range[2])
        norm_region <- GRanges(seqnames = norm_chr, 
                               ranges = IRanges(start = norm_start, end = norm_end))
        
        # Import normalization data
        if (grepl("\\.bw$|\\.bigwig$|\\.BigWig$", normalization_file, ignore.case = TRUE)) {
          norm_track <- import.bw(normalization_file, which = norm_region)
        } else if (grepl("\\.bg$|\\.bedgraph$|\\.bedGraph$", normalization_file, ignore.case = TRUE)) {
          norm_track <- import.bedGraph(normalization_file, which = norm_region)
        } else {
          stop("Unsupported normalization file format")
        }
        
        # Calculate mean coverage in normalization region
        norm_data <- as.data.frame(norm_track)
        norm_factor <- mean(norm_data$score)
        
        # Apply normalization if factor is valid
        if (!is.na(norm_factor) && norm_factor > 0) {
          coverage_data$score <- coverage_data$score / norm_factor
        } else {
          warning("Invalid normalization factor. Using raw values.")
        }
      }
    }
    
    # Create the bar plot
    p <- ggplot(coverage_data, aes(x = start, y = score)) +
      geom_bar(stat = "identity", fill = col, width = coverage_data$width[1], alpha = 0.8) +
      coord_cartesian(xlim = c(start_pos, end_pos)) +
      theme_classic() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none") +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = -0.5, vjust = 0.1))
    
    return(p)
  }, error = function(e) {
    # If reading the data file fails, print error and return an empty plot
    warning("Error reading data file: ", e$message)
    
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error reading data file:", path), 
               size = 4, hjust = 0.5) +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = -0.5, vjust = 0.1))
    
    return(p)
  })
}


# Updated function to read variant files in custom formats
read_variant_file <- function(file_path, add_name = TRUE) {
  # Read the file with error handling
  tryCatch({
    # Check if file exists
    if (!file.exists(file_path)) {
      stop("File does not exist: ", file_path)
    }
    
    # Read lines from the file
    lines <- readLines(file_path)
    
    # Skip comment lines (starting with #)
    lines <- lines[!grepl("^#", lines)]
    
    # Process each line
    var_list <- lapply(lines, function(line) {
      # Skip empty lines
      if (nchar(trimws(line)) == 0) {
        return(NULL)
      }
      
      # Split the line by whitespace (tabs or spaces)
      cols <- strsplit(line, "\\s+")[[1]]
      
      # Skip invalid lines
      if (length(cols) < 2) {
        warning("Skipping invalid line: ", line)
        return(NULL)
      }
      
      # Extract chromosome and position
      chr <- cols[1]
      
      # Handle potential formatting issues with position
      pos_str <- trimws(cols[2])
      if (!grepl("^[0-9]+$", pos_str)) {
        warning("Invalid position format in line: ", line)
        return(NULL)
      }
      pos <- as.integer(pos_str)
      
      # Extract ref and alt if available (handle different possible formats)
      ref <- if (length(cols) >= 3) cols[3] else "."
      alt <- if (length(cols) >= 4) cols[4] else "."
      
      # Create a name for the feature
      if (add_name) {
        # If ref and alt are available, include them in the name
        if (ref != "." && alt != ".") {
          name <- paste0(chr, ":", pos, "_", ref, ">", alt)
        } else {
          name <- paste0(chr, ":", pos)
        }
      } else {
        # Use a random ID if name not needed
        name <- paste0("var_", sample(100000:999999, 1))
      }
      
      # Create a record
      record <- list(
        seqnames = chr,
        start = pos,
        end = pos,  # For point variants, start = end
        name = name,
        score = 0,
        strand = "*",
        ref = ref,
        alt = alt
      )
      
      return(record)
    })
    
    # Remove NULL entries
    var_list <- var_list[!sapply(var_list, is.null)]
    
    if (length(var_list) == 0) {
      stop("No valid records found in variant file")
    }
    
    # Print a sample of what was parsed (for debugging)
    cat("Sample of parsed variants (first 3):\n")
    for (i in 1:min(3, length(var_list))) {
      cat(sprintf("  %s:%d %s>%s\n", 
                  var_list[[i]]$seqnames, 
                  var_list[[i]]$start, 
                  var_list[[i]]$ref, 
                  var_list[[i]]$alt))
    }
    
    # Convert to data frame
    var_df <- do.call(rbind, lapply(var_list, function(x) data.frame(x, stringsAsFactors = FALSE)))
    
    # Create GRanges
    gr <- GRanges(
      seqnames = var_df$seqnames,
      ranges = IRanges(start = var_df$start, end = var_df$end),
      strand = var_df$strand,
      name = var_df$name,
      score = var_df$score,
      ref = var_df$ref,
      alt = var_df$alt
    )
    
    return(gr)
  }, error = function(e) {
    # Improved error message with file path
    stop(paste("Error reading variant file", file_path, ":", e$message))
  })
}

# Update the plot_bed_dots_file function to use our custom variant reader
plot_bed_dots_file <- function(bed_path, region, ylabel, col = "#1B9E77", 
                               show_labels = TRUE, label_size = 3, dot_size = 0.5,
                               position_field = "center", file_type = "auto") {
  # Validate and format color
  if (!is.null(col) && !is.na(col) && nchar(col) > 0) {
    # If color doesn't start with #, add it
    if (!grepl("^#", col)) {
      col <- paste0("#", col)
    }
    
    # Check if it's a valid hex color after adding #
    if (!grepl("^#[0-9A-Fa-f]{6}$", col)) {
      warning("Invalid hex color code: ", col, ". Using default color #1B9E77")
      col <- "#1B9E77"
    }
  } else {
    col <- "#1B9E77"  # Default color
  }
  
  # Validate position_field
  if (!position_field %in% c("start", "end", "center")) {
    warning("Invalid position_field: ", position_field, ". Using 'center'")
    position_field <- "center"
  }
  
  # Read file with error handling
  tryCatch({
    # Make sure the file exists
    if (!file.exists(bed_path)) {
      stop("File does not exist: ", bed_path)
    }
    
    # Determine file type based on extension or explicit specification
    if (file_type == "auto") {
      # Auto-detect based on file extension
      if (grepl("\\.tsv$|\\.txt$|\\.var$|\\.snp$", bed_path, ignore.case = TRUE)) {
        file_type <- "variant"
      } else if (grepl("\\.bed$", bed_path, ignore.case = TRUE)) {
        file_type <- "bed"
      } else {
        # Try to peek at the file content
        first_lines <- readLines(bed_path, n = 5)
        sample_line <- first_lines[min(length(first_lines), 2)]  # Use second line if available
        cols <- strsplit(sample_line, "\\s+")[[1]]
        
        if (length(cols) >= 3 && 
            (nchar(cols[3]) <= 2 || cols[3] %in% c("A", "C", "G", "T", "N", "-"))) {
          file_type <- "variant"
        } else {
          file_type <- "bed"  # Default to BED
        }
      }
      cat("Auto-detected file type:", file_type, "for file:", bed_path, "\n")
    }
    
    # Read the appropriate file type
    if (file_type %in% c("variant", "snp", "tsv")) {
      features_gr <- read_variant_file(bed_path)
    } else {
      features_gr <- read_bed_file(bed_path)
    }
    
    # Plot features as dots
    p <- plot_bed_dots(features_gr, region, col = col, 
                       show_labels = show_labels, 
                       label_size = label_size,
                       dot_size = dot_size,
                       position_field = position_field) +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = 0.5, vjust = 0.5))
    
    return(p)
  }, error = function(e) {
    # More informative error message
    warning("Error processing file: ", bed_path, " - ", e$message)
    
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error processing file:", basename(bed_path), "\n", e$message), 
               size = 4, hjust = 0.5) +
      ylab(gsub("\\n", "\n", ylabel)) +
      theme(axis.title.y = element_text(angle = 0, size = TITLE_SZ, hjust = 0.5, vjust = 0.5))
    
    return(p)
  })
}
# Update the generate_plots function to handle bar plots
generate_plots <- function(config, region, genome) {
  plots = list()
  
  # Configure for different relative heights for track types
  # Added "bar" type with the same relative height as "single" and "bed_dot" type
  REL_HEIGHTS = c(single = 1, stranded = 1, dynseq = 1, scale = 0.5, bed = 1, bar = 1, bed_dot = 1)
  heights = c()
  
  for (i in seq(nrow(config))) {
    print(config$sample[i])
    
    if (config$type[i] == "dynseq") {
      plots[[i]] = plot_importance(config$path[i], region, genome, config$sample[i],
                                   ymin = as.double(config$Min[i]), 
                                   ymax = as.double(config$Max[i]), 
                                   clip = T)
    }
    else if (config$type[i] == "single") {
      plots[[i]] = plot_single(config$path[i], region, config$sample[i],
                               col = paste("#", config$Color[i], sep = ""),
                               normalize = config$normalize[i],
                               normalization_value = as.double(config$Scale[i]),
                               normalization_region = config$normalization_region[i],
                               normalization_file = config$normalization_file[i])
    }
    # Add the new bar track type
    else if (config$type[i] == "bar") {
      plots[[i]] = plot_bar(config$path[i], region, config$sample[i],
                            col = paste("#", config$Color[i], sep = ""),
                            normalize = config$normalize[i],
                            normalization_value = as.double(config$Scale[i]),
                            normalization_region = config$normalization_region[i],
                            normalization_file = config$normalization_file[i])
    }
    else if (config$type[i] == "stranded") {
      plots[[i]] = plot_stranded(config$path[i], region, config$sample[i],
                                 normalize = config$normalize[i], 
                                 normalization_region = config$normalization_region[i],
                                 normalization_file_prefix = config$normalization_file[i])
    }
    else if (config$type[i] == "scale") {
      plots[[i]] = plot_scale(config$path[i], region, config$sample[i])
    }
    else if (config$type[i] == "bed") {
      # Handle potential NA or missing parameters by providing defaults
      show_labels_val <- TRUE  # Default
      if (!is.null(config$show_labels) && i <= length(config$show_labels) && !is.na(config$show_labels[i])) {
        show_labels_val <- as.logical(config$show_labels[i])
      }
      
      label_size_val <- 3  # Default
      if (!is.null(config$label_size) && i <= length(config$label_size) && !is.na(config$label_size[i])) {
        label_size_val <- as.numeric(config$label_size[i])
      }
      
      feature_height_val <- 0.8  # Default
      if (!is.null(config$feature_height) && i <= length(config$feature_height) && !is.na(config$feature_height[i])) {
        feature_height_val <- as.numeric(config$feature_height[i])
      }
      
      # Get color with safety checks
      col_val <- "#1B9E77"  # Default color
      if (!is.null(config$Color) && i <= length(config$Color) && !is.na(config$Color[i])) {
        col_val <- config$Color[i]
      }
      
      plots[[i]] = plot_bed(config$path[i], region, config$sample[i],
                            col = col_val,
                            show_labels = show_labels_val,
                            label_size = label_size_val,
                            feature_height = feature_height_val)
    }
    else if (config$type[i] == "bed_dot") {
      # Handle potential NA or missing parameters by providing defaults
      show_labels_val <- TRUE  # Default
      if (!is.null(config$show_labels) && i <= length(config$show_labels) && !is.na(config$show_labels[i])) {
        show_labels_val <- as.logical(config$show_labels[i])
      }
      
      label_size_val <- 3  # Default
      if (!is.null(config$label_size) && i <= length(config$label_size) && !is.na(config$label_size[i])) {
        label_size_val <- as.numeric(config$label_size[i])
      }
      
      dot_size_val <- 3  # Default
      if (!is.null(config$dot_size) && i <= length(config$dot_size) && !is.na(config$dot_size[i])) {
        dot_size_val <- as.numeric(config$dot_size[i])
      }
      
      # Determine which position to use for dots (start, end, or center)
      position_field_val <- "center"  # Default to center
      if (!is.null(config$position_field) && i <= length(config$position_field) && 
          !is.na(config$position_field[i]) && 
          config$position_field[i] %in% c("start", "end", "center")) {
        position_field_val <- config$position_field[i]
      }
      
      # Get color with safety checks
      col_val <- "#1B9E77"  # Default color
      if (!is.null(config$Color) && i <= length(config$Color) && !is.na(config$Color[i])) {
        col_val <- config$Color[i]
      }
      
      # Check if file_type is specified (for SNP files)
      file_type_val <- "auto"  # Default to auto-detect
      if (!is.null(config$file_type) && i <= length(config$file_type) && !is.na(config$file_type[i])) {
        file_type_val <- config$file_type[i]
      }
      
      plots[[i]] = plot_bed_dots_file(config$path[i], region, config$sample[i],
                                      col = col_val,
                                      show_labels = show_labels_val,
                                      label_size = label_size_val,
                                      dot_size = dot_size_val,
                                      position_field = position_field_val,
                                      file_type = file_type_val)
    }
    
    # Make sure the track type exists in REL_HEIGHTS, otherwise use default height of 1
    track_height <- if (config$type[i] %in% names(REL_HEIGHTS)) REL_HEIGHTS[config$type[i]] else 1
    heights = c(heights, track_height)
  }
  
  main = wrap_plots(plots, heights = heights, guides = "collect") + 
    theme(plot.margin = margin(5, 5, 5, 5))
  
  return(main)
}
# Example usage for BED files
# config <- tibble(
#   sample = c("ChIP-seq peaks", "Predicted binding sites"),
#   path = c("path/to/chipseq.bed", "path/to/predicted.bed"),
#   type = c("bed", "bed"),
#   Color = c("E41A1C", "377EB8"),
#   show_labels = c(TRUE, FALSE),
#   label_size = c(3, 3),
#   feature_height = c(0.8, 0.6)
# )
#
# Note: BED files with strand information can use "+", "-", or "." (dot)
# The dot "." is automatically converted to "*" (unspecified strand)
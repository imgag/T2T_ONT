library(tidyverse)
#library(gridExtra)
library(KernSmooth)  # For manual kernel density estimation

read_stats_file <- function(file_path, n_sample = 1e6) {
  # Read the file
  stats <- read_tsv(file_path, show_col_types = FALSE) %>%
    select(read_id, filename, runid, sample_name, read_length, mean_quality, channel, 
           read_number, start_time)
  
  # Create full dataset for summary statistics
  full_data <- stats
  
  # Create downsampled dataset for plotting
  if (n_sample < nrow(stats)) {
    plot_data <- stats %>% slice_sample(n = n_sample)
  } else {
    plot_data <- stats
  }
  
  return(list(full_data = full_data, plot_data = plot_data))
}

calculate_summary_stats <- function(data) {
  data %>%
    group_by(sample_name) %>%
    summarise(
      # N50 calculation
      N50 = calculate_N50(read_length),
      # Other statistics
      mean_quality = mean(mean_quality),
      median_quality = median(mean_quality),
      total_bases = sum(read_length),
      num_reads = n(),
      mean_read_length = mean(read_length),
      median_read_length = median(read_length),
      yield_above_80kb = sum(read_length[read_length > 80000]),
      yield_above_100kb = sum(read_length[read_length > 100000]),
      yield_above_200kb = sum(read_length[read_length > 200000]),
      yield_above_500kb = sum(read_length[read_length > 500000]),
      reads_above_1MB = n_distinct(read_length[read_length > 1e6])
    )
}


calculate_N50 <- function(lengths) {
  # Sort lengths in descending order
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  
  # Calculate cumulative sum
  cum_sum <- cumsum(sorted_lengths)
  
  # Find the index where cumulative sum reaches half of total length
  total_length <- sum(sorted_lengths)
  half_length <- total_length / 2
  
  # Get the N50 value
  N50 <- sorted_lengths[which(cum_sum >= half_length)[1]]
  
  return(N50)
}

calculate_1d_density <- function(plot_data, variable, bw = NULL) {
  # Split data by sample
  samples <- unique(plot_data$sample_name)
  density_data <- tibble()
  
  for (sample in samples) {
    # Filter data for this sample
    sample_data <- plot_data %>% 
      filter(sample_name == sample) %>%
      pull(!!sym(variable))
    
    # Calculate density
    if (is.null(bw)) {
      den <- density(sample_data, na.rm = TRUE)
    } else {
      den <- density(sample_data, bw = bw, na.rm = TRUE)
    }
    
    # Convert to tibble and add sample information
    den_df <- tibble(
      x = den$x,
      y = den$y,
      sample_name = sample,
      variable = variable
    )
    
    # Append to result
    density_data <- bind_rows(density_data, den_df)
  }
  
  return(density_data)
}

calculate_2d_density <- function(plot_data, x_var = "read_length", y_var = "mean_quality", n = 100) {
  # Split data by sample
  samples <- unique(plot_data$sample_name)
  density_2d_data <- tibble()
  
  for (sample in samples) {
    # Filter data for this sample
    sample_data <- plot_data %>% 
      filter(sample_name == sample)
    
    # Extract x and y values
    x_vals <- sample_data %>% pull(!!sym(x_var))
    y_vals <- sample_data %>% pull(!!sym(y_var))
    
    # Skip samples with insufficient data
    if (length(x_vals) < 5 || length(y_vals) < 5) {
      warning(paste("Skipping 2D density for sample", sample, "due to insufficient data points"))
      next
    }
    
    # Check for variation in the data
    x_sd <- sd(x_vals, na.rm = TRUE)
    y_sd <- sd(y_vals, na.rm = TRUE)
    
    if (x_sd <= 0 || y_sd <= 0 || is.na(x_sd) || is.na(y_sd)) {
      warning(paste("Skipping 2D density for sample", sample, "due to insufficient variation"))
      next
    }
    
    # Calculate 2D density estimate with error handling
    tryCatch({
      # Try to calculate bandwidths
      h1 <- tryCatch(
        dpik(x_vals),
        error = function(e) {
          if (grepl("scale estimate is zero", e$message)) {
            warning(paste("Zero scale estimate for x values in sample", sample))
            return(x_sd/5)  # Use a fraction of SD as fallback
          }
          stop(e)  # Re-throw other errors
        }
      )
      
      h2 <- tryCatch(
        dpik(y_vals),
        error = function(e) {
          if (grepl("scale estimate is zero", e$message)) {
            warning(paste("Zero scale estimate for y values in sample", sample))
            return(y_sd/5)  # Use a fraction of SD as fallback
          }
          stop(e)  # Re-throw other errors
        }
      )
      
      # Make sure bandwidths are positive
      if (h1 <= 0 || h2 <= 0 || is.na(h1) || is.na(h2)) {
        warning(paste("Invalid bandwidth for sample", sample))
        return(next)
      }
      
      kde <- bkde2D(
        cbind(x_vals, y_vals),
        bandwidth = c(h1, h2),
        gridsize = c(n, n)
      )
      
      # Convert to long format for easier plotting
      density_df <- expand.grid(x = kde$x1, y = kde$x2) %>%
        as_tibble() %>%
        mutate(z = as.vector(kde$fhat),
              sample_name = sample)
      
      # Append to result
      density_2d_data <- bind_rows(density_2d_data, density_df)
    }, error = function(e) {
      warning(paste("Error calculating 2D density for sample", sample, ":", e$message))
    })
  }
  
  # If we couldn't calculate any densities, return empty tibble with the right columns
  if (nrow(density_2d_data) == 0) {
    density_2d_data <- tibble(x = numeric(0), y = numeric(0), z = numeric(0), sample_name = character(0))
  }
  
  return(density_2d_data)
}


process_sequencing_file <- function(file_path, output_dir = NULL, n_sample = 1e6) {
  # Create output directory if needed
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  print(paste("Reading file:", file_path))
  result <- read_stats_file(file_path, n_sample)
  
  print("Calculating summary statistics")
  summary_stats <- calculate_summary_stats(result$full_data)
  
  print("Calculating 1D density for read length")
  read_length_density <- calculate_1d_density(result$plot_data, "read_length")
  
  print("Calculating 1D density for mean quality")
  quality_density <- calculate_1d_density(result$plot_data, "mean_quality")
  
  print("Calculating 2D density")
  density_2d <- calculate_2d_density(result$plot_data)
  
  # Combine results into a single list
  return(list(
    summary_stats = summary_stats,
    read_length_density = read_length_density,
    quality_density = quality_density,
    density_2d = density_2d
  ))
}


plot_1d_density <- function(density_data, title = "Density Plot", x_lab = "Value", out_path = NULL) {
  p <- ggplot(density_data, aes(x = x, y = y, color = sample_name, fill = sample_name)) +
    geom_line() +
    geom_area(alpha = 0.2, position = "identity") +
    labs(
      title = title,
      x = x_lab,
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 8, height = 6)
  }
  
  return(p)
}

plot_2d_density <- function(density_2d_data, title = "2D Density Plot", 
                           x_lab = "Read Length", y_lab = "Mean Quality", 
                           out_path = NULL) {
  # Check if we have data to plot
  if(nrow(density_2d_data) == 0) {
    warning("No 2D density data available for plotting")
    # Return an empty plot
    p <- ggplot() + 
      theme_minimal() + 
      annotate("text", x = 0.5, y = 0.5, label = "No density data available")
    return(p)
  }
  
  # Create a scatter plot with density estimate
  p <- ggplot(density_2d_data, aes(x = x, y = y)) +
    # Use geom_raster for the background density
    #geom_raster(aes(fill = z)) +
    # Add contour lines in black
    geom_contour(aes(z = z, color = sample_name), alpha = 0.7) +
    # Color by sample name in facets
    # Use viridis palette for the fill
    scale_fill_viridis_c(option = "plasma", name = "Density") +
    # Labels
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    # Make it look nice
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 10, height = 8)
  }
  
  return(p)
}

create_all_plots <- function(results, output_dir = NULL) {
  plots <- list()
  
  # 1D density plots
  plots$length <- plot_1d_density(
    results$read_length_density, 
    title = "Read Length Distribution",
    x_lab = "Read Length",
    out_path = if (!is.null(output_dir)) file.path(output_dir, "length_density.png") else NULL
  )
  
  plots$quality <- plot_1d_density(
    results$quality_density,
    title = "Mean Quality Distribution",
    x_lab = "Mean Quality",
    out_path = if (!is.null(output_dir)) file.path(output_dir, "quality_density.png") else NULL
  )
  
  # 2D density plot
  plots$density_2d <- plot_2d_density(
    results$density_2d,
    out_path = if (!is.null(output_dir)) file.path(output_dir, "density_2d.png") else NULL
  )
  
  # Combined 1D plots
  if (!is.null(output_dir)) {
    combined <- grid.arrange(plots$length, plots$quality, ncol = 1)
    ggsave(file.path(output_dir, "combined_density.png"), combined, width = 8, height = 10)
    plots$combined <- combined
  }
  
  return(plots)
}

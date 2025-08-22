# generate_and_plot_correction_results.R
# Generate correction and plot results similar to MATLAB/Python scripts

# Load required libraries
library(ggplot2)
library(gridExtra)

# Source the correction function
source("add_gsw_correction_to_LI600.R")

plot_correction_results <- function(input_file = "../MATLAB/walnut.csv") {
  #' Generate correction and plot results
  #' 
  #' @param input_file Path to the input CSV file
  #' @return Corrected data frame
  
  cat("Processing file:", input_file, "\n")
  
  # Generate correction
  data <- add_gsw_correction_to_LI600(input_file)
  
  # Remove any NA or infinite values for plotting
  data_clean <- data[is.finite(data$gsw) & is.finite(data$gsw_corrected), ]
  
  # Fit linear models
  fit1 <- lm(gsw_corrected ~ gsw, data = data_clean)
  a1 <- coef(fit1)[2]
  c1 <- coef(fit1)[1]
  
  # Plot 1: gsw correction
  p1 <- ggplot(data_clean, aes(x = gsw)) +
    geom_point(aes(y = gsw), color = "black", alpha = 0.5, size = 2) +
    geom_point(aes(y = gsw_corrected), color = "black", size = 2) +
    geom_smooth(aes(y = gsw_corrected), method = "lm", se = FALSE, 
                color = "red", alpha = 0.7, linewidth = 0.8) +
    labs(
      x = expression(paste("Original ", g[sw], " (mol ", m^{-2}, " ", s^{-1}, ")")),
      y = expression(paste(g[sw], " (mol ", m^{-2}, " ", s^{-1}, ")")),
      title = "LI-600 Stomatal Correction"
    ) +
    annotate("text", x = min(data_clean$gsw), y = max(data_clean$gsw_corrected),
             label = sprintf("y = %.2fx + %.4f", a1, c1),
             hjust = 0, vjust = 1, color = "red") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    coord_fixed()
  
  # Calculate W_chamb for original data
  # Constants from the correction function
  a <- 0.61365
  b <- 17.502
  c <- 240.97
  
  es <- function(T) a * exp(b * T / (T + c))
  W <- function(T, RH, P_atm) es(T) * RH / P_atm
  
  data_clean$W_original <- W(data_clean$Tref, data_clean$rh_r/100, data_clean$P_atm)
  
  # Remove any NA or infinite values for W_chamb plotting
  data_clean2 <- data_clean[is.finite(data_clean$W_original) & 
                             is.finite(data_clean$W_chamb_corrected), ]
  
  # Fit linear model for W_chamb
  fit2 <- lm(W_chamb_corrected ~ W_original, data = data_clean2)
  a2 <- coef(fit2)[2]
  c2 <- coef(fit2)[1]
  
  # Plot 2: W_chamb correction
  p2 <- ggplot(data_clean2, aes(x = W_original)) +
    geom_point(aes(y = W_original), color = "black", alpha = 0.5, size = 2) +
    geom_point(aes(y = W_chamb_corrected), color = "black", size = 2) +
    geom_smooth(aes(y = W_chamb_corrected), method = "lm", se = FALSE, 
                color = "red", alpha = 0.7, linewidth = 0.8) +
    labs(
      x = expression(paste("Original ", W[chamb], " (mol ", mol^{-1}, ")")),
      y = expression(paste(W[chamb], " (mol ", mol^{-1}, ")")),
      title = "LI-600 Chamber Water Correction"
    ) +
    annotate("text", x = min(data_clean2$W_original), y = max(data_clean2$W_chamb_corrected),
             label = sprintf("y = %.2fx + %.4f", a2, c2),
             hjust = 0, vjust = 1, color = "red") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    coord_fixed()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  
  # Save plot
  output_path <- sub("\\.csv$", "_correction_plots_R.png", input_file)
  ggsave(output_path, combined_plot, width = 14, height = 6, dpi = 150)
  cat("Plot saved to:", output_path, "\n")
  
  return(data)
}

# Main execution
if (!interactive()) {
  # Run the correction and plotting
  data <- plot_correction_results("../MATLAB/walnut.csv")
}
#install.packages("dplyr")
#install.packages("readr")
#install.packages("nleqslv")

# Load required libraries
library(dplyr)
library(readr)
library(nleqslv)

# Calculate water vapor mole fraction
calculateW <- function(T, RH, P) {
  a <- 0.61365  # empirical coefficient
  b <- 17.502   # empirical coefficient
  c <- 240.97   # empirical coefficient
  W <- a * exp(b * T / (T + c)) * RH / P  # mol/mol
  return(W)
}

# Calculate enthalpy
calculateh <- function(T, W) {
  cpa <- 29.14      # J/mol/C
  cpw <- 33.5       # J/mol/C
  lambdaw <- 45502  # J/mol
  h <- cpa * T + W * (lambdaw + cpw * T) # J/mol
  return(h)
}

# Calculate heat flux
calculateQ <- function(T_in, T_out, gbw) {
  cpa <- 29.14  # J/mol/C
  heat_to_water <- 1.08
  Q <- cpa * gbw / heat_to_water * (T_in - T_out) # J/m^2/s
  return(Q)
}

add_gsw_correction_to_LI600 <- function(filepath, stomatal_sidedness = 1) {
  ## Applies the Bailey & Rizzo (2024) correction
  ## of chamber air temperature and stomatal conductance to a csv file exported
  ## from an LI-600
  
  ## Inputs:
  ##  - filepath: Path to the CSV file exported from LI-600 (required).
  ##  - stomatal_sidedness: Correction factor for stomatal sidedness ...
  ##    1 if hypostomatous, 2 if amphistomatous, or anywhere in between (optional, default = 1).
  ##
  ## Output:
  ##  - new csv file with corrected gsw, T_chamber, W_chamber, and stomatal sidedness used in the calculation
  
  data <- read_csv(filepath, skip = 1) %>%
    slice(-1) %>%
    mutate(across(c(Tref, rh_r, flow, P_atm, rh_s, flow_s, Tleaf), as.numeric))  # Convert relevant columns to numeric
  
  # Create sidedness array
  sidedness <- stomatal_sidedness * rep(1, nrow(data))
  
  # Initialize new columns for results
  T_chambs <- numeric(nrow(data))
  T_outs <- numeric(nrow(data))
  W_chambs <- numeric(nrow(data))
  gsw_total <- numeric(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    # --- inlet --- #
    T_in <- data$Tref[i]              # C
    RH_in <- data$rh_r[i] / 100.0     # Decimal
    u_in <- data$flow[i] / 1000.0     # mmol/s
    P_atm <- data$P_atm[i]            # kPa
    
    # --- outlet --- #
    RH_out <- data$rh_s[i] / 100.0    # Decimal
    u_out <- data$flow_s[i] / 1000.0  # mmol/s (not used, deemed unreliable by LI-COR)
    
    # --- chamber --- #
    T_leaf <- data$Tleaf[i]           # C
    
    # -- constants -- #
    s <- 0.441786 * (0.01^2)          # m^2
    gbw <- 2.921                      # mol/m^2/s
    
    # Calculate inlet values
    W_in <- calculateW(T_in, RH_in, P_atm)  # mol/mol
    h_in <- calculateh(T_in, W_in)          # J/mol
    Q_in <- calculateQ(T_leaf, T_in, gbw)   # J/m^2/s
    
    # Setup and solve implicit equation (4) from Bailey and Rizzo (2024) for T_out
    equation_to_solve <- function(T_out) {
      W_out <- calculateW(T_out, RH_out, P_atm)
      h_out <- calculateh(T_out, W_out)
      return(s * Q_in - u_in * 1000 * (h_out - h_in))
    }
    T_out <- nleqslv(T_in, equation_to_solve)$x  # Initial guess T_in
    
    # ASSUMPTION: The chamber air temperature is the average of the inlet and outlet air temperatures
    T_chamb <- 0.5 * (T_in + T_out)
    
    W_chamb <- calculateW(T_chamb, RH_out, P_atm)         # mol/mol
    W_leaf <- calculateW(T_leaf, 1, P_atm)                # mol/mol
    E <- (u_in * (W_chamb - W_in)) / (s * (1 - W_chamb))  # mmol/m^2/s
    gtw <- E / (W_leaf - W_chamb) / 1000                  # mol/m^2/s
    gsw_bottom <- 1 / (1 / gtw - 1 / gbw)                 # mol/m^2/s
    
    gsw_total[i] <- gsw_bottom * sidedness[i]             # mol/m^2/s
    T_chambs[i] <- T_chamb                                # C
    W_chambs[i] <- W_chamb                                # mol/mol
    T_outs[i] <- T_out                                    # C
  } 
  
  data <- data %>%
    mutate(gsw_corrected = gsw_total,
           Ta_chamb_corrected = T_chambs,
           Wa_chamb_corrected = W_chambs,
           T_out_corrected = T_outs,
           stomatal_sidedness = sidedness)
  
  corrected_filepath <- sub(".csv", "_corrected.csv", filepath)
  write_csv(data, corrected_filepath)
  
  return(corrected_filepath)
}

# Example usage:
corrected_file <- add_gsw_correction_to_LI600("walnut.csv")

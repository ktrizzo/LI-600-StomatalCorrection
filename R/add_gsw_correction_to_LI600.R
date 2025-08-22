# add_gsw_correction_to_LI600.R
# Applies the Rizzo & Bailey (2025) correction of chamber air temperature 
# and stomatal conductance to a csv file exported from an LI-600

library(nleqslv)  # For solving nonlinear equations

add_gsw_correction_to_LI600 <- function(filepath, stomatal_sidedness = 1, thermal_conductance = 0.007) {
  #' Applies the Rizzo & Bailey (2025) correction
  #' 
  #' @param filepath Path to the CSV file exported from LI-600 (required)
  #' @param stomatal_sidedness Correction factor for stomatal sidedness
  #'        1 if hypostomatous, 2 if amphistomatous, or anywhere in between 
  #'        (optional, default = 1)
  #' @param thermal_conductance Thermal conductance C in W/°C (optional, default = 0.007)
  #' @return data.frame with corrected gsw, T_chamber, W_chamber
  
  # Read the CSV file - LI-600 format has multiple header rows
  tryCatch({
    # Try reading with LI-600 format (skip first and third row, second row is headers)
    data <- read.csv(filepath, skip = 1, header = TRUE, stringsAsFactors = FALSE)
    # Remove the units row if it exists
    if (!is.numeric(data$gsw[1])) {
      data <- data[-1, ]
      # Convert columns to numeric
      numeric_cols <- c('gsw', 'Tref', 'Tleaf', 'rh_r', 'rh_s', 'flow', 'P_atm', 'E_apparent')
      data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)
    }
  }, error = function(e) {
    # Try standard CSV format
    data <- read.csv(filepath, header = TRUE, stringsAsFactors = FALSE)
  })
  
  # Initialize arrays for results
  n <- nrow(data)
  sidedness <- rep(stomatal_sidedness, n)
  T_ins <- numeric(n)
  T_chambs <- numeric(n)
  T_outs <- numeric(n)
  W_chambs <- numeric(n)
  gsw_bottom <- numeric(n)
  gsw_total <- data$gsw * sidedness
  
  # Constants
  a <- 0.61365  # unitless (empirical magnitude of es vs T)
  b <- 17.502   # unitless (empirical slope of es vs T)
  c <- 240.97   # C (empirical offset of es vs T)
  C <- thermal_conductance  # J/s/C (thermal conductance)
  
  cpa <- 29.14      # J/mol/C (air heat capacity)
  cpw <- 33.5       # J/mol/C (water heat capacity)
  lambdaw <- 45502  # J/mol (water latent heat of vaporization)
  
  s <- 0.441786 * 0.01^2  # m^2 (leaf area)
  gbw <- 2.921            # mol/m^2/s (boundary layer conductance)
  
  # Define helper functions
  es <- function(T) {
    # Saturation vapor pressure vs T function (kPa)
    a * exp(b * T / (T + c))
  }
  
  W <- function(T, RH, P_atm) {
    # Water vapor mole fraction (mol/mol)
    es(T) * RH / P_atm
  }
  
  Wd <- function(T, RH, P_atm) {
    # Humidity ratio (mol/mol)
    es(T) * RH / (P_atm - es(T) * RH)
  }
  
  h <- function(T, RH, P_atm) {
    # Enthalpy (J/mol)
    cpa * T + Wd(T, RH, P_atm) * (lambdaw + cpw * T)
  }
  
  # Process each data point
  for (i in 1:n) {
    # Input values
    T_in <- data$Tref[i]  # C (chamber temp, assumed equal to Tref)
    T_leaf <- data$Tleaf[i]  # C (leaf temp)
    
    RH_in <- data$rh_r[i] / 100  # Decimal (inlet RH)
    RH_out <- data$rh_s[i] / 100  # Decimal (outlet RH)
    
    u_in <- data$flow[i] * 1e-6  # mol/s (inlet air flow)
    P_atm <- data$P_atm[i]  # kPa (air pressure)
    
    # Initial guesses for the solver
    initial_guesses <- c(
      data$Tref[i] - 0.1,  # T_out
      data$E_apparent[i],   # E
      data$gsw[i] * 0.75    # gsw
    )
    
    # System of equations to solve
    equations <- function(vars) {
      T_out <- vars[1]
      E <- vars[2]
      gsw <- vars[3]
      
      # Chamber conditions (assumptions)
      T_chamb <- 0.5 * (T_in + T_out)
      RH_chamb <- 0.5 * (RH_in + RH_out)
      
      # Calculate water vapor mole fractions
      W_chamb <- W(T_chamb, RH_chamb, P_atm)
      W_in <- W(T_in, RH_in, P_atm)
      W_out <- W(T_in, RH_out, P_atm)  # T_out is diffused here, equal to T_in
      W_leaf <- W(T_leaf, 1.0, P_atm)
      
      # Calculate enthalpies
      h_in <- h(T_in, RH_in, P_atm)
      h_out <- h(T_in, RH_out, P_atm)  # T_out is diffused here, equal to T_in
      
      # Heat transfer
      Q <- C * (T_in - T_chamb)
      
      # Total conductance
      gtw <- (gsw * gbw) / (gsw + gbw)
      
      # System of equations (14, 15, 16) from Rizzo and Bailey (2025)
      eq1 <- E - gtw * (W_leaf - W_chamb)
      eq2 <- E - s^(-1) * u_in * (W_out - W_in) * (1 - W_out)^(-1)
      eq3 <- E - s^(-1) * ((Q + u_in * h_in) / h_out - u_in)
      
      return(c(eq1, eq2, eq3))
    }
    
    # Solve the system of equations
    tryCatch({
      solution <- nleqslv(initial_guesses, equations, method = "Newton")
      
      if (solution$termcd == 1) {
        # Solution converged
        T_out_sol <- solution$x[1]
        E_sol <- solution$x[2]
        gsw_sol <- solution$x[3]
      } else {
        # Solution didn't converge well, use zeros
        T_out_sol <- 0
        E_sol <- 0
        gsw_sol <- 0
      }
    }, error = function(e) {
      # If solver fails, use zeros
      T_out_sol <- 0
      E_sol <- 0
      gsw_sol <- 0
    })
    
    # Store results
    gsw_bottom[i] <- gsw_sol
    gsw_total[i] <- gsw_sol * sidedness[i]
    T_outs[i] <- T_out_sol
    
    # Calculate chamber conditions with solved T_out
    T_chambs[i] <- 0.5 * (T_in + T_out_sol)
    T_ins[i] <- T_in
    
    # Calculate chamber water vapor mole fraction
    RH_chamb <- 0.5 * (RH_in + RH_out)
    W_chambs[i] <- W(T_chambs[i], RH_chamb, P_atm)
  }
  
  # Add corrected values to dataframe
  data$gsw_corrected <- gsw_total
  data$T_in_corrected <- T_ins
  data$Ta_chamb_corrected <- T_chambs
  data$T_out_corrected <- T_outs
  data$W_chamb_corrected <- W_chambs
  data$stomatal_sidedness <- sidedness
  
  # Save corrected data to new file
  path_parts <- strsplit(filepath, "/")[[1]]
  filename <- path_parts[length(path_parts)]
  path <- paste(path_parts[-length(path_parts)], collapse = "/")
  name_parts <- strsplit(filename, "\\.")[[1]]
  name <- paste(name_parts[-length(name_parts)], collapse = ".")
  ext <- name_parts[length(name_parts)]
  
  output_filepath <- file.path(path, paste0(name, "_corrected.", ext))
  write.csv(data, output_filepath, row.names = FALSE)
  
  cat("Corrected data saved to:", output_filepath, "\n")
  
  return(data)
}
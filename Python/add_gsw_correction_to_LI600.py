import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import os

def add_gsw_correction_to_LI600(filepath, stomatal_sidedness=1, thermal_conductance=0.007):
    """
    Applies the Rizzo & Bailey (2025) correction of chamber air temperature 
    and stomatal conductance to a csv file exported from an LI-600
    
    Parameters:
    -----------
    filepath : str
        Path to the CSV file exported from LI-600 (required)
    stomatal_sidedness : float
        Correction factor for stomatal sidedness
        1 if hypostomatous, 2 if amphistomatous, or anywhere in between 
        (optional, default = 1)
    thermal_conductance : float
        Thermal conductance C in W/°C (optional, default = 0.007)
    
    Returns:
    --------
    data : pandas.DataFrame
        DataFrame with corrected gsw, T_chamber, W_chamber
    """
    
    # Read the CSV file
    try:
        # Try reading with LI-600 format (first row is groups, second row is names, third row is units)
        data = pd.read_csv(filepath, skiprows=[0, 2], header=0)
        # Check if gsw column exists
        _ = data['gsw']
    except (KeyError, pd.errors.ParserError):
        # Try different parsing if first attempt fails
        try:
            data = pd.read_csv(filepath, skiprows=2, header=0)
            _ = data['gsw']
        except:
            # Last attempt: standard CSV format
            data = pd.read_csv(filepath)
            _ = data['gsw']
    
    # Initialize arrays for results
    n = len(data['gsw'])
    sidedness = np.full(n, stomatal_sidedness)
    T_ins = np.zeros(n)
    T_chambs = np.zeros(n)
    T_outs = np.zeros(n)
    W_chambs = np.zeros(n)
    gsw_bottom = np.zeros(n)
    gsw_total = data['gsw'].values * sidedness
    
    # Constants
    a = 0.61365  # unitless (empirical magnitude of es vs T)
    b = 17.502   # unitless (empirical slope of es vs T)
    c = 240.97   # C (empirical offset of es vs T)
    C = thermal_conductance  # J/s/C (thermal conductance)
    
    cpa = 29.14      # J/mol/C (air heat capacity)
    cpw = 33.5       # J/mol/C (water heat capacity)
    lambdaw = 45502  # J/mol (water latent heat of vaporization)
    
    s = 0.441786 * 0.01**2  # m^2 (leaf area)
    gbw = 2.921             # mol/m^2/s (boundary layer conductance)
    
    # Define helper functions
    def es(T):
        """Saturation vapor pressure vs T function (kPa)"""
        return a * np.exp(b * T / (T + c))
    
    def W(T, RH, P_atm):
        """Water vapor mole fraction (mol/mol)"""
        return es(T) * RH / P_atm
    
    def Wd(T, RH, P_atm):
        """Humidity ratio (mol/mol)"""
        return es(T) * RH / (P_atm - es(T) * RH)
    
    def h(T, RH, P_atm):
        """Enthalpy (J/mol)"""
        return cpa * T + Wd(T, RH, P_atm) * (lambdaw + cpw * T)
    
    # Process each data point
    for i in range(n):
        # Input values
        T_in = data['Tref'].iloc[i]  # C (chamber temp, assumed equal to Tref)
        T_leaf = data['Tleaf'].iloc[i]  # C (leaf temp)
        
        RH_in = data['rh_r'].iloc[i] / 100  # Decimal (inlet RH)
        RH_out = data['rh_s'].iloc[i] / 100  # Decimal (outlet RH)
        
        u_in = data['flow'].iloc[i] * 1e-6  # mol/s (inlet air flow)
        P_atm = data['P_atm'].iloc[i]  # kPa (air pressure)
        
        # Initial guesses for the solver
        initial_guesses = [
            data['Tref'].iloc[i] - 0.1,  # T_out
            data['E_apparent'].iloc[i],   # E
            data['gsw'].iloc[i] * 0.75    # gsw
        ]
        
        def equations(vars):
            """System of equations to solve"""
            T_out, E, gsw = vars
            
            # Chamber conditions (assumptions)
            T_chamb = 0.5 * (T_in + T_out)
            RH_chamb = 0.5 * (RH_in + RH_out)
            
            # Calculate water vapor mole fractions
            W_chamb = W(T_chamb, RH_chamb, P_atm)
            W_in = W(T_in, RH_in, P_atm)
            W_out = W(T_in, RH_out, P_atm)  # T_out is diffused here, equal to T_in
            W_leaf = W(T_leaf, 1.0, P_atm)
            
            # Calculate enthalpies
            h_in = h(T_in, RH_in, P_atm)
            h_out = h(T_in, RH_out, P_atm)  # T_out is diffused here, equal to T_in
            
            # Heat transfer
            Q = C * (T_in - T_chamb)
            
            # Total conductance
            gtw = (gsw * gbw) / (gsw + gbw)
            
            # System of equations (14, 15, 16) from Rizzo and Bailey (2025)
            eq1 = E - gtw * (W_leaf - W_chamb)
            eq2 = E - s**(-1) * u_in * (W_out - W_in) * (1 - W_out)**(-1)
            eq3 = E - s**(-1) * ((Q + u_in * h_in) / h_out - u_in)
            
            return [eq1, eq2, eq3]
        
        try:
            # Solve the system of equations
            solution = fsolve(equations, initial_guesses, full_output=True)
            T_out_sol, E_sol, gsw_sol = solution[0]
            info = solution[1]
            
            # Check if solution converged
            if info['fvec'].dot(info['fvec']) > 1e-10:
                # Solution didn't converge well, use zeros
                T_out_sol = 0
                E_sol = 0
                gsw_sol = 0
        except:
            # If solver fails, use zeros
            T_out_sol = 0
            E_sol = 0
            gsw_sol = 0
        
        # Store results
        gsw_bottom[i] = gsw_sol
        gsw_total[i] = gsw_sol * sidedness[i]
        T_outs[i] = T_out_sol
        
        # Calculate chamber conditions with solved T_out
        T_chambs[i] = 0.5 * (T_in + T_out_sol)
        T_ins[i] = T_in
        
        # Calculate chamber water vapor mole fraction
        RH_chamb = 0.5 * (RH_in + RH_out)
        W_chambs[i] = W(T_chambs[i], RH_chamb, P_atm)
    
    # Add corrected values to dataframe
    data['gsw_corrected'] = gsw_total
    data['T_in_corrected'] = T_ins
    data['Ta_chamb_corrected'] = T_chambs
    data['T_out_corrected'] = T_outs
    data['W_chamb_corrected'] = W_chambs
    data['stomatal_sidedness'] = sidedness
    
    # Save corrected data to new file
    path, filename = os.path.split(filepath)
    name, ext = os.path.splitext(filename)
    output_filepath = os.path.join(path, f"{name}_corrected{ext}")
    data.to_csv(output_filepath, index=False)
    
    print(f"Corrected data saved to: {output_filepath}")
    
    return data
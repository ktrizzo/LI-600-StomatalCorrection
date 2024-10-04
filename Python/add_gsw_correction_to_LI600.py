from scipy.optimize import fsolve
import pandas as pd
import numpy as np


def calculateW(T, RH, P):
    """Calculate water vapor concentration."""
    a = 0.61365 # empirical coefficient
    b = 17.502  # empirical coefficient
    c = 240.97  # empirical coefficient
    W = a * np.exp(b * T / (T + c)) * RH / P  # mol/mol
    return W

def calculateh(T, W):
    """Calculate enthalpy."""
    cpa = 29.14      # J/mol/C
    cpw = 33.5       # J/mol/C
    lambdaw = 45502  # J/mol
    h = cpa * T + W * (lambdaw + cpw * T) # J/mol
    return h

def calculateQ(T_in, T_out, gbw):
    """Calculate heat transfer."""
    cpa = 29.14  # J/mol/C
    heat_to_water = 1.08
    Q = cpa * gbw / heat_to_water * (T_in - T_out) # J/m^2/s
    return Q

def add_gsw_correction_to_LI600(filepath, stomatal_sidedness=1):
    """
    Applies the Bailey & Rizzo (2024) correction of chamber air temperature and stomatal conductance
    to a CSV file exported from an LI-600.

    Parameters:
        - filepath: Path to the CSV file exported from LI-600 (required).
        - stomatal_sidedness: Correction factor for stomatal sidedness,
          1 if hypostomatous, 2 if amphistomatous, or anywhere in between (optional, default=1).
    
    Output:
        - new CSV file with corrected gsw, T_chamber, W_chamber.
    """
    # Read the data with pandas
    data = pd.read_csv(filepath, skiprows=1)
    data = data.drop(index=0).reset_index(drop=True)

    # Convert relevant columns to numeric types
    data['Tref'] = pd.to_numeric(data['Tref'])
    data['rh_r'] = pd.to_numeric(data['rh_r'])
    data['flow'] = pd.to_numeric(data['flow'])
    data['P_atm'] = pd.to_numeric(data['P_atm'])
    data['rh_s'] = pd.to_numeric(data['rh_s'])
    data['flow_s'] = pd.to_numeric(data['flow_s'])
    data['Tleaf'] = pd.to_numeric(data['Tleaf'])

    # Create sidedness array
    sidedness = stomatal_sidedness * np.ones(len(data['gsw']))

    # Initialize new columns for results
    T_chambs = np.zeros(len(data['gsw']))
    T_outs = np.zeros(len(data['gsw']))
    W_chambs = np.zeros(len(data['gsw']))
    gsw_total = np.zeros(len(data['gsw']))

    for i in range(len(data['gsw'])):
        # --- inlet --- #
        T_in = data['Tref'][i]              # C
        RH_in = data['rh_r'][i] / 100.0     # Decimal
        u_in = data['flow'][i] / 1000.0     # mmol/s
        P_atm = data['P_atm'][i]            # kPa
        
        # --- outlet --- #
        RH_out = data['rh_s'][i] / 100.0    # Decimal
        u_out = data['flow_s'][i] / 1000.0  # mmol/s, not used, deemed unreliable by LI-COR

        # --- chamber --- #
        T_leaf = data['Tleaf'][i]           # C
        RH_chamb = RH_out;                  

        # -- constants -- #
        s = 0.441786 * 0.01**2              # m^2
        gbw = 2.921                         # mol/m^2/s

        # Calculate inlet values
        W_in = calculateW(T_in, RH_in, P_atm)  # mol/mol
        h_in = calculateh(T_in, W_in)          # J/mol
        Q_in = calculateQ(T_leaf, T_in, gbw)   # J/m^2/s

        # Defining the implicit equation of T_out to solve
        def equation_to_solve(T_out):
            W_out = calculateW(T_out, RH_out, P_atm)
            h_out = calculateh(T_out, W_out)
            return s * Q_in - u_in * 1000 * (h_out - h_in)

        # Use fsolve to find the root
        T_out = fsolve(equation_to_solve, T_in)  # Initial guess T_in

        # ASSUMPTION: The chamber air temperature is the average of the inlet and outlet air temperatures
        T_chamb = 0.5 * (T_in + T_out)

        W_chamb = calculateW(T_chamb, RH_chamb, P_atm)
        W_leaf = calculateW(T_leaf, 1, P_atm)
        E = (u_in * (W_chamb - W_in)) / (s * (1 - W_chamb))  # mmol/m^2/s
        gtw = E / (W_leaf - W_chamb) / 1000
        gsw_bottom = 1 / (1 / gtw - 1 / gbw)

        gsw_total[i] = gsw_bottom * sidedness[i]
        T_chambs[i] = T_chamb
        W_chambs[i] = W_chamb
        T_outs[i] = T_out
        




    # Add new columns to data
    data['gsw_corrected'] = gsw_total
    data['Ta_chamb_corrected'] = T_chambs
    data['Wa_chamb_corrected'] = W_chambs
    data['T_out_corrected'] = T_outs
    data['stomatal_sidedness'] = sidedness

    # Write out the corrected data
    corrected_filepath = filepath.replace('.csv', '_corrected.csv')
    data.to_csv(corrected_filepath, index=False)

    return corrected_filepath

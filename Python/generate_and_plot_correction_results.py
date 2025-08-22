import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from add_gsw_correction_to_LI600 import add_gsw_correction_to_LI600
import os

def linear_fit(x, a, c):
    """Linear function for fitting: y = a*x + c"""
    return a * x + c

def plot_correction_results(input_file="walnut.csv"):
    """
    Generate correction and plot results similar to MATLAB script
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file
    """
    
    # Generate correction
    print(f"Processing file: {input_file}")
    data = add_gsw_correction_to_LI600(input_file)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: gsw correction
    x1 = data['gsw'].values
    y1 = data['gsw_corrected'].values
    
    # Remove any NaN or infinite values
    mask1 = np.isfinite(x1) & np.isfinite(y1)
    x1_clean = x1[mask1]
    y1_clean = y1[mask1]
    
    # Plot original vs original (1:1 line)
    ax1.scatter(x1_clean, x1_clean, c='black', alpha=0.5, label='Original $g_{sw}$')
    # Plot original vs corrected
    ax1.scatter(x1_clean, y1_clean, c='black', s=50, label='Corrected $g_{sw}$')
    
    # Fit linear model
    if len(x1_clean) > 1:
        popt1, _ = curve_fit(linear_fit, x1_clean, y1_clean)
        a1, c1 = popt1
        x_fit1 = np.linspace(x1_clean.min(), x1_clean.max(), 100)
        y_fit1 = linear_fit(x_fit1, a1, c1)
        ax1.plot(x_fit1, y_fit1, 'r-', alpha=0.7, 
                label=f'y = {a1:.2f}x + {c1:.4f}')
    
    ax1.set_xlabel('Original $g_{sw}$ (mol m$^{-2}$ s$^{-1}$)', fontsize=12)
    ax1.set_ylabel('$g_{sw}$ (mol m$^{-2}$ s$^{-1}$)', fontsize=12)
    ax1.set_title('LI-600 Stomatal Correction', fontsize=14)
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='box')
    
    # Plot 2: W_chamb correction
    # Calculate original W_chamb
    a = 0.61365  # unitless (empirical magnitude of es vs T)
    b = 17.502   # unitless (empirical slope of es vs T)
    c = 240.97   # C (empirical offset of es vs T)
    
    def es(T):
        """Saturation vapor pressure vs T function (kPa)"""
        return a * np.exp(b * T / (T + c))
    
    def W(T, RH, P_atm):
        """Water vapor mole fraction (mol/mol)"""
        return es(T) * RH / P_atm
    
    x2 = W(data['Tref'].values, data['rh_r'].values/100, data['P_atm'].values)
    y2 = data['W_chamb_corrected'].values
    
    # Remove any NaN or infinite values
    mask2 = np.isfinite(x2) & np.isfinite(y2)
    x2_clean = x2[mask2]
    y2_clean = y2[mask2]
    
    # Plot original vs original (1:1 line)
    ax2.scatter(x2_clean, x2_clean, c='black', alpha=0.5, label='Original $W_{chamb}$')
    # Plot original vs corrected
    ax2.scatter(x2_clean, y2_clean, c='black', s=50, label='Corrected $W_{chamb}$')
    
    # Fit linear model
    if len(x2_clean) > 1:
        popt2, _ = curve_fit(linear_fit, x2_clean, y2_clean)
        a2, c2 = popt2
        x_fit2 = np.linspace(x2_clean.min(), x2_clean.max(), 100)
        y_fit2 = linear_fit(x_fit2, a2, c2)
        ax2.plot(x_fit2, y_fit2, 'r-', alpha=0.7, 
                label=f'y = {a2:.2f}x + {c2:.4f}')
    
    ax2.set_xlabel('Original $W_{chamb}$ (mol mol$^{-1}$)', fontsize=12)
    ax2.set_ylabel('$W_{chamb}$ (mol mol$^{-1}$)', fontsize=12)
    ax2.set_title('LI-600 Chamber Water Correction', fontsize=14)
    ax2.legend(loc='lower right')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal', adjustable='box')
    
    # Adjust layout and save
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.splitext(input_file)[0] + '_correction_plots.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    
    # Show plot
    plt.show()
    
    return data

if __name__ == "__main__":
    # Run the correction and plotting
    # Use the walnut.csv file from the MATLAB folder
    data = plot_correction_results("walnut.csv")
"""
LI-600 Stomatal Correction Tool - Graphical User Interface
A Streamlit-based GUI for applying Rizzo & Bailey (2025) corrections to LI-600 data
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
import io
import os
from datetime import datetime

# Page configuration
st.set_page_config(
    page_title="LI-600 Stomatal Correction",
    layout="wide"
)

# Add custom CSS for better styling
st.markdown("""
    <style>
    .stApp {
        max-width: 1200px;
        margin: 0 auto;
    }
    </style>
    """, unsafe_allow_html=True)

@st.cache_data
def add_gsw_correction_to_LI600(file_buffer, stomatal_sidedness=1, thermal_conductance=0.007):
    """
    Applies the Rizzo & Bailey (2025) correction of chamber air temperature 
    and stomatal conductance to a csv file exported from an LI-600
    """
    
    # Read the CSV file
    try:
        # Try reading with LI-600 format (first row is groups, second row is names, third row is units)
        data = pd.read_csv(file_buffer, skiprows=[0, 2], header=0)
        # Reset buffer for potential re-reading
        file_buffer.seek(0)
    except:
        file_buffer.seek(0)
        try:
            data = pd.read_csv(file_buffer, skiprows=2, header=0)
            file_buffer.seek(0)
        except:
            file_buffer.seek(0)
            data = pd.read_csv(file_buffer)
            file_buffer.seek(0)
    
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
    
    # Progress bar
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Process each data point
    for i in range(n):
        # Update progress
        progress = (i + 1) / n
        progress_bar.progress(progress)
        status_text.text(f'Processing row {i+1} of {n}...')
        
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
    
    # Clear progress indicators
    progress_bar.empty()
    status_text.empty()
    
    # Add corrected values to dataframe
    data['gsw_corrected'] = gsw_total
    data['T_in_corrected'] = T_ins
    data['Ta_chamb_corrected'] = T_chambs
    data['T_out_corrected'] = T_outs
    data['W_chamb_corrected'] = W_chambs
    data['stomatal_sidedness'] = sidedness
    
    return data

def create_plots(data):
    """Create comparison plots for the corrected data"""
    
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
    ax1.scatter(x1_clean, x1_clean, c='black', alpha=0.3, label='Original $g_{sw}$', s=30)
    # Plot original vs corrected
    ax1.scatter(x1_clean, y1_clean, c='black', s=30, label='Corrected $g_{sw}$')
    
    # Fit linear model
    if len(x1_clean) > 1:
        def linear_fit(x, a, c):
            return a * x + c
        
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
    ax2.scatter(x2_clean, x2_clean, c='black', alpha=0.3, label='Original $W_{chamb}$', s=30)
    # Plot original vs corrected
    ax2.scatter(x2_clean, y2_clean, c='black', s=30, label='Corrected $W_{chamb}$')
    
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
    
    # Adjust layout
    plt.tight_layout()
    
    return fig

# Main app
def main():
    # Title and description
    st.title("LI-600 Stomatal Correction Tool")
    st.markdown("""
    This tool applies the **Rizzo & Bailey (2025)** psychrometric temperature correction 
    to address the systematic positive bias in stomatal conductance measurements from 
    LI-600 porometer CSV exports. The correction is particularly important for measurements 
    with high stomatal conductance (>0.25 mol m⁻² s⁻¹) or low relative humidity.
    """)
    
    # Sidebar for settings
    with st.sidebar:
        st.header(":material/settings: Settings")
        
        stomatal_sidedness = st.slider(
            "Stomatal Sidedness",
            min_value=1.0,
            max_value=2.0,
            value=1.0,
            step=0.1,
            help="1.0 = hypostomatous (stomata on one side), 2.0 = amphistomatous (stomata on both sides)"
        )
        
        st.markdown("---")
        
        # Initialize thermal_conductance with default
        thermal_conductance = 0.007
        
        with st.expander("Advanced Settings"):
            thermal_conductance = st.number_input(
                "Thermal Conductance C (W/°C)",
                min_value=0.001,
                max_value=0.1,
                value=0.007,
                step=0.001,
                format="%.4f",
                help="Thermal conductance between inlet air and chamber air. Default value of 0.007 was empirically calibrated by Rizzo & Bailey (2025)."
            )
        
        st.markdown("---")
        st.markdown("""
        ### About
        
        This tool implements a psychrometric temperature correction for the 
        persistent positive bias observed in stomatal conductance measured by 
        the LI-600 open flow-through porometer.
        
        The correction addresses the instrument's assumption of constant air 
        temperature throughout the flow stream by:
        - Applying psychrometric principles to estimate temperature changes
        - Solving an augmented system of equations with energy balance constraints
        - Accounting for heat addition from instrument electronics and fans
        
        The correction reduces measurement bias from ~34% to ~11% and is most 
        critical at high stomatal conductance (>0.25 mol m⁻² s⁻¹) and low 
        relative humidity conditions.
        
        **Reference:** 
        Rizzo, K.T. & Bailey, B.N. (2025). A psychrometric temperature 
        correction for the persistent positive bias observed in stomatal 
        conductance measured by the open flow-through LI-600 porometer.
        """)
    
    # Main content area - single column layout
    st.header(":material/upload_file: Upload Data")
    
    uploaded_file = st.file_uploader(
        "Choose an LI-600 CSV file",
        type=['csv'],
        help="Drag and drop or click to browse for your LI-600 exported CSV file"
    )
    
    if uploaded_file is not None:
        st.success(f":material/check_circle: File uploaded: {uploaded_file.name}")
        
        # Show preview of original data
        st.subheader("Data Preview (first 5 rows)")
        try:
            preview_df = pd.read_csv(uploaded_file, skiprows=[0, 2], header=0, nrows=5)
            uploaded_file.seek(0)  # Reset file pointer
            st.dataframe(preview_df[['gsw', 'Tref', 'Tleaf', 'rh_r', 'rh_s', 'flow', 'P_atm']].round(3))
        except:
            st.warning("Preview not available - file will still be processed")
            uploaded_file.seek(0)
        
        # Process section
        st.markdown("---")
        st.header(":material/autorenew: Process Data")
        
        if st.button(":material/play_arrow: Run Correction", type="primary", use_container_width=False):
            
            with st.spinner("Processing your data..."):
                # Run correction
                corrected_data = add_gsw_correction_to_LI600(uploaded_file, stomatal_sidedness, thermal_conductance)
                
                # Create plots
                fig = create_plots(corrected_data)
            
            st.success(":material/check_circle: Correction completed successfully!")
            
            # Store in session state for download
            st.session_state['corrected_data'] = corrected_data
            st.session_state['figure'] = fig
            st.session_state['filename'] = uploaded_file.name
    
    # Results section
    if 'corrected_data' in st.session_state:
        st.markdown("---")
        st.header(":material/insights: Results")
        
        # Display plots
        st.pyplot(st.session_state['figure'])
        
        # Summary statistics
        col1, col2, col3 = st.columns(3)
        
        data = st.session_state['corrected_data']
        
        with col1:
            st.metric(
                "Mean gsw correction",
                f"{(data['gsw_corrected'].mean() - data['gsw'].mean()):.4f}",
                f"{((data['gsw_corrected'].mean() / data['gsw'].mean() - 1) * 100):.1f}%"
            )
        
        with col2:
            st.metric(
                "Mean T_chamber change",
                f"{(data['Ta_chamb_corrected'].mean() - data['Tref'].mean()):.2f}°C"
            )
        
        with col3:
            st.metric(
                "Data points processed",
                len(data)
            )
        
        # Download section
        st.markdown("---")
        st.header(":material/download: Download Results")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Prepare CSV for download
            csv = data.to_csv(index=False)
            original_name = st.session_state['filename']
            corrected_name = original_name.replace('.csv', '_corrected.csv')
            
            st.download_button(
                label=":material/download: Download Corrected CSV",
                data=csv,
                file_name=corrected_name,
                mime='text/csv',
                use_container_width=True
            )
        
        with col2:
            # Save plot for download
            img_buffer = io.BytesIO()
            st.session_state['figure'].savefig(img_buffer, format='png', dpi=150, bbox_inches='tight')
            img_buffer.seek(0)
            
            plot_name = original_name.replace('.csv', '_plots.png')
            
            st.download_button(
                label=":material/image: Download Plots (PNG)",
                data=img_buffer,
                file_name=plot_name,
                mime='image/png',
                use_container_width=True
            )
        
        # Show corrected data preview
        with st.expander("View Corrected Data (first 10 rows)"):
            display_cols = ['gsw', 'gsw_corrected', 'Ta_chamb_corrected', 'T_out_corrected', 'W_chamb_corrected']
            st.dataframe(data[display_cols].head(10).round(4))

if __name__ == "__main__":
    main()
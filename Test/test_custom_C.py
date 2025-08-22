#!/usr/bin/env python3
"""Test all implementations with custom thermal conductance value"""

import sys
import os
import subprocess
import pandas as pd
import numpy as np

# Add Python implementation to path
sys.path.append('../Python')
from add_gsw_correction_to_LI600 import add_gsw_correction_to_LI600

def test_custom_C():
    """Test all implementations with custom thermal conductance C=0.01"""
    
    input_file = 'walnut.csv'
    custom_C = 0.01  # Different from default 0.007
    
    print(f"Testing all implementations with custom C = {custom_C}")
    print("=" * 60)
    
    # Test Python implementation
    print("\n1. Testing Python implementation...")
    try:
        data_python = add_gsw_correction_to_LI600(input_file, 
                                                  stomatal_sidedness=1.0, 
                                                  thermal_conductance=custom_C)
        gsw_python = data_python['gsw_corrected'].values
        print(f"   Python: First 5 gsw_corrected values: {gsw_python[:5]}")
        print(f"   Python: Mean gsw_corrected = {np.mean(gsw_python):.6f}")
    except Exception as e:
        print(f"   Python ERROR: {e}")
        gsw_python = None
    
    # Test MATLAB implementation
    print("\n2. Testing MATLAB implementation...")
    matlab_cmd = f"""
    cd ../MATLAB;
    data = add_gsw_correction_to_LI600('walnut.csv', 1.0, {custom_C});
    gsw = data.gsw_corrected;
    fprintf('MATLAB: First 5 values: [%f, %f, %f, %f, %f]\\n', gsw(1:5));
    fprintf('MATLAB: Mean = %f\\n', mean(gsw));
    exit;
    """
    
    try:
        result = subprocess.run(['matlab', '-batch', matlab_cmd], 
                              capture_output=True, text=True, timeout=30)
        print(f"   {result.stdout}")
    except Exception as e:
        print(f"   MATLAB ERROR: {e}")
    
    # Test R implementation
    print("\n3. Testing R implementation...")
    r_script = f"""
    source('../R/add_gsw_correction_to_LI600.R')
    data <- add_gsw_correction_to_LI600('../R/walnut.csv', 
                                        stomatal_sidedness = 1.0, 
                                        thermal_conductance = {custom_C})
    cat('R: First 5 values:', head(data$gsw_corrected, 5), '\\n')
    cat('R: Mean =', mean(data$gsw_corrected), '\\n')
    """
    
    try:
        result = subprocess.run(['Rscript', '-e', r_script], 
                              capture_output=True, text=True, timeout=30)
        print(f"   {result.stdout}")
    except Exception as e:
        print(f"   R ERROR: {e}")
    
    # Compare with default C value
    print("\n4. Comparing with default C = 0.007...")
    try:
        data_default = add_gsw_correction_to_LI600(input_file, 
                                                   stomatal_sidedness=1.0, 
                                                   thermal_conductance=0.007)
        gsw_default = data_default['gsw_corrected'].values
        
        if gsw_python is not None:
            diff = np.mean(np.abs(gsw_python - gsw_default))
            print(f"   Mean absolute difference: {diff:.6f}")
            print(f"   Max absolute difference: {np.max(np.abs(gsw_python - gsw_default)):.6f}")
            
            # Check that results are different (confirming C parameter is being used)
            if diff > 1e-6:
                print("   ✓ Confirmed: Custom C value is being used (results differ from default)")
            else:
                print("   ⚠ Warning: Results same as default - C parameter may not be working")
    except Exception as e:
        print(f"   Comparison ERROR: {e}")
    
    print("\n" + "=" * 60)
    print("Testing complete!")

if __name__ == "__main__":
    test_custom_C()
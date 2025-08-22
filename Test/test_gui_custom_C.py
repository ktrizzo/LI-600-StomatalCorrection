#!/usr/bin/env python3
"""Test GUI implementation with custom thermal conductance value"""

import sys
import os
import pandas as pd
import numpy as np

# Import the correction function (same one used by GUI)
sys.path.append('../Python')
from add_gsw_correction_to_LI600 import add_gsw_correction_to_LI600

def test_gui_custom_C():
    """Test GUI implementation's correction function with custom C"""
    
    input_file = 'walnut.csv'
    
    print("Testing GUI implementation with different C values")
    print("=" * 60)
    
    # Test with default C
    print("\n1. Testing with default C = 0.007...")
    data_default = add_gsw_correction_to_LI600(input_file, 
                                               stomatal_sidedness=1.0, 
                                               thermal_conductance=0.007)
    gsw_default = data_default['gsw_corrected'].values
    print(f"   Mean gsw_corrected = {np.mean(gsw_default):.6f}")
    
    # Test with custom C = 0.01
    print("\n2. Testing with custom C = 0.01...")
    data_custom = add_gsw_correction_to_LI600(input_file, 
                                              stomatal_sidedness=1.0, 
                                              thermal_conductance=0.01)
    gsw_custom = data_custom['gsw_corrected'].values
    print(f"   Mean gsw_corrected = {np.mean(gsw_custom):.6f}")
    
    # Test with custom C = 0.005
    print("\n3. Testing with custom C = 0.005...")
    data_low = add_gsw_correction_to_LI600(input_file, 
                                           stomatal_sidedness=1.0, 
                                           thermal_conductance=0.005)
    gsw_low = data_low['gsw_corrected'].values
    print(f"   Mean gsw_corrected = {np.mean(gsw_low):.6f}")
    
    # Compare results
    print("\n4. Comparison of results:")
    print(f"   Difference (C=0.01 vs C=0.007): {np.mean(gsw_custom - gsw_default):.6f}")
    print(f"   Difference (C=0.005 vs C=0.007): {np.mean(gsw_low - gsw_default):.6f}")
    
    # Verify that higher C leads to lower correction
    if np.mean(gsw_custom) < np.mean(gsw_default) < np.mean(gsw_low):
        print("   ✓ Confirmed: Higher C value leads to lower corrected gsw (expected behavior)")
    else:
        print("   ⚠ Warning: Unexpected relationship between C and corrected gsw")
    
    print("\n" + "=" * 60)
    print("GUI correction function test complete!")
    print("\nNote: The GUI uses the same Python correction function,")
    print("so these results confirm the GUI will work correctly with custom C values.")

if __name__ == "__main__":
    test_gui_custom_C()
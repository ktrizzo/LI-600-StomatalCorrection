import pandas as pd
import numpy as np
import sys
import os

def compare_csv_files(matlab_file, python_file):
    """
    Compare the corrected values from MATLAB and Python implementations
    """
    
    # Read both CSV files
    print(f"Reading MATLAB output: {matlab_file}")
    matlab_data = pd.read_csv(matlab_file)
    
    print(f"Reading Python output: {python_file}")
    python_data = pd.read_csv(python_file)
    
    # Key columns to compare
    columns_to_compare = [
        'gsw_corrected',
        'T_in_corrected', 
        'Ta_chamb_corrected',
        'T_out_corrected',
        'W_chamb_corrected'
    ]
    
    print("\n" + "="*80)
    print("COMPARISON OF CORRECTED VALUES")
    print("="*80)
    
    for col in columns_to_compare:
        if col in matlab_data.columns and col in python_data.columns:
            matlab_values = matlab_data[col].values
            python_values = python_data[col].values
            
            # Calculate differences
            differences = matlab_values - python_values
            abs_differences = np.abs(differences)
            
            # Remove NaN values for statistics
            valid_mask = ~np.isnan(differences)
            valid_diffs = differences[valid_mask]
            valid_abs_diffs = abs_differences[valid_mask]
            
            # Calculate statistics
            if len(valid_diffs) > 0:
                mean_diff = np.mean(valid_diffs)
                std_diff = np.std(valid_diffs)
                max_abs_diff = np.max(valid_abs_diffs)
                mean_abs_diff = np.mean(valid_abs_diffs)
                
                # Calculate relative error where possible
                matlab_nonzero = matlab_values[valid_mask]
                nonzero_mask = matlab_nonzero != 0
                if np.any(nonzero_mask):
                    relative_errors = valid_abs_diffs[nonzero_mask] / np.abs(matlab_nonzero[nonzero_mask]) * 100
                    mean_rel_error = np.mean(relative_errors)
                    max_rel_error = np.max(relative_errors)
                else:
                    mean_rel_error = 0
                    max_rel_error = 0
                
                print(f"\n{col}:")
                print(f"  Number of valid comparisons: {len(valid_diffs)}")
                print(f"  Mean difference (MATLAB - Python): {mean_diff:.6e}")
                print(f"  Std dev of differences: {std_diff:.6e}")
                print(f"  Mean absolute difference: {mean_abs_diff:.6e}")
                print(f"  Max absolute difference: {max_abs_diff:.6e}")
                print(f"  Mean relative error: {mean_rel_error:.3f}%")
                print(f"  Max relative error: {max_rel_error:.3f}%")
                
                # Show first few differences for inspection
                print(f"  First 5 values comparison:")
                for i in range(min(5, len(matlab_values))):
                    print(f"    Row {i+1}: MATLAB={matlab_values[i]:.6f}, Python={python_values[i]:.6f}, Diff={differences[i]:.6e}")
        else:
            print(f"\n{col}: Column not found in one or both files")
    
    # Check if the original columns are the same (they should be)
    print("\n" + "="*80)
    print("ORIGINAL DATA CONSISTENCY CHECK")
    print("="*80)
    
    original_cols = ['gsw', 'Tref', 'Tleaf', 'rh_r', 'rh_s', 'flow', 'P_atm']
    
    for col in original_cols:
        if col in matlab_data.columns and col in python_data.columns:
            matlab_vals = matlab_data[col].values
            python_vals = python_data[col].values
            
            if np.allclose(matlab_vals, python_vals, rtol=1e-10, equal_nan=True):
                print(f"{col}: ✓ Identical")
            else:
                max_diff = np.max(np.abs(matlab_vals - python_vals))
                print(f"{col}: ✗ Different (max diff: {max_diff:.6e})")
    
    return matlab_data, python_data

if __name__ == "__main__":
    # File paths
    matlab_file = "walnut_corrected_matlab.csv"
    python_file = "walnut_corrected_python.csv"
    
    # Check if files exist
    if not os.path.exists(matlab_file):
        print(f"Error: MATLAB file not found: {matlab_file}")
        sys.exit(1)
    
    if not os.path.exists(python_file):
        print(f"Error: Python file not found: {python_file}")
        sys.exit(1)
    
    # Run comparison
    matlab_data, python_data = compare_csv_files(matlab_file, python_file)
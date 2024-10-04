# Import necessary libraries and the function
from add_gsw_correction_to_LI600 import *

def main():
    # Specify the path to the CSV file
    csv_file_path = 'walnut.csv'
    
    # Run correction
    new_file_default = add_gsw_correction_to_LI600(csv_file_path)
    print(f"Corrected file created: {new_file_default}")


if __name__ == "__main__":
    main()

# LI-600 Stomatal Correction Tool - Graphical User Interface

A user-friendly web-based interface for applying Rizzo & Bailey (2025) corrections to LI-600 porometer data.

## Features

-  **Drag-and-drop file upload** - Simply drag your CSV file into the browser
-  **Real-time visualization** - See correction plots immediately after processing
- **Adjustable parameters** - Control stomatal sidedness with an intuitive slider
-  **Easy downloads** - Download corrected CSV and plots with one click

## Quick Start

### Prerequisites

- Python 3.7 or higher
- Git (to clone the repository)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/ktrizzo/LI-600-StomatalCorrection.git
cd LI-600-StomatalCorrection/Graphical\ User\ Interface
```

2. Create a virtual environment (recommended):
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

### Running the Application

1. Start the Streamlit app:
```bash
streamlit run app.py
```

2. Your browser will automatically open to `http://localhost:8501`

3. If the browser doesn't open automatically, manually navigate to: `http://localhost:8501`

## Using the Tool

1. **Upload your file**: Drag and drop or click to browse for your LI-600 CSV export

2. **Adjust settings**: Use the sidebar slider to set stomatal sidedness:
   - 1.0 = Hypostomatous (stomata on one side only)
   - 2.0 = Amphistomatous (stomata on both sides)
   - Values between 1.0-2.0 for partial coverage

3. **Process data**: Click the "🚀 Run Correction" button

4. **View results**: 
   - See comparison plots of original vs corrected values
   - Review summary statistics
   - Preview corrected data

5. **Download results**:
   - Download corrected CSV file
   - Download plots as PNG image


## Input File Format

The tool accepts CSV files exported directly from the LI-600 porometer. Required columns include:
- `gsw` - Stomatal conductance
- `Tref` - Reference temperature
- `Tleaf` - Leaf temperature
- `rh_r` - Reference humidity
- `rh_s` - Sample humidity
- `flow` - Air flow rate
- `P_atm` - Atmospheric pressure
- `E_apparent` - Apparent transpiration

## Output Files

The tool generates:
1. **Corrected CSV** - Original data plus correction columns:
   - `gsw_corrected` - Corrected stomatal conductance
   - `Ta_chamb_corrected` - Corrected chamber temperature
   - `T_out_corrected` - Corrected outlet temperature
   - `W_chamb_corrected` - Corrected chamber water mole fraction

2. **Plots (PNG)** - Visual comparison of corrections

## Troubleshooting

### Port already in use
If you see an error about port 8501 being in use:
```bash
streamlit run app.py --server.port 8502
```

### Package installation issues
If you encounter issues installing packages:
```bash
pip install --upgrade pip
pip install -r requirements.txt --no-cache-dir
```

### File upload issues
- Ensure your CSV file is properly formatted from LI-600 export
- Check that the file isn't corrupted or incomplete
- Try opening the file in Excel/Numbers to verify it's readable

## Citation

If you use this tool in your research, please cite:

Rizzo, K. T. & Bailey, B. N. (2025). *A psychrometric temperature correction for the persistent positive bias observed in stomatal conductance measured by the open flow-through LI-600 porometer.*.

## Support

For issues or questions, please open an issue on the [GitHub repository](https://github.com/ktrizzo/LI-600-StomatalCorrection/issues).

## License

This project is licensed under the MIT License - see the LICENSE file for details.
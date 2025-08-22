# LI-600 Stomatal Conductance Correction

## Overview

This repository provides implementations of the **Rizzo & Bailey (2025)** psychrometric temperature correction for addressing the systematic positive bias observed in stomatal conductance measurements from the LI-600 open flow-through porometer.

## Background

The LI-600 porometer, while offering high-throughput measurement capabilities, has been shown to produce systematically biased stomatal conductance values compared to reference infrared gas analyzer (IRGA) systems like the LI-6800. This bias:

- Increases exponentially with stomatal conductance
- Can exceed 100% at high conductance values (>0.8 mol m⁻² s⁻¹)
- Is most pronounced under low relative humidity conditions
- Results from the instrument's assumption of constant air temperature throughout the flow stream

### The Correction

The psychrometric correction addresses this bias by:
- Relaxing the isothermal assumption in the flow path
- Applying energy balance constraints to estimate temperature variation
- Solving an augmented system of equations that accounts for:
  - Psychrometric cooling from water evaporation
  - Heat addition from instrument electronics and fans
  - Temperature gradients between inlet and outlet air streams

When applied, this correction:
- Reduces mean bias
- Improves correlation between porometer and IRGA measurements
- Is particularly critical for measurements with:
  - High stomatal conductance (>0.25 mol m⁻² s⁻¹)
  - Low relative humidity (<40%)

## Available Implementations

This repository provides **four different ways** to apply the correction to your LI-600 data:

### 1. MATLAB Implementation (`/MATLAB`)
- Original implementation
- Function: `add_gsw_correction_to_LI600.m`
- Includes plotting capabilities
- Requires MATLAB with Symbolic Math Toolbox

### 2. Python Implementation (`/Python`)
- Full-featured Python version
- Function: `add_gsw_correction_to_LI600.py`
- Includes visualization scripts
- Requirements: numpy, pandas, scipy, matplotlib

### 3. R Implementation (`/R`)
- Statistical computing implementation
- Function: `add_gsw_correction_to_LI600.R`
- Includes ggplot2-based visualization
- Requirements: nleqslv, ggplot2, gridExtra

### 4. Graphical User Interface (`/Graphical User Interface`)
- **User-friendly web-based interface**
- No coding required
- Drag-and-drop file upload
- Real-time visualization
- Built with Streamlit (Python-based)
- **Recommended for users unfamiliar with programming**

## Quick Start

### For Non-Programmers (GUI)
```bash
# Clone the repository
git clone https://github.com/ktrizzo/LI-600-StomatalCorrection.git
cd "LI-600-StomatalCorrection/Graphical User Interface"

# Install requirements
pip install -r requirements.txt

# Run the GUI
streamlit run app.py
```
Then open your browser to http://localhost:8501 and drag-drop your LI-600 CSV file.

### For MATLAB Users
```matlab
% Add to MATLAB path
addpath('LI-600-StomatalCorrection/MATLAB')

% Apply correction
data = add_gsw_correction_to_LI600('your_file.csv', stomatal_sidedness);
```

### For Python Users
```python
from add_gsw_correction_to_LI600 import add_gsw_correction_to_LI600

# Apply correction
data = add_gsw_correction_to_LI600('your_file.csv', stomatal_sidedness=1.0)
```

### For R Users
```r
source("add_gsw_correction_to_LI600.R")

# Apply correction
data <- add_gsw_correction_to_LI600("your_file.csv", stomatal_sidedness = 1.0)
```

## Input Requirements

All implementations accept CSV files exported directly from the LI-600 with required columns:
- `gsw` - Stomatal conductance
- `Tref` - Reference temperature  
- `Tleaf` - Leaf temperature
- `rh_r` - Reference relative humidity
- `rh_s` - Sample relative humidity
- `flow` - Air flow rate
- `P_atm` - Atmospheric pressure
- `E_apparent` - Apparent transpiration rate

## Output

Each implementation generates:
- **Corrected CSV file** with additional columns:
  - `gsw_corrected` - Corrected stomatal conductance
  - `Ta_chamb_corrected` - Corrected chamber air temperature
  - `T_out_corrected` - Corrected outlet temperature
  - `W_chamb_corrected` - Corrected chamber water mole fraction
- **Visualization plots** comparing original vs corrected values

## Validation

The correction has been validated on:
- 610 paired measurements across 25 angiosperm species
- Both field and controlled conditions
- Wide range of stomatal conductance values (0-2.0 mol m⁻² s⁻¹)
- Various relative humidity conditions (30-70%)

## Citation

If you use this correction in your research, please cite:

**Rizzo, K.T. & Bailey, B.N. (2025).** A psychrometric temperature correction for the persistent positive bias observed in stomatal conductance measured by the open flow-through LI-600 porometer. *[Journal details to be added upon publication]*

## Support

For issues, questions, or contributions, please:
- Open an issue on the [GitHub repository](https://github.com/ktrizzo/LI-600-StomatalCorrection/issues)
- Contact the authors (see paper for contact information)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This work was supported by [funding information from paper].

---

**Note:** The thermal conductance parameter (C = 0.007 W/°C) used in the correction was empirically calibrated using paired LI-600/LI-6800 measurements. While this value performed well across diverse species and conditions, users should be aware that instrument-specific variations may exist. You may recalibrate this coefficient yourself with controlled LI-600 and LI-6800 measurements on the same leaf patches.
=======
# LI-600-StomatalCorrection
Routines for correcting stomatal conductance reported by the LI-600 Porometer according to the method of Rizzo and Bailey (2025)
>>>>>>> a6144b4e0418bbf673535d5b4875ba9d471b1489

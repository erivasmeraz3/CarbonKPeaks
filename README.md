# CarbonKPeaks - C1s XAS Peak Fitting Software

**Version 1.0**

CarbonKPeaks is a graphical user interface (GUI) application for fitting and analyzing Carbon K-edge (C1s) X-ray Absorption Spectroscopy (XAS) data. It provides automated peak fitting with a double arctangent baseline and Gaussian peaks, along with comprehensive visualization and export capabilities.

## Features

- **Automated Peak Fitting**: Fits C1s XAS spectra using 8 Gaussian peaks with a double arctangent baseline
- **Interactive Visualization**: Real-time display of spectra, fitted peaks, and residuals
- **Multi-Sample Comparison**: Compare and analyze multiple samples simultaneously
- **Publication-Quality Exports**: Generate 300 DPI plots suitable for publications
- **Comprehensive Data Export**: Export fit parameters, raw data, and metadata in multiple formats

## Installation

### Option 1: Run from Executable (Windows)
1. Download `CarbonKPeaks.exe` from the `dist` folder
2. Double-click to run - no installation required

### Option 2: Run from Source
```bash
# Required Python packages
pip install numpy pandas matplotlib lmfit tkinter

# Run the application
python c1s_peak_viewer_gui_final.py
```

## Input Data Format

The software accepts CSV files with two columns:
- **Column 1**: Energy (eV) - typically 280-320 eV range
- **Column 2**: Normalized intensity

Example:
```csv
energy,y
279.81,0.00331595
280.21,-0.03049639
280.61,0.02211656
...
```

## Peak Assignments

The software fits 8 peaks corresponding to different carbon functional groups:

| Peak | Energy (eV) | Assignment |
|------|-------------|------------|
| 1 | 284.0 | Quinone C=O |
| 2 | 285.0 | Aromatic C=C (1s -> pi*) |
| 3 | 286.2 | Phenolic/Ketone C-OH/C=O |
| 4 | 287.4 | Aliphatic C-H |
| 5 | 288.6 | Carboxyl COOH |
| 6 | 289.3 | O-alkyl C-O |
| 7 | 290.3 | Carbonate CO3 |
| 8 | 292.0 | Sigma* C-C |

## User Interface

### Tab 1: Peak Fitting
- **Sample List**: Browse and select samples from a folder
- **Peak Parameters**: View fitted peak centers, FWHM, and area percentages
- **Fit/Quant Checkboxes**:
  - **F (Fit)**: Include peak in fitting
  - **Q (Quant)**: Include peak in quantification
- **Baseline Controls**: Manual adjustment of baseline parameters (when enabled)
- **Plot Area**: Shows spectrum with fit, individual peaks, and residuals

### Tab 2: Analysis & Comparison
- **Sample Selection**: Select multiple samples for comparison
- **Generate Plots**:
  - **Statistics**: Summary grid with fit quality for each sample
  - **Comparison Charts**: Stacked bar charts and composition tables
  - **Spectral Overlay**: Overlay multiple spectra with adjustable y-offset
- **Export Options**: Export figures, data, or complete analysis packages

### Tab 3: Peak Configuration
- **Custom Peak Centers**: Enable/disable custom peak positions
- **Peak Ranges**: Adjust allowed fitting ranges for each peak

## How to Use

### Basic Workflow

1. **Load Spectra**
   - Click "Browse" to select a folder containing CSV files
   - Or use command line: `python c1s_peak_viewer_gui_final.py --spectra_dir /path/to/files`

2. **Fit Spectra**
   - Select a sample from the list
   - Click "Fit" to perform automatic fitting
   - View results in the plot area

3. **Adjust Fitting Options**
   - Use F checkboxes to include/exclude peaks from fitting
   - Use Q checkboxes to include/exclude peaks from quantification
   - Enable "Manual Baseline" for fine-tuning

4. **Compare Multiple Samples**
   - Go to "Analysis & Comparison" tab
   - Select samples and click "Add Selected"
   - Generate different plot types for comparison

5. **Export Results**
   - Use export buttons in the main tab for single samples
   - Use "Export Complete Analysis" for comprehensive multi-sample exports

## Export Options

### Single Sample Export
- **Export Plot**: Save current plot as PNG/PDF
- **Export Fit Data**: Save peak parameters as CSV
- **Export All (Complete)**: Full data package including:
  - Fit plots (PDF & PNG)
  - Peak parameters CSV
  - Spectrum data CSV (raw, fit, baseline, residual, individual peaks)
  - Metadata JSON
  - Fit report TXT

### Multi-Sample Analysis Export
- **Export Figure**: Save current analysis plot
- **Export Data (CSV)**: Save comparison data
- **Export Complete Analysis**: Full analysis package including:
  - All 3 plot types (Statistics, Comparison, Spectral Overlay)
  - Individual sample data folders
  - Combined summary CSV
  - Analysis report

## Output Files Description

### Peak Parameters CSV
| Column | Description |
|--------|-------------|
| Peak | Display name |
| Center_eV | Peak center position |
| Height | Peak height |
| FWHM_eV | Full width at half maximum |
| Area | Integrated peak area |
| Percentage | Percent of total quantified area |

### Metadata JSON
Contains:
- Fit quality metrics (R-squared, chi-squared, AIC, BIC)
- Baseline parameters
- FWHM parameters
- Peak summary
- Settings used

### Spectrum Data CSV
| Column | Description |
|--------|-------------|
| Energy_eV | Energy values |
| Raw_Intensity | Original data |
| Total_Fit | Total fitted curve |
| Baseline | Double arctangent baseline |
| Residual | Fit residual |
| Peak_[name] | Individual peak contributions |

## Tips for Best Results

1. **Data Quality**: Ensure spectra are properly normalized before fitting
2. **Energy Range**: Data should cover approximately 280-320 eV
3. **Peak Selection**: Disable peaks that are not expected in your sample
4. **Baseline**: Use manual baseline adjustment for challenging spectra
5. **Comparison**: Use spectral overlay to visually compare samples before quantitative analysis

## Troubleshooting

**Problem**: Fit does not converge
- Try adjusting peak center ranges in the Peak Configuration tab
- Disable peaks that are not present in your sample

**Problem**: Poor baseline fit
- Enable manual baseline and adjust step heights
- Check that your data is properly normalized

**Problem**: Peaks appear at wrong positions
- Enable "Custom Peak Centers" and adjust ranges
- Verify your data energy calibration

## License

MIT License

## Contact

For questions, bug reports, or feature requests, contact the development team.

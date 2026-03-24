# CarbonKPeaks Tutorial
## A Step-by-Step Guide to C1s XAS Peak Fitting

---

## Introduction

This tutorial demonstrates how to use CarbonKPeaks to analyze Carbon K-edge XAS data. We will use cluster analysis data from sample P7S_199 as our example dataset.

**Sample Files Used:**
- C_P7S_199_cluster0.csv
- C_P7S_199_cluster1.csv
- C_P7S_199_cluster2.csv
- C_P7S_199_cluster3.csv
- C_P7S_199_cluster4.csv

These represent 5 different spectral clusters from a single sample, allowing us to compare carbon speciation across different regions.

---

## Step 1: Loading Your Data

### 1.1 Launch the Application
- Double-click `CarbonKPeaks.exe` or run from Python
- The main window will appear with three tabs

### 1.2 Browse to Your Data Folder
1. Click the **"Browse"** button in the left panel
2. Navigate to the folder containing your normalized CSV files
3. Select the folder and click "Select Folder"

### 1.3 Verify Files Loaded
- Your samples will appear in the "Samples" list
- The first sample will be automatically selected and displayed

**What You Should See:**
```
Samples:
  C_P7S_199_cluster0
  C_P7S_199_cluster1
  C_P7S_199_cluster2
  C_P7S_199_cluster3
  C_P7S_199_cluster4
```

---

## Step 2: Fitting Individual Spectra

### 2.1 Select a Sample
- Click on `C_P7S_199_cluster0` in the sample list
- The raw spectrum will display in the plot area

### 2.2 Perform the Fit
1. Click the **"Fit"** button
2. The software will automatically:
   - Fit the double arctangent baseline
   - Fit all 8 Gaussian peaks
   - Calculate peak areas and percentages

### 2.3 Understanding the Display

**Top Plot - Spectrum and Fit:**
- Black line: Raw data
- Red dashed line: Total fit
- Blue dashed line: Baseline
- Colored peaks: Individual Gaussian components

**Bottom Plot - Residual:**
- Shows the difference between data and fit
- A good fit has residuals randomly distributed around zero

### 2.4 Reviewing Peak Parameters

The peak parameters table shows:

| Column | Meaning |
|--------|---------|
| F | Checkbox - Include in Fit |
| Q | Checkbox - Include in Quantification |
| Peak | Carbon functional group |
| Center | Peak position (eV) |
| FWHM | Peak width |
| %Area | Percentage of total area |

**Example Results for Cluster 0:**
```
Peak          Center    FWHM    %Area
Quinone       284.0     0.82     5.2%
Aromatic      285.0     0.82    32.5%
Phenolic      286.2     0.82    15.2%
Aliphatic     287.4     0.82    12.8%
Carboxyl      288.6     0.82    11.3%
O-alkyl       289.3     0.82    14.2%
Carbonate     290.3     1.20     3.5%
Sigma*        292.0     2.50     5.3%
```

---

## Step 3: Customizing the Fit

### 3.1 Selecting Peaks for Fitting

If certain peaks are not present in your sample:
1. Uncheck the **"F"** box for that peak
2. Click **"Fit"** again
3. The fit will exclude that peak

**Common Scenarios:**
- Organic matter: May not have carbonate peak
- Carbonates: Strong 290.3 eV peak
- Aromatic-rich: Dominant 285.0 eV peak

### 3.2 Selecting Peaks for Quantification

To focus quantification on specific peaks:
1. Uncheck the **"Q"** box for peaks to exclude
2. Percentages will recalculate automatically

**Example:** To quantify only the main carbon species:
- Keep Q checked for: Aromatic, Aliphatic, Carboxyl, O-alkyl
- Uncheck Q for: Quinone, Phenolic, Carbonate, Sigma*

### 3.3 Using Toggle All Buttons
- **"All Fit"**: Select/deselect all peaks for fitting
- **"All Quant"**: Select/deselect all peaks for quantification

---

## Step 4: Manual Baseline Adjustment

For challenging spectra, manual baseline control can improve fits.

### 4.1 Enable Manual Baseline
1. Check the **"Manual Baseline"** checkbox
2. Baseline control sliders become active

### 4.2 Adjust Baseline Parameters

**Step 1 (Pre-edge):**
- **Center**: Position of first step (~284 eV)
- **Height**: Height of first step
- **Width**: Steepness of transition

**Step 2 (Post-edge):**
- **Center**: Position of second step (~292 eV)
- **Height**: Height of second step
- **Width**: Steepness of transition

### 4.3 Update the Fit
- Click **"Fit"** after adjusting baseline
- The fit will use your baseline as a starting point

---

## Step 5: Multi-Sample Comparison

### 5.1 Navigate to Analysis Tab
- Click on the **"Analysis & Comparison"** tab

### 5.2 Select Samples for Comparison
1. In the sample list, Ctrl+Click to select multiple samples:
   - C_P7S_199_cluster0
   - C_P7S_199_cluster1
   - C_P7S_199_cluster2
   - C_P7S_199_cluster3
   - C_P7S_199_cluster4
2. Click **"Add Selected"**
3. Selected samples appear in "Selected for Analysis"

### 5.3 Generate Statistics Plot
1. Click **"Statistics"**
2. View a grid showing fit quality for each sample:
   - Sample name
   - R-squared value
   - Reduced chi-squared
   - Visual quality indicator

### 5.4 Generate Comparison Charts
1. Click **"Comparison Charts"**
2. Four plots appear:
   - **Stacked Bar Chart**: Composition comparison
   - **Grouped Bar Chart**: Side-by-side peak comparison
   - **Line Plot**: Trends across samples
   - **Data Table**: Numerical values

### 5.5 Generate Spectral Overlay
1. Click **"Spectral Overlay"**
2. All selected spectra display with their fits
3. Adjust **"Y-Axis Offset"** slider to separate spectra
4. Check/uncheck **"Show Individual Peaks"** as needed

### 5.6 Zoom and Navigate
- Use **mouse scroll wheel** to zoom in/out
- Use the **toolbar** at the bottom for pan, zoom, and home

---

## Step 6: Exporting Results

### 6.1 Single Sample Export

From the main fitting tab:

**Export Plot:**
- Click "Export Plot"
- Choose format (PNG, PDF, SVG, TIFF)
- 300 DPI publication quality

**Export Fit Data:**
- Click "Export Fit Data"
- Saves peak parameters as CSV

**Export All (Complete):**
- Click "Export All (Complete)"
- Select output folder
- Creates comprehensive package:
  ```
  sample_name/
  ├── sample_name_fit_plot.pdf
  ├── sample_name_fit_plot.png
  ├── sample_name_peak_parameters.csv
  ├── sample_name_spectrum_data.csv
  ├── sample_name_metadata.json
  └── sample_name_fit_report.txt
  ```

### 6.2 Multi-Sample Analysis Export

From the Analysis & Comparison tab:

**Export Complete Analysis:**
1. Click "Export Complete Analysis"
2. Select output folder
3. Creates comprehensive package:
   ```
   analysis_export_YYYYMMDD_HHMMSS/
   ├── plots/
   │   ├── statistics_summary.pdf
   │   ├── statistics_summary.png
   │   ├── comparison_charts.pdf
   │   ├── comparison_charts.png
   │   ├── spectral_overlay.pdf
   │   └── spectral_overlay.png
   ├── sample_data/
   │   ├── C_P7S_199_cluster0/
   │   ├── C_P7S_199_cluster1/
   │   ├── C_P7S_199_cluster2/
   │   ├── C_P7S_199_cluster3/
   │   └── C_P7S_199_cluster4/
   ├── analysis_summary.csv
   ├── all_samples_peak_parameters.csv
   └── analysis_report.txt
   ```

---

## Step 7: Interpreting Results

### 7.1 Understanding Peak Percentages

Peak percentages represent the relative abundance of each carbon functional group:

| High Aromatic (>30%) | Suggests |
|---------------------|----------|
| Black carbon, char, or graphitic material | Combustion products |

| High Aliphatic (>20%) | Suggests |
|----------------------|----------|
| Organic matter with long carbon chains | Lipids, waxes |

| High Carboxyl (>15%) | Suggests |
|---------------------|----------|
| Oxidized organic matter | Degraded material |

| High O-alkyl (>15%) | Suggests |
|---------------------|----------|
| Polysaccharides, cellulose | Plant material |

### 7.2 Comparing Clusters

For the P7S_199 clusters, look for:
- Which cluster has the highest aromatic content?
- Are there spatial variations in carbon speciation?
- Do certain clusters show oxidation signatures?

### 7.3 Quality Control

**Good Fit Indicators:**
- R-squared > 0.99
- Residuals randomly scattered
- Peak positions match expected values

**Warning Signs:**
- Systematic residual patterns
- Peaks shifted far from expected positions
- Very high chi-squared values

---

## Quick Reference

### Keyboard Shortcuts
- **Mouse Wheel**: Zoom in/out on plots
- **Ctrl+Click**: Select multiple samples

### Common Tasks

| Task | Action |
|------|--------|
| Fit spectrum | Click "Fit" |
| Compare samples | Tab 2 -> Add Selected -> Generate plots |
| Export everything | "Export Complete Analysis" |
| Adjust baseline | Enable "Manual Baseline" |
| Change peak selection | Toggle F/Q checkboxes |

### File Formats Required
- Input: CSV with "energy" and "y" columns
- Energy range: ~280-320 eV
- Normalized intensity values

---

## Conclusion

CarbonKPeaks provides a streamlined workflow for:
1. Loading normalized C1s XAS data
2. Fitting peaks with physical constraints
3. Comparing multiple samples
4. Exporting publication-ready results

For questions or issues, refer to the README.md or contact the development team.

---

*CarbonKPeaks - Making C1s XAS Analysis Accessible*

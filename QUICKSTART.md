# CarbonKPeaks Quick Start Guide

## 5-Minute Getting Started

### 1. Launch
Double-click `CarbonKPeaks.exe`

### 2. Load Data
Click **Browse** -> Navigate to folder with CSV files -> **Select Folder**

### 3. Fit a Spectrum
Select a sample -> Click **Fit**

### 4. Compare Samples
- Go to **Analysis & Comparison** tab
- Select samples (Ctrl+Click for multiple)
- Click **Add Selected**
- Click **Comparison Charts**

### 5. Export
Click **Export Complete Analysis** -> Choose folder -> Done!

---

## Data Requirements

Your CSV files must have:
```
energy,y
280.0,0.001
280.5,0.002
...
```

- Energy in eV (280-320 range)
- Normalized intensity values
- Header row with "energy" and "y"

---

## Peak Assignments

| Peak | Energy (eV) | Assignment |
|------|-------------|------------|
| 1 | 284.0 | Quinone |
| 2 | 285.0 | Aromatic |
| 3 | 286.2 | Phenolic |
| 4 | 287.4 | Aliphatic |
| 5 | 288.6 | Carboxyl |
| 6 | 289.3 | O-alkyl |
| 7 | 290.3 | Carbonate |
| 8 | 292.0 | Sigma* |

---

## Key Buttons

| Button | What it does |
|--------|--------------|
| Browse | Load data folder |
| Fit | Fit current spectrum |
| Export All | Full data package |
| Statistics | Multi-sample quality grid |
| Comparison Charts | Bar charts & tables |
| Spectral Overlay | Stacked spectra plot |

---

## Checkboxes

- **F** = Include in **F**itting
- **Q** = Include in **Q**uantification

Uncheck peaks that aren't in your sample!

---

## Output Files

After "Export Complete Analysis":
```
your_export_folder/
├── plots/           <- Publication figures
├── sample_data/     <- Individual CSVs
├── analysis_summary.csv
└── analysis_report.txt
```

---

*That's it! You're ready to analyze C1s XAS data.*

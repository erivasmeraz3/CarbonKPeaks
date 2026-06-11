# CarbonKPeaks — C1s XAS Peak Fitting Software

**Version 1.2**

CarbonKPeaks is a desktop application for fitting and quantifying Carbon K-edge (C 1s) X-ray absorption spectra. It models each spectrum as eight Gaussian peaks on a double-arctangent baseline, assigns the peaks to carbon functional groups, and reports their relative abundances — with tools for testing how many peaks the data actually support and for comparing many samples at once.

For a step-by-step walkthrough with a worked example, see [TUTORIAL.md](TUTORIAL.md).

## Features

- Automated fitting of eight Gaussian peaks on a double-arctangent baseline (lmfit)
- Immediate fit on sample selection, cached per sample
- Fit-quality and uncertainty tools: Box Plot Test, Model Selection, and Monte Carlo (plus an Auto-Reduce tool, beta)
- FWHM coupling modes (per-group or locked per-group) for consistent multi-sample fitting
- Multi-sample comparison figures and C 1s diagnostic ratios
- Session save/load, and a one-row-per-sample Quick Summary CSV
- 300 DPI publication figures

## Installation

### Run the executable (Windows)
Download `CarbonKPeaks.exe` from the [Releases page](https://github.com/erivasmeraz3/CarbonKPeaks/releases) and double-click it — no installation required.

### Run from source
```bash
pip install -r requirements.txt   # tkinter ships with Python
python c1s_peak_viewer_gui_final.py
```

## Quick Start

1. **Load** — click **Add Files** and select spectra, or **Import Athena (.prj)** for an Athena project.
2. **Fit** — click a sample in the list; it fits immediately and the result displays. Fits are cached per sample.
3. **Adjust** — uncheck a peak's Include box to drop it from the fit; set custom centers or an FWHM mode in the Peak Centers tab.
4. **Compare** — on the Analysis & Comparison tab, select samples, click **Add Selected**, and pick a figure.
5. **Export** — **Quick Summary CSV** writes one row per sample with statistics, percentages, and diagnostic ratios.

## Input Data Format

Two-column spectra: energy in eV (≈280–320) and normalized intensity. Accepted formats are `.csv` (with an `energy,y` header), `.dat`, `.txt`, `.xmu`, `.xdi`, `.nor`, plus Athena `.prj` projects. Plain-text files are read as whitespace-separated columns with `#` comment lines ignored.

```csv
energy,y
279.81,0.00331595
280.21,-0.03049639
...
```

## Peak Assignments

The eight-peak scheme follows Solomon et al. (2005). The centers below are the defaults; each is adjustable in the Peak Centers tab.

| Peak | Center (eV) | Assignment |
|------|-------------|------------|
| 1 | 284.4 | Quinone / aromatic C=O |
| 2 | 285.3 | Aromatic C=C (1s→π*) |
| 3 | 286.8 | Phenolic / ketone C–OH, C=O |
| 4 | 287.6 | Aliphatic C–H |
| 5 | 288.2 | Carboxyl COOH (1s→π*) |
| 6 | 289.3 | O-alkyl C–O |
| 7 | 290.3 | Carbonate CO₃ (1s→π*) / organic C σ* — see note |
| 8 | 291.5 | C–C (1s→σ*) |

The ~290.3 eV peak is assignment-ambiguous. In carbonate-bearing samples it is the carbonate 1s→π* resonance, the diagnostic feature for inorganic carbonate; in organic-matter spectra the same 290–292 eV window is instead dominated by 1s→σ* transitions of organic carbon (the σ* peaks in Solomon's scheme). Interpret peak 7 according to your sample type. Peaks 1–6 share one fitted width; the carbonate and σ* peaks have their own.

## Interface

The application has three tabs:

- **Fit Viewer** — load and fit spectra, toggle peaks, adjust the baseline, and run the fit-quality and uncertainty tools. The controls, sample list, and plot are separated by draggable dividers.
- **Analysis & Comparison** — build comparison figures (stacked and grouped bars, mean ± std, ratio heatmap, pie, summary table, spectral overlay, per-group trends) and read the diagnostic ratios.
- **Peak Centers** — set custom peak positions and ranges, and the FWHM values, bounds, and coupling mode.

## Output Files

- **Peak parameters CSV** — one row per peak: name, center (eV), height, FWHM (eV), area, and percentage of quantified area.
- **Spectrum data CSV** — energy, raw intensity, total fit, baseline, residual, and each individual peak.
- **Metadata JSON** — fit quality (R², χ², AIC, BIC), baseline and FWHM parameters, peak summary, and settings.

The full set of export options is described in the [tutorial](TUTORIAL.md#step-8-exporting-results).

## Troubleshooting

- **Fit does not converge** — disable peaks that are not present, or widen their center ranges in the Peak Centers tab.
- **Poor baseline** — enable manual baseline adjustment, and confirm the data are normalized.
- **Peaks at the wrong position** — enable Custom Peak Centers and check the energy calibration.

## Contact

Bug reports and feature requests please reach out.
email: erivasmeraz@lbl.gov

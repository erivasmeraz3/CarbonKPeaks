# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CarbonKPeaks is a Python desktop application (Tkinter + Matplotlib) for fitting and analyzing Carbon K-edge (C1s) X-ray Absorption Spectroscopy (XAS) data. It uses lmfit for non-linear least-squares fitting with a double arctangent baseline and up to 8 Gaussian peaks.

## Commands

```bash
# Run the GUI application
python c1s_peak_viewer_gui_final.py
python c1s_peak_viewer_gui_final.py --spectra_dir /path/to/csv/files

# Run command-line fitting on a single spectrum
python c1s_fitter_optimized.py <spectrum_file> [--output DIR] [--plot]

# Build standalone Windows executable (installs PyInstaller if needed)
python build.py
# Output: dist/CarbonKPeaks.exe
```

No test suite, linter, or formatter is configured.

## Architecture

The application has two core modules:

### `c1s_fitter_optimized.py` (~475 lines) — Fitting Engine
Pure computation, no GUI. Key functions:
- `double_arctangent()` / `gaussian()` — baseline and peak shape primitives
- `total_model()` — complete spectral model (baseline + 8 peaks)
- `estimate_baseline_parameters()` — data-adaptive initial guesses
- `setup_parameters_optimized()` — lmfit parameter initialization with constraints
- `fit_spectrum()` — executes fitting, returns lmfit result
- `load_spectrum()` — reads .csv, .xmu, .xdi, .nor formats
- `calculate_peak_areas()` — integrates fitted peaks for quantification
- `validate_baseline()` — checks baseline quality

The model uses three FWHM groups: `main_fwhm` (peaks 1-6), `carb_fwhm` (peak 7, carbonate), `sigma_fwhm` (peak 8, sigma*).

### `c1s_peak_viewer_gui_final.py` (~4,600 lines) — GUI Application
Single class `C1sPeakViewerFinal` containing the entire GUI. Imports fitting functions from the fitter module. Organized into three tabs:
- **Tab 1 (Peak Fitting)**: Load spectra, fit individual samples, F/Q checkboxes per peak, manual baseline controls with lock options, matplotlib canvas with residuals
- **Tab 2 (Analysis & Comparison)**: Multi-sample comparison with statistics grids, stacked bar charts, spectral overlays
- **Tab 3 (Peak Configuration)**: Custom peak center positions and allowed fitting ranges

Key patterns:
- `self.file_state_cache` — dict keyed by file path storing per-file F/Q checkbox state, fit results, energy, and intensity so switching between files preserves state
- Fitting runs in a background thread (`threading`) to keep the GUI responsive
- Matplotlib embedded via `FigureCanvasTkAgg` with `TkAgg` backend (set before pyplot import)
- The GUI falls back to importing `c1s_fitter_final` if `c1s_fitter_optimized` is unavailable

### `build.py` — PyInstaller Build Script
Builds a single-file Windows executable. Handles splash screen, icon embedding, hidden imports for scipy/lmfit, and data file bundling.

## Fitting Model Details

8 Gaussian peaks on a double arctangent baseline for C1s XAS:

| Peak | ~eV   | Assignment        |
|------|-------|-------------------|
| 1    | 284.0 | Quinone C=O       |
| 2    | 285.0 | Aromatic C=C      |
| 3    | 286.2 | Phenolic C-OH     |
| 4    | 287.4 | Aliphatic C-H     |
| 5    | 288.6 | Carboxyl COOH     |
| 6    | 289.3 | O-alkyl C-O       |
| 7    | 290.3 | Carbonate CO3     |
| 8    | 292.0 | Sigma* C-C        |

Baseline: two arctangent step functions at ~290.5 eV and ~292 eV. Constraints are data-adaptive — baseline bounds are computed from the actual data range.

## Dependencies

Python 3.13+ with: numpy, scipy, pandas, matplotlib, seaborn, lmfit, Pillow, scikit-learn. Managed via `.venv` virtual environment (no requirements.txt).

## Directory Layout

- `archive/` — old script versions, development docs, test data (not part of the active application)
- `dist/` — built executable output
- `carbonpeaks.ico`, `carbonpeaks_*.png` — application icons

# CarbonKPeaks Tutorial
## A Step-by-Step Guide to C1s XAS Peak Fitting

---

## Introduction

This tutorial demonstrates how to use CarbonKPeaks to analyze Carbon K-edge XAS data. The example dataset is a set of cluster spectra from sample P7S_199:

- C_P7S_199_cluster0.csv
- C_P7S_199_cluster1.csv
- C_P7S_199_cluster2.csv
- C_P7S_199_cluster3.csv
- C_P7S_199_cluster4.csv

These represent five spectral clusters from a single sample, which makes them a convenient set for comparing carbon speciation across regions of one specimen.

---

## Step 1: Loading Your Data

### 1.1 Launch the application
Double-click `CarbonKPeaks.exe` or run `python c1s_peak_viewer_gui_final.py`. The main window opens with three tabs: Fit Viewer, Analysis & Comparison, and Peak Centers.

### 1.2 Add your spectra
1. Click **Add Files** in the File Control panel.
2. Select the spectrum files. The dialog accepts `.csv`, `.dat`, `.txt`, `.xmu`, `.xdi`, and `.nor` files, and several files can be selected at once.
  *Note* You may need to be sure to select the option to view all file extensions.
3. The application remembers the folder between launches, so subsequent sessions open the dialog where you left off.

If your data live in an Athena project, use **Import Athena (.prj)** instead. Each group in the project is extracted, normalized using the project's stored parameters, and written out as a CSV alongside the project file.

### 1.3 Verify the files loaded
The samples appear in the Samples list, and the first sample is selected and fit automatically:

```
Samples:
  C_P7S_199_cluster0
  C_P7S_199_cluster1
  C_P7S_199_cluster2
  C_P7S_199_cluster3
  C_P7S_199_cluster4
```

The order of the list can be changed with the **▲ Up** / **▼ Down** buttons or **Sort Naturally (A-Z / 0-9)**, which sorts embedded numbers numerically rather than alphabetically.

The three panels of the Fit Viewer — controls, sample list, and plot — are separated by draggable dividers. If long sample names are cut off, drag the divider to the right of the list to widen it.

---

## Step 2: Fitting Individual Spectra

### 2.1 Select a sample
Click `C_P7S_199_cluster0` in the sample list. There is no separate Fit button: selection triggers the fit, which runs the double-arctangent baseline and all enabled Gaussian peaks against the data. Each sample's fit is cached, so returning to a sample you have already fit is immediate.

If a saved fit state (`<sample>_fit.json`) exists next to the file, the application offers to load it instead of refitting.

### 2.2 Understanding the display

**Top plot — spectrum and fit:**
- Black points: raw data
- Red line: total fit
- Dashed line: baseline
- Colored fills: individual Gaussian components

**Bottom plot — residual:**
The difference between data and fit. A good fit has residuals scattered randomly around zero; structure in the residual indicates a missing or misplaced component.

### 2.3 Reviewing peak parameters

The peak parameters table shows one row per peak:

| Column | Meaning |
|--------|---------|
| Include | Y/- — whether the peak is in the fit and the quantification |
| Peak | Carbon functional group |
| Center | Fitted peak position (eV) |
| FWHM | Fitted peak width (eV) |
| %Area | Percentage of the total quantified area |
| ± | Relative area uncertainty, propagated from the fit's height and width standard errors; shown as an em dash when the fit cannot estimate it |

**Example results for cluster 0:**
```
Peak          Center    FWHM    %Area
Quinone       284.4     0.82     5.2%
Aromatic      285.3     0.82    32.5%
Phenolic      286.8     0.82    15.2%
Aliphatic     287.6     0.82    12.8%
Carboxyl      288.2     0.82    11.3%
O-alkyl       289.3     0.82    14.2%
Carbonate     290.3     0.85     3.5%
Sigma*        291.5     1.60     5.3%
```

Note that the six main peaks share one fitted FWHM; the carbonate and sigma* peaks have their own. The Fit Quality panel reports R², reduced χ², and the AIC and BIC information criteria for the current fit.

The ~290.3 eV "carbonate" peak is assignment-ambiguous: it is the carbonate 1s→π* resonance in carbonate-bearing samples, but in organic-matter spectra the same 290–292 eV region reflects 1s→σ* transitions of organic carbon. Read peak 7 in light of what your sample is likely to contain.

---

## Step 3: Customizing the Fit

### 3.1 Selecting peaks

If a peak is not present in your sample, uncheck its **Include** box. The spectrum refits immediately without that peak, restarting the optimizer from default starting values so the solution is not biased toward the previous one. The **All** checkbox toggles every peak at once.

Excluded peaks are removed from both the fit and the percentage calculation, so the reported percentages always refer to the set of peaks you consider present.

### 3.2 Custom peak centers and ranges

The defaults suit most natural organic matter, but the **Peak Centers** tab allows full control:

1. Check **Use Custom Peak Centers**.
2. Edit the center position and allowed range (±) for any peak. The min/max labels update as you type.
3. Return to the Fit Viewer; subsequent fits use the custom configuration.

### 3.3 FWHM settings and coupling modes

The same tab configures peak widths, independently of the custom-centers checkbox. Each group — main peaks (1–6), carbonate (7), and sigma* (8) — has a starting value and allowed bounds. Beneath the table, the **FWHM Coupling Mode** offers the two arrangements supported by the C K-edge literature, both of which keep the broad sigma* peak at its own width rather than forcing it to the narrow π* width:

- **Per-group** (default, after Solomon): peaks 1–6 share one fitted width; carbonate and sigma* vary independently within their bounds.
- **Locked per-group**: each group is fixed at its configured value and does not vary during the fit.

Locked per-group is useful when comparing many samples: fixing the widths removes a source of sample-to-sample variation, so differences in the percentages reflect composition rather than fitting freedom. To pin a group to an exact width while still letting the fit run, set that group's **Min FWHM = Max FWHM** in the table above. The button **Use current sample's FWHM as session default** copies the converged widths of the current fit into the group fields and switches to locked per-group mode; follow it with **File → Refit All Samples** to apply the constraint across the whole set.

---

## Step 4: Manual Baseline Adjustment

For challenging spectra the automatic baseline can be overridden.

1. Click **Show Advanced Baseline Controls** in the Fit Viewer to expand the panel.
2. Check **Manual Adjustment**. The baseline sliders become active.
3. Adjust the two arctangent steps. Step 1 is the main absorption edge; step 2 accounts for the continuum rise near 292 eV. Each step has a center, height, and width (steepness).
4. **Lock Heights** and **Lock Widths** tie the two steps' parameters together, which halves the number of choices when a single constraint suffices.

With manual adjustment active, the baseline is held fixed at your values and only the peaks are fit.

---

## Step 5: How Many Peaks Does the Data Support?

A fit with all eight peaks is not necessarily the right fit; overlapping Gaussians can trade area without changing the total curve. The Actions panel provides tools to check whether a fit is well constrained and to attach uncertainties to it.

- **Box Plot Test** repeats the fit from many random starting conditions and shows the spread of each peak's percentage. A tight spread indicates the fit converges to a single, well-defined minimum; a wide spread means the peak set has more freedom than the data can constrain, and should be reconsidered.
- **Model Selection** compares candidate peak subsets by information criteria (AIC/BIC) and reports the ranking, which is useful when two chemical interpretations are plausible.
- **Monte Carlo** propagates spectral noise through the fit to give uncertainty estimates on the percentages.

A defensible workflow is: fit with a physically motivated peak set (disabling peaks not expected in the sample), confirm the fit is well constrained with the **Box Plot Test**, and then run **Monte Carlo** to attach uncertainties to the reported percentages. Use **Model Selection** when you need to choose between competing peak sets.

> **Auto-Reduce (Beta).** The Auto-Reduce and Batch Auto-Reduce buttons remove peaks automatically by ΔAIC. This feature is still under testing; confirm any reduced peak set with the Box Plot Test before relying on it.

---

## Step 6: Multi-Sample Comparison

### 6.1 Add samples to the comparison
1. Open the **Analysis & Comparison** tab.
2. Select the five cluster samples in the list (Ctrl+click for several, Shift+click for a range).
3. Click **Add Selected**. The samples appear under "Selected for Analysis"; any that have not been fit yet are fit on the way in.

### 6.2 The Figure Viewer
Choose a figure from the dropdown; it renders immediately:

- **Stacked Bar** — composition of each sample in one bar
- **Grouped Bar** — functional groups side by side across samples
- **Mean ± Std** — average composition with variability
- **Ratio Heatmap** — each sample relative to the first
- **Pie Chart** — average composition
- **Summary Table** — the numbers behind the figures
- **Spectral Overlay** — spectra with fits, offset vertically; the Overlay Options panel sets the offset and toggles individual peaks
- **Trend: <group>** — one functional group followed across the sample series

The ◀ / ▶ buttons cycle through figure types. The Summary (Average) panel shows mean percentages for the selection, and the Ratios panel reports three diagnostics: the aromaticity index Ar/(Ar+Al), the carboxyl-to-aromatic ratio, and the oxygenated-carbon fraction (phenolic + carboxyl + O-alkyl + carbonate over total).

### 6.3 Zoom and navigate
The scroll wheel zooms the plot; the toolbar below it provides pan, zoom-to-rectangle, and reset.

---

## Step 7: Saving Your Work

The File menu separates two kinds of saving:

- **Save Fit State** (Ctrl+S) writes the current sample's fit — parameters, included peaks, and data — to a JSON file. Saving it as `<sample>_fit.json` next to the spectrum makes the application offer to reload it whenever that file is opened.
- **Save Session** (Ctrl+Shift+S) writes the entire analysis: every loaded sample, every cached fit, the peak configuration, and the FWHM mode. **Load Session** (Ctrl+Shift+L) restores it.

Session files store fits as parameter snapshots. When a restored sample is first displayed, it is refit under the current constraints so that everything on screen reflects one consistent model; run **File → Refit All Samples** after loading to bring the whole set up to date at once.

---

## Step 8: Exporting Results

### 8.1 Single-sample exports (Fit Viewer tab)

- **Export Plot** — the current figure as PNG, PDF, SVG, or TIFF at 300 DPI.
- **Export Fit Data** — peak parameters as CSV.
- **Export All (Complete)** — a folder per sample:
  ```
  sample_name/
  ├── sample_name_fit_plot.pdf
  ├── sample_name_fit_plot.png
  ├── sample_name_peak_parameters.csv
  ├── sample_name_spectrum_data.csv
  ├── sample_name_metadata.json
  └── sample_name_fit_report.txt
  ```
- **Batch Export All Samples** — the same package for every loaded sample plus a combined summary CSV.
- **Quick Summary CSV** — one row per sample with R², reduced χ², AIC, BIC, peak count, all percentages, and the three diagnostic ratios. Cached fits are reused, so this is fast even for large sample sets, and it is usually the most convenient route into further analysis elsewhere.

### 8.2 Comparison exports (Analysis & Comparison tab)

- **Export Current Figure** — the displayed figure as an image file.
- **Export Data (CSV)** — the comparison percentages as a table.
- **Export Report (PDF)** — a multi-page PDF of the figure types you select.
- **Export All Figures** — each selected figure type as PDF and PNG files in a folder.

## Quick Reference

### Keyboard shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+S / Ctrl+L | Save / load the current sample's fit state |
| Ctrl+Shift+S / Ctrl+Shift+L | Save / load the full session |
| Mouse wheel | Zoom plots; scroll the control panel |
| Ctrl+click, Shift+click | Multi-select samples |

### Common tasks

| Task | Action |
|------|--------|
| Fit a spectrum | Select it in the Samples list |
| Remove a peak from the model | Uncheck its Include box |
| Find the minimal peak set | Auto-Reduce |
| Apply one set of constraints to all samples | File → Refit All Samples |
| Compare samples | Analysis tab → Add Selected → choose a figure |
| Get all numbers as one table | Quick Summary CSV |

### Input formats

Two-column energy/intensity data in `.csv` (with header), `.dat`, `.txt`, `.xmu`, `.xdi`, or `.nor`; Athena `.prj` projects via the import button. Energy in eV over roughly 280–320, intensity normalized.

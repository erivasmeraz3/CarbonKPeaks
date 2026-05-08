"""
C1s XAS Peak Fitting - OPTIMIZED VERSION

Improvements over c1s_fitter_final.py:
- Data-adaptive baseline constraints
- Better initial guesses based on data
- Ensures baseline doesn't exceed data range
- Improved convergence for difficult spectra
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters
import pandas as pd
from pathlib import Path
import argparse
import gzip
import re


def double_arctangent(x, center1, height1, width1, center2, height2, width2):
    """Double arctangent baseline."""
    step1 = height1 * (0.5 + (1/np.pi) * np.arctan((x - center1) / abs(width1)))
    step2 = height2 * (0.5 + (1/np.pi) * np.arctan((x - center2) / abs(width2)))
    return step1 + step2


def gaussian(x, center, height, fwhm):
    """Gaussian peak."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def total_model(x, arc1_center, arc1_height, arc1_width, arc2_center, arc2_height, arc2_width,
                c1, c2, c3, c4, c5, c6, c7, c8,
                h1, h2, h3, h4, h5, h6, h7, h8,
                main_fwhm, carb_fwhm, sigma_fwhm):
    """Complete model: baseline + 8 peaks."""
    baseline = double_arctangent(x, arc1_center, arc1_height, arc1_width,
                                 arc2_center, arc2_height, arc2_width)

    # Main peaks (shared FWHM)
    peaks = gaussian(x, c1, h1, main_fwhm)
    peaks += gaussian(x, c2, h2, main_fwhm)
    peaks += gaussian(x, c3, h3, main_fwhm)
    peaks += gaussian(x, c4, h4, main_fwhm)
    peaks += gaussian(x, c5, h5, main_fwhm)
    peaks += gaussian(x, c6, h6, main_fwhm)

    # Carbonate
    peaks += gaussian(x, c7, h7, carb_fwhm)

    # Sigma*
    peaks += gaussian(x, c8, h8, sigma_fwhm)

    return baseline + peaks


def estimate_baseline_parameters(energy, intensity):
    """
    Estimate baseline parameters from data characteristics.

    Returns smart initial guesses and bounds based on actual data.
    """
    # Data characteristics
    data_min = np.min(intensity)
    data_max = np.max(intensity)
    data_range = data_max - data_min

    # Find intensity at key energy points
    # Pre-edge region (before 285 eV)
    pre_edge_mask = energy < 285
    if np.any(pre_edge_mask):
        pre_edge_level = np.median(intensity[pre_edge_mask])
    else:
        pre_edge_level = data_min

    # Post-edge region (around 292-294 eV)
    post_edge_mask = (energy >= 292) & (energy <= 294)
    if np.any(post_edge_mask):
        post_edge_level = np.median(intensity[post_edge_mask])
    else:
        post_edge_level = data_max

    # Mid-edge region (around 290 eV)
    mid_edge_mask = (energy >= 289) & (energy <= 291)
    if np.any(mid_edge_mask):
        mid_edge_level = np.median(intensity[mid_edge_mask])
    else:
        mid_edge_level = (pre_edge_level + post_edge_level) / 2

    # Estimate step heights
    # First step: from pre-edge to mid-edge
    step1_height = mid_edge_level - pre_edge_level

    # Second step: from mid-edge to post-edge
    step2_height = post_edge_level - mid_edge_level

    # Ensure heights are positive and reasonable
    step1_height = max(0.1, min(step1_height, data_range * 0.8))
    step2_height = max(0.1, min(step2_height, data_range * 0.8))

    return {
        'data_min': data_min,
        'data_max': data_max,
        'data_range': data_range,
        'pre_edge_level': pre_edge_level,
        'post_edge_level': post_edge_level,
        'mid_edge_level': mid_edge_level,
        'step1_height': step1_height,
        'step2_height': step2_height
    }


def setup_parameters_optimized(energy, intensity):
    """
    Set up parameters with data-adaptive constraints.

    Key improvements:
    - Baseline heights bounded by actual data range
    - Initial guesses based on data characteristics
    - Constraint: total baseline cannot exceed data max
    """
    params = Parameters()

    # Estimate data characteristics
    est = estimate_baseline_parameters(energy, intensity)

    # Baseline - Step 1
    params.add('arc1_center', value=290.5, min=290.2, max=290.8)
    params.add('arc1_height',
              value=est['step1_height'],
              min=est['data_range'] * 0.1,  # At least 10% of range
              max=est['data_range'] * 0.9)  # At most 90% of range
    params.add('arc1_width', value=0.4, min=0.2, max=0.6)

    # Baseline - Step 2
    params.add('arc2_center', value=292.0, min=291.7, max=292.3)
    params.add('arc2_height',
              value=est['step2_height'],
              min=est['data_range'] * 0.1,
              max=est['data_range'] * 0.9)
    params.add('arc2_width', value=0.4, min=0.2, max=0.6)

    # Constraint: Combined baseline at post-edge should not exceed data max
    # At x >> arc2_center: baseline ≈ arc1_height + arc2_height
    # This should be less than or equal to data_max
    params.add('baseline_total',
              expr='arc1_height + arc2_height',
              min=0,
              max=est['data_max'] * 1.1)  # Allow 10% overshoot for fitting tolerance

    # Main peaks
    peak_positions = [284.4, 285.3, 286.8, 287.6, 288.2, 289.3]

    # Estimate reasonable peak height based on data range
    max_peak_height = est['data_range'] * 1.5  # Peaks can be 1.5x the data range

    for i, pos in enumerate(peak_positions, 1):
        params.add(f'c{i}', value=pos, min=pos-0.3, max=pos+0.3)
        params.add(f'h{i}', value=est['data_range'] * 0.3, min=0, max=max_peak_height)

    params.add('main_fwhm', value=1.5, min=0.8, max=2.0)

    # Carbonate
    params.add('c7', value=290.3, min=290.0, max=290.6)
    params.add('h7', value=est['data_range'] * 0.2, min=0, max=max_peak_height * 0.8)
    params.add('carb_fwhm', value=0.8, min=0.5, max=1.0)

    # Sigma*
    params.add('c8', value=291.5, min=291.3, max=294.0)
    params.add('h8', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
    params.add('sigma_fwhm', value=0.6, min=0.4, max=2.0)

    return params


def parse_athena_prj(filepath):
    """Parse an Athena .prj file and return list of spectrum dicts.

    Athena project files are gzip-compressed Perl Data::Dumper output.
    Each record has: $old_group, @args (metadata), @x (energy), @y (raw mu).
    """
    with gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as f:
        content = f.read()

    spectra = []
    blocks = re.split(r"\$old_group\s*=\s*'([^']+)';", content)

    for i in range(1, len(blocks), 2):
        group_name = blocks[i]
        block = blocks[i + 1] if i + 1 < len(blocks) else ""

        label_match = re.search(r"'label','([^']*)'", block)
        label = label_match.group(1) if label_match else group_name

        args_match = re.search(r"@args\s*=\s*\((.*?)\);", block, re.DOTALL)
        params = {}
        if args_match:
            args_str = args_match.group(1)
            tokens = re.findall(r"'([^']*)'|([^,'\s]+)", args_str)
            flat = [t[0] if t[0] else t[1] for t in tokens]
            for j in range(0, len(flat) - 1, 2):
                params[flat[j]] = flat[j + 1]

        x_match = re.search(r"@x\s*=\s*\(([^)]+)\)", block)
        if not x_match:
            continue
        x_vals = [float(v.strip("' ")) for v in x_match.group(1).split(",") if v.strip("' ")]

        y_match = re.search(r"@y\s*=\s*\(([^)]+)\)", block)
        if not y_match:
            continue
        y_vals = [float(v.strip("' ")) for v in y_match.group(1).split(",") if v.strip("' ")]

        energy = np.array(x_vals)
        mu = np.array(y_vals)

        if len(energy) != len(mu):
            min_len = min(len(energy), len(mu))
            energy = energy[:min_len]
            mu = mu[:min_len]

        spectra.append({
            'group': group_name,
            'label': label,
            'energy': energy,
            'mu': mu,
            'params': params,
        })

    return spectra


def normalize_athena_spectrum(energy, mu, params):
    """Normalize a spectrum using Athena's stored parameters.

    Uses the stored pre-edge line (bkg_slope, bkg_int) and edge step
    (bkg_step / bkg_fitted_step) for standard XANES normalization:
        norm(E) = (mu(E) - pre_edge_line(E)) / edge_step
    """
    e0 = float(params.get('bkg_e0', 284.0))
    eshift = float(params.get('bkg_eshift', 0))
    e0 += eshift

    # Pre-edge line: mu_pre = slope * E + intercept
    slope = float(params.get('bkg_slope', 0))
    intercept = float(params.get('bkg_int', 0))
    pre_edge = slope * energy + intercept

    mu_sub = mu - pre_edge

    # Edge step for normalization
    step = float(params.get('bkg_step', 0))
    if step == 0 or abs(step) < 1e-12:
        step = float(params.get('bkg_fitted_step', 1.0))

    if abs(step) > 1e-12:
        normalized = mu_sub / step
    else:
        pre1 = float(params.get('bkg_pre1', -20))
        pre2 = float(params.get('bkg_pre2', -10))
        nor1 = float(params.get('bkg_nor1', 25))
        nor2 = float(params.get('bkg_nor2', 50))

        pre_mask = (energy >= (e0 + pre1)) & (energy <= (e0 + pre2))
        post_mask = (energy >= (e0 + nor1)) & (energy <= (e0 + nor2))

        if np.sum(pre_mask) >= 2 and np.sum(post_mask) >= 2:
            pre_val = np.mean(mu[pre_mask])
            post_val = np.mean(mu[post_mask])
            fallback_step = post_val - pre_val
            if abs(fallback_step) > 1e-12:
                normalized = (mu - pre_val) / fallback_step
            else:
                normalized = mu_sub
        else:
            normalized = mu_sub

    return normalized


def load_spectrum(file_path):
    """Load spectrum from various formats."""
    file_path = Path(file_path)
    suffix = file_path.suffix.lower()

    if suffix == '.csv':
        df = pd.read_csv(file_path)
        # Try standard column names
        if 'energy' in df.columns and 'y' in df.columns:
            energy = df['energy'].values
            intensity = df['y'].values
        else:
            energy = df.iloc[:, 0].values
            intensity = df.iloc[:, 1].values

    elif suffix in ['.xmu', '.xdi', '.nor']:
        data = np.loadtxt(file_path, comments='#')
        energy = data[:, 0]
        intensity = data[:, 1]

    else:
        data = np.loadtxt(file_path)
        energy = data[:, 0]
        intensity = data[:, 1]

    # Filter to ≤294 eV (matching v5_fixed)
    mask = energy <= 294
    return energy[mask], intensity[mask]


def fit_spectrum(energy, intensity):
    """Fit spectrum using lmfit with optimized parameters."""
    model = Model(total_model)
    params = setup_parameters_optimized(energy, intensity)
    result = model.fit(intensity, params, x=energy, method='leastsq', max_nfev=10000)
    return result


def calculate_peak_areas(energy, result):
    """Calculate peak areas using the analytical Gaussian integral.

    Area = height * sigma * sqrt(2*pi) = height * fwhm * sqrt(pi / (4*ln(2)))
    This is exact and avoids truncation errors at the energy range boundaries.
    """
    areas = {}
    names = ['Quinone', 'Aromatic', 'Phenolic', 'Aliphatic',
             'Carboxyl', 'O_alkyl', 'Carbonate', 'Sigma_star']

    # Analytical area factor: sqrt(pi / (4 * ln(2)))
    area_factor = np.sqrt(np.pi / (4 * np.log(2)))

    for i, name in enumerate(names, 1):
        height = result.params[f'h{i}'].value

        if i <= 6:
            fwhm = result.params['main_fwhm'].value
        elif i == 7:
            fwhm = result.params['carb_fwhm'].value
        else:
            fwhm = result.params['sigma_fwhm'].value

        areas[name] = height * fwhm * area_factor

    return areas


def validate_baseline(energy, result, intensity):
    """
    Validate that baseline is reasonable.

    Returns (is_valid, issues) where issues is a list of problems found.
    """
    issues = []

    # Calculate baseline
    baseline = double_arctangent(energy,
                                result.params['arc1_center'].value,
                                result.params['arc1_height'].value,
                                result.params['arc1_width'].value,
                                result.params['arc2_center'].value,
                                result.params['arc2_height'].value,
                                result.params['arc2_width'].value)

    # Check 1: Baseline should not exceed data max significantly
    data_max = np.max(intensity)
    data_min = np.min(intensity)
    data_range = data_max - data_min
    baseline_max = np.max(baseline)
    baseline_min = np.min(baseline)

    if baseline_max > data_max + data_range * 0.2:
        issues.append(f"Baseline exceeds data max: {baseline_max:.3f} > {data_max:.3f}")

    # Check 2: Baseline should not go below data min significantly
    if baseline_min < data_min - data_range * 0.2:
        issues.append(f"Baseline below data min: {baseline_min:.3f} < {data_min:.3f}")

    # Check 3: Post-edge baseline should be reasonable
    post_edge_mask = energy >= 292
    if np.any(post_edge_mask):
        post_edge_baseline = np.mean(baseline[post_edge_mask])
        post_edge_data = np.mean(intensity[post_edge_mask])

        if post_edge_baseline > post_edge_data * 1.3:
            issues.append(f"Post-edge baseline too high: {post_edge_baseline:.3f} vs data {post_edge_data:.3f}")

    # Check 4: Second step should not be unreasonably large
    arc1_height = result.params['arc1_height'].value
    arc2_height = result.params['arc2_height'].value

    if arc2_height > arc1_height * 3:
        issues.append(f"Second step much larger than first: {arc2_height:.3f} vs {arc1_height:.3f}")

    is_valid = len(issues) == 0
    return is_valid, issues


def print_summary(result, energy, intensity):
    """Print fit summary with baseline validation."""
    r_squared = 1 - result.residual.var() / np.var(result.data)

    print("\n" + "="*80)
    print("C1s XAS PEAK FITTING SUMMARY (Optimized)")
    print("="*80)

    print(f"\nFit Quality:")
    print(f"  R^2 = {r_squared:.4f}")
    print(f"  Reduced chi^2 = {result.redchi:.4e}")
    print(f"  AIC = {result.aic:.2f}")
    print(f"  BIC = {result.bic:.2f}")

    # Validate baseline
    is_valid, issues = validate_baseline(energy, result, intensity)
    print(f"\nBaseline Validation:")
    if is_valid:
        print("  Status: VALID")
    else:
        print("  Status: WARNING - Issues detected:")
        for issue in issues:
            print(f"    - {issue}")

    areas = calculate_peak_areas(energy, result)
    total_area = sum(areas.values())
    normalized = {k: v/total_area*100 for k, v in areas.items()}

    print(f"\nPeak Parameters:")
    print(f"{'Peak':<14} {'Center (eV)':<12} {'FWHM (eV)':<12} {'% Total':<10}")
    print("-"*80)

    names = list(areas.keys())
    for i, name in enumerate(names, 1):
        center = result.params[f'c{i}'].value

        if i <= 6:
            fwhm = result.params['main_fwhm'].value
        elif i == 7:
            fwhm = result.params['carb_fwhm'].value
        else:
            fwhm = result.params['sigma_fwhm'].value

        print(f"{name:<14} {center:<12.3f} {fwhm:<12.3f} {normalized[name]:<10.1f}")

    print(f"\nBaseline Parameters:")
    print(f"  Step 1: center={result.params['arc1_center'].value:.2f} eV, "
          f"height={result.params['arc1_height'].value:.3f}, "
          f"width={result.params['arc1_width'].value:.3f}")
    print(f"  Step 2: center={result.params['arc2_center'].value:.2f} eV, "
          f"height={result.params['arc2_height'].value:.3f}, "
          f"width={result.params['arc2_width'].value:.3f}")
    print(f"  Total baseline height: {result.params['arc1_height'].value + result.params['arc2_height'].value:.3f}")

    print(f"\nAromaticity: {normalized['Aromatic']:.1f}%")
    print("="*80 + "\n")


def plot_fit(energy, intensity, result, title='C1s XAS Fit', save_path=None):
    """Create plot with baseline validation."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Data and fit
    ax1.plot(energy, intensity, 'ko', markersize=3, alpha=0.6, label='Data')
    ax1.plot(energy, result.best_fit, 'r-', linewidth=2, label='Total Fit')

    # Baseline
    baseline = double_arctangent(energy,
                                 result.params['arc1_center'].value,
                                 result.params['arc1_height'].value,
                                 result.params['arc1_width'].value,
                                 result.params['arc2_center'].value,
                                 result.params['arc2_height'].value,
                                 result.params['arc2_width'].value)
    ax1.plot(energy, baseline, 'k--', linewidth=1.5, alpha=0.5, label='Baseline')

    # Individual peaks
    colors = plt.cm.tab10(np.linspace(0, 1, 8))
    names = ['Quinone', 'Aromatic', 'Phenolic', 'Aliphatic',
             'Carboxyl', 'O-alkyl', 'Carbonate', 'Sigma*']

    for i, (name, color) in enumerate(zip(names, colors), 1):
        center = result.params[f'c{i}'].value
        height = result.params[f'h{i}'].value

        if i <= 6:
            fwhm = result.params['main_fwhm'].value
        elif i == 7:
            fwhm = result.params['carb_fwhm'].value
        else:
            fwhm = result.params['sigma_fwhm'].value

        peak = gaussian(energy, center, height, fwhm)
        ax1.plot(energy, peak + baseline, '--', color=color,
                linewidth=1, alpha=0.7, label=name)

    ax1.set_xlabel('Energy (eV)', fontsize=11)
    ax1.set_ylabel('Normalized Absorption', fontsize=11)

    # Add baseline validation to title
    is_valid, _ = validate_baseline(energy, result, intensity)
    validation_status = "✓" if is_valid else "⚠"
    ax1.set_title(f'{title} {validation_status}', fontsize=13, fontweight='bold')

    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(energy.min(), 294)

    # Residuals
    r_squared = 1 - result.residual.var() / np.var(result.data)
    ax2.plot(energy, result.residual, 'b-', linewidth=1)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Energy (eV)', fontsize=11)
    ax2.set_ylabel('Residuals', fontsize=11)
    ax2.set_title(f'Residuals (R² = {r_squared:.4f})', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(energy.min(), 294)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")

    return fig


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(description='C1s XAS fitter (optimized with adaptive baseline)')
    parser.add_argument('input_file', help='Spectrum file (.xmu, .csv, etc.)')
    parser.add_argument('--output', default=None, help='Output directory')
    parser.add_argument('--plot', action='store_true', help='Show plot')

    args = parser.parse_args()

    # Load
    print(f"Loading: {args.input_file}")
    energy, intensity = load_spectrum(args.input_file)
    print(f"Loaded {len(energy)} points ({energy.min():.2f} - {energy.max():.2f} eV)")

    # Show data characteristics
    est = estimate_baseline_parameters(energy, intensity)
    print(f"\nData characteristics:")
    print(f"  Min: {est['data_min']:.3f}")
    print(f"  Max: {est['data_max']:.3f}")
    print(f"  Range: {est['data_range']:.3f}")
    print(f"  Estimated step 1 height: {est['step1_height']:.3f}")
    print(f"  Estimated step 2 height: {est['step2_height']:.3f}")

    # Fit
    print("\nFitting with lmfit (optimized constraints)...")
    result = fit_spectrum(energy, intensity)

    # Results
    print_summary(result, energy, intensity)

    # Save
    input_path = Path(args.input_file)
    output_dir = Path(args.output) if args.output else input_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    base_name = input_path.stem
    plot_path = output_dir / f"{base_name}_fit_optimized.pdf"
    report_path = output_dir / f"{base_name}_fit_report_optimized.txt"

    with open(report_path, 'w') as f:
        f.write(result.fit_report())
        f.write("\n\nBaseline Validation:\n")
        is_valid, issues = validate_baseline(energy, result, intensity)
        if is_valid:
            f.write("Status: VALID\n")
        else:
            f.write("Status: WARNING\n")
            for issue in issues:
                f.write(f"  - {issue}\n")
    print(f"Fit report: {report_path}")

    plot_fit(energy, intensity, result, title=f'C1s XAS: {base_name}', save_path=plot_path)

    if args.plot:
        plt.show()
    else:
        plt.close()


if __name__ == '__main__':
    main()

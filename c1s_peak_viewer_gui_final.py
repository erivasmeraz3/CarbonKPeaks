#!/usr/bin/env python3
"""
CarbonKPeaks - C1s XAS Peak Fitting Interactive GUI

Refined features:
- Sample list with search
- Manual baseline adjustment with lock option
- Statistics tab for single sample analysis
- Multi-sample comparison tab
- Proper baseline/fit display options
"""

import os
import sys
import json
from pathlib import Path
import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

# Set matplotlib backend before importing pyplot
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import numpy as np
import pandas as pd
import seaborn as sns

# Import C1s fitter
from c1s_fitter_optimized import (load_spectrum, fit_spectrum, calculate_peak_areas,
                                   double_arctangent, gaussian, validate_baseline)


class C1sPeakViewerFinal:
    def __init__(self, root, spectra_dir=None):
        self.root = root
        self.root.title("CarbonKPeaks - C1s XAS Analysis")
        self.root.geometry("1800x1000")

        # Set application icon (try PNG first for better quality, fallback to ICO)
        self._set_app_icon()

        # Data storage
        self.spectra_dir = Path(spectra_dir) if spectra_dir else Path.home()
        self.spectra_files = []
        self.current_index = 0
        self.current_result = None
        self.current_energy = None
        self.current_intensity = None
        self.file_list_modified = False  # Track if files were added manually
        self._fit_in_progress = False  # Track if a fit is currently running

        # Comparison samples (list of tuples: (name, result, energy, peak_sel))
        self.comparison_samples = []

        # Per-file state cache: {file_path_str: {fit_peaks, quant_peaks, result, energy, intensity}}
        self.file_state_cache = {}

        # Display options
        self.show_baseline = tk.BooleanVar(value=True)
        self.show_peaks = tk.BooleanVar(value=True)

        # Manual baseline adjustment
        self.manual_baseline = tk.BooleanVar(value=False)
        self.lock_heights = tk.BooleanVar(value=False)
        self.lock_widths = tk.BooleanVar(value=False)

        self.baseline_params = {
            'arc1_center': tk.DoubleVar(value=290.5),
            'arc1_height': tk.DoubleVar(value=0.5),
            'arc1_width': tk.DoubleVar(value=0.4),
            'arc2_center': tk.DoubleVar(value=292.0),
            'arc2_height': tk.DoubleVar(value=0.5),
            'arc2_width': tk.DoubleVar(value=0.4)
        }

        # Trace height and width changes for bidirectional locking
        self.baseline_params['arc1_height'].trace('w', self.on_height1_change)
        self.baseline_params['arc2_height'].trace('w', self.on_height2_change)
        self.baseline_params['arc1_width'].trace('w', self.on_width1_change)
        self.baseline_params['arc2_width'].trace('w', self.on_width2_change)

        # Peak names
        self.peak_names = ['Quinone', 'Aromatic', 'Phenolic', 'Aliphatic',
                          'Carboxyl', 'O_alkyl', 'Carbonate', 'Sigma_star']
        self.peak_display_names = ['Quinone', 'Aromatic', 'Phenolic', 'Aliphatic',
                                   'Carboxyl', 'O-alkyl', 'Carbonate', 'Sigma*']
        self.peak_colors = plt.cm.tab10(np.linspace(0, 1, 8))
        # Two separate checkboxes: one for fitting, one for quantitation
        self.peak_in_fit = {name: tk.BooleanVar(value=True) for name in self.peak_names}
        self.peak_in_quant = {name: tk.BooleanVar(value=True) for name in self.peak_names}
        # Peak center parameters (default values and ranges from c1s_fitter_optimized.py)
        self.default_peak_centers = {
            'Quinone': 284.4,
            'Aromatic': 285.3,
            'Phenolic': 286.8,
            'Aliphatic': 287.6,
            'Carboxyl': 288.2,
            'O_alkyl': 289.3,
            'Carbonate': 290.3,
            'Sigma_star': 291.5
        }

        self.peak_center_vars = {name: tk.DoubleVar(value=self.default_peak_centers[name])
                                for name in self.peak_names}
        self.peak_range_vars = {name: tk.DoubleVar(value=0.3) for name in self.peak_names}
        # Sigma* has larger range
        self.peak_range_vars['Sigma_star'].set(2.7)  # 291.3 to 294.0

        # FWHM parameters (3 groups)
        self.fwhm_main_var = tk.DoubleVar(value=1.5)  # Peaks 1-6 (default: 1.5, range: 0.8-2.0)
        self.fwhm_main_min_var = tk.DoubleVar(value=0.8)
        self.fwhm_main_max_var = tk.DoubleVar(value=2.0)

        self.fwhm_carb_var = tk.DoubleVar(value=0.8)  # Peak 7 Carbonate (default: 0.8, range: 0.5-1.0)
        self.fwhm_carb_min_var = tk.DoubleVar(value=0.5)
        self.fwhm_carb_max_var = tk.DoubleVar(value=1.0)

        self.fwhm_sigma_var = tk.DoubleVar(value=0.6)  # Peak 8 Sigma* (default: 0.6, range: 0.4-2.0)
        self.fwhm_sigma_min_var = tk.DoubleVar(value=0.4)
        self.fwhm_sigma_max_var = tk.DoubleVar(value=2.0)

        # Peak editing mode
        self.custom_peak_centers = tk.BooleanVar(value=False)

        # Auto-reduce threshold (delta-AIC below which a peak is expendable)
        self.auto_reduce_threshold = tk.DoubleVar(value=2.0)

        # Setup GUI
        self.setup_gui()

        # Create menu bar (File > Save/Load Fit State, Session)
        self._create_menu_bar()

        # Only auto-load if a directory was specified via command line
        if spectra_dir is not None:
            self.load_spectra_list()
            if self.spectra_files:
                self.load_current_spectrum()

    def _get_file_key(self, index=None):
        """Return a string key for the file at the given index (default: current)."""
        if index is None:
            index = self.current_index
        return str(self.spectra_files[index])

    def _save_file_state(self):
        """Snapshot current F/Q checkbox states and fit result into the cache."""
        if not self.spectra_files:
            return
        key = self._get_file_key()
        self.file_state_cache[key] = {
            'fit_peaks': {name: self.peak_in_fit[name].get() for name in self.peak_names},
            'quant_peaks': {name: self.peak_in_quant[name].get() for name in self.peak_names},
            'result': self.current_result,
            'energy': self.current_energy,
            'intensity': self.current_intensity,
        }

    def _restore_file_state(self):
        """Restore cached F/Q states and fit result for the current file.

        Returns True on cache hit, False on miss.
        """
        if not self.spectra_files:
            return False
        key = self._get_file_key()
        cached = self.file_state_cache.get(key)
        if cached is None:
            return False
        for name in self.peak_names:
            self.peak_in_fit[name].set(cached['fit_peaks'][name])
            self.peak_in_quant[name].set(cached['quant_peaks'][name])
        self.current_result = cached['result']
        self.current_energy = cached['energy']
        self.current_intensity = cached['intensity']
        return True

    # ------------------------------------------------------------------
    # Menu Bar  (File > Save/Load Fit State & Session)
    # ------------------------------------------------------------------

    def _create_menu_bar(self):
        """Create the application menu bar."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)

        file_menu.add_command(label="Save Fit State", accelerator="Ctrl+S",
                              command=self.save_fit_state)
        file_menu.add_command(label="Load Fit State", accelerator="Ctrl+L",
                              command=self.load_fit_state)
        file_menu.add_separator()
        file_menu.add_command(label="Save Session", accelerator="Ctrl+Shift+S",
                              command=self.save_session)
        file_menu.add_command(label="Load Session", accelerator="Ctrl+Shift+L",
                              command=self.load_session)

        self.root.bind_all('<Control-s>', lambda e: self.save_fit_state())
        self.root.bind_all('<Control-l>', lambda e: self.load_fit_state())
        self.root.bind_all('<Control-Shift-S>', lambda e: self.save_session())
        self.root.bind_all('<Control-Shift-L>', lambda e: self.load_session())

    # ------------------------------------------------------------------
    # Serialize / Deserialize helpers
    # ------------------------------------------------------------------

    def _serialize_fit_state(self, file_path_str=None):
        """Serialize the fit state for one spectrum to a JSON-safe dict.

        Args:
            file_path_str: Optional file path string. If None, uses current file.

        Returns:
            dict suitable for json.dump()
        """
        if file_path_str is not None:
            cached = self.file_state_cache.get(file_path_str)
            if cached is None:
                return None
            fit_peaks = cached['fit_peaks']
            quant_peaks = cached['quant_peaks']
            result = cached['result']
            energy = cached['energy']
        else:
            fit_peaks = {name: self.peak_in_fit[name].get() for name in self.peak_names}
            quant_peaks = {name: self.peak_in_quant[name].get() for name in self.peak_names}
            result = self.current_result
            energy = self.current_energy

        if result is None:
            return None

        # Serialize lmfit parameters
        param_dict = {}
        for pname, par in result.params.items():
            param_dict[pname] = {
                'value': par.value,
                'vary': par.vary,
                'min': par.min if np.isfinite(par.min) else None,
                'max': par.max if np.isfinite(par.max) else None,
                'expr': par.expr,
            }

        # Fit statistics
        r_squared = 1 - result.residual.var() / np.var(result.data)
        try:
            areas = calculate_peak_areas(energy, result)
        except Exception:
            areas = {name: 0.0 for name in self.peak_names}
        selected_areas = {k: v for k, v in areas.items() if quant_peaks.get(k, False)}
        total_sel = sum(selected_areas.values())
        percentages = {}
        for k in areas:
            if quant_peaks.get(k, False) and total_sel > 0:
                percentages[k] = areas[k] / total_sel * 100
            else:
                percentages[k] = 0.0

        state = {
            'version': 1,
            'fit_peaks': fit_peaks,
            'quant_peaks': quant_peaks,
            'parameters': param_dict,
            'fit_statistics': {
                'r_squared': float(r_squared),
                'redchi': float(result.redchi) if hasattr(result, 'redchi') else None,
                'aic': float(result.aic) if hasattr(result, 'aic') else None,
                'bic': float(result.bic) if hasattr(result, 'bic') else None,
                'n_data': int(result.ndata) if hasattr(result, 'ndata') else len(energy),
                'n_variables': int(result.nvarys) if hasattr(result, 'nvarys') else 0,
            },
            'peak_areas': {k: float(v) for k, v in areas.items()},
            'peak_percentages': {k: float(v) for k, v in percentages.items()},
        }
        return state

    def _reconstruct_result(self, energy, intensity, param_dict):
        """Rebuild an lmfit ModelResult from saved parameter values.

        Args:
            energy: numpy array of energy values
            intensity: numpy array of intensity values
            param_dict: dict of parameter dicts from _serialize_fit_state

        Returns:
            lmfit ModelResult with the saved parameter values
        """
        from lmfit import Model, Parameters
        from c1s_fitter_optimized import total_model

        params = Parameters()

        # First pass: add non-expression parameters
        for pname, pinfo in param_dict.items():
            if pinfo.get('expr') is not None:
                continue
            pmin = pinfo['min'] if pinfo['min'] is not None else -np.inf
            pmax = pinfo['max'] if pinfo['max'] is not None else np.inf
            params.add(pname, value=pinfo['value'], vary=False,
                      min=pmin, max=pmax)

        # Second pass: add expression parameters
        for pname, pinfo in param_dict.items():
            if pinfo.get('expr') is None:
                continue
            params.add(pname, expr=pinfo['expr'])

        # Run a minimal fit to get a valid ModelResult object
        model = Model(total_model)
        result = model.fit(intensity, params, x=energy, method='leastsq', max_nfev=1)
        return result

    # ------------------------------------------------------------------
    # Per-file Save / Load
    # ------------------------------------------------------------------

    def save_fit_state(self):
        """Save current fit state as a JSON sidecar file next to the spectrum."""
        if not self.spectra_files or self.current_result is None:
            messagebox.showwarning("No Fit", "Load and fit a spectrum first.")
            return

        state = self._serialize_fit_state()
        if state is None:
            messagebox.showwarning("No Fit", "No fit result to save.")
            return

        spectrum_path = self.spectra_files[self.current_index]
        default_name = f"{spectrum_path.stem}_fit.json"
        default_dir = str(spectrum_path.parent)

        save_path = filedialog.asksaveasfilename(
            title="Save Fit State",
            initialdir=default_dir,
            initialfile=default_name,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json")]
        )
        if not save_path:
            return

        try:
            with open(save_path, 'w') as f:
                json.dump(state, f, indent=2)
            messagebox.showinfo("Saved", f"Fit state saved to:\n{Path(save_path).name}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save fit state:\n{e}")

    def load_fit_state(self):
        """Load fit state from a JSON sidecar file."""
        if not self.spectra_files:
            messagebox.showwarning("No Spectrum", "Load a spectrum first.")
            return

        spectrum_path = self.spectra_files[self.current_index]

        # Try auto-detecting sidecar next to current file
        auto_path = spectrum_path.parent / f"{spectrum_path.stem}_fit.json"
        if auto_path.exists():
            initial_dir = str(spectrum_path.parent)
            initial_file = auto_path.name
        else:
            initial_dir = str(spectrum_path.parent)
            initial_file = ""

        json_path = filedialog.askopenfilename(
            title="Load Fit State",
            initialdir=initial_dir,
            initialfile=initial_file,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not json_path:
            return

        self._apply_fit_state(json_path)

    def _apply_fit_state(self, json_path):
        """Read a fit-state JSON and apply it to the current spectrum."""
        try:
            with open(json_path, 'r') as f:
                state = json.load(f)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read JSON:\n{e}")
            return

        if self.current_energy is None or self.current_intensity is None:
            messagebox.showwarning("No Data", "Load spectrum data first.")
            return

        # Restore F/Q checkboxes
        for name in self.peak_names:
            if name in state.get('fit_peaks', {}):
                self.peak_in_fit[name].set(state['fit_peaks'][name])
            if name in state.get('quant_peaks', {}):
                self.peak_in_quant[name].set(state['quant_peaks'][name])

        # Reconstruct result
        param_dict = state.get('parameters', {})
        if param_dict:
            try:
                self.current_result = self._reconstruct_result(
                    self.current_energy, self.current_intensity, param_dict)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to reconstruct fit:\n{e}")
                return

        # Update all displays
        self.update_peak_tree_display()
        self.update_plot()
        self.update_statistics()
        self.update_peak_parameters()
        self.update_statistics_tab()
        self._save_file_state()

        messagebox.showinfo("Loaded", f"Fit state loaded from:\n{Path(json_path).name}")

    # ------------------------------------------------------------------
    # Session Save / Load
    # ------------------------------------------------------------------

    def save_session(self):
        """Save the entire session (global config + all per-file states) to JSON."""
        save_path = filedialog.asksaveasfilename(
            title="Save Session",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json")],
            initialfile="carbonkpeaks_session.json"
        )
        if not save_path:
            return

        # Make sure current file state is cached
        self._save_file_state()

        # Global configuration
        global_config = {
            'custom_peak_centers': self.custom_peak_centers.get(),
            'peak_centers': {name: self.peak_center_vars[name].get()
                            for name in self.peak_names},
            'peak_ranges': {name: self.peak_range_vars[name].get()
                           for name in self.peak_names},
            'fwhm_main': self.fwhm_main_var.get(),
            'fwhm_main_min': self.fwhm_main_min_var.get(),
            'fwhm_main_max': self.fwhm_main_max_var.get(),
            'fwhm_carb': self.fwhm_carb_var.get(),
            'fwhm_carb_min': self.fwhm_carb_min_var.get(),
            'fwhm_carb_max': self.fwhm_carb_max_var.get(),
            'fwhm_sigma': self.fwhm_sigma_var.get(),
            'fwhm_sigma_min': self.fwhm_sigma_min_var.get(),
            'fwhm_sigma_max': self.fwhm_sigma_max_var.get(),
            'manual_baseline': self.manual_baseline.get(),
            'lock_heights': self.lock_heights.get(),
            'lock_widths': self.lock_widths.get(),
            'baseline_params': {k: v.get() for k, v in self.baseline_params.items()},
            'auto_reduce_threshold': self.auto_reduce_threshold.get(),
        }

        # Per-file states
        per_file = {}
        for file_key, cached in self.file_state_cache.items():
            state = self._serialize_fit_state(file_key)
            if state is not None:
                per_file[file_key] = state

        session = {
            'version': 1,
            'global_config': global_config,
            'spectra_files': [str(f) for f in self.spectra_files],
            'current_index': self.current_index,
            'per_file_states': per_file,
        }

        try:
            with open(save_path, 'w') as f:
                json.dump(session, f, indent=2)
            messagebox.showinfo("Saved",
                                f"Session saved ({len(per_file)} spectra with fits).")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save session:\n{e}")

    def load_session(self):
        """Load a session from JSON, restoring global config and all per-file states."""
        json_path = filedialog.askopenfilename(
            title="Load Session",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not json_path:
            return

        try:
            with open(json_path, 'r') as f:
                session = json.load(f)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read session file:\n{e}")
            return

        # Restore global config
        gc = session.get('global_config', {})
        if 'custom_peak_centers' in gc:
            self.custom_peak_centers.set(gc['custom_peak_centers'])
        for name in self.peak_names:
            if name in gc.get('peak_centers', {}):
                self.peak_center_vars[name].set(gc['peak_centers'][name])
            if name in gc.get('peak_ranges', {}):
                self.peak_range_vars[name].set(gc['peak_ranges'][name])
        if 'fwhm_main' in gc:
            self.fwhm_main_var.set(gc['fwhm_main'])
        if 'fwhm_main_min' in gc:
            self.fwhm_main_min_var.set(gc['fwhm_main_min'])
        if 'fwhm_main_max' in gc:
            self.fwhm_main_max_var.set(gc['fwhm_main_max'])
        if 'fwhm_carb' in gc:
            self.fwhm_carb_var.set(gc['fwhm_carb'])
        if 'fwhm_carb_min' in gc:
            self.fwhm_carb_min_var.set(gc['fwhm_carb_min'])
        if 'fwhm_carb_max' in gc:
            self.fwhm_carb_max_var.set(gc['fwhm_carb_max'])
        if 'fwhm_sigma' in gc:
            self.fwhm_sigma_var.set(gc['fwhm_sigma'])
        if 'fwhm_sigma_min' in gc:
            self.fwhm_sigma_min_var.set(gc['fwhm_sigma_min'])
        if 'fwhm_sigma_max' in gc:
            self.fwhm_sigma_max_var.set(gc['fwhm_sigma_max'])
        if 'manual_baseline' in gc:
            self.manual_baseline.set(gc['manual_baseline'])
        if 'lock_heights' in gc:
            self.lock_heights.set(gc['lock_heights'])
        if 'lock_widths' in gc:
            self.lock_widths.set(gc['lock_widths'])
        for k, v in gc.get('baseline_params', {}).items():
            if k in self.baseline_params:
                self.baseline_params[k].set(v)
        if 'auto_reduce_threshold' in gc:
            self.auto_reduce_threshold.set(gc['auto_reduce_threshold'])

        # Load spectra files
        file_paths = session.get('spectra_files', [])
        valid_paths = []
        missing = []
        for fp in file_paths:
            p = Path(fp)
            if p.exists():
                valid_paths.append(p)
            else:
                missing.append(fp)

        if missing:
            messagebox.showwarning("Missing Files",
                                   f"{len(missing)} spectrum file(s) not found and will be skipped.")

        self.spectra_files = valid_paths
        self.file_state_cache.clear()

        # Populate sample listbox
        self.sample_listbox.delete(0, tk.END)
        for fp in self.spectra_files:
            self.sample_listbox.insert(tk.END, fp.stem)

        # Restore per-file states
        per_file = session.get('per_file_states', {})
        restored_count = 0
        for file_key, state in per_file.items():
            # Check the file still exists in our loaded list
            if not Path(file_key).exists():
                continue
            # Load data for this file
            try:
                energy, intensity = load_spectrum(Path(file_key))
            except Exception:
                continue

            # Restore checkboxes
            fit_peaks = state.get('fit_peaks', {name: True for name in self.peak_names})
            quant_peaks = state.get('quant_peaks', {name: True for name in self.peak_names})

            # Reconstruct result
            param_dict = state.get('parameters', {})
            result = None
            if param_dict:
                try:
                    result = self._reconstruct_result(energy, intensity, param_dict)
                except Exception:
                    pass

            self.file_state_cache[file_key] = {
                'fit_peaks': fit_peaks,
                'quant_peaks': quant_peaks,
                'result': result,
                'energy': energy,
                'intensity': intensity,
            }
            restored_count += 1

        # Navigate to saved index
        saved_index = session.get('current_index', 0)
        if self.spectra_files:
            self.current_index = min(saved_index, len(self.spectra_files) - 1)
            self.load_current_spectrum()

        messagebox.showinfo("Loaded",
                            f"Session loaded: {len(valid_paths)} spectra, "
                            f"{restored_count} with restored fits.")

    def on_height1_change(self, *args):
        """Lock step 2 height to step 1 if locked."""
        if self.lock_heights.get() and not getattr(self, '_updating_heights', False):
            self._updating_heights = True
            self.baseline_params['arc2_height'].set(self.baseline_params['arc1_height'].get())
            self._updating_heights = False
            self.update_plot()

    def on_height2_change(self, *args):
        """Lock step 1 height to step 2 if locked (bidirectional)."""
        if self.lock_heights.get() and not getattr(self, '_updating_heights', False):
            self._updating_heights = True
            self.baseline_params['arc1_height'].set(self.baseline_params['arc2_height'].get())
            self._updating_heights = False
            self.update_plot()

    def on_width1_change(self, *args):
        """Lock step 2 width to step 1 if locked."""
        if self.lock_widths.get() and not getattr(self, '_updating_widths', False):
            self._updating_widths = True
            self.baseline_params['arc2_width'].set(self.baseline_params['arc1_width'].get())
            self._updating_widths = False
            self.update_plot()

    def on_width2_change(self, *args):
        """Lock step 1 width to step 2 if locked (bidirectional)."""
        if self.lock_widths.get() and not getattr(self, '_updating_widths', False):
            self._updating_widths = True
            self.baseline_params['arc1_width'].set(self.baseline_params['arc2_width'].get())
            self._updating_widths = False
            self.update_plot()

    def setup_gui(self):
        """Setup the GUI layout with notebook tabs."""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Tab 1: Fit Viewer
        self.fit_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.fit_tab, text="Fit Viewer")
        self.setup_fit_tab(self.fit_tab)

        # Tab 2: Analysis & Comparison (merged statistics and comparison)
        self.analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.analysis_tab, text="Analysis & Comparison")
        self.setup_analysis_tab(self.analysis_tab)

        # Tab 3: Peak Centers Configuration
        self.peak_config_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.peak_config_tab, text="Peak Centers")
        self.setup_peak_config_tab(self.peak_config_tab)

    def setup_fit_tab(self, parent):
        """Setup fit viewer tab."""
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True)

        # Left panel container with scrollbar (wider to show all columns)
        left_container = ttk.Frame(main_container, width=320)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        left_container.pack_propagate(False)

        # Create canvas and scrollbar for left panel
        left_canvas = tk.Canvas(left_container)
        left_scrollbar = ttk.Scrollbar(left_container, orient=tk.VERTICAL, command=left_canvas.yview)
        left_panel = ttk.Frame(left_canvas)

        left_panel.bind("<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))

        left_canvas.create_window((0, 0), window=left_panel, anchor="nw", width=300)
        left_canvas.configure(yscrollcommand=left_scrollbar.set)

        # Bind mouse wheel to scroll the left panel only when mouse is over it
        def on_left_scroll(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        left_canvas.bind("<Enter>", lambda e: left_canvas.bind_all("<MouseWheel>", on_left_scroll))
        left_canvas.bind("<Leave>", lambda e: left_canvas.unbind_all("<MouseWheel>"))

        left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Middle panel (sample list)
        middle_panel = ttk.Frame(main_container, width=200)
        middle_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        middle_panel.pack_propagate(False)

        # Right panel (plots)
        right_panel = ttk.Frame(main_container)
        right_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.setup_controls(left_panel)
        self.setup_sample_list(middle_panel)
        self.setup_plot_area(right_panel)

    def setup_controls(self, parent):
        """Setup control panel."""
        # File controls
        file_frame = ttk.LabelFrame(parent, text="File Control", padding=5)
        file_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Button(file_frame, text="Add Files",
                  command=self.add_files).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="Remove Selected",
                  command=self.remove_file).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="Clear All",
                  command=self.clear_files).pack(fill=tk.X, pady=2)

        # Fit statistics
        stats_frame = ttk.LabelFrame(parent, text="Fit Quality", padding=5)
        stats_frame.pack(fill=tk.X, pady=(0, 5))
        self.stats_text = tk.Text(stats_frame, height=10, width=28, font=('Courier', 8))
        self.stats_text.pack(fill=tk.BOTH, expand=True)
        self.stats_text.config(state=tk.DISABLED)

        # Display options
        display_frame = ttk.LabelFrame(parent, text="Display Options", padding=5)
        display_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Checkbutton(display_frame, text="Show Baseline",
                       variable=self.show_baseline,
                       command=self.update_plot).pack(anchor=tk.W)
        ttk.Checkbutton(display_frame, text="Show Peaks",
                       variable=self.show_peaks,
                       command=self.update_plot).pack(anchor=tk.W)

        # Peak parameters with checkboxes
        self.setup_peak_parameters(parent)

        # Baseline adjustment (collapsible)
        baseline_frame = ttk.LabelFrame(parent, text="Baseline Control (Advanced)", padding=5)
        baseline_frame.pack(fill=tk.X, pady=(0, 5))

        # Expand/collapse button
        self.baseline_expanded = tk.BooleanVar(value=True)  # Expanded by default
        expand_btn = ttk.Checkbutton(baseline_frame, text="Show Advanced Baseline Controls",
                                     variable=self.baseline_expanded,
                                     command=self.toggle_baseline_visibility)
        expand_btn.pack(anchor=tk.W)

        # Container for baseline controls (visible by default now)
        self.baseline_container = ttk.Frame(baseline_frame)
        self.baseline_container.pack(fill=tk.BOTH, expand=True, pady=(5,0))  # Pack immediately since expanded

        ttk.Checkbutton(self.baseline_container, text="Manual Adjustment",
                       variable=self.manual_baseline,
                       command=self.toggle_manual_baseline).pack(anchor=tk.W, pady=(5,0))

        # Lock options
        lock_frame = ttk.Frame(self.baseline_container)
        lock_frame.pack(fill=tk.X, pady=(5, 5))
        ttk.Checkbutton(lock_frame, text="Lock Heights",
                       variable=self.lock_heights).pack(side=tk.LEFT)
        ttk.Checkbutton(lock_frame, text="Lock Widths",
                       variable=self.lock_widths).pack(side=tk.LEFT, padx=(10, 0))

        # Scrollable frame for baseline parameters (taller to show all 6 controls)
        canvas = tk.Canvas(self.baseline_container, height=200)
        scrollbar = ttk.Scrollbar(self.baseline_container, orient="vertical", command=canvas.yview)
        self.baseline_controls_frame = ttk.Frame(canvas)

        self.baseline_controls_frame.bind("<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        canvas.create_window((0, 0), window=self.baseline_controls_frame, anchor="nw", width=280)
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.setup_baseline_controls()

        # Actions
        action_frame = ttk.LabelFrame(parent, text="Actions", padding=3)
        action_frame.pack(fill=tk.X)
        ttk.Button(action_frame, text="Refit",
                  command=self.refit_spectrum).pack(fill=tk.X, pady=2)
        ttk.Button(action_frame, text="Box Plot Test",
                  command=self.run_box_plot_test).pack(fill=tk.X, pady=2)
        ttk.Button(action_frame, text="Model Selection",
                  command=self.run_model_selection).pack(fill=tk.X, pady=2)
        ttk.Button(action_frame, text="Monte Carlo",
                  command=self.run_monte_carlo).pack(fill=tk.X, pady=2)

        # Auto-reduce frame
        reduce_frame = ttk.LabelFrame(parent, text="Auto-Reduce", padding=3)
        reduce_frame.pack(fill=tk.X, pady=(5, 0))

        threshold_row = ttk.Frame(reduce_frame)
        threshold_row.pack(fill=tk.X, pady=2)
        ttk.Label(threshold_row, text="\u0394AIC:").pack(side=tk.LEFT)
        ttk.Spinbox(threshold_row, from_=0, to=10, increment=0.5, width=5,
                    textvariable=self.auto_reduce_threshold).pack(side=tk.LEFT, padx=5)

        ttk.Button(reduce_frame, text="Auto-Reduce",
                  command=self.run_auto_reduce).pack(fill=tk.X, pady=2)
        ttk.Button(reduce_frame, text="Batch Auto-Reduce",
                  command=self.run_batch_auto_reduce).pack(fill=tk.X, pady=2)

        # Export options
        export_frame = ttk.LabelFrame(parent, text="Export", padding=3)
        export_frame.pack(fill=tk.X, pady=(5, 0))
        ttk.Button(export_frame, text="Export Plot",
                  command=self.export_plot).pack(fill=tk.X, pady=1)
        ttk.Button(export_frame, text="Export Fit Data",
                  command=self.export_fit_data).pack(fill=tk.X, pady=1)
        ttk.Button(export_frame, text="Export All (Complete)",
                  command=self.export_all_data).pack(fill=tk.X, pady=1)
        ttk.Button(export_frame, text="Batch Export All Samples",
                  command=self.batch_export_all).pack(fill=tk.X, pady=1)

    def setup_peak_parameters(self, parent):
        """Setup peak parameters table with resizable columns using Treeview."""
        frame = ttk.LabelFrame(parent, text="Peak Parameters", padding=3)
        frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))

        # Checkbox frame above the treeview
        cb_frame = ttk.Frame(frame)
        cb_frame.pack(fill=tk.X, pady=(0, 3))

        ttk.Label(cb_frame, text="Toggle:", font=('Arial', 8)).pack(side=tk.LEFT)
        self.select_all_fit = tk.BooleanVar(value=True)
        self.select_all_quant = tk.BooleanVar(value=True)
        ttk.Checkbutton(cb_frame, text="All Fit", variable=self.select_all_fit,
                       command=self.toggle_all_fit).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(cb_frame, text="All Quant", variable=self.select_all_quant,
                       command=self.toggle_all_quant).pack(side=tk.LEFT)

        # Create Treeview with columns
        columns = ('peak', 'center', 'fwhm', 'area')
        self.peak_tree = ttk.Treeview(frame, columns=columns, show='tree headings',
                                       height=8, selectmode='none')

        # Define column headings and widths (user can resize by dragging)
        self.peak_tree.heading('#0', text='F/Q', anchor=tk.W)
        self.peak_tree.heading('peak', text='Peak', anchor=tk.W)
        self.peak_tree.heading('center', text='Center', anchor=tk.CENTER)
        self.peak_tree.heading('fwhm', text='FWHM', anchor=tk.CENTER)
        self.peak_tree.heading('area', text='%Area', anchor=tk.CENTER)

        # Set initial column widths (resizable by user)
        self.peak_tree.column('#0', width=50, minwidth=40, stretch=False)
        self.peak_tree.column('peak', width=70, minwidth=50)
        self.peak_tree.column('center', width=55, minwidth=40)
        self.peak_tree.column('fwhm', width=50, minwidth=40)
        self.peak_tree.column('area', width=50, minwidth=40)

        # Scrollbar
        tree_scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.peak_tree.yview)
        self.peak_tree.configure(yscrollcommand=tree_scroll.set)

        self.peak_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Bind click events for toggling checkboxes
        self.peak_tree.bind('<Button-1>', self.on_peak_tree_click)

        # Store tree item IDs
        self.peak_tree_items = {}
        self.peak_labels = {}
        self.peak_checkboxes = {}

    def toggle_all_fit(self):
        """Toggle all fit checkboxes."""
        state = self.select_all_fit.get()
        for name in self.peak_names:
            self.peak_in_fit[name].set(state)
        self.on_fit_checkbox_changed()
        self.update_peak_tree_display()

    def toggle_all_quant(self):
        """Toggle all quantification checkboxes."""
        state = self.select_all_quant.get()
        for name in self.peak_names:
            self.peak_in_quant[name].set(state)
        self.on_quant_checkbox_changed()
        self.update_peak_tree_display()

    def on_peak_tree_click(self, event):
        """Handle clicks on the peak tree to toggle F/Q checkboxes."""
        region = self.peak_tree.identify_region(event.x, event.y)
        if region == 'tree':  # Clicked on the tree column (#0)
            item = self.peak_tree.identify_row(event.y)
            if item and item in self.peak_tree_items.values():
                # Find which peak this is
                for peak_name, item_id in self.peak_tree_items.items():
                    if item_id == item:
                        # Determine if click was on F or Q based on x position
                        col_width = self.peak_tree.column('#0', 'width')
                        if event.x < col_width / 2:
                            # Toggle Fit
                            current = self.peak_in_fit[peak_name].get()
                            self.peak_in_fit[peak_name].set(not current)
                            self.on_fit_checkbox_changed()
                        else:
                            # Toggle Quant
                            current = self.peak_in_quant[peak_name].get()
                            self.peak_in_quant[peak_name].set(not current)
                            self.on_quant_checkbox_changed()
                        self.update_peak_tree_display()
                        break

    def update_peak_tree_display(self):
        """Update the checkbox display in the tree."""
        for peak_name, item_id in self.peak_tree_items.items():
            fit_char = 'Y' if self.peak_in_fit[peak_name].get() else '-'
            quant_char = 'Y' if self.peak_in_quant[peak_name].get() else '-'
            self.peak_tree.item(item_id, text=f"{fit_char} {quant_char}")

    def toggle_baseline_visibility(self):
        """Toggle visibility of baseline controls."""
        if self.baseline_expanded.get():
            self.baseline_container.pack(fill=tk.BOTH, expand=True, pady=(5,0))
        else:
            self.baseline_container.pack_forget()

    def setup_baseline_controls(self):
        """Setup baseline parameter controls."""
        params = [
            ('Step 1 Center', 'arc1_center', 290.2, 290.8, 0.1),
            ('Step 1 Height', 'arc1_height', 0.0, 2.0, 0.05),
            ('Step 1 Width', 'arc1_width', 0.1, 1.0, 0.05),
            ('Step 2 Center', 'arc2_center', 291.5, 292.5, 0.1),
            ('Step 2 Height', 'arc2_height', 0.0, 2.0, 0.05),
            ('Step 2 Width', 'arc2_width', 0.1, 1.0, 0.05)
        ]

        self.baseline_scales = {}

        for i, (label, param, min_val, max_val, resolution) in enumerate(params):
            frame = ttk.Frame(self.baseline_controls_frame)
            frame.grid(row=i, column=0, sticky=tk.EW, pady=2)

            ttk.Label(frame, text=label, width=12, font=('Arial', 7)).pack(side=tk.LEFT)

            value_label = ttk.Label(frame, text=f"{self.baseline_params[param].get():.3f}",
                                   width=6, font=('Arial', 7))
            value_label.pack(side=tk.RIGHT)

            scale = ttk.Scale(frame, from_=min_val, to=max_val,
                            variable=self.baseline_params[param],
                            orient=tk.HORIZONTAL,
                            command=lambda v, p=param, l=value_label: self.update_baseline_value(p, v, l))
            scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

            self.baseline_scales[param] = (scale, value_label)
            scale.state(['disabled'])

    def setup_sample_list(self, parent):
        """Setup sample list with scrollbar."""
        list_frame = ttk.LabelFrame(parent, text="Samples", padding=5)
        list_frame.pack(fill=tk.BOTH, expand=True)

        # Search box
        search_frame = ttk.Frame(list_frame)
        search_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.search_var = tk.StringVar()
        self.search_var.trace('w', self.filter_sample_list)
        search_entry = ttk.Entry(search_frame, textvariable=self.search_var)
        search_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        # Listbox with scrollbar (multi-select enabled)
        listbox_frame = ttk.Frame(list_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(listbox_frame, orient=tk.VERTICAL)
        self.sample_listbox = tk.Listbox(listbox_frame,
                                        yscrollcommand=scrollbar.set,
                                        font=('Courier', 9),
                                        selectmode=tk.EXTENDED)  # Single click selects, Ctrl+click adds, Shift+click ranges
        scrollbar.config(command=self.sample_listbox.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.sample_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.sample_listbox.bind('<<ListboxSelect>>', self.on_sample_select)

        # Navigation
        nav_frame = ttk.Frame(list_frame)
        nav_frame.pack(fill=tk.X, pady=(5, 0))
        ttk.Button(nav_frame, text="◄ Prev",
                  command=self.previous_spectrum).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(nav_frame, text="Next ►",
                  command=self.next_spectrum).pack(side=tk.LEFT, expand=True, fill=tk.X)

    def setup_plot_area(self, parent):
        """Setup matplotlib plot area."""
        self.fig = Figure(figsize=(12, 8), dpi=100)
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, parent)
        toolbar.update()

    def setup_analysis_tab(self, parent):
        """Setup combined analysis and comparison tab with zoomable plots."""
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Use PanedWindow for resizable sections
        paned = ttk.PanedWindow(main_container, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        # Left panel: Sample selection and controls (with scrollable content)
        left_panel = ttk.Frame(paned, width=280)
        paned.add(left_panel, weight=0)

        # Create a canvas with scrollbar for the left panel content
        left_canvas = tk.Canvas(left_panel, width=270, highlightthickness=0)
        left_scrollbar = ttk.Scrollbar(left_panel, orient=tk.VERTICAL, command=left_canvas.yview)
        left_scrollable = ttk.Frame(left_canvas)

        left_scrollable.bind("<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))

        left_canvas.create_window((0, 0), window=left_scrollable, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar.set)

        # Enable mouse wheel scrolling on the left panel
        def on_left_panel_scroll(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        left_canvas.bind_all("<MouseWheel>", on_left_panel_scroll)
        left_canvas.bind("<Enter>", lambda e: left_canvas.bind_all("<MouseWheel>", on_left_panel_scroll))
        left_canvas.bind("<Leave>", lambda e: left_canvas.unbind_all("<MouseWheel>"))

        left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Sample selection
        sample_frame = ttk.LabelFrame(left_scrollable, text="Select Samples", padding=5)
        sample_frame.pack(fill=tk.X, pady=(0, 5))

        # Search box
        search_frame = ttk.Frame(sample_frame)
        search_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.stats_search_var = tk.StringVar()
        self.stats_search_var.trace('w', self.filter_stats_samples)
        ttk.Entry(search_frame, textvariable=self.stats_search_var).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        # Sample listbox (multi-select with EXTENDED mode)
        list_frame = ttk.Frame(sample_frame)
        list_frame.pack(fill=tk.BOTH, expand=True)

        list_scroll = ttk.Scrollbar(list_frame, orient=tk.VERTICAL)
        self.stats_samples_listbox = tk.Listbox(list_frame, height=10,
                                                 selectmode=tk.EXTENDED,
                                                 yscrollcommand=list_scroll.set,
                                                 font=('Courier', 9))
        list_scroll.config(command=self.stats_samples_listbox.yview)
        list_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.stats_samples_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Also create available_samples_listbox as alias for compatibility
        self.available_samples_listbox = self.stats_samples_listbox

        # Selected samples display
        selected_frame = ttk.LabelFrame(sample_frame, text="Selected for Analysis", padding=5)
        selected_frame.pack(fill=tk.X, pady=(5, 0))

        self.comparison_listbox = tk.Listbox(selected_frame, height=4, font=('Courier', 8))
        self.comparison_listbox.pack(fill=tk.X)

        # Action buttons
        btn_frame = ttk.Frame(sample_frame)
        btn_frame.pack(fill=tk.X, pady=(5, 0))
        ttk.Button(btn_frame, text="Add Selected",
                  command=self.add_selected_to_comparison).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Remove",
                  command=self.remove_from_comparison).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Clear",
                  command=self.clear_comparison).pack(side=tk.LEFT, padx=2)

        # Plot type selection
        plot_frame = ttk.LabelFrame(left_scrollable, text="Generate Plots", padding=5)
        plot_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Button(plot_frame, text="Statistics",
                  command=self.generate_multi_sample_stats).pack(fill=tk.X, pady=2)
        ttk.Button(plot_frame, text="Comparison Charts",
                  command=self.generate_comparison).pack(fill=tk.X, pady=2)
        ttk.Button(plot_frame, text="Spectral Overlay",
                  command=self.generate_spectral_overlay).pack(fill=tk.X, pady=2)

        # Spectral overlay options
        overlay_frame = ttk.LabelFrame(left_scrollable, text="Overlay Options", padding=5)
        overlay_frame.pack(fill=tk.X, pady=(0, 5))

        # Y-axis offset control (increased range: 0-2.0)
        ttk.Label(overlay_frame, text="Y-Axis Offset:").pack(anchor=tk.W)
        offset_container = ttk.Frame(overlay_frame)
        offset_container.pack(fill=tk.X)
        self.y_offset_var = tk.DoubleVar(value=0.5)
        self.offset_scale = ttk.Scale(offset_container, from_=0.0, to=2.0,
                                 variable=self.y_offset_var,
                                 orient=tk.HORIZONTAL)
        self.offset_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.offset_value_label = ttk.Label(offset_container, text="0.50", width=5)
        self.offset_value_label.pack(side=tk.LEFT)
        self.y_offset_var.trace('w', lambda *args: self.offset_value_label.config(
            text=f"{self.y_offset_var.get():.2f}"))

        # Show peaks checkbox
        self.show_overlay_peaks = tk.BooleanVar(value=True)
        ttk.Checkbutton(overlay_frame, text="Show Individual Peaks",
                       variable=self.show_overlay_peaks).pack(anchor=tk.W)

        # Summary statistics display
        stats_frame = ttk.LabelFrame(left_scrollable, text="Summary (Average)", padding=5)
        stats_frame.pack(fill=tk.X, pady=(0, 5))

        metrics = [
            ('Aromatic:', 'aromaticity'),
            ('Aliphatic:', 'aliphatic'),
            ('Carboxyl:', 'carboxyl'),
            ('O-alkyl:', 'oalkyl'),
            ('Carbonate:', 'carbonate'),
        ]

        self.stat_labels = {}
        for label, key in metrics:
            row_frame = ttk.Frame(stats_frame)
            row_frame.pack(fill=tk.X)
            ttk.Label(row_frame, text=label, font=('Arial', 9)).pack(side=tk.LEFT)
            self.stat_labels[key] = ttk.Label(row_frame, text="--", font=('Arial', 9, 'bold'))
            self.stat_labels[key].pack(side=tk.RIGHT)

        self.stat_labels['total'] = ttk.Label(stats_frame, text="0 samples", font=('Arial', 8, 'italic'))
        self.stat_labels['total'].pack(anchor=tk.W, pady=(5, 0))

        # Status label
        self.stats_status_label = self.stat_labels['total']

        # Export buttons
        export_frame = ttk.LabelFrame(left_scrollable, text="Export", padding=5)
        export_frame.pack(fill=tk.X)

        ttk.Button(export_frame, text="Export Figure",
                  command=self.export_analysis_figure).pack(fill=tk.X, pady=2)
        ttk.Button(export_frame, text="Export Data (CSV)",
                  command=self.export_comparison_data).pack(fill=tk.X, pady=2)
        ttk.Button(export_frame, text="Export Complete Analysis",
                  command=self.export_complete_analysis).pack(fill=tk.X, pady=2)

        # Right panel: Zoomable plot area
        right_panel = ttk.Frame(paned)
        paned.add(right_panel, weight=1)

        # Create figure with higher DPI for publication quality
        self.analysis_fig = Figure(figsize=(12, 10), dpi=100)
        self.analysis_canvas = FigureCanvasTkAgg(self.analysis_fig, master=right_panel)
        self.analysis_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add navigation toolbar for zoom/pan
        toolbar_frame = ttk.Frame(right_panel)
        toolbar_frame.pack(fill=tk.X)
        self.analysis_toolbar = NavigationToolbar2Tk(self.analysis_canvas, toolbar_frame)
        self.analysis_toolbar.update()

        # Bind mouse wheel for zoom
        self.analysis_canvas.get_tk_widget().bind('<MouseWheel>', self.on_analysis_scroll)
        self.analysis_canvas.get_tk_widget().bind('<Button-4>', self.on_analysis_scroll)  # Linux scroll up
        self.analysis_canvas.get_tk_widget().bind('<Button-5>', self.on_analysis_scroll)  # Linux scroll down

        # Store current zoom level
        self._analysis_zoom_level = 1.0

        # Create aliases for compatibility with existing methods
        self.stats_fig = self.analysis_fig
        self.stats_canvas = self.analysis_canvas
        self.comparison_fig = self.analysis_fig
        self.comparison_canvas = self.analysis_canvas

    def on_analysis_scroll(self, event):
        """Handle mouse wheel scroll for zooming the analysis plot."""
        # Get the current axes
        if not self.analysis_fig.axes:
            return

        # Find which axes the mouse is over
        ax = None
        for a in self.analysis_fig.axes:
            if a.in_axes(event):
                ax = a
                break

        # If not over any axes, zoom all of them
        if ax is None:
            axes_to_zoom = [a for a in self.analysis_fig.axes if a.get_visible()]
        else:
            axes_to_zoom = [ax]

        # Determine zoom direction
        if event.num == 5 or (hasattr(event, 'delta') and event.delta < 0):
            # Scroll down - zoom out
            scale_factor = 1.15
        else:
            # Scroll up - zoom in
            scale_factor = 0.87

        for ax in axes_to_zoom:
            # Skip axes that are turned off (like tables)
            if not ax.axison:
                continue

            # Get current limits
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            # Calculate center
            x_center = (xlim[0] + xlim[1]) / 2
            y_center = (ylim[0] + ylim[1]) / 2

            # Calculate new limits centered on center
            new_width = (xlim[1] - xlim[0]) * scale_factor
            new_height = (ylim[1] - ylim[0]) * scale_factor

            new_xlim = [x_center - new_width / 2, x_center + new_width / 2]
            new_ylim = [y_center - new_height / 2, y_center + new_height / 2]

            ax.set_xlim(new_xlim)
            ax.set_ylim(new_ylim)

        self.analysis_canvas.draw_idle()

    def export_analysis_figure(self):
        """Export analysis figure as publication-quality image."""
        filepath = filedialog.asksaveasfilename(
            title="Export Analysis Figure",
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("TIFF files", "*.tiff")
            ]
        )
        if filepath:
            self.analysis_fig.savefig(filepath, dpi=300, bbox_inches='tight',
                                      facecolor='white', edgecolor='none')
            messagebox.showinfo("Success", f"Figure exported to:\n{filepath}")

    def setup_peak_config_tab(self, parent):
        """Setup peak centers configuration tab."""
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Instructions
        instr_frame = ttk.Frame(main_container)
        instr_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(instr_frame, text="Peak Center Configuration",
                 font=('Arial', 12, 'bold')).pack(anchor=tk.W)
        ttk.Label(instr_frame,
                 text="Adjust peak center positions and allowed ranges for fitting. Changes apply to new fits only.",
                 wraplength=700).pack(anchor=tk.W, pady=(5, 0))

        # Enable custom centers
        ttk.Checkbutton(instr_frame, text="Use Custom Peak Centers",
                       variable=self.custom_peak_centers).pack(anchor=tk.W, pady=(10, 0))

        # Peak configuration table
        table_frame = ttk.LabelFrame(main_container, text="Peak Centers and Ranges", padding=10)
        table_frame.pack(fill=tk.BOTH, expand=True)

        # Scrollable frame
        canvas = tk.Canvas(table_frame)
        scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=canvas.yview)
        peaks_frame = ttk.Frame(canvas)

        peaks_frame.bind("<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        canvas.create_window((0, 0), window=peaks_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Header row
        headers = ["Peak", "Default Center (eV)", "Center (eV)", "±Range (eV)", "Min (eV)", "Max (eV)"]
        for col, header in enumerate(headers):
            ttk.Label(peaks_frame, text=header, font=('Arial', 9, 'bold')).grid(
                row=0, column=col, padx=5, pady=5, sticky=tk.W)

        # Peak rows
        self.peak_center_entries = {}
        self.peak_range_entries = {}
        self.peak_min_labels = {}
        self.peak_max_labels = {}

        for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
            row = i

            # Peak name
            ttk.Label(peaks_frame, text=display_name,
                     foreground=matplotlib.colors.rgb2hex(self.peak_colors[i-1])).grid(
                row=row, column=0, padx=5, pady=5, sticky=tk.W)

            # Default center
            default = self.default_peak_centers[name]
            ttk.Label(peaks_frame, text=f"{default:.2f}").grid(
                row=row, column=1, padx=5, pady=5)

            # Current center (editable)
            center_entry = ttk.Entry(peaks_frame, width=8, textvariable=self.peak_center_vars[name])
            center_entry.grid(row=row, column=2, padx=5, pady=5)
            self.peak_center_entries[name] = center_entry

            # Range (editable)
            range_entry = ttk.Entry(peaks_frame, width=8, textvariable=self.peak_range_vars[name])
            range_entry.grid(row=row, column=3, padx=5, pady=5)
            self.peak_range_entries[name] = range_entry

            # Min (calculated)
            min_label = ttk.Label(peaks_frame, text=f"{default - self.peak_range_vars[name].get():.2f}")
            min_label.grid(row=row, column=4, padx=5, pady=5)
            self.peak_min_labels[name] = min_label

            # Max (calculated)
            max_label = ttk.Label(peaks_frame, text=f"{default + self.peak_range_vars[name].get():.2f}")
            max_label.grid(row=row, column=5, padx=5, pady=5)
            self.peak_max_labels[name] = max_label

            # Trace changes to update min/max
            self.peak_center_vars[name].trace('w', lambda *args, n=name: self.update_peak_range_labels(n))
            self.peak_range_vars[name].trace('w', lambda *args, n=name: self.update_peak_range_labels(n))

        # FWHM configuration section
        fwhm_frame = ttk.LabelFrame(main_container, text="FWHM (Full Width at Half Maximum) Settings", padding=10)
        fwhm_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Label(fwhm_frame, text="Configure FWHM values and ranges for peak groups",
                 font=('Arial', 9)).pack(anchor=tk.W, pady=(0, 10))

        # FWHM table
        fwhm_table = ttk.Frame(fwhm_frame)
        fwhm_table.pack(fill=tk.X)

        # Headers
        headers_fwhm = ["Peak Group", "Default FWHM", "FWHM (eV)", "Min FWHM", "Max FWHM"]
        for col, header in enumerate(headers_fwhm):
            ttk.Label(fwhm_table, text=header, font=('Arial', 9, 'bold')).grid(
                row=0, column=col, padx=5, pady=5, sticky=tk.W)

        # Main FWHM (peaks 1-6)
        ttk.Label(fwhm_table, text="Main Peaks (1-6)").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        ttk.Label(fwhm_table, text="1.5").grid(row=1, column=1, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_main_var).grid(row=1, column=2, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_main_min_var).grid(row=1, column=3, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_main_max_var).grid(row=1, column=4, padx=5, pady=5)

        # Carbonate FWHM (peak 7)
        ttk.Label(fwhm_table, text="Carbonate (7)").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        ttk.Label(fwhm_table, text="0.8").grid(row=2, column=1, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_carb_var).grid(row=2, column=2, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_carb_min_var).grid(row=2, column=3, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_carb_max_var).grid(row=2, column=4, padx=5, pady=5)

        # Sigma* FWHM (peak 8)
        ttk.Label(fwhm_table, text="Sigma* (8)").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
        ttk.Label(fwhm_table, text="0.6").grid(row=3, column=1, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_sigma_var).grid(row=3, column=2, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_sigma_min_var).grid(row=3, column=3, padx=5, pady=5)
        ttk.Entry(fwhm_table, width=8, textvariable=self.fwhm_sigma_max_var).grid(row=3, column=4, padx=5, pady=5)

        # Reset button
        btn_frame = ttk.Frame(main_container)
        btn_frame.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(btn_frame, text="Reset to Defaults",
                  command=self.reset_peak_centers).pack(side=tk.LEFT, padx=5)

    def update_peak_range_labels(self, peak_name):
        """Update min/max labels when center or range changes."""
        try:
            center = self.peak_center_vars[peak_name].get()
            range_val = self.peak_range_vars[peak_name].get()
        except (tk.TclError, ValueError):
            return

        self.peak_min_labels[peak_name].config(text=f"{center - range_val:.2f}")
        self.peak_max_labels[peak_name].config(text=f"{center + range_val:.2f}")

    def reset_peak_centers(self):
        """Reset all peak centers and FWHM to defaults."""
        for name in self.peak_names:
            self.peak_center_vars[name].set(self.default_peak_centers[name])
            if name == 'Sigma_star':
                self.peak_range_vars[name].set(2.7)
            else:
                self.peak_range_vars[name].set(0.3)

        # Reset FWHM values
        self.fwhm_main_var.set(1.5)
        self.fwhm_main_min_var.set(0.8)
        self.fwhm_main_max_var.set(2.0)

        self.fwhm_carb_var.set(0.8)
        self.fwhm_carb_min_var.set(0.5)
        self.fwhm_carb_max_var.set(1.0)

        self.fwhm_sigma_var.set(0.6)
        self.fwhm_sigma_min_var.set(0.4)
        self.fwhm_sigma_max_var.set(2.0)

        messagebox.showinfo("Success", "Peak centers and FWHM reset to default values")

    def toggle_manual_baseline(self):
        """Enable/disable manual baseline controls."""
        state = 'normal' if self.manual_baseline.get() else 'disabled'

        for scale, _ in self.baseline_scales.values():
            scale.state(['!disabled'] if state == 'normal' else ['disabled'])

        if self.manual_baseline.get() and self.current_result:
            for param in self.baseline_params.keys():
                value = self.current_result.params[param].value
                self.baseline_params[param].set(value)

        self.update_plot()

    def update_baseline_value(self, param, value, label):
        """Update baseline parameter value label and plot."""
        label.config(text=f"{float(value):.3f}")
        if self.manual_baseline.get():
            self.update_plot()

    def load_spectra_list(self):
        """Load list of spectra files."""
        if not self.spectra_dir.exists():
            messagebox.showerror("Error", f"Directory not found:\n{self.spectra_dir}")
            return

        self.spectra_files = sorted(list(self.spectra_dir.glob('*.csv')))
        self.file_state_cache.clear()

        if not self.spectra_files:
            messagebox.showwarning("Warning", f"No CSV files found")
            return

        self.update_sample_listbox()

        # Populate available samples listbox for comparison
        if hasattr(self, 'available_samples_listbox'):
            self.available_samples_listbox.delete(0, tk.END)
            for file_path in self.spectra_files:
                self.available_samples_listbox.insert(tk.END, file_path.stem)

        # Populate stats samples listbox
        if hasattr(self, 'stats_samples_listbox'):
            self.stats_samples_listbox.delete(0, tk.END)
            for file_path in self.spectra_files:
                self.stats_samples_listbox.insert(tk.END, file_path.stem)

    def update_sample_listbox(self):
        """Update sample listbox with optional filtering."""
        self.sample_listbox.delete(0, tk.END)
        search_term = self.search_var.get().lower()

        for i, file_path in enumerate(self.spectra_files):
            name = file_path.stem
            if search_term == '' or search_term in name.lower():
                self.sample_listbox.insert(tk.END, name)

                if i == self.current_index:
                    self.sample_listbox.selection_set(self.sample_listbox.size() - 1)
                    self.sample_listbox.see(self.sample_listbox.size() - 1)

    def filter_sample_list(self, *args):
        """Filter sample list based on search."""
        self.update_sample_listbox()

    def on_sample_select(self, event):
        """Handle sample selection from listbox (supports multi-select)."""
        selection = self.sample_listbox.curselection()
        if not selection:
            return

        if len(selection) == 1:
            # Single selection - load normally
            selected_name = self.sample_listbox.get(selection[0])
            for i, file_path in enumerate(self.spectra_files):
                if file_path.stem == selected_name:
                    self._save_file_state()
                    self.current_index = i
                    self.load_current_spectrum()
                    break
        else:
            # Multiple selections - overlay them
            self.plot_multiple_samples(selection)

    def add_files(self):
        """Add individual CSV files to the list."""
        files = filedialog.askopenfilenames(
            title="Select CSV Files",
            initialdir=self.spectra_dir if hasattr(self, 'spectra_dir') else Path.home(),
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if files:
            # Convert to Path objects and add to list
            for file_path in files:
                file_path = Path(file_path)
                if file_path not in self.spectra_files:
                    self.spectra_files.append(file_path)

            # Sort files by name
            self.spectra_files.sort()
            self.file_list_modified = True

            # Update the listboxes
            self.update_sample_listbox()
            if hasattr(self, 'available_samples_listbox'):
                self.available_samples_listbox.delete(0, tk.END)
                for file_path in self.spectra_files:
                    self.available_samples_listbox.insert(tk.END, file_path.stem)

            # Update stats samples listbox
            if hasattr(self, 'stats_samples_listbox'):
                self.stats_samples_listbox.delete(0, tk.END)
                for file_path in self.spectra_files:
                    self.stats_samples_listbox.insert(tk.END, file_path.stem)

            # If this is the first file(s), load the first one
            if len(self.spectra_files) > 0 and self.current_result is None:
                self.current_index = 0
                self.load_current_spectrum()

    def remove_file(self):
        """Remove selected file from the list."""
        selection = self.sample_listbox.curselection()
        if not selection:
            messagebox.showwarning("Warning", "No file selected")
            return

        selected_name = self.sample_listbox.get(selection[0])

        # Find and remove the file
        for i, file_path in enumerate(self.spectra_files):
            if file_path.stem == selected_name:
                del self.spectra_files[i]
                self.file_list_modified = True

                # Update the listboxes
                self.update_sample_listbox()
                if hasattr(self, 'available_samples_listbox'):
                    self.available_samples_listbox.delete(0, tk.END)
                    for file_path in self.spectra_files:
                        self.available_samples_listbox.insert(tk.END, file_path.stem)

                # Adjust current index if needed
                if len(self.spectra_files) == 0:
                    self.current_index = 0
                    self.current_result = None
                    self.current_energy = None
                    self.current_intensity = None
                    self.ax1.clear()
                    self.ax2.clear()
                    self.canvas.draw()
                else:
                    if self.current_index >= len(self.spectra_files):
                        self.current_index = len(self.spectra_files) - 1
                    self.load_current_spectrum()
                break

    def clear_files(self):
        """Clear all files from the list."""
        if not self.spectra_files:
            return

        if messagebox.askyesno("Confirm", "Remove all files from the list?"):
            self.spectra_files = []
            self.file_state_cache.clear()
            self.file_list_modified = True
            self.current_index = 0
            self.current_result = None
            self.current_energy = None
            self.current_intensity = None

            # Update listboxes
            self.sample_listbox.delete(0, tk.END)
            if hasattr(self, 'available_samples_listbox'):
                self.available_samples_listbox.delete(0, tk.END)

            # Clear plots
            self.ax1.clear()
            self.ax2.clear()
            self.canvas.draw()

    def fit_spectrum_with_custom_centers(self, energy, intensity, manual_baseline=False,
                                          baseline_params=None, lock_heights=False, lock_widths=False):
        """Fit spectrum with custom peak centers from GUI."""
        from lmfit import Parameters, Model
        from c1s_fitter_optimized import estimate_baseline_parameters, total_model, double_arctangent, gaussian

        params = Parameters()

        # Estimate data characteristics for baseline
        est = estimate_baseline_parameters(energy, intensity)

        # Baseline - use manual values if specified
        if manual_baseline and baseline_params:
            arc1_center = baseline_params['arc1_center'].get()
            arc1_height = baseline_params['arc1_height'].get()
            arc1_width = baseline_params['arc1_width'].get()
            arc2_center = baseline_params['arc2_center'].get()
            arc2_height = baseline_params['arc2_height'].get()
            arc2_width = baseline_params['arc2_width'].get()

            # Fixed baseline parameters
            params.add('arc1_center', value=arc1_center, vary=False)
            params.add('arc1_height', value=arc1_height, vary=False)
            params.add('arc1_width', value=arc1_width, vary=False)
            params.add('arc2_center', value=arc2_center, vary=False)

            # If heights are locked, constrain arc2_height to arc1_height
            if lock_heights:
                params.add('arc2_height', expr='arc1_height')
            else:
                params.add('arc2_height', value=arc2_height, vary=False)

            # If widths are locked, constrain arc2_width to arc1_width
            if lock_widths:
                params.add('arc2_width', expr='arc1_width')
            else:
                params.add('arc2_width', value=arc2_width, vary=False)

            params.add('baseline_total', expr='arc1_height + arc2_height')
        else:
            # Default baseline parameters
            params.add('arc1_center', value=290.5, min=290.2, max=290.8)
            params.add('arc1_height',
                      value=est['step1_height'],
                      min=est['data_range'] * 0.1,
                      max=est['data_range'] * 0.9)
            params.add('arc1_width', value=0.4, min=0.2, max=0.6)

            params.add('arc2_center', value=292.0, min=291.7, max=292.3)
            params.add('arc2_height',
                      value=est['step2_height'],
                      min=est['data_range'] * 0.1,
                      max=est['data_range'] * 0.9)
            params.add('arc2_width', value=0.4, min=0.2, max=0.6)

            params.add('baseline_total',
                      expr='arc1_height + arc2_height',
                      min=0,
                      max=est['data_max'] * 1.1)

        # Main peaks with CUSTOM centers (respect peak_in_fit)
        max_peak_height = est['data_range'] * 1.5

        for i, name in enumerate(self.peak_names[:6], 1):  # First 6 peaks
            center = self.peak_center_vars[name].get()
            range_val = self.peak_range_vars[name].get()
            params.add(f'c{i}', value=center, min=center-range_val, max=center+range_val)
            if self.peak_in_fit[name].get():
                params.add(f'h{i}', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
            else:
                params.add(f'h{i}', value=0, vary=False)

        # Main FWHM from GUI
        params.add('main_fwhm',
                  value=self.fwhm_main_var.get(),
                  min=self.fwhm_main_min_var.get(),
                  max=self.fwhm_main_max_var.get())

        # Carbonate (peak 7) with custom center
        center = self.peak_center_vars['Carbonate'].get()
        range_val = self.peak_range_vars['Carbonate'].get()
        params.add('c7', value=center, min=center-range_val, max=center+range_val)
        if self.peak_in_fit['Carbonate'].get():
            params.add('h7', value=est['data_range'] * 0.2, min=0, max=max_peak_height * 0.8)
        else:
            params.add('h7', value=0, vary=False)
        # Carbonate FWHM from GUI
        params.add('carb_fwhm',
                  value=self.fwhm_carb_var.get(),
                  min=self.fwhm_carb_min_var.get(),
                  max=self.fwhm_carb_max_var.get())

        # Sigma* (peak 8) with custom center
        center = self.peak_center_vars['Sigma_star'].get()
        range_val = self.peak_range_vars['Sigma_star'].get()
        params.add('c8', value=center, min=center-range_val, max=center+range_val)
        if self.peak_in_fit['Sigma_star'].get():
            params.add('h8', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
        else:
            params.add('h8', value=0, vary=False)
        # Sigma* FWHM from GUI
        params.add('sigma_fwhm',
                  value=self.fwhm_sigma_var.get(),
                  min=self.fwhm_sigma_min_var.get(),
                  max=self.fwhm_sigma_max_var.get())

        # Perform fit using Model approach (more reliable)
        from lmfit import Model
        result = Model(total_model).fit(intensity, params, x=energy, method='leastsq', max_nfev=10000)

        return result

    def refit_with_enabled_peaks(self):
        """Refit spectrum with only enabled peaks (disabled peaks have h=0 fixed).

        Respects custom peak centers, manual baseline, and peak_in_fit settings.
        """
        if self.current_energy is None or self.current_intensity is None:
            return

        from lmfit import Model, Parameters
        from c1s_fitter_optimized import estimate_baseline_parameters, total_model

        energy = self.current_energy
        intensity = self.current_intensity

        params = Parameters()

        # Estimate data characteristics for baseline
        est = estimate_baseline_parameters(energy, intensity)

        # Baseline parameters - respect manual baseline if active
        if self.manual_baseline.get():
            params.add('arc1_center', value=self.baseline_params['arc1_center'].get(), vary=False)
            params.add('arc1_height', value=self.baseline_params['arc1_height'].get(), vary=False)
            params.add('arc1_width', value=self.baseline_params['arc1_width'].get(), vary=False)
            params.add('arc2_center', value=self.baseline_params['arc2_center'].get(), vary=False)
            if self.lock_heights.get():
                params.add('arc2_height', expr='arc1_height')
            else:
                params.add('arc2_height', value=self.baseline_params['arc2_height'].get(), vary=False)
            if self.lock_widths.get():
                params.add('arc2_width', expr='arc1_width')
            else:
                params.add('arc2_width', value=self.baseline_params['arc2_width'].get(), vary=False)
            params.add('baseline_total', expr='arc1_height + arc2_height')
        else:
            # Use current fit result baseline values as starting points if available
            if self.current_result is not None:
                bl_vals = self.current_result.params
                arc1_c = bl_vals['arc1_center'].value
                arc1_h = bl_vals['arc1_height'].value
                arc1_w = bl_vals['arc1_width'].value
                arc2_c = bl_vals['arc2_center'].value
                arc2_h = bl_vals['arc2_height'].value
                arc2_w = bl_vals['arc2_width'].value
            else:
                arc1_c = 290.5
                arc1_h = est['step1_height']
                arc1_w = 0.4
                arc2_c = 292.0
                arc2_h = est['step2_height']
                arc2_w = 0.4

            params.add('arc1_center', value=arc1_c, min=290.2, max=290.8)
            params.add('arc1_height',
                      value=arc1_h,
                      min=est['data_range'] * 0.01,
                      max=est['data_range'] * 1.2)
            params.add('arc1_width', value=arc1_w, min=0.2, max=0.6)

            params.add('arc2_center', value=arc2_c, min=291.7, max=292.3)
            params.add('arc2_height',
                      value=arc2_h,
                      min=est['data_range'] * 0.01,
                      max=est['data_range'] * 1.2)
            params.add('arc2_width', value=arc2_w, min=0.2, max=0.6)

            params.add('baseline_total',
                      expr='arc1_height + arc2_height',
                      min=0,
                      max=est['data_max'] * 1.5)

        # Peak parameters - respect enabled/disabled state
        max_peak_height = est['data_range'] * 1.5

        # Map peak names to indices
        peak_indices = {name: i for i, name in enumerate(self.peak_names, 1)}

        for i, name in enumerate(self.peak_names[:6], 1):  # First 6 main peaks
            # Use custom centers if enabled, otherwise default
            if self.custom_peak_centers.get():
                center = self.peak_center_vars[name].get()
                range_val = self.peak_range_vars[name].get()
            else:
                center = self.default_peak_centers[name]
                range_val = 0.3

            params.add(f'c{i}', value=center, min=center-range_val, max=center+range_val)

            # If peak is disabled, fix height to 0
            if self.peak_in_fit[name].get():
                params.add(f'h{i}', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
            else:
                params.add(f'h{i}', value=0, vary=False)

        # Main FWHM
        if self.custom_peak_centers.get():
            params.add('main_fwhm',
                      value=self.fwhm_main_var.get(),
                      min=self.fwhm_main_min_var.get(),
                      max=self.fwhm_main_max_var.get())
        else:
            params.add('main_fwhm', value=1.5, min=0.8, max=2.0)

        # Carbonate (peak 7)
        name = 'Carbonate'
        if self.custom_peak_centers.get():
            center = self.peak_center_vars[name].get()
            range_val = self.peak_range_vars[name].get()
        else:
            center = self.default_peak_centers[name]
            range_val = 0.3

        params.add('c7', value=center, min=center-range_val, max=center+range_val)
        if self.peak_in_fit[name].get():
            params.add('h7', value=est['data_range'] * 0.2, min=0, max=max_peak_height * 0.8)
        else:
            params.add('h7', value=0, vary=False)

        if self.custom_peak_centers.get():
            params.add('carb_fwhm',
                      value=self.fwhm_carb_var.get(),
                      min=self.fwhm_carb_min_var.get(),
                      max=self.fwhm_carb_max_var.get())
        else:
            params.add('carb_fwhm', value=0.8, min=0.5, max=1.0)

        # Sigma* (peak 8)
        name = 'Sigma_star'
        if self.custom_peak_centers.get():
            center = self.peak_center_vars[name].get()
            range_val = self.peak_range_vars[name].get()
        else:
            center = self.default_peak_centers[name]
            range_val = 0.5

        params.add('c8', value=center, min=center-range_val, max=center+range_val)
        if self.peak_in_fit[name].get():
            params.add('h8', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
        else:
            params.add('h8', value=0, vary=False)

        if self.custom_peak_centers.get():
            params.add('sigma_fwhm',
                      value=self.fwhm_sigma_var.get(),
                      min=self.fwhm_sigma_min_var.get(),
                      max=self.fwhm_sigma_max_var.get())
        else:
            params.add('sigma_fwhm', value=0.6, min=0.4, max=2.0)

        # Perform fit
        try:
            result = Model(total_model).fit(intensity, params, x=energy, method='leastsq', max_nfev=10000)
            self.current_result = result
        except Exception as e:
            print(f"Refit failed: {e}")

    def _build_fit_params(self, est, peak_set=None, baseline_values=None):
        """Build lmfit Parameters for enabled peaks.

        Args:
            est: dict from estimate_baseline_parameters()
            peak_set: set of peak names to enable (None = use current checkboxes)
            baseline_values: optional dict with keys arc1_center, arc1_height,
                arc1_width, arc2_center, arc2_height, arc2_width to use as
                starting points (e.g. from a prior fit)

        Returns:
            lmfit.Parameters object ready for fitting
        """
        from lmfit import Parameters

        params = Parameters()
        max_peak_height = est['data_range'] * 1.5

        # --- Baseline parameters ---
        # Use wider bounds matching the original run_model_selection
        if baseline_values is not None:
            bl = baseline_values
            params.add('arc1_center', value=bl['arc1_center'], min=290.2, max=290.8)
            params.add('arc1_height', value=bl['arc1_height'],
                      min=est['data_range'] * 0.01, max=est['data_range'] * 1.2)
            params.add('arc1_width', value=bl['arc1_width'], min=0.2, max=0.6)
            params.add('arc2_center', value=bl['arc2_center'], min=291.7, max=292.3)
            params.add('arc2_height', value=bl['arc2_height'],
                      min=est['data_range'] * 0.01, max=est['data_range'] * 1.2)
            params.add('arc2_width', value=bl['arc2_width'], min=0.2, max=0.6)
        else:
            params.add('arc1_center', value=290.5, min=290.2, max=290.8)
            params.add('arc1_height', value=est['step1_height'],
                      min=est['data_range'] * 0.01, max=est['data_range'] * 1.2)
            params.add('arc1_width', value=0.4, min=0.2, max=0.6)
            params.add('arc2_center', value=292.0, min=291.7, max=292.3)
            params.add('arc2_height', value=est['step2_height'],
                      min=est['data_range'] * 0.01, max=est['data_range'] * 1.2)
            params.add('arc2_width', value=0.4, min=0.2, max=0.6)

        params.add('baseline_total', expr='arc1_height + arc2_height',
                  min=0, max=est['data_max'] * 1.5)

        # --- Peak parameters (peaks 1-6, shared main_fwhm) ---
        any_main_enabled = False
        for i, name in enumerate(self.peak_names[:6], 1):
            if self.custom_peak_centers.get():
                center = self.peak_center_vars[name].get()
                range_val = self.peak_range_vars[name].get()
            else:
                center = self.default_peak_centers[name]
                range_val = 0.3

            enabled = (name in peak_set) if peak_set is not None else self.peak_in_fit[name].get()
            if enabled:
                params.add(f'c{i}', value=center, min=center - range_val, max=center + range_val)
                params.add(f'h{i}', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
                any_main_enabled = True
            else:
                params.add(f'c{i}', value=center, min=center - range_val, max=center + range_val,
                          vary=False)
                params.add(f'h{i}', value=0, vary=False)

        # Main FWHM — only free if at least one of peaks 1-6 is enabled
        if self.custom_peak_centers.get():
            params.add('main_fwhm', value=self.fwhm_main_var.get(),
                      min=self.fwhm_main_min_var.get(), max=self.fwhm_main_max_var.get(),
                      vary=any_main_enabled)
        else:
            params.add('main_fwhm', value=1.5, min=0.8, max=2.0, vary=any_main_enabled)

        # --- Carbonate (peak 7) ---
        cname = 'Carbonate'
        if self.custom_peak_centers.get():
            center = self.peak_center_vars[cname].get()
            range_val = self.peak_range_vars[cname].get()
        else:
            center = self.default_peak_centers[cname]
            range_val = 0.3

        carb_enabled = (cname in peak_set) if peak_set is not None else self.peak_in_fit[cname].get()
        if carb_enabled:
            params.add('c7', value=center, min=center - range_val, max=center + range_val)
            params.add('h7', value=est['data_range'] * 0.2, min=0, max=max_peak_height * 0.8)
            if self.custom_peak_centers.get():
                params.add('carb_fwhm', value=self.fwhm_carb_var.get(),
                          min=self.fwhm_carb_min_var.get(), max=self.fwhm_carb_max_var.get())
            else:
                params.add('carb_fwhm', value=0.8, min=0.5, max=1.0)
        else:
            params.add('c7', value=center, vary=False)
            params.add('h7', value=0, vary=False)
            params.add('carb_fwhm', value=0.8, vary=False)

        # --- Sigma* (peak 8) ---
        sname = 'Sigma_star'
        if self.custom_peak_centers.get():
            center = self.peak_center_vars[sname].get()
            range_val = self.peak_range_vars[sname].get()
        else:
            center = self.default_peak_centers[sname]
            range_val = 0.5

        sigma_enabled = (sname in peak_set) if peak_set is not None else self.peak_in_fit[sname].get()
        if sigma_enabled:
            params.add('c8', value=center, min=center - range_val, max=center + range_val)
            params.add('h8', value=est['data_range'] * 0.3, min=0, max=max_peak_height)
            if self.custom_peak_centers.get():
                params.add('sigma_fwhm', value=self.fwhm_sigma_var.get(),
                          min=self.fwhm_sigma_min_var.get(), max=self.fwhm_sigma_max_var.get())
            else:
                params.add('sigma_fwhm', value=0.6, min=0.4, max=2.0)
        else:
            params.add('c8', value=center, vary=False)
            params.add('h8', value=0, vary=False)
            params.add('sigma_fwhm', value=0.6, vary=False)

        return params

    # ------------------------------------------------------------------
    # Core helpers (GUI-free) for model selection & box plot
    # ------------------------------------------------------------------

    def _run_model_selection_core(self, energy, intensity, peak_set, progress_cb=None):
        """Run leave-one-out AIC model selection (no GUI).

        Args:
            energy: numpy array
            intensity: numpy array
            peak_set: set of peak name strings currently in the model
            progress_cb: optional callable(step, total) for progress

        Returns:
            list of result dicts, each with keys:
            label, peaks, n_peaks, n_free, aic, bic, r_squared, abbe, resid_std
        """
        from lmfit import Model
        from c1s_fitter_optimized import estimate_baseline_parameters, total_model

        est = estimate_baseline_parameters(energy, intensity)
        n_data = len(energy)

        # Step 1: fit the full model first to establish good baseline values
        full_params = self._build_fit_params(est, peak_set=peak_set)
        baseline_vals = None
        try:
            full_result = Model(total_model).fit(intensity, full_params, x=energy,
                                                 method='leastsq', max_nfev=10000)
            baseline_vals = {
                'arc1_center': full_result.params['arc1_center'].value,
                'arc1_height': full_result.params['arc1_height'].value,
                'arc1_width': full_result.params['arc1_width'].value,
                'arc2_center': full_result.params['arc2_center'].value,
                'arc2_height': full_result.params['arc2_height'].value,
                'arc2_width': full_result.params['arc2_width'].value,
            }
        except Exception:
            pass

        # Step 2: build configurations — full model + leave-one-out for each peak
        configs = [("Full model", set(peak_set))]
        for name in sorted(peak_set):
            idx = self.peak_names.index(name)
            display = self.peak_display_names[idx]
            configs.append((f"Without {display}", peak_set - {name}))

        results = []
        for cfg_idx, (label, test_set) in enumerate(configs):
            params = self._build_fit_params(est, peak_set=test_set,
                                            baseline_values=baseline_vals)

            # Count free parameters
            n_free = sum(1 for p in params.values()
                        if p.vary and p.expr is None)

            try:
                result = Model(total_model).fit(intensity, params, x=energy,
                                                method='leastsq', max_nfev=10000)
                residuals = result.residual
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((intensity - np.mean(intensity))**2)
                r_squared = 1 - ss_res / ss_tot
                resid_std = np.std(residuals)
                aic, bic = self._compute_aic_bic(residuals, n_free)
                abbe = self._abbe_criterion(residuals)

                results.append({
                    'label': label, 'peaks': test_set.copy(),
                    'n_peaks': len(test_set), 'n_free': n_free,
                    'aic': aic, 'bic': bic, 'r_squared': r_squared,
                    'abbe': abbe, 'resid_std': resid_std,
                })
            except Exception:
                results.append({
                    'label': label, 'peaks': test_set.copy(),
                    'n_peaks': len(test_set), 'n_free': n_free,
                    'aic': np.inf, 'bic': np.inf, 'r_squared': 0.0,
                    'abbe': 0.0, 'resid_std': np.inf,
                })

            if progress_cb:
                progress_cb(cfg_idx + 1, len(configs))

        return results

    def _run_box_plot_core(self, energy, intensity, peak_names_set,
                           n_iterations=32, progress_cb=None):
        """Run box plot overfitting test (no GUI).

        Args:
            energy: numpy array
            intensity: numpy array
            peak_names_set: set of peak name strings to include
            n_iterations: number of random-start fits
            progress_cb: optional callable(step, total)

        Returns:
            (all_areas, all_r_squared) where all_areas is
            {peak_name: [area_values]} and all_r_squared is a list
        """
        from lmfit import Model, Parameters
        from c1s_fitter_optimized import estimate_baseline_parameters, total_model

        est = estimate_baseline_parameters(energy, intensity)
        max_peak_height = est['data_range'] * 1.5
        rng = np.random.default_rng()

        enabled_peaks = [(self.peak_names.index(n) + 1, n)
                         for n in self.peak_names if n in peak_names_set]

        all_areas = {name: [] for _, name in enabled_peaks}
        all_r_squared = []

        for iteration in range(n_iterations):
            # Build params with random starting heights
            params = self._build_fit_params(est, peak_set=peak_names_set)

            # Randomize heights for enabled peaks
            for i, name in enabled_peaks:
                pname = f'h{i}'
                if params[pname].vary:
                    hi = params[pname].max if np.isfinite(params[pname].max) else max_peak_height
                    params[pname].value = rng.uniform(0, hi)

            try:
                result = Model(total_model).fit(intensity, params, x=energy,
                                                method='leastsq', max_nfev=10000)
                areas = calculate_peak_areas(energy, result)
                for _, pname in enabled_peaks:
                    all_areas[pname].append(areas.get(pname, 0.0))
                r_sq = 1 - result.residual.var() / np.var(result.data)
                all_r_squared.append(r_sq)
            except Exception:
                for _, pname in enabled_peaks:
                    all_areas[pname].append(np.nan)
                all_r_squared.append(np.nan)

            if progress_cb:
                progress_cb(iteration + 1, n_iterations)

        return all_areas, all_r_squared

    def _run_monte_carlo_core(self, energy, intensity, peak_names_set,
                              fit_result, n_iterations=32, progress_cb=None):
        """Run Monte Carlo noise perturbation test (no GUI).

        Generates synthetic spectra by adding Gaussian noise to the best-fit
        curve, refits each, and returns the distribution of peak areas.

        Args:
            energy: numpy array
            intensity: numpy array
            peak_names_set: set of peak name strings to include
            fit_result: lmfit result to extract best-fit curve and noise level
            n_iterations: number of synthetic fits
            progress_cb: optional callable(step, total)

        Returns:
            (all_areas, all_r_squared) same format as _run_box_plot_core
        """
        from lmfit import Model
        from c1s_fitter_optimized import estimate_baseline_parameters, total_model

        est = estimate_baseline_parameters(energy, intensity)
        rng = np.random.default_rng()

        noise_level = np.std(fit_result.residual)
        best_fit_curve = fit_result.best_fit.copy()

        # Extract baseline values from fit result
        bl = fit_result.params
        bl_vals = {
            'arc1_center': bl['arc1_center'].value,
            'arc1_height': bl['arc1_height'].value,
            'arc1_width': bl['arc1_width'].value,
            'arc2_center': bl['arc2_center'].value,
            'arc2_height': bl['arc2_height'].value,
            'arc2_width': bl['arc2_width'].value,
        }

        enabled_peaks = [(self.peak_names.index(n) + 1, n)
                         for n in self.peak_names if n in peak_names_set]

        all_areas = {name: [] for _, name in enabled_peaks}
        all_r_squared = []

        for iteration in range(n_iterations):
            # Synthetic spectrum: best-fit + random noise
            synthetic = best_fit_curve + rng.normal(0, noise_level, len(energy))

            # Build params using helper (current baseline as starting point)
            params = self._build_fit_params(est, baseline_values=bl_vals,
                                            peak_set=peak_names_set)

            # Seed peak heights from current fit (test noise sensitivity, not starting conditions)
            for pi, pname in enabled_peaks:
                hkey = f'h{pi}'
                if hkey in fit_result.params:
                    params[hkey].set(value=fit_result.params[hkey].value)

            try:
                result = Model(total_model).fit(synthetic, params, x=energy,
                                                method='leastsq', max_nfev=10000)
                areas = calculate_peak_areas(energy, result)
                for _, pname in enabled_peaks:
                    all_areas[pname].append(areas.get(pname, 0.0))
                r_sq = 1 - result.residual.var() / np.var(result.data)
                all_r_squared.append(r_sq)
            except Exception:
                for _, pname in enabled_peaks:
                    all_areas[pname].append(np.nan)
                all_r_squared.append(np.nan)

            if progress_cb:
                progress_cb(iteration + 1, n_iterations)

        return all_areas, all_r_squared

    # ------------------------------------------------------------------
    # Automated iterative model reduction
    # ------------------------------------------------------------------

    def _auto_reduce_core(self, energy, intensity, initial_peak_set,
                          threshold, progress_cb=None):
        """Iteratively remove expendable peaks using AIC+BIC consensus, box plot CV,
        and Monte Carlo noise sensitivity.

        Strategy:
        1. Primary: remove peaks where BOTH ΔAIC and ΔBIC < threshold
        2. Secondary: if Phase 1 keeps all peaks, run box plot test. If any peak
           has CV > 20% (overfitting), remove the peak with the highest CV
           that also has ΔAIC < relaxed and ΔBIC < relaxed.
        3. Tertiary: if Phase 2 finds no removals, run Monte Carlo noise test.
           If any peak has CV > 20% and ΔAIC+ΔBIC < relaxed, remove it.
        4. Repeat until the model is stable (no removals in any phase).

        Args:
            energy, intensity: spectrum data
            initial_peak_set: set of peak name strings
            threshold: delta-AIC/BIC below which a peak is considered expendable
            progress_cb: callable(message_str) for logging

        Returns:
            (removal_log, final_peak_set, final_model_results, box_plot_data)
            where box_plot_data is (all_areas, all_r_squared) or None
        """
        CV_OVERFIT_THRESHOLD = 20.0   # CV% above which a peak is poorly determined
        AIC_RELAXED_THRESHOLD = 7.0   # max ΔAIC/ΔBIC for CV-driven removal

        current_peaks = set(initial_peak_set)
        removal_log = []
        iteration = 0

        while len(current_peaks) >= 2:
            iteration += 1
            if progress_cb:
                progress_cb(f"\nIteration {iteration}: testing {len(current_peaks)} peaks...")

            # --- Phase 1: AIC+BIC consensus model selection ---
            ms_results = self._run_model_selection_core(energy, intensity, current_peaks)

            full_aic = None
            full_bic = None
            for r in ms_results:
                if r['label'] == 'Full model':
                    full_aic = r['aic']
                    full_bic = r['bic']
                    break
            if full_aic is None or not np.isfinite(full_aic):
                if progress_cb:
                    progress_cb("  Full model fit failed. Stopping.")
                break

            # Build ΔAIC and ΔBIC lookups: peak_name -> delta
            delta_aic_map = {}
            delta_bic_map = {}
            candidates = []
            for r in ms_results:
                if r['label'] == 'Full model':
                    continue
                d_aic = r['aic'] - full_aic
                d_bic = r['bic'] - full_bic
                candidates.append((d_aic, d_bic, r['label'], r['peaks']))
                # Extract removed peak name
                removed = current_peaks - r['peaks']
                if removed:
                    rname = removed.pop()
                    delta_aic_map[rname] = d_aic
                    delta_bic_map[rname] = d_bic

            if progress_cb:
                for d_aic, d_bic, label, _ in sorted(candidates, key=lambda x: x[0]):
                    progress_cb(f"  {label}: \u0394AIC = {d_aic:.2f}, \u0394BIC = {d_bic:.2f}")

            # Check primary criterion: BOTH ΔAIC < threshold AND ΔBIC < threshold
            expendable = [(d_aic, d_bic, lbl, ps) for d_aic, d_bic, lbl, ps in candidates
                          if d_aic < threshold and d_bic < threshold]

            if expendable:
                # Remove the most expendable (lowest ΔAIC)
                expendable.sort(key=lambda x: x[0])
                best_d_aic, best_d_bic, best_label, best_peaks = expendable[0]
                removed = current_peaks - best_peaks
                removed_name = removed.pop() if removed else "?"
                removed_display = self.peak_display_names[self.peak_names.index(removed_name)] \
                    if removed_name in self.peak_names else removed_name

                removal_log.append({
                    'iteration': iteration,
                    'removed': removed_display,
                    'removed_name': removed_name,
                    'delta_aic': best_d_aic,
                    'delta_bic': best_d_bic,
                    'remaining': len(best_peaks),
                    'criterion': 'AIC+BIC',
                })
                if progress_cb:
                    progress_cb(f"  >> Removing {removed_display} "
                                f"(\u0394AIC = {best_d_aic:.2f}, \u0394BIC = {best_d_bic:.2f}, "
                                f"criterion: AIC+BIC)")
                current_peaks = best_peaks
                continue  # re-enter loop with reduced model

            # --- Phase 2: Box plot overfitting check (random starting conditions) ---
            if progress_cb:
                progress_cb(f"  No peaks with \u0394AIC AND \u0394BIC < {threshold:.1f}.")
                progress_cb(f"  Running box plot test (random starting conditions)...")

            box_areas, box_r2 = self._run_box_plot_core(
                energy, intensity, current_peaks, n_iterations=32)

            # Compute per-peak CV (std / mean)
            peak_cvs = {}
            for name in current_peaks:
                vals = np.array(box_areas.get(name, []))
                valid = vals[~np.isnan(vals)]
                if len(valid) > 0:
                    mean_val = np.mean(valid)
                    cv = (np.std(valid) / mean_val * 100) if mean_val > 0 else 0.0
                else:
                    cv = 0.0
                peak_cvs[name] = cv

            max_cv = max(peak_cvs.values()) if peak_cvs else 0
            if progress_cb:
                for name in sorted(peak_cvs, key=lambda n: peak_cvs[n], reverse=True):
                    disp = self.peak_display_names[self.peak_names.index(name)]
                    progress_cb(f"    {disp}: CV = {peak_cvs[name]:.1f}%")

            if max_cv < CV_OVERFIT_THRESHOLD:
                if progress_cb:
                    progress_cb(f"  Box plot OK (max CV = {max_cv:.1f}% < {CV_OVERFIT_THRESHOLD}%). "
                                f"Model is stable.")
                break

            # Find the peak with highest CV that has ΔAIC and ΔBIC < relaxed threshold
            overfit_candidates = [
                (peak_cvs[name], name) for name in current_peaks
                if peak_cvs[name] >= CV_OVERFIT_THRESHOLD
                and delta_aic_map.get(name, np.inf) < AIC_RELAXED_THRESHOLD
                and delta_bic_map.get(name, np.inf) < AIC_RELAXED_THRESHOLD
            ]

            if not overfit_candidates:
                if progress_cb:
                    progress_cb(f"  Overfitting detected (max CV = {max_cv:.1f}%) but "
                                f"all overfit peaks have \u0394AIC or \u0394BIC >= {AIC_RELAXED_THRESHOLD:.1f}.")

                # --- Phase 3: Monte Carlo noise-sensitivity check ---
                if progress_cb:
                    progress_cb(f"  Running Monte Carlo test (noise sensitivity)...")

                # Fit current peak set to get baseline result for Monte Carlo
                from lmfit import Model as _Model
                from c1s_fitter_optimized import estimate_baseline_parameters as _est_bl, total_model as _tm
                _est = _est_bl(energy, intensity)
                _params = self._build_fit_params(_est, peak_set=current_peaks)
                try:
                    _fit_result = _Model(_tm).fit(intensity, _params, x=energy,
                                                   method='leastsq', max_nfev=10000)
                except Exception:
                    if progress_cb:
                        progress_cb("  Monte Carlo fit failed. Stopping.")
                    break

                mc_areas, mc_r2 = self._run_monte_carlo_core(
                    energy, intensity, current_peaks, _fit_result, n_iterations=32)

                # Compute per-peak CV from Monte Carlo
                mc_cvs = {}
                for name in current_peaks:
                    vals = np.array(mc_areas.get(name, []))
                    valid = vals[~np.isnan(vals)]
                    if len(valid) > 0:
                        mean_val = np.mean(valid)
                        cv = (np.std(valid) / mean_val * 100) if mean_val > 0 else 0.0
                    else:
                        cv = 0.0
                    mc_cvs[name] = cv

                mc_max_cv = max(mc_cvs.values()) if mc_cvs else 0
                if progress_cb:
                    for name in sorted(mc_cvs, key=lambda n: mc_cvs[n], reverse=True):
                        disp = self.peak_display_names[self.peak_names.index(name)]
                        progress_cb(f"    {disp}: MC CV = {mc_cvs[name]:.1f}%")

                if mc_max_cv < CV_OVERFIT_THRESHOLD:
                    if progress_cb:
                        progress_cb(f"  Monte Carlo OK (max CV = {mc_max_cv:.1f}% "
                                    f"< {CV_OVERFIT_THRESHOLD}%). Model is stable.")
                    break

                # Find peak with highest MC CV that has ΔAIC+ΔBIC < relaxed
                mc_overfit = [
                    (mc_cvs[name], name) for name in current_peaks
                    if mc_cvs[name] >= CV_OVERFIT_THRESHOLD
                    and delta_aic_map.get(name, np.inf) < AIC_RELAXED_THRESHOLD
                    and delta_bic_map.get(name, np.inf) < AIC_RELAXED_THRESHOLD
                ]

                if not mc_overfit:
                    if progress_cb:
                        progress_cb(f"  Monte Carlo sensitivity detected (max CV = {mc_max_cv:.1f}%) "
                                    f"but all sensitive peaks have \u0394AIC or \u0394BIC >= "
                                    f"{AIC_RELAXED_THRESHOLD:.1f}. Stopping.")
                    break

                mc_overfit.sort(key=lambda x: x[0], reverse=True)
                mc_worst_cv, mc_worst_name = mc_overfit[0]
                mc_worst_display = self.peak_display_names[self.peak_names.index(mc_worst_name)]
                mc_worst_d_aic = delta_aic_map.get(mc_worst_name, 0)
                mc_worst_d_bic = delta_bic_map.get(mc_worst_name, 0)

                removal_log.append({
                    'iteration': iteration,
                    'removed': mc_worst_display,
                    'removed_name': mc_worst_name,
                    'delta_aic': mc_worst_d_aic,
                    'delta_bic': mc_worst_d_bic,
                    'remaining': len(current_peaks) - 1,
                    'criterion': f'Monte Carlo (CV={mc_worst_cv:.0f}%)',
                })
                if progress_cb:
                    progress_cb(f"  >> Removing {mc_worst_display} "
                                f"(MC CV = {mc_worst_cv:.0f}%, \u0394AIC = {mc_worst_d_aic:.2f}, "
                                f"\u0394BIC = {mc_worst_d_bic:.2f}, criterion: Monte Carlo)")
                current_peaks = current_peaks - {mc_worst_name}
                continue  # loop back to Phase 1

            # Remove the peak with the highest CV (Phase 2)
            overfit_candidates.sort(key=lambda x: x[0], reverse=True)
            worst_cv, worst_name = overfit_candidates[0]
            worst_display = self.peak_display_names[self.peak_names.index(worst_name)]
            worst_d_aic = delta_aic_map.get(worst_name, 0)
            worst_d_bic = delta_bic_map.get(worst_name, 0)

            removal_log.append({
                'iteration': iteration,
                'removed': worst_display,
                'removed_name': worst_name,
                'delta_aic': worst_d_aic,
                'delta_bic': worst_d_bic,
                'remaining': len(current_peaks) - 1,
                'criterion': f'Box plot (CV={worst_cv:.0f}%)',
            })
            if progress_cb:
                progress_cb(f"  >> Removing {worst_display} "
                            f"(CV = {worst_cv:.0f}%, \u0394AIC = {worst_d_aic:.2f}, "
                            f"\u0394BIC = {worst_d_bic:.2f}, criterion: box plot)")
            current_peaks = current_peaks - {worst_name}

        # Final validation
        final_results = self._run_model_selection_core(energy, intensity, current_peaks)

        if progress_cb:
            progress_cb(f"\nRunning final box plot validation on {len(current_peaks)}-peak model...")
        box_data = None
        if len(current_peaks) >= 2:
            box_data = self._run_box_plot_core(energy, intensity, current_peaks,
                                               n_iterations=32, progress_cb=None)

        return removal_log, current_peaks, final_results, box_data

    def run_auto_reduce(self):
        """Run auto-reduce on the current spectrum (GUI wrapper)."""
        if self.current_energy is None or self.current_intensity is None:
            messagebox.showwarning("No Data", "Load a spectrum first.")
            return

        enabled = {name for name in self.peak_names if self.peak_in_fit[name].get()}
        if len(enabled) < 2:
            messagebox.showwarning("Too Few Peaks",
                                   "At least 2 enabled peaks are needed.")
            return

        threshold = self.auto_reduce_threshold.get()

        # Progress window with scrolling text log
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Auto-Reduce")
        progress_win.geometry("520x400")
        progress_win.transient(self.root)
        progress_win.grab_set()
        progress_win.update_idletasks()
        px = self.root.winfo_x() + (self.root.winfo_width() - 520) // 2
        py = self.root.winfo_y() + (self.root.winfo_height() - 400) // 2
        progress_win.geometry(f"+{px}+{py}")

        sample_name = self.spectra_files[self.current_index].stem if self.spectra_files else "Spectrum"
        ttk.Label(progress_win, text=f"Auto-Reduce: {sample_name}",
                  font=('Helvetica', 11, 'bold')).pack(pady=(5, 2))
        ttk.Label(progress_win, text=f"\u0394AIC threshold: {threshold:.1f}").pack()

        log_text = tk.Text(progress_win, wrap=tk.WORD, font=('Courier', 9),
                          state=tk.DISABLED, height=20)
        log_scroll = ttk.Scrollbar(progress_win, orient=tk.VERTICAL, command=log_text.yview)
        log_text.configure(yscrollcommand=log_scroll.set)
        log_scroll.pack(side=tk.RIGHT, fill=tk.Y, padx=(0, 5), pady=5)
        log_text.pack(fill=tk.BOTH, expand=True, padx=(5, 0), pady=5)

        def append_log(msg):
            def _do():
                if not progress_win.winfo_exists():
                    return
                log_text.config(state=tk.NORMAL)
                log_text.insert(tk.END, msg + "\n")
                log_text.see(tk.END)
                log_text.config(state=tk.DISABLED)
            self.root.after(0, _do)

        def do_work():
            energy = self.current_energy.copy()
            intensity = self.current_intensity.copy()
            removal_log, final_peaks, final_results, box_data = \
                self._auto_reduce_core(energy, intensity, enabled, threshold,
                                       progress_cb=append_log)
            self.root.after(0, lambda: _on_complete(
                removal_log, final_peaks, final_results, box_data))

        def _on_complete(removal_log, final_peaks, final_results, box_data):
            progress_win.destroy()

            # Update checkboxes
            for name in self.peak_names:
                self.peak_in_fit[name].set(name in final_peaks)

            # Refit with the reduced model
            self.refit_with_enabled_peaks()
            self.update_peak_tree_display()
            self.update_plot()
            self.update_statistics()
            self.update_peak_parameters()
            self.update_statistics_tab()
            self._save_file_state()

            # Show summary window
            self._show_auto_reduce_summary(removal_log, final_peaks,
                                           final_results, box_data)

        thread = threading.Thread(target=do_work, daemon=True)
        thread.start()

    def _show_auto_reduce_summary(self, removal_log, final_peaks,
                                  final_results, box_data):
        """Display auto-reduce summary window with log, stats, and box plot."""
        win = tk.Toplevel(self.root)
        win.title("Auto-Reduce Results")
        win.geometry("900x700")
        win.transient(self.root)

        sample_name = self.spectra_files[self.current_index].stem if self.spectra_files else "Spectrum"
        ttk.Label(win, text=f"Auto-Reduce: {sample_name}",
                  font=('Helvetica', 13, 'bold')).pack(pady=(10, 5))

        # Removal log
        log_frame = ttk.LabelFrame(win, text="Iteration Log", padding=5)
        log_frame.pack(fill=tk.X, padx=10, pady=5)

        log_text = tk.Text(log_frame, height=6, wrap=tk.WORD, font=('Courier', 9))
        if removal_log:
            for entry in removal_log:
                criterion = entry.get('criterion', 'AIC+BIC')
                d_bic = entry.get('delta_bic', None)
                bic_str = f", \u0394BIC = {d_bic:.2f}" if d_bic is not None else ""
                log_text.insert(tk.END,
                    f"  Iter {entry['iteration']}: removed {entry['removed']} "
                    f"(\u0394AIC = {entry['delta_aic']:.2f}{bic_str}, {criterion}), "
                    f"{entry['remaining']} peaks remain\n")
        else:
            log_text.insert(tk.END, "  No peaks removed — all are needed.\n")
        log_text.config(state=tk.DISABLED)
        log_text.pack(fill=tk.X)

        # Final model stats
        final_display = [self.peak_display_names[self.peak_names.index(n)]
                        for n in self.peak_names if n in final_peaks]
        stats_frame = ttk.LabelFrame(win, text="Final Model", padding=5)
        stats_frame.pack(fill=tk.X, padx=10, pady=5)

        stats_text = tk.Text(stats_frame, height=4, wrap=tk.WORD, font=('Courier', 9))
        stats_text.insert(tk.END, f"  Peaks ({len(final_peaks)}): {', '.join(final_display)}\n")
        if final_results:
            full = next((r for r in final_results if r['label'] == 'Full model'), None)
            if full:
                stats_text.insert(tk.END,
                    f"  AIC = {full['aic']:.1f}  |  R\u00b2 = {full['r_squared']:.5f}  |  "
                    f"Abbe = {full['abbe']:.3f}\n")
        stats_text.config(state=tk.DISABLED)
        stats_text.pack(fill=tk.X)

        # Box plot figure
        if box_data is not None:
            all_areas, all_r_squared = box_data
            enabled_peaks = [(self.peak_names.index(n) + 1, n)
                            for n in self.peak_names if n in final_peaks]

            peak_names_display = []
            area_data = []
            for _, name in enabled_peaks:
                vals = np.array(all_areas.get(name, []))
                valid = vals[~np.isnan(vals)]
                if len(valid) > 0:
                    area_data.append(valid)
                    idx = self.peak_names.index(name)
                    peak_names_display.append(self.peak_display_names[idx])

            if area_data:
                fig = Figure(figsize=(9, 4.5), dpi=90)
                ax = fig.add_subplot(111)

                cvs = []
                for vals in area_data:
                    mean_val = np.mean(vals)
                    cv = (np.std(vals) / mean_val * 100) if mean_val > 0 else 0.0
                    cvs.append(cv)

                bp = ax.boxplot(area_data, patch_artist=True, widths=0.5,
                                medianprops=dict(color='black', linewidth=1.5),
                                showmeans=True,
                                meanprops=dict(marker='D', markerfacecolor='blue', markersize=5))

                for patch, cv in zip(bp['boxes'], cvs):
                    if cv < 10:
                        color = '#4CAF50'
                    elif cv < 20:
                        color = '#FFC107'
                    else:
                        color = '#F44336'
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

                ax.set_xticklabels(peak_names_display, fontsize=10, rotation=30, ha='right')
                ax.set_ylabel('Peak Area (a.u.)', fontsize=11)

                max_cv = max(cvs) if cvs else 0
                if max_cv < 10:
                    verdict = "ROBUST"
                    vc = '#4CAF50'
                elif max_cv < 20:
                    verdict = "CAUTION"
                    vc = '#FFC107'
                else:
                    verdict = "OVERFITTING"
                    vc = '#F44336'

                ax.set_title(f"Box Plot Validation: {verdict}", fontsize=12,
                            fontweight='bold', color=vc)
                for i, (cv, _) in enumerate(zip(cvs, peak_names_display)):
                    y_top = np.max(area_data[i])
                    ax.text(i + 1, y_top, f"CV={cv:.0f}%", ha='center',
                            va='bottom', fontsize=8, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='y')
                fig.tight_layout()

                canvas = FigureCanvasTkAgg(fig, master=win)
                canvas.draw()
                canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

    def run_batch_auto_reduce(self):
        """Run auto-reduce on all loaded spectra."""
        if not self.spectra_files:
            messagebox.showwarning("No Data", "Load spectra first.")
            return

        threshold = self.auto_reduce_threshold.get()
        n_files = len(self.spectra_files)

        # Progress window
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Batch Auto-Reduce")
        progress_win.geometry("550x400")
        progress_win.transient(self.root)
        progress_win.grab_set()
        progress_win.update_idletasks()
        px = self.root.winfo_x() + (self.root.winfo_width() - 550) // 2
        py = self.root.winfo_y() + (self.root.winfo_height() - 400) // 2
        progress_win.geometry(f"+{px}+{py}")

        ttk.Label(progress_win, text=f"Batch Auto-Reduce ({n_files} spectra)",
                  font=('Helvetica', 11, 'bold')).pack(pady=(5, 2))
        progress_var = tk.DoubleVar(value=0)
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var,
                                       maximum=n_files, length=480)
        progress_bar.pack(padx=10, pady=5)
        progress_label = ttk.Label(progress_win, text=f"0 / {n_files}")
        progress_label.pack()
        log_text = tk.Text(progress_win, wrap=tk.WORD, font=('Courier', 9),
                          state=tk.DISABLED)
        log_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        batch_results = []

        def do_batch():
            for idx, fpath in enumerate(self.spectra_files):
                file_key = str(fpath)
                try:
                    energy, intensity = load_spectrum(fpath)
                except Exception as e:
                    batch_results.append({
                        'name': fpath.stem, 'error': str(e),
                        'removed': [], 'final_count': 0,
                        'r_squared': 0, 'aic': np.inf, 'verdict': 'ERROR'
                    })
                    _update(idx + 1, f"{fpath.stem}: load error")
                    continue

                # Determine initial peak set from cache or default (all on)
                cached = self.file_state_cache.get(file_key)
                if cached:
                    initial_peaks = {n for n, v in cached['fit_peaks'].items() if v}
                else:
                    initial_peaks = set(self.peak_names)

                if len(initial_peaks) < 2:
                    initial_peaks = set(self.peak_names)

                removal_log, final_peaks, final_results, box_data = \
                    self._auto_reduce_core(energy, intensity, initial_peaks, threshold)

                # Extract stats
                full_r = next((r for r in (final_results or [])
                              if r['label'] == 'Full model'), None)
                r_sq = full_r['r_squared'] if full_r else 0
                aic = full_r['aic'] if full_r else np.inf

                # Box plot verdict
                verdict = 'N/A'
                if box_data:
                    all_areas, _ = box_data
                    max_cv = 0
                    for name in final_peaks:
                        vals = np.array(all_areas.get(name, []))
                        valid = vals[~np.isnan(vals)]
                        if len(valid) > 0:
                            mean_val = np.mean(valid)
                            cv = (np.std(valid) / mean_val * 100) if mean_val > 0 else 0
                            max_cv = max(max_cv, cv)
                    if max_cv < 10:
                        verdict = 'ROBUST'
                    elif max_cv < 20:
                        verdict = 'CAUTION'
                    else:
                        verdict = 'OVERFITTING'

                removed_names = [e['removed'] for e in removal_log]
                batch_results.append({
                    'name': fpath.stem, 'removed': removed_names,
                    'final_count': len(final_peaks),
                    'r_squared': r_sq, 'aic': aic, 'verdict': verdict,
                })

                # Update file_state_cache
                fit_peaks = {n: (n in final_peaks) for n in self.peak_names}
                quant_peaks = {n: (n in final_peaks) for n in self.peak_names}

                # Refit with final peaks to get result
                from lmfit import Model
                from c1s_fitter_optimized import estimate_baseline_parameters, total_model
                est = estimate_baseline_parameters(energy, intensity)
                params = self._build_fit_params(est, peak_set=final_peaks)
                try:
                    res = Model(total_model).fit(intensity, params, x=energy,
                                                 method='leastsq', max_nfev=10000)
                except Exception:
                    res = None

                self.file_state_cache[file_key] = {
                    'fit_peaks': fit_peaks,
                    'quant_peaks': quant_peaks,
                    'result': res,
                    'energy': energy,
                    'intensity': intensity,
                }

                rm_str = ', '.join(removed_names) if removed_names else 'none'
                _update(idx + 1, f"{fpath.stem}: removed [{rm_str}] -> "
                                 f"{len(final_peaks)} peaks, {verdict}")

            self.root.after(0, _on_complete)

        def _update(count, msg):
            def _do():
                if not progress_win.winfo_exists():
                    return
                progress_var.set(count)
                progress_label.config(text=f"{count} / {n_files}")
                log_text.config(state=tk.NORMAL)
                log_text.insert(tk.END, msg + "\n")
                log_text.see(tk.END)
                log_text.config(state=tk.DISABLED)
            self.root.after(0, _do)

        def _on_complete():
            progress_win.destroy()
            # Reload current spectrum from cache
            self.load_current_spectrum()
            self._show_batch_auto_reduce_summary(batch_results)

        thread = threading.Thread(target=do_batch, daemon=True)
        thread.start()

    def _show_batch_auto_reduce_summary(self, batch_results):
        """Display batch auto-reduce results in a summary table."""
        win = tk.Toplevel(self.root)
        win.title("Batch Auto-Reduce Results")
        win.geometry("900x500")
        win.transient(self.root)

        ttk.Label(win, text=f"Batch Auto-Reduce: {len(batch_results)} spectra",
                  font=('Helvetica', 13, 'bold')).pack(pady=(10, 5))

        columns = ('sample', 'removed', 'final', 'r_squared', 'aic', 'verdict')
        tree = ttk.Treeview(win, columns=columns, show='headings', height=20)
        tree.heading('sample', text='Sample')
        tree.heading('removed', text='Peaks Removed')
        tree.heading('final', text='# Final')
        tree.heading('r_squared', text='R\u00b2')
        tree.heading('aic', text='AIC')
        tree.heading('verdict', text='Box Plot')

        tree.column('sample', width=180, anchor='w')
        tree.column('removed', width=250, anchor='w')
        tree.column('final', width=60, anchor='center')
        tree.column('r_squared', width=90, anchor='e')
        tree.column('aic', width=90, anchor='e')
        tree.column('verdict', width=100, anchor='center')

        tree.tag_configure('robust', background='#C8E6C9')
        tree.tag_configure('caution', background='#FFF9C4')
        tree.tag_configure('overfitting', background='#FFCDD2')
        tree.tag_configure('error', background='#E0E0E0')

        for r in batch_results:
            if 'error' in r:
                tag = 'error'
                tree.insert('', 'end', values=(
                    r['name'], f"ERROR: {r['error']}", '', '', '', ''
                ), tags=(tag,))
                continue

            rm = ', '.join(r['removed']) if r['removed'] else 'None'
            tag = r['verdict'].lower() if r['verdict'] in ('ROBUST', 'CAUTION', 'OVERFITTING') else ''

            tree.insert('', 'end', values=(
                r['name'], rm, r['final_count'],
                f"{r['r_squared']:.5f}",
                f"{r['aic']:.1f}" if np.isfinite(r['aic']) else 'N/A',
                r['verdict']
            ), tags=(tag,))

        scrollbar = ttk.Scrollbar(win, orient='vertical', command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        tree_frame = ttk.Frame(win)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        tree.pack(in_=tree_frame, side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(in_=tree_frame, side=tk.RIGHT, fill=tk.Y)

    def run_box_plot_test(self):
        """Run box plot overfitting test (Moeini et al. 2021, adapted for C1s XAS).

        Performs N fits with randomized starting peak heights to test whether
        the fit converges to a single global minimum. Narrow box plots indicate
        a robust fit; wide box plots suggest overfitting / multiple local minima.
        """
        if self.current_energy is None or self.current_intensity is None:
            messagebox.showwarning("No Data", "Load a spectrum first.")
            return

        N_ITERATIONS = 32

        # Collect enabled peak info
        peak_names_set = {name for name in self.peak_names if self.peak_in_fit[name].get()}
        enabled_peaks = [(self.peak_names.index(n) + 1, n)
                         for n in self.peak_names if n in peak_names_set]

        if len(enabled_peaks) < 2:
            messagebox.showwarning("Too Few Peaks",
                                   "At least 2 enabled peaks are needed for the box plot test.")
            return

        # Create progress window
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Box Plot Test")
        progress_win.geometry("320x100")
        progress_win.resizable(False, False)
        progress_win.transient(self.root)
        progress_win.grab_set()

        progress_win.update_idletasks()
        px = self.root.winfo_x() + (self.root.winfo_width() - 320) // 2
        py = self.root.winfo_y() + (self.root.winfo_height() - 100) // 2
        progress_win.geometry(f"+{px}+{py}")

        ttk.Label(progress_win, text=f"Running {N_ITERATIONS} fits with random starting conditions...",
                  wraplength=280).pack(padx=10, pady=(10, 5))
        progress_var = tk.DoubleVar(value=0)
        ttk.Progressbar(progress_win, variable=progress_var,
                        maximum=N_ITERATIONS, length=280).pack(padx=10, pady=5)
        progress_label = ttk.Label(progress_win, text=f"0 / {N_ITERATIONS}")
        progress_label.pack()

        def do_fits():
            energy = self.current_energy.copy()
            intensity = self.current_intensity.copy()

            def on_progress(step, total):
                self.root.after(0, lambda s=step: _update_progress(s))

            all_areas, all_r_squared = self._run_box_plot_core(
                energy, intensity, peak_names_set,
                n_iterations=N_ITERATIONS, progress_cb=on_progress)
            self.root.after(0, lambda: _show_results(all_areas, all_r_squared))

        def _update_progress(count):
            if not progress_win.winfo_exists():
                return
            progress_var.set(count)
            progress_label.config(text=f"{count} / {N_ITERATIONS}")

        def _show_results(all_areas, all_r_squared):
            if progress_win.winfo_exists():
                progress_win.destroy()
            self._show_box_plot_results(enabled_peaks, all_areas, all_r_squared)

        thread = threading.Thread(target=do_fits, daemon=True)
        thread.start()

    def _show_box_plot_results(self, enabled_peaks, all_areas, all_r_squared):
        """Display box plot test results in a new window."""
        # Build data arrays (drop NaN rows)
        peak_names_display = []
        area_data = []
        for _, name in enabled_peaks:
            vals = np.array(all_areas[name])
            valid = vals[~np.isnan(vals)]
            if len(valid) > 0:
                area_data.append(valid)
                idx = self.peak_names.index(name)
                peak_names_display.append(self.peak_display_names[idx])

        if not area_data:
            messagebox.showerror("Error", "All fits failed.")
            return

        # Compute coefficient of variation (CV = std/mean) for each peak
        cvs = []
        for vals in area_data:
            mean_val = np.mean(vals)
            cv = (np.std(vals) / mean_val * 100) if mean_val > 0 else 0.0
            cvs.append(cv)

        # Create results window
        result_win = tk.Toplevel(self.root)
        result_win.title("Box Plot Overfitting Test")
        result_win.geometry("900x650")
        result_win.transient(self.root)

        fig = Figure(figsize=(10, 7), dpi=90)
        ax = fig.add_subplot(111)

        # Draw box plot
        bp = ax.boxplot(area_data, patch_artist=True, widths=0.5,
                        medianprops=dict(color='black', linewidth=1.5))

        # Color boxes: green if CV < 10%, yellow if 10-20%, red if > 20%
        for i, (patch, cv) in enumerate(zip(bp['boxes'], cvs)):
            if cv < 10:
                color = '#4CAF50'  # green
            elif cv < 20:
                color = '#FFC107'  # amber
            else:
                color = '#F44336'  # red
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_xticklabels(peak_names_display, fontsize=10, rotation=30, ha='right')
        ax.set_ylabel('Peak Area (a.u.)', fontsize=11)

        # Build title with overall assessment
        r_sq_arr = np.array(all_r_squared)
        r_sq_valid = r_sq_arr[~np.isnan(r_sq_arr)]
        max_cv = max(cvs) if cvs else 0
        if max_cv < 10:
            verdict = "ROBUST - single global minimum"
            verdict_color = '#4CAF50'
        elif max_cv < 20:
            verdict = "CAUTION - some parameter sensitivity"
            verdict_color = '#FFC107'
        else:
            verdict = "OVERFITTING - multiple local minima detected"
            verdict_color = '#F44336'

        sample_name = self.spectra_files[self.current_index].stem if self.spectra_files else "Unknown"
        ax.set_title(f"Box Plot Test: {sample_name}\n"
                     f"{len(area_data[0])} fits with random starting conditions",
                     fontsize=12, fontweight='bold')

        # Add CV annotations above each box
        for i, (cv, name) in enumerate(zip(cvs, peak_names_display)):
            y_top = np.max(area_data[i])
            ax.text(i + 1, y_top, f"CV={cv:.0f}%", ha='center', va='bottom',
                    fontsize=8, fontweight='bold')

        # Add verdict text
        ax.text(0.5, 0.02, verdict, transform=ax.transAxes,
                fontsize=13, fontweight='bold', color=verdict_color,
                ha='center', va='bottom',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor=verdict_color, alpha=0.9))

        # Add R² summary
        if len(r_sq_valid) > 0:
            ax.text(0.98, 0.98,
                    f"R² range: {r_sq_valid.min():.4f} - {r_sq_valid.max():.4f}\n"
                    f"R² mean: {np.mean(r_sq_valid):.4f}",
                    transform=ax.transAxes, fontsize=9, va='top', ha='right',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

        # Legend for colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#4CAF50', alpha=0.7, label='CV < 10%: Robust'),
            Patch(facecolor='#FFC107', alpha=0.7, label='CV 10-20%: Caution'),
            Patch(facecolor='#F44336', alpha=0.7, label='CV > 20%: Overfitting'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

        ax.grid(True, alpha=0.3, axis='y')
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=result_win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, result_win)
        toolbar.update()

    def run_model_selection(self):
        """Run AIC/BIC model selection by leave-one-out peak removal.

        Tests removing each enabled peak one at a time and compares models
        using AIC, BIC, R², Abbe criterion, and residual std dev.
        """
        if self.current_energy is None or self.current_intensity is None:
            messagebox.showwarning("No Data", "Load a spectrum first.")
            return

        peak_names_set = {name for name in self.peak_names if self.peak_in_fit[name].get()}

        if len(peak_names_set) < 2:
            messagebox.showwarning("Too Few Peaks",
                                   "At least 2 enabled peaks are needed for model selection.")
            return

        n_configs = len(peak_names_set) + 1  # full model + leave-one-out

        # Create progress window
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Model Selection")
        progress_win.geometry("320x100")
        progress_win.resizable(False, False)
        progress_win.transient(self.root)
        progress_win.grab_set()

        progress_win.update_idletasks()
        px = self.root.winfo_x() + (self.root.winfo_width() - 320) // 2
        py = self.root.winfo_y() + (self.root.winfo_height() - 100) // 2
        progress_win.geometry(f"+{px}+{py}")

        ttk.Label(progress_win, text=f"Testing {n_configs} model configurations...",
                  wraplength=280).pack(padx=10, pady=(10, 5))
        progress_var = tk.DoubleVar(value=0)
        ttk.Progressbar(progress_win, variable=progress_var,
                        maximum=n_configs, length=280).pack(padx=10, pady=5)
        progress_label = ttk.Label(progress_win, text=f"0 / {n_configs}")
        progress_label.pack()

        def do_model_fits():
            energy = self.current_energy.copy()
            intensity = self.current_intensity.copy()

            def on_progress(step, total):
                self.root.after(0, lambda s=step: _update_progress(s))

            results = self._run_model_selection_core(
                energy, intensity, peak_names_set, progress_cb=on_progress)
            self.root.after(0, lambda: _show_results(results))

        def _update_progress(count):
            if not progress_win.winfo_exists():
                return
            progress_var.set(count)
            progress_label.config(text=f"{count} / {n_configs}")

        def _show_results(results):
            if progress_win.winfo_exists():
                progress_win.destroy()
            self._show_model_selection_results(results)

        thread = threading.Thread(target=do_model_fits, daemon=True)
        thread.start()

    def _show_model_selection_results(self, results):
        """Display model selection results in a comparison table."""
        if not results:
            messagebox.showerror("Error", "No model selection results.")
            return

        # Sort by AIC
        results.sort(key=lambda r: r['aic'])

        best_aic = results[0]['aic']
        best_bic = min(r['bic'] for r in results)

        # Create results window
        result_win = tk.Toplevel(self.root)
        result_win.title("Model Selection Results")
        result_win.geometry("850x500")
        result_win.transient(self.root)

        sample_name = self.spectra_files[self.current_index].stem if self.spectra_files else "Unknown"
        header = ttk.Label(result_win,
                          text=f"Model Selection: {sample_name}",
                          font=('Helvetica', 13, 'bold'))
        header.pack(pady=(10, 5))

        # Table frame
        table_frame = ttk.Frame(result_win)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Treeview for results table
        columns = ('rank', 'model', 'peaks', 'params', 'aic', 'delta_aic',
                   'bic', 'delta_bic', 'r_squared', 'abbe', 'resid_std')
        tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=12)

        tree.heading('rank', text='Rank')
        tree.heading('model', text='Model')
        tree.heading('peaks', text='# Peaks')
        tree.heading('params', text='# Params')
        tree.heading('aic', text='AIC')
        tree.heading('delta_aic', text='\u0394AIC')
        tree.heading('bic', text='BIC')
        tree.heading('delta_bic', text='\u0394BIC')
        tree.heading('r_squared', text='R\u00b2')
        tree.heading('abbe', text='Abbe')
        tree.heading('resid_std', text='Resid SD')

        tree.column('rank', width=40, anchor='center')
        tree.column('model', width=150, anchor='w')
        tree.column('peaks', width=55, anchor='center')
        tree.column('params', width=60, anchor='center')
        tree.column('aic', width=75, anchor='e')
        tree.column('delta_aic', width=65, anchor='e')
        tree.column('bic', width=75, anchor='e')
        tree.column('delta_bic', width=65, anchor='e')
        tree.column('r_squared', width=70, anchor='e')
        tree.column('abbe', width=60, anchor='e')
        tree.column('resid_std', width=70, anchor='e')

        # Color tags
        tree.tag_configure('green', background='#C8E6C9')
        tree.tag_configure('yellow', background='#FFF9C4')
        tree.tag_configure('red', background='#FFCDD2')

        for rank, r in enumerate(results, 1):
            delta_aic = r['aic'] - best_aic
            delta_bic = r['bic'] - best_bic

            # Color coding based on delta-AIC
            if delta_aic < 2:
                tag = 'green'
            elif delta_aic < 7:
                tag = 'yellow'
            else:
                tag = 'red'

            tree.insert('', 'end', values=(
                rank,
                r['label'],
                r['n_peaks'],
                r['n_free'],
                f"{r['aic']:.1f}",
                f"{delta_aic:.1f}",
                f"{r['bic']:.1f}",
                f"{delta_bic:.1f}",
                f"{r['r_squared']:.5f}",
                f"{r['abbe']:.3f}",
                f"{r['resid_std']:.4f}",
            ), tags=(tag,))

        # Scrollbar
        scrollbar = ttk.Scrollbar(table_frame, orient='vertical', command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Recommendation text
        rec_frame = ttk.LabelFrame(result_win, text="Interpretation", padding=8)
        rec_frame.pack(fill=tk.X, padx=10, pady=(5, 10))

        best = results[0]
        supported = [r for r in results if (r['aic'] - best_aic) < 2]
        plausible = [r for r in results if 2 <= (r['aic'] - best_aic) < 7]

        lines = []
        lines.append(f"Best model (lowest AIC): {best['label']}")
        lines.append(f"  AIC = {best['aic']:.1f}, BIC = {best['bic']:.1f}, "
                     f"R\u00b2 = {best['r_squared']:.5f}")
        lines.append("")

        if len(supported) > 1:
            others = [r['label'] for r in supported if r is not best]
            lines.append(f"Substantial support (\u0394AIC < 2): {', '.join(others)}")
        else:
            lines.append("No other models have substantial support (\u0394AIC < 2).")

        if plausible:
            names = [r['label'] for r in plausible]
            lines.append(f"Less support (\u0394AIC 2-7): {', '.join(names)}")

        not_supported = [r for r in results if (r['aic'] - best_aic) >= 7]
        if not_supported:
            names = [r['label'] for r in not_supported]
            lines.append(f"Essentially no support (\u0394AIC > 7): {', '.join(names)}")

        lines.append("")
        lines.append("Legend:  Green = substantial support  |  "
                     "Yellow = less support  |  Red = no support")

        rec_text = tk.Text(rec_frame, height=8, wrap=tk.WORD, font=('Courier', 9))
        rec_text.insert('1.0', '\n'.join(lines))
        rec_text.config(state=tk.DISABLED)
        rec_text.pack(fill=tk.X)

    # ------------------------------------------------------------------
    # Monte Carlo Analysis (Moeini et al. 2021, Section 2.5)
    # ------------------------------------------------------------------

    def run_monte_carlo(self):
        """Run Monte Carlo analysis to assess parameter uncertainty due to noise.

        Generates N synthetic spectra by adding Gaussian noise to the current
        best-fit, refits each, and reports the distribution of peak areas.
        """
        if self.current_energy is None or self.current_intensity is None:
            messagebox.showwarning("No Data", "Load a spectrum first.")
            return
        if self.current_result is None:
            messagebox.showwarning("No Fit", "Fit the spectrum first.")
            return

        N_ITERATIONS = 32

        # Collect enabled peaks
        enabled_peaks = []
        for i, name in enumerate(self.peak_names, 1):
            if self.peak_in_fit[name].get():
                enabled_peaks.append((i, name))

        if len(enabled_peaks) < 1:
            messagebox.showwarning("No Peaks", "Enable at least one peak.")
            return

        # Noise level from current fit residuals
        noise_level = np.std(self.current_result.residual)
        best_fit_curve = self.current_result.best_fit.copy()

        # Create progress window
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Monte Carlo Analysis")
        progress_win.geometry("320x100")
        progress_win.resizable(False, False)
        progress_win.transient(self.root)
        progress_win.grab_set()

        progress_win.update_idletasks()
        px = self.root.winfo_x() + (self.root.winfo_width() - 320) // 2
        py = self.root.winfo_y() + (self.root.winfo_height() - 100) // 2
        progress_win.geometry(f"+{px}+{py}")

        ttk.Label(progress_win,
                  text=f"Running {N_ITERATIONS} fits on synthetic noisy spectra...",
                  wraplength=280).pack(padx=10, pady=(10, 5))
        progress_var = tk.DoubleVar(value=0)
        ttk.Progressbar(progress_win, variable=progress_var,
                        maximum=N_ITERATIONS, length=280).pack(padx=10, pady=5)
        progress_label = ttk.Label(progress_win, text=f"0 / {N_ITERATIONS}")
        progress_label.pack()

        all_areas = {name: [] for _, name in enabled_peaks}
        all_r_squared = []
        all_synthetic = []
        all_fit_curves = []
        # Per-iteration peak curves: list of dicts {peak_name: curve_array}
        all_peak_curves = []

        def do_fits():
            from lmfit import Model
            from c1s_fitter_optimized import (estimate_baseline_parameters,
                                              total_model, gaussian)

            energy = self.current_energy
            est = estimate_baseline_parameters(energy, self.current_intensity)
            rng = np.random.default_rng()

            # Extract baseline values once (outside loop)
            bl = self.current_result.params
            bl_vals = {
                'arc1_center': bl['arc1_center'].value,
                'arc1_height': bl['arc1_height'].value,
                'arc1_width': bl['arc1_width'].value,
                'arc2_center': bl['arc2_center'].value,
                'arc2_height': bl['arc2_height'].value,
                'arc2_width': bl['arc2_width'].value,
            }

            peak_names_set = {name for _, name in enabled_peaks}

            for iteration in range(N_ITERATIONS):
                # Synthetic spectrum: best-fit + random noise
                synthetic = best_fit_curve + rng.normal(0, noise_level, len(energy))
                all_synthetic.append(synthetic.copy())

                # Build params using helper (current baseline as starting point)
                params = self._build_fit_params(est, baseline_values=bl_vals,
                                                peak_set=peak_names_set)

                # Seed peak heights from current fit so we purely test
                # noise sensitivity, not starting-condition effects
                cur = self.current_result.params
                for pi, pname in enabled_peaks:
                    params[f'h{pi}'].set(value=cur[f'h{pi}'].value)

                try:
                    result = Model(total_model).fit(synthetic, params, x=energy,
                                                    method='leastsq', max_nfev=10000)
                    areas = calculate_peak_areas(energy, result)
                    for _, pname in enabled_peaks:
                        all_areas[pname].append(areas.get(pname, 0.0))
                    r_sq = 1 - result.residual.var() / np.var(result.data)
                    all_r_squared.append(r_sq)
                    all_fit_curves.append(result.best_fit.copy())

                    # Extract individual peak curves
                    peaks = {}
                    for pi, pname in enabled_peaks:
                        center = result.params[f'c{pi}'].value
                        height = result.params[f'h{pi}'].value
                        if pi <= 6:
                            fwhm = result.params['main_fwhm'].value
                        elif pi == 7:
                            fwhm = result.params['carb_fwhm'].value
                        else:
                            fwhm = result.params['sigma_fwhm'].value
                        peaks[pname] = gaussian(energy, center, height, fwhm)
                    all_peak_curves.append(peaks)
                except Exception:
                    for _, pname in enabled_peaks:
                        all_areas[pname].append(np.nan)
                    all_r_squared.append(np.nan)
                    all_fit_curves.append(None)
                    all_peak_curves.append(None)

                self.root.after(0, lambda it=iteration + 1: _update_progress(it))

            self.root.after(0, _show_results)

        def _update_progress(count):
            if not progress_win.winfo_exists():
                return
            progress_var.set(count)
            progress_label.config(text=f"{count} / {N_ITERATIONS}")

        def _show_results():
            if progress_win.winfo_exists():
                progress_win.destroy()
            self._show_monte_carlo_results(enabled_peaks, all_areas, all_r_squared,
                                           noise_level, all_synthetic,
                                           all_fit_curves, best_fit_curve,
                                           all_peak_curves)

        thread = threading.Thread(target=do_fits, daemon=True)
        thread.start()

    def _show_monte_carlo_results(self, enabled_peaks, all_areas, all_r_squared,
                                  noise_level, all_synthetic, all_fit_curves,
                                  best_fit_curve, all_peak_curves):
        """Display Monte Carlo analysis results: fit overlay and box plots."""
        from matplotlib.patches import Patch
        from c1s_fitter_optimized import gaussian

        peak_names_display = []
        area_data = []
        for _, name in enabled_peaks:
            vals = np.array(all_areas[name])
            valid = vals[~np.isnan(vals)]
            if len(valid) > 0:
                area_data.append(valid)
                idx = self.peak_names.index(name)
                peak_names_display.append(self.peak_display_names[idx])

        if not area_data:
            messagebox.showerror("Error", "All Monte Carlo fits failed.")
            return

        # CV for each peak
        cvs = []
        stds = []
        for vals in area_data:
            mean_val = np.mean(vals)
            cv = (np.std(vals) / mean_val * 100) if mean_val > 0 else 0.0
            cvs.append(cv)
            stds.append(np.std(vals))

        sample_name = (self.spectra_files[self.current_index].stem
                       if self.spectra_files else "Unknown")
        energy = self.current_energy

        # Build peak color map: peak_name -> matplotlib color
        peak_color_map = {}
        for pi, pname in enabled_peaks:
            peak_color_map[pname] = self.peak_colors[pi - 1]

        # Create results window with two subplots
        result_win = tk.Toplevel(self.root)
        result_win.title("Monte Carlo Analysis")
        result_win.geometry("1100x750")
        result_win.transient(self.root)

        fig = Figure(figsize=(12, 7), dpi=90)

        # --- Left panel: fit overlay ---
        ax_fits = fig.add_subplot(121)

        # Plot each synthetic spectrum and its fit + individual peaks
        for i, (syn, fit_curve) in enumerate(zip(all_synthetic, all_fit_curves)):
            ax_fits.plot(energy, syn, '-', color='#BBDEFB', linewidth=0.5,
                         alpha=0.4)
            if fit_curve is not None:
                ax_fits.plot(energy, fit_curve, '-', color='#E57373',
                             linewidth=0.5, alpha=0.4)
            # Individual peak curves
            if i < len(all_peak_curves) and all_peak_curves[i] is not None:
                for pname, curve in all_peak_curves[i].items():
                    ax_fits.plot(energy, curve, '-',
                                 color=peak_color_map[pname],
                                 linewidth=0.4, alpha=0.25)

        # Overlay the original best-fit prominently
        ax_fits.plot(energy, best_fit_curve, '-', color='#B71C1C', linewidth=2,
                     label='Original best fit')

        # Overlay original peaks prominently
        if self.current_result is not None:
            for pi, pname in enabled_peaks:
                center = self.current_result.params[f'c{pi}'].value
                height = self.current_result.params[f'h{pi}'].value
                if pi <= 6:
                    fwhm = self.current_result.params['main_fwhm'].value
                elif pi == 7:
                    fwhm = self.current_result.params['carb_fwhm'].value
                else:
                    fwhm = self.current_result.params['sigma_fwhm'].value
                peak_curve = gaussian(energy, center, height, fwhm)
                didx = self.peak_names.index(pname)
                ax_fits.plot(energy, peak_curve, '--',
                             color=peak_color_map[pname], linewidth=1.5,
                             alpha=0.8, label=self.peak_display_names[didx])

        ax_fits.plot(energy, self.current_intensity, 'o', color='#0D47A1',
                     markersize=1.5, alpha=0.5, label='Measured data')

        ax_fits.set_xlabel('Energy (eV)', fontsize=10)
        ax_fits.set_ylabel('Intensity (a.u.)', fontsize=10)
        ax_fits.set_title(f"Monte Carlo Fits (N={len(all_synthetic)})",
                          fontsize=11, fontweight='bold')
        ax_fits.legend(fontsize=7, loc='upper left')
        ax_fits.grid(True, alpha=0.3)

        # --- Right panel: box plots ---
        ax_box = fig.add_subplot(122)

        bp = ax_box.boxplot(area_data, patch_artist=True, widths=0.5,
                            medianprops=dict(color='black', linewidth=1.5),
                            showmeans=True,
                            meanprops=dict(marker='D', markerfacecolor='blue', markersize=5))

        for i, (patch, cv) in enumerate(zip(bp['boxes'], cvs)):
            if cv < 10:
                color = '#4CAF50'
            elif cv < 20:
                color = '#FFC107'
            else:
                color = '#F44336'
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax_box.set_xticklabels(peak_names_display, fontsize=9, rotation=30,
                               ha='right')
        ax_box.set_ylabel('Peak Area (a.u.)', fontsize=10)
        ax_box.set_title("Peak Area Distributions", fontsize=11, fontweight='bold')

        # CV annotations
        for i, (cv, name) in enumerate(zip(cvs, peak_names_display)):
            y_top = np.max(area_data[i])
            ax_box.text(i + 1, y_top, f"CV={cv:.0f}%", ha='center', va='bottom',
                        fontsize=8, fontweight='bold')

        # Verdict
        max_cv = max(cvs) if cvs else 0
        if max_cv < 10:
            verdict = "ROBUST - parameters well-constrained by data"
            verdict_color = '#4CAF50'
        elif max_cv < 20:
            verdict = "CAUTION - some noise sensitivity"
            verdict_color = '#FFC107'
        else:
            verdict = "SENSITIVE - high parameter uncertainty from noise"
            verdict_color = '#F44336'

        ax_box.text(0.5, 0.02, verdict, transform=ax_box.transAxes,
                    fontsize=10, fontweight='bold', color=verdict_color,
                    ha='center', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                              edgecolor=verdict_color, alpha=0.9))

        # Info text box
        r_sq_arr = np.array(all_r_squared)
        r_sq_valid = r_sq_arr[~np.isnan(r_sq_arr)]

        info_lines = [f"Noise \u03c3: {noise_level:.4f}"]
        if len(r_sq_valid) > 0:
            info_lines.append(f"Mean R\u00b2: {np.mean(r_sq_valid):.4f}")
        info_lines.append("")
        info_lines.append("Peak uncertainty (\u03c3 of area):")
        for name, sd in zip(peak_names_display, stds):
            info_lines.append(f"  {name}: {sd:.3f}")

        ax_box.text(0.98, 0.98, '\n'.join(info_lines),
                    transform=ax_box.transAxes, fontsize=7, va='top', ha='right',
                    family='monospace',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                              alpha=0.8))

        # Legend
        legend_elements = [
            Patch(facecolor='#4CAF50', alpha=0.7, label='CV < 10%: Robust'),
            Patch(facecolor='#FFC107', alpha=0.7, label='CV 10-20%: Caution'),
            Patch(facecolor='#F44336', alpha=0.7, label='CV > 20%: Sensitive'),
        ]
        ax_box.legend(handles=legend_elements, loc='upper left', fontsize=8)

        ax_box.grid(True, alpha=0.3, axis='y')

        fig.suptitle(f"Monte Carlo Analysis: {sample_name}", fontsize=13,
                     fontweight='bold')
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        canvas = FigureCanvasTkAgg(fig, master=result_win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, result_win)
        toolbar.update()

    def plot_multiple_samples(self, selection_indices):
        """Plot multiple selected samples overlaid on the same plot."""
        self.ax1.clear()
        self.ax2.clear()

        # Color palette for multiple samples
        colors = plt.cm.tab10(np.linspace(0, 0.9, len(selection_indices)))

        all_residuals = []
        sample_names = []

        for idx, color in zip(selection_indices, colors):
            selected_name = self.sample_listbox.get(idx)

            # Find the file
            for i, file_path in enumerate(self.spectra_files):
                if file_path.stem == selected_name:
                    try:
                        energy, intensity = load_spectrum(file_path)

                        # Use custom centers if enabled
                        if self.custom_peak_centers.get():
                            result = self.fit_spectrum_with_custom_centers(energy, intensity)
                        else:
                            result = fit_spectrum(energy, intensity)

                        # Plot data and fit
                        self.ax1.plot(energy, intensity, 'o', color=color, markersize=2,
                                     alpha=0.4, label=f'{selected_name} (data)')
                        self.ax1.plot(energy, result.best_fit, '-', color=color,
                                     linewidth=2, label=f'{selected_name} (fit)')

                        # Collect residuals
                        all_residuals.append((energy, result.residual, color, selected_name))
                        sample_names.append(selected_name)

                    except Exception as e:
                        print(f"Error loading {selected_name}: {e}")
                    break

        # Plot residuals for all samples
        for energy, residual, color, name in all_residuals:
            self.ax2.plot(energy, residual, '-', color=color, linewidth=1,
                         alpha=0.7, label=name)

        self.ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)

        # Formatting
        self.ax1.set_xlabel('Energy (eV)', fontsize=11)
        self.ax1.set_ylabel('Normalized Absorption', fontsize=11)
        self.ax1.set_title(f'C1s XAS: Multiple Samples ({len(selection_indices)} selected)',
                          fontsize=12, fontweight='bold')
        self.ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        self.ax1.grid(True, alpha=0.3)
        self.ax1.set_xlim(280, 294)

        self.ax2.set_xlabel('Energy (eV)', fontsize=11)
        self.ax2.set_ylabel('Residuals', fontsize=11)
        self.ax2.set_title('Fit Residuals', fontsize=11, fontweight='bold')
        self.ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        self.ax2.grid(True, alpha=0.3)
        self.ax2.set_xlim(280, 294)

        self.fig.tight_layout()
        self.canvas.draw()

        # Disable statistics update for multi-select
        self.stats_text.config(state=tk.NORMAL)
        self.stats_text.delete('1.0', tk.END)
        self.stats_text.insert('1.0', f"Multi-select mode\n{len(selection_indices)} samples\n\n" +
                               "Select a single sample\nto view statistics")
        self.stats_text.config(state=tk.DISABLED)

    def load_current_spectrum(self):
        """Load and fit current spectrum, using cached state when available."""
        if not self.spectra_files:
            return

        file_path = self.spectra_files[self.current_index]

        try:
            if self._restore_file_state():
                # Cache hit — skip loading/fitting
                pass
            else:
                # Cache miss — load data, fit, then cache
                self.current_energy, self.current_intensity = load_spectrum(file_path)

                # Check for a saved fit-state sidecar file
                sidecar = file_path.parent / f"{file_path.stem}_fit.json"
                if sidecar.exists():
                    load_it = messagebox.askyesno(
                        "Saved Fit Found",
                        f"A saved fit state was found:\n{sidecar.name}\n\n"
                        "Load it instead of refitting?")
                    if load_it:
                        self._apply_fit_state(str(sidecar))
                        return
                # Fit using refit_with_enabled_peaks which respects peak_in_fit,
                # custom centers, and manual baseline settings
                self.refit_with_enabled_peaks()
                self._save_file_state()

            self.update_peak_tree_display()
            self.update_plot()
            self.update_statistics()
            self.update_peak_parameters()
            self.update_statistics_tab()

            # Update listbox selection
            self.sample_listbox.selection_clear(0, tk.END)
            for i in range(self.sample_listbox.size()):
                if self.sample_listbox.get(i) == file_path.stem:
                    self.sample_listbox.selection_set(i)
                    self.sample_listbox.see(i)
                    break

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load spectrum:\n{str(e)}")

    def update_plot(self):
        """Update the plot."""
        if self.current_result is None:
            return

        self.ax1.clear()
        self.ax2.clear()

        energy = self.current_energy
        intensity = self.current_intensity
        result = self.current_result

        # Get baseline
        if self.manual_baseline.get():
            baseline = double_arctangent(energy,
                                        self.baseline_params['arc1_center'].get(),
                                        self.baseline_params['arc1_height'].get(),
                                        self.baseline_params['arc1_width'].get(),
                                        self.baseline_params['arc2_center'].get(),
                                        self.baseline_params['arc2_height'].get(),
                                        self.baseline_params['arc2_width'].get())
        else:
            baseline = double_arctangent(energy,
                                        result.params['arc1_center'].value,
                                        result.params['arc1_height'].value,
                                        result.params['arc1_width'].value,
                                        result.params['arc2_center'].value,
                                        result.params['arc2_height'].value,
                                        result.params['arc2_width'].value)

        # Determine what to plot
        if self.show_baseline.get():
            # Show everything with baseline
            self.ax1.plot(energy, intensity, 'ko', markersize=3, alpha=0.6, label='Data')
            self.ax1.plot(energy, result.best_fit, 'r-', linewidth=2, label='Total Fit')
            self.ax1.plot(energy, baseline, 'k--', linewidth=1.5, alpha=0.5, label='Baseline')
            y_offset = baseline
            ylabel = 'Normalized Absorption'
        else:
            # Show baseline-subtracted but keep fit
            self.ax1.plot(energy, intensity - baseline, 'ko', markersize=3, alpha=0.6, label='Data - Baseline')
            self.ax1.plot(energy, result.best_fit - baseline, 'r-', linewidth=2, label='Fit - Baseline')
            y_offset = np.zeros_like(baseline)
            ylabel = 'Baseline-Subtracted Absorption'

        # Plot peaks
        if self.show_peaks.get():
            for i, (name, display_name, color) in enumerate(zip(self.peak_names, self.peak_display_names, self.peak_colors), 1):
                if not self.peak_in_fit[name].get():
                    continue

                center = result.params[f'c{i}'].value
                height = result.params[f'h{i}'].value

                if i <= 6:
                    fwhm = result.params['main_fwhm'].value
                elif i == 7:
                    fwhm = result.params['carb_fwhm'].value
                else:
                    fwhm = result.params['sigma_fwhm'].value

                peak = gaussian(energy, center, height, fwhm)
                self.ax1.plot(energy, peak + y_offset, '--', color=color,
                             linewidth=1.5, alpha=0.7, label=display_name)

        self.ax1.set_xlabel('Energy (eV)', fontsize=11)
        self.ax1.set_ylabel(ylabel, fontsize=11)

        title = f'C1s XAS: {self.spectra_files[self.current_index].stem}'
        if validate_baseline and not self.manual_baseline.get():
            is_valid, _ = validate_baseline(energy, result, intensity)
            title += " [OK]" if is_valid else " [!]"

        self.ax1.set_title(title, fontsize=12, fontweight='bold')
        self.ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        self.ax1.grid(True, alpha=0.3)
        self.ax1.set_xlim(energy.min(), 294)

        # Residuals
        r_squared = 1 - result.residual.var() / np.var(result.data)
        self.ax2.plot(energy, result.residual, 'b-', linewidth=1)
        self.ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
        self.ax2.set_xlabel('Energy (eV)', fontsize=11)
        self.ax2.set_ylabel('Residuals', fontsize=11)
        self.ax2.set_title(f'Residuals (R² = {r_squared:.4f})', fontsize=11, fontweight='bold')
        self.ax2.grid(True, alpha=0.3)
        self.ax2.set_xlim(energy.min(), 294)

        self.fig.tight_layout()
        self.canvas.draw()

    @staticmethod
    def _abbe_criterion(residuals):
        """Compute the Abbe criterion for residual autocorrelation.

        Returns a value near 1.0 for randomly distributed residuals (good).
        Values significantly < 1 indicate structured/autocorrelated residuals.
        """
        n = len(residuals)
        if n < 3:
            return 1.0
        ss_diff = np.sum(np.diff(residuals) ** 2)
        ss_total = np.sum((residuals - np.mean(residuals)) ** 2)
        if ss_total == 0:
            return 1.0
        return (ss_diff / (2 * (n - 1))) / (ss_total / (n - 1))

    @staticmethod
    def _compute_aic_bic(residuals, n_free):
        """Compute RSS-based AIC and BIC (consistent with model selection)."""
        n = len(residuals)
        ss_res = np.sum(residuals ** 2)
        aic = n * np.log(ss_res / n) + 2 * n_free
        bic = n * np.log(ss_res / n) + n_free * np.log(n)
        return aic, bic

    def update_statistics(self):
        """Update fit statistics."""
        if self.current_result is None:
            return

        result = self.current_result
        r_squared = 1 - result.residual.var() / np.var(result.data)
        resid_std = np.std(result.residual)
        abbe = self._abbe_criterion(result.residual)
        n_free = sum(1 for p in result.params.values()
                     if p.vary and p.expr is None)
        aic, bic = self._compute_aic_bic(result.residual, n_free)

        stats_text = f"R² = {r_squared:.4f}\n"
        stats_text += f"χ²ᵣ = {result.redchi:.2e}\n"
        stats_text += f"σᵣ  = {resid_std:.2e}\n"
        stats_text += f"AIC = {aic:.1f}\n"
        stats_text += f"BIC = {bic:.1f}\n"

        # Abbe criterion: ~1.0 = random residuals, <<1 = structured
        abbe_flag = "[OK]" if abbe > 0.5 else "[!]"
        stats_text += f"Abbe = {abbe:.3f} {abbe_flag}\n"
        n_peaks = sum(1 for v in self.peak_in_fit.values() if v.get())
        stats_text += f"Peaks = {n_peaks}\n"
        stats_text += f"Points = {len(self.current_energy)}\n"

        if self.manual_baseline.get():
            stats_text += f"\nBaseline (manual):\n"
            stats_text += f"S1: {self.baseline_params['arc1_center'].get():.2f} eV\n"
            stats_text += f"    H={self.baseline_params['arc1_height'].get():.3f}\n"
            stats_text += f"S2: {self.baseline_params['arc2_center'].get():.2f} eV\n"
            stats_text += f"    H={self.baseline_params['arc2_height'].get():.3f}"
        else:
            stats_text += f"\nBaseline:\n"
            stats_text += f"S1: {result.params['arc1_center'].value:.2f} eV\n"
            stats_text += f"    H={result.params['arc1_height'].value:.3f}\n"
            stats_text += f"S2: {result.params['arc2_center'].value:.2f} eV\n"
            stats_text += f"    H={result.params['arc2_height'].value:.3f}"

        self.stats_text.config(state=tk.NORMAL)
        self.stats_text.delete('1.0', tk.END)
        self.stats_text.insert('1.0', stats_text)
        self.stats_text.config(state=tk.DISABLED)

    def update_peak_parameters(self):
        """Update peak parameters table using Treeview."""
        if self.current_result is None:
            return

        # Clear existing tree items
        for item in self.peak_tree.get_children():
            self.peak_tree.delete(item)
        self.peak_tree_items.clear()

        try:
            # Calculate areas
            areas = calculate_peak_areas(self.current_energy, self.current_result)

            # Validate that all expected peaks are present
            for expected_name in self.peak_names:
                if expected_name not in areas:
                    areas[expected_name] = 0.0

            # Calculate percentages for peaks in QUANTITATION only
            selected_areas = {k: v for k, v in areas.items() if self.peak_in_quant[k].get()}
            total_selected_area = sum(selected_areas.values())

            if total_selected_area == 0:
                percentages = {k: 0.0 for k in areas.keys()}
            else:
                percentages = {k: v/total_selected_area*100 for k, v in selected_areas.items()}
                # Non-selected peaks get 0%
                for k in areas.keys():
                    if k not in percentages:
                        percentages[k] = 0.0

        except Exception as e:
            print(f"Error calculating peak areas: {e}")
            areas = {name: 0.0 for name in self.peak_names}
            percentages = {name: 0.0 for name in self.peak_names}

        result = self.current_result

        # Add peak rows to Treeview
        for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
            # Get checkbox states
            in_fit = self.peak_in_fit[name].get()
            fit_char = 'Y' if in_fit else '-'
            quant_char = 'Y' if self.peak_in_quant[name].get() else '-'

            if in_fit:
                # Active peak: show full parameters
                center = result.params[f'c{i}'].value
                if i <= 6:
                    fwhm = result.params['main_fwhm'].value
                elif i == 7:
                    fwhm = result.params['carb_fwhm'].value
                else:
                    fwhm = result.params['sigma_fwhm'].value
                pct = percentages[name]
                values = (display_name, f"{center:.2f}", f"{fwhm:.2f}", f"{pct:.1f}")
            else:
                # Excluded peak: show dashes
                values = (display_name, "---", "---", "---")

            # Insert into tree
            item_id = self.peak_tree.insert('', 'end',
                                            text=f"{fit_char} {quant_char}",
                                            values=values,
                                            tags=(name,))

            self.peak_tree_items[name] = item_id

            # Apply color tag: full color if in fit, gray if excluded
            if in_fit:
                color_hex = matplotlib.colors.rgb2hex(self.peak_colors[i-1])
            else:
                color_hex = '#999999'
            self.peak_tree.tag_configure(name, foreground=color_hex)

    def on_fit_checkbox_changed(self):
        """Handle fit checkbox change - requires refitting."""
        if self._fit_in_progress:
            return  # Ignore if already fitting

        def do_refit():
            self._fit_in_progress = True
            try:
                # Run the actual fitting
                self.refit_with_enabled_peaks()
                # Schedule UI updates on main thread
                self.root.after(0, self._finish_refit)
            except Exception as e:
                print(f"Refit error: {e}")
                self._fit_in_progress = False

        # Show busy cursor
        self.root.config(cursor="wait")
        self.root.update()

        # Run fit in background thread
        thread = threading.Thread(target=do_refit, daemon=True)
        thread.start()

    def _finish_refit(self):
        """Called on main thread after refit completes."""
        self._fit_in_progress = False
        self.root.config(cursor="")
        self.update_peak_parameters()
        self.update_plot()
        self.update_statistics()
        self.update_statistics_tab()
        self._save_file_state()

    def on_quant_checkbox_changed(self):
        """Handle quantitation checkbox change - just recalculate percentages, no refit."""
        self.update_peak_parameters()
        self.update_statistics_tab()
        self._save_file_state()

    def update_statistics_tab(self):
        """Update statistics tab for current sample."""
        if self.current_result is None:
            return

        try:
            areas = calculate_peak_areas(self.current_energy, self.current_result)
            # Use peak_in_quant for quantitation (not peak_in_fit)
            selected_areas = {k: v for k, v in areas.items() if self.peak_in_quant[k].get()}
            total_area = sum(selected_areas.values())

            if total_area > 0:
                percentages = {k: v/total_area*100 for k, v in selected_areas.items()}
            else:
                percentages = {k: 0.0 for k in selected_areas.keys()}

            # Update summary statistics
            self.stat_labels['aromaticity'].config(text=f"{percentages.get('Aromatic', 0.0):.1f}%")
            self.stat_labels['aliphatic'].config(text=f"{percentages.get('Aliphatic', 0.0):.1f}%")
            self.stat_labels['carboxyl'].config(text=f"{percentages.get('Carboxyl', 0.0):.1f}%")
            self.stat_labels['oalkyl'].config(text=f"{percentages.get('O_alkyl', 0.0):.1f}%")
            self.stat_labels['carbonate'].config(text=f"{percentages.get('Carbonate', 0.0):.1f}%")

            num_selected = sum(1 for enabled in self.peak_in_quant.values() if enabled.get())
            self.stat_labels['total'].config(text=f"{num_selected} peaks")

            # Update visualizations
            self.update_statistics_plots(percentages, selected_areas)

        except Exception as e:
            print(f"Error updating statistics tab: {e}")

    def update_statistics_plots(self, percentages, selected_areas):
        """Update statistics visualization plots with publication quality."""
        self.stats_fig.clear()

        selected_names = [self.peak_display_names[i] for i, name in enumerate(self.peak_names)
                         if self.peak_in_quant[name].get() and percentages.get(name, 0) > 0]
        selected_values = [percentages[name] for name in self.peak_names
                          if self.peak_in_quant[name].get() and percentages.get(name, 0) > 0]
        selected_colors = [self.peak_colors[i] for i, name in enumerate(self.peak_names)
                          if self.peak_in_quant[name].get() and percentages.get(name, 0) > 0]

        if not selected_values:
            ax = self.stats_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'No peaks selected for quantification',
                   ha='center', va='center', fontsize=14)
            ax.axis('off')
            self.stats_canvas.draw()
            return

        # Publication quality settings
        plt.rcParams.update({
            'font.family': 'Arial',
            'font.size': 12,
            'axes.linewidth': 1.5,
            'axes.labelweight': 'bold',
            'xtick.major.width': 1.5,
            'ytick.major.width': 1.5,
        })

        # Pie chart
        ax1 = self.stats_fig.add_subplot(121)
        wedges, texts, autotexts = ax1.pie(selected_values, labels=selected_names,
                                           colors=selected_colors, autopct='%1.1f%%',
                                           startangle=90, textprops={'fontsize': 11})
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(11)
        for text in texts:
            text.set_fontsize(11)
            text.set_fontweight('bold')
        ax1.set_title('Peak Distribution', fontsize=14, fontweight='bold', pad=15)

        # Bar chart
        ax2 = self.stats_fig.add_subplot(122)
        x_pos = np.arange(len(selected_names))
        bars = ax2.bar(x_pos, selected_values, color=selected_colors, alpha=0.85,
                       edgecolor='black', linewidth=1.2)
        ax2.set_ylabel('Percentage (%)', fontsize=13, fontweight='bold')
        ax2.set_title('Peak Percentages', fontsize=14, fontweight='bold', pad=15)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(selected_names, rotation=45, ha='right', fontsize=11, fontweight='bold')
        ax2.tick_params(axis='y', labelsize=11)
        ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax2.set_axisbelow(True)  # Grid behind bars
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        for bar, value in zip(bars, selected_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{value:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        self.stats_fig.tight_layout(pad=2.0)
        self.stats_canvas.draw()

    def filter_stats_samples(self, *args):
        """Filter stats samples list based on search."""
        if not hasattr(self, 'stats_samples_listbox'):
            return
        search_term = self.stats_search_var.get().lower()
        self.stats_samples_listbox.delete(0, tk.END)

        for file_path in self.spectra_files:
            name = file_path.stem
            if search_term == '' or search_term in name.lower():
                self.stats_samples_listbox.insert(tk.END, name)

    def generate_multi_sample_stats(self):
        """Generate statistics for multiple selected samples."""
        selections = self.stats_samples_listbox.curselection()
        if not selections:
            messagebox.showwarning("Warning", "No samples selected. Please select samples from the list.")
            return

        # Collect data from selected samples
        all_percentages = {}
        sample_names = []

        for idx in selections:
            sample_name = self.stats_samples_listbox.get(idx)
            sample_names.append(sample_name)

            # Find the sample file
            sample_path = None
            for f in self.spectra_files:
                if f.stem == sample_name:
                    sample_path = f
                    break

            if sample_path:
                try:
                    energy, intensity = load_spectrum(sample_path)
                    result = self.fit_spectrum_with_custom_centers(energy, intensity) if self.custom_peak_centers.get() else fit_spectrum(energy, intensity)
                    areas = calculate_peak_areas(energy, result)
                    # Use peak_in_quant for quantitation
                    selected_areas = {k: v for k, v in areas.items() if self.peak_in_quant[k].get()}
                    total = sum(selected_areas.values())
                    if total > 0:
                        all_percentages[sample_name] = {k: v/total*100 for k, v in selected_areas.items()}
                    else:
                        all_percentages[sample_name] = {k: 0.0 for k in selected_areas.keys()}
                except Exception as e:
                    print(f"Error processing {sample_name}: {e}")

        if not all_percentages:
            messagebox.showerror("Error", "Could not process any samples.")
            return

        # Update status label
        self.stats_status_label.config(text=f"Showing statistics for {len(sample_names)} samples")

        # Calculate average statistics
        avg_percentages = {}
        for peak_name in self.peak_names:
            if self.peak_in_quant[peak_name].get():
                values = [all_percentages[s].get(peak_name, 0) for s in sample_names]
                avg_percentages[peak_name] = np.mean(values)

        # Update summary labels with averages
        self.stat_labels['aromaticity'].config(text=f"{avg_percentages.get('Aromatic', 0.0):.1f}% (avg)")
        self.stat_labels['aliphatic'].config(text=f"{avg_percentages.get('Aliphatic', 0.0):.1f}% (avg)")
        self.stat_labels['carboxyl'].config(text=f"{avg_percentages.get('Carboxyl', 0.0):.1f}% (avg)")
        self.stat_labels['oalkyl'].config(text=f"{avg_percentages.get('O_alkyl', 0.0):.1f}% (avg)")
        self.stat_labels['carbonate'].config(text=f"{avg_percentages.get('Carbonate', 0.0):.1f}% (avg)")
        self.stat_labels['total'].config(text=f"{len(sample_names)} samples")

        # Generate publication quality plots
        self.generate_multi_sample_plots(all_percentages, sample_names)

    def generate_multi_sample_plots(self, all_percentages, sample_names):
        """Generate publication quality plots for multiple samples."""
        self.stats_fig.clear()

        # Publication quality settings
        plt.rcParams.update({
            'font.family': 'Arial',
            'font.size': 11,
            'axes.linewidth': 1.5,
            'axes.labelweight': 'bold',
            'xtick.major.width': 1.5,
            'ytick.major.width': 1.5,
        })

        # Get enabled peaks
        enabled_peaks = [name for name in self.peak_names if self.peak_in_quant[name].get()]
        enabled_display = [self.peak_display_names[i] for i, name in enumerate(self.peak_names)
                          if self.peak_in_quant[name].get()]
        enabled_colors = [self.peak_colors[i] for i, name in enumerate(self.peak_names)
                         if self.peak_in_quant[name].get()]

        if len(sample_names) == 1:
            # Single sample - pie and bar
            percentages = all_percentages[sample_names[0]]
            selected_values = [percentages.get(name, 0) for name in enabled_peaks if percentages.get(name, 0) > 0]
            selected_names = [enabled_display[i] for i, name in enumerate(enabled_peaks) if percentages.get(name, 0) > 0]
            selected_colors = [enabled_colors[i] for i, name in enumerate(enabled_peaks) if percentages.get(name, 0) > 0]

            ax1 = self.stats_fig.add_subplot(121)
            wedges, texts, autotexts = ax1.pie(selected_values, labels=selected_names,
                                               colors=selected_colors, autopct='%1.1f%%',
                                               startangle=90, textprops={'fontsize': 11})
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
            ax1.set_title(f'{sample_names[0]}', fontsize=12, fontweight='bold')

            ax2 = self.stats_fig.add_subplot(122)
            x_pos = np.arange(len(selected_names))
            bars = ax2.bar(x_pos, selected_values, color=selected_colors, alpha=0.85,
                           edgecolor='black', linewidth=1.2)
            ax2.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
            ax2.set_xticks(x_pos)
            ax2.set_xticklabels(selected_names, rotation=45, ha='right', fontsize=10)
            ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

        else:
            # Multiple samples - grouped bar with error bars
            ax1 = self.stats_fig.add_subplot(121)

            # Calculate means and stds
            means = []
            stds = []
            for peak_name in enabled_peaks:
                values = [all_percentages[s].get(peak_name, 0) for s in sample_names]
                means.append(np.mean(values))
                stds.append(np.std(values))

            x_pos = np.arange(len(enabled_peaks))
            bars = ax1.bar(x_pos, means, yerr=stds, color=enabled_colors, alpha=0.85,
                           edgecolor='black', linewidth=1.2, capsize=5, error_kw={'linewidth': 1.5})
            ax1.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
            ax1.set_title(f'Average Peak Distribution (n={len(sample_names)})',
                         fontsize=13, fontweight='bold')
            ax1.set_xticks(x_pos)
            ax1.set_xticklabels(enabled_display, rotation=45, ha='right', fontsize=10, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

            # Grouped bar chart for each sample
            ax2 = self.stats_fig.add_subplot(122)
            x = np.arange(len(enabled_peaks))
            width = 0.8 / len(sample_names)
            colors = plt.cm.tab10(np.linspace(0, 1, len(sample_names)))

            for i, sample_name in enumerate(sample_names):
                vals = [all_percentages[sample_name].get(pname, 0) for pname in enabled_peaks]
                offset = (i - len(sample_names)/2 + 0.5) * width
                ax2.bar(x + offset, vals, width, label=sample_name[:20], alpha=0.85,
                       edgecolor='black', linewidth=0.8, color=colors[i])

            ax2.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
            ax2.set_title('Sample Comparison', fontsize=13, fontweight='bold')
            ax2.set_xticks(x)
            ax2.set_xticklabels(enabled_display, rotation=45, ha='right', fontsize=10, fontweight='bold')
            ax2.legend(fontsize=8, loc='upper right')
            ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

        self.stats_fig.tight_layout(pad=2.0)
        self.stats_canvas.draw()

    def _current_peak_sel(self):
        """Return a snapshot of current F/Q checkbox states."""
        return {
            'fit_peaks': {name: self.peak_in_fit[name].get() for name in self.peak_names},
            'quant_peaks': {name: self.peak_in_quant[name].get() for name in self.peak_names},
        }

    def add_selected_to_comparison(self):
        """Add selected samples from available list to comparison."""
        selections = self.available_samples_listbox.curselection()
        if not selections:
            messagebox.showwarning("Warning", "No samples selected")
            return

        added_count = 0
        for idx in selections:
            sample_name = self.available_samples_listbox.get(idx)

            # Check if already in comparison
            already_added = any(comp_name == sample_name for comp_name, _, _, _ in self.comparison_samples)
            if already_added:
                continue

            # Find the sample file
            sample_path = None
            for f in self.spectra_files:
                if f.stem == sample_name:
                    sample_path = f
                    break

            if sample_path:
                try:
                    # Use cached state if available, else fit fresh
                    cached = self.file_state_cache.get(str(sample_path))
                    if cached and cached['result'] is not None:
                        result = cached['result']
                        energy = cached['energy']
                        peak_sel = {
                            'fit_peaks': dict(cached['fit_peaks']),
                            'quant_peaks': dict(cached['quant_peaks']),
                        }
                    else:
                        energy, intensity = load_spectrum(sample_path)
                        result = self.fit_spectrum_with_custom_centers(energy, intensity) if self.custom_peak_centers.get() else fit_spectrum(energy, intensity)
                        peak_sel = self._current_peak_sel()

                    # Add to comparison
                    self.comparison_samples.append((sample_name, result, energy, peak_sel))
                    self.comparison_listbox.insert(tk.END, sample_name)
                    added_count += 1
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to load {sample_name}:\n{str(e)}")

        if added_count > 0:
            messagebox.showinfo("Success", f"Added {added_count} sample(s) to comparison")

    def remove_from_comparison(self):
        """Remove selected sample from comparison."""
        selection = self.comparison_listbox.curselection()
        if not selection:
            return

        index = selection[0]
        self.comparison_listbox.delete(index)
        del self.comparison_samples[index]

    def clear_comparison(self):
        """Clear all comparison samples."""
        self.comparison_listbox.delete(0, tk.END)
        self.comparison_samples = []

    def generate_comparison(self):
        """Generate multi-sample comparison using per-sample peak selections."""
        if len(self.comparison_samples) < 2:
            messagebox.showwarning("Warning", "Add at least 2 samples to compare")
            return

        self.comparison_fig.clear()

        # Union of all per-sample quant peaks for chart layout
        union_quant = set()
        for _name, _result, _energy, peak_sel in self.comparison_samples:
            for pn in self.peak_names:
                if peak_sel['quant_peaks'].get(pn, False):
                    union_quant.add(pn)

        enabled_peaks = [name for name in self.peak_names if name in union_quant]
        enabled_display = [self.peak_display_names[i] for i, name in enumerate(self.peak_names)
                          if name in union_quant]
        enabled_colors = [self.peak_colors[i] for i, name in enumerate(self.peak_names)
                         if name in union_quant]

        if not enabled_peaks:
            ax = self.comparison_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'No peaks enabled for quantification',
                   ha='center', va='center', fontsize=14)
            ax.axis('off')
            self.comparison_canvas.draw()
            return

        # Calculate areas per sample using each sample's own quant peak set
        all_percentages = {}
        for name, result, energy, peak_sel in self.comparison_samples:
            areas = calculate_peak_areas(energy, result)
            sample_quant = [pn for pn in self.peak_names if peak_sel['quant_peaks'].get(pn, False)]
            enabled_areas = {k: v for k, v in areas.items() if k in sample_quant}
            total = sum(enabled_areas.values())
            if total > 0:
                pcts = {k: v/total*100 for k, v in enabled_areas.items()}
            else:
                pcts = {}
            # Fill 0 for union peaks not in this sample's quant set
            all_percentages[name] = {k: pcts.get(k, 0.0) for k in enabled_peaks}

        sample_names = [name for name, _, _, _ in self.comparison_samples]

        # Plot 1: Grouped bar chart
        ax1 = self.comparison_fig.add_subplot(221)
        x = np.arange(len(enabled_peaks))
        width = 0.8 / len(self.comparison_samples)

        for i, (name, _, _, _) in enumerate(self.comparison_samples):
            vals = [all_percentages[name].get(pname, 0) for pname in enabled_peaks]
            offset = (i - len(self.comparison_samples)/2 + 0.5) * width
            ax1.bar(x + offset, vals, width, label=name, alpha=0.8)

        ax1.set_ylabel('Percentage (%)')
        ax1.set_title('Peak Distribution Comparison')
        ax1.set_xticks(x)
        ax1.set_xticklabels(enabled_display, rotation=45, ha='right')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3, axis='y')

        # Plot 2: Ratio heatmap (only enabled peaks)
        ax2 = self.comparison_fig.add_subplot(222)

        n_samples = len(self.comparison_samples)
        n_peaks = len(enabled_peaks)

        # Create ratio matrix: peaks x samples (relative to first sample)
        ratio_matrix = np.zeros((n_peaks, n_samples))
        reference_name = sample_names[0]

        for i, peak_name in enumerate(enabled_peaks):
            for j, sample_name in enumerate(sample_names):
                val_sample = all_percentages[sample_name].get(peak_name, 0)
                val_ref = all_percentages[reference_name].get(peak_name, 0)

                if val_ref > 0.1:
                    ratio_matrix[i, j] = val_sample / val_ref
                else:
                    ratio_matrix[i, j] = np.nan

        sns.heatmap(ratio_matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                   center=1.0, vmin=0, vmax=3,
                   xticklabels=sample_names,
                   yticklabels=enabled_display,
                   cbar_kws={'label': 'Ratio'},
                   ax=ax2)
        ax2.set_title(f'Peak Ratio (vs {reference_name})')

        # Plot 3: Stacked bar chart (only enabled peaks)
        ax3 = self.comparison_fig.add_subplot(223)

        bottoms = np.zeros(len(self.comparison_samples))
        for i, (pname, disp_name, color) in enumerate(zip(enabled_peaks, enabled_display, enabled_colors)):
            vals = [all_percentages[name].get(pname, 0) for name, _, _, _ in self.comparison_samples]
            ax3.bar(range(len(self.comparison_samples)), vals, bottom=bottoms,
                   label=disp_name, color=color, alpha=0.8)
            bottoms += vals

        ax3.set_ylabel('Percentage (%)')
        ax3.set_title('Stacked Peak Distribution')
        ax3.set_xticks(range(len(self.comparison_samples)))
        ax3.set_xticklabels(sample_names, rotation=45, ha='right')
        ax3.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)

        # Plot 4: Summary table as grid (only enabled peaks)
        ax4 = self.comparison_fig.add_subplot(224)
        ax4.axis('off')

        # Create a table with samples as columns and peaks as rows
        n_samples = len(sample_names)
        n_peaks_to_show = min(len(enabled_peaks), 6)  # Show up to 6 peaks

        # Build table data
        cell_text = []
        for pname, dname in zip(enabled_peaks[:n_peaks_to_show], enabled_display[:n_peaks_to_show]):
            row = [f"{all_percentages[name].get(pname, 0):.1f}%" for name in sample_names]
            cell_text.append(row)

        # Truncate sample names for display
        col_labels = [name[:12] + '...' if len(name) > 12 else name for name in sample_names]
        row_labels = enabled_display[:n_peaks_to_show]

        # Create table
        table = ax4.table(cellText=cell_text,
                         rowLabels=row_labels,
                         colLabels=col_labels,
                         cellLoc='center',
                         loc='center',
                         colColours=['lightblue'] * n_samples,
                         rowColours=[self.peak_colors[self.peak_names.index(p)] for p in enabled_peaks[:n_peaks_to_show]])

        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        ax4.set_title('Peak Percentages by Sample', fontsize=11, fontweight='bold')

        self.comparison_fig.tight_layout()
        self.comparison_canvas.draw()

    def previous_spectrum(self):
        """Load previous spectrum."""
        if not self.spectra_files:
            return
        self._save_file_state()
        self.current_index = (self.current_index - 1) % len(self.spectra_files)
        self.load_current_spectrum()

    def next_spectrum(self):
        """Load next spectrum."""
        if not self.spectra_files:
            return
        self._save_file_state()
        self.current_index = (self.current_index + 1) % len(self.spectra_files)
        self.load_current_spectrum()

    def refit_spectrum(self):
        """Refit current spectrum, respecting peak_in_fit checkboxes."""
        if self.current_energy is None:
            return

        try:
            self.refit_with_enabled_peaks()

            self.update_plot()
            self.update_statistics()
            self.update_peak_parameters()
            self.update_statistics_tab()
            self._save_file_state()
            messagebox.showinfo("Success", "Spectrum refitted")
        except Exception as e:
            messagebox.showerror("Error", f"Refit failed:\n{str(e)}")

    def export_plot(self):
        """Export current plot."""
        if self.current_result is None:
            messagebox.showwarning("Warning", "No fit to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("PNG files", "*.png")],
            initialfile=f"{self.spectra_files[self.current_index].stem}_fit.pdf"
        )

        if filename:
            try:
                self.fig.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot exported")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed:\n{str(e)}")

    def export_fit_data(self):
        """Export fit data for current spectrum as CSV."""
        if self.current_result is None:
            messagebox.showwarning("Warning", "No fit data to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile=f"{self.spectra_files[self.current_index].stem}_fit_data.csv"
        )

        if filename:
            try:
                result = self.current_result
                energy = self.current_energy

                # Calculate areas
                areas = calculate_peak_areas(energy, result)
                selected_areas = {k: v for k, v in areas.items() if self.peak_in_fit[k].get()}
                total = sum(selected_areas.values())

                # Build data structure
                export_data = []

                # Peak information
                for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
                    if self.peak_in_fit[name].get():
                        center = result.params[f'c{i}'].value
                        height = result.params[f'h{i}'].value

                        if i <= 6:
                            fwhm = result.params['main_fwhm'].value
                        elif i == 7:
                            fwhm = result.params['carb_fwhm'].value
                        else:
                            fwhm = result.params['sigma_fwhm'].value

                        area = selected_areas[name]
                        percentage = (area / total * 100) if total > 0 else 0

                        export_data.append({
                            'Peak': display_name,
                            'Center_eV': center,
                            'Height': height,
                            'FWHM_eV': fwhm,
                            'Area': area,
                            'Percentage': percentage
                        })

                # Create DataFrame and export
                df = pd.DataFrame(export_data)

                # Add metadata header
                with open(filename, 'w') as f:
                    f.write(f"# Sample: {self.spectra_files[self.current_index].stem}\n")
                    r_squared = 1 - result.residual.var() / np.var(result.data)
                    f.write(f"# R-squared: {r_squared:.6f}\n")
                    f.write(f"# Reduced chi-square: {result.redchi:.6e}\n")
                    f.write(f"# AIC: {result.aic:.2f}\n")
                    f.write(f"# Baseline arc1_center: {result.params['arc1_center'].value:.4f}\n")
                    f.write(f"# Baseline arc1_height: {result.params['arc1_height'].value:.4f}\n")
                    f.write(f"# Baseline arc2_center: {result.params['arc2_center'].value:.4f}\n")
                    f.write(f"# Baseline arc2_height: {result.params['arc2_height'].value:.4f}\n")
                    f.write("#\n")

                    # Write DataFrame
                    df.to_csv(f, index=False)

                messagebox.showinfo("Success", f"Fit data exported to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed:\n{str(e)}")

    def export_comparison_data(self):
        """Export comparison data as CSV."""
        if len(self.comparison_samples) < 2:
            messagebox.showwarning("Warning", "No comparison data available")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="comparison_data.csv"
        )

        if filename:
            try:
                # Calculate areas for all samples using per-sample quant peaks
                all_data = []
                for name, result, energy, peak_sel in self.comparison_samples:
                    areas = calculate_peak_areas(energy, result)
                    selected_areas = {k: v for k, v in areas.items()
                                      if peak_sel['quant_peaks'].get(k, False)}
                    total = sum(selected_areas.values())

                    row = {'Sample': name}
                    for peak_name in self.peak_names:
                        if peak_sel['quant_peaks'].get(peak_name, False):
                            area = selected_areas[peak_name]
                            percentage = (area / total * 100) if total > 0 else 0
                            row[f'{peak_name}_area'] = area
                            row[f'{peak_name}_percent'] = percentage

                    all_data.append(row)

                # Create DataFrame and export
                df = pd.DataFrame(all_data)
                df.to_csv(filename, index=False)
                messagebox.showinfo("Success", f"Comparison data exported to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed:\n{str(e)}")

    def export_complete_analysis(self):
        """Export complete analysis package for all selected samples.

        Includes: statistics plots, fit data, raw data, and metadata for each sample.
        """
        if len(self.comparison_samples) < 1:
            messagebox.showwarning("Warning", "Add at least 1 sample for analysis export")
            return

        # Ask for output directory
        output_dir = filedialog.askdirectory(
            title="Select Output Directory for Complete Analysis Export"
        )

        if not output_dir:
            return

        output_dir = Path(output_dir)
        from datetime import datetime
        import json
        from c1s_fitter_optimized import total_model, double_arctangent, gaussian

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        export_name = f"analysis_export_{timestamp}"
        export_dir = output_dir / export_name
        export_dir.mkdir(exist_ok=True)

        try:
            # 1. Export all three plot types
            plots_dir = export_dir / "plots"
            plots_dir.mkdir(exist_ok=True)

            # Store current figure state
            original_fig_state = self.analysis_fig.axes.copy() if self.analysis_fig.axes else []

            # Generate and export Statistics plot
            self.generate_multi_sample_stats()
            self.analysis_fig.savefig(plots_dir / "statistics_summary.pdf", dpi=300, bbox_inches='tight')
            self.analysis_fig.savefig(plots_dir / "statistics_summary.png", dpi=300, bbox_inches='tight')

            # Generate and export Comparison Charts
            self.generate_comparison()
            self.analysis_fig.savefig(plots_dir / "comparison_charts.pdf", dpi=300, bbox_inches='tight')
            self.analysis_fig.savefig(plots_dir / "comparison_charts.png", dpi=300, bbox_inches='tight')

            # Generate and export Spectral Overlay
            self.generate_spectral_overlay()
            self.analysis_fig.savefig(plots_dir / "spectral_overlay.pdf", dpi=300, bbox_inches='tight')
            self.analysis_fig.savefig(plots_dir / "spectral_overlay.png", dpi=300, bbox_inches='tight')

            # 2. Export comprehensive data for each sample
            samples_dir = export_dir / "sample_data"
            samples_dir.mkdir(exist_ok=True)

            summary_data = []
            all_peak_data = []

            for name, result, energy, peak_sel in self.comparison_samples:
                sample_dir = samples_dir / name
                sample_dir.mkdir(exist_ok=True)

                # Calculate areas using this sample's own quant peak set
                areas = calculate_peak_areas(energy, result)
                selected_areas_quant = {k: v for k, v in areas.items()
                                        if peak_sel['quant_peaks'].get(k, False)}
                total_quant = sum(selected_areas_quant.values())

                # Fit quality metrics
                r_squared = 1 - result.residual.var() / np.var(result.data)

                # Peak parameters for this sample
                peak_data = []
                for i, (peak_name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
                    center = result.params[f'c{i}'].value
                    height = result.params[f'h{i}'].value

                    if i <= 6:
                        fwhm = result.params['main_fwhm'].value
                    elif i == 7:
                        fwhm = result.params['carb_fwhm'].value
                    else:
                        fwhm = result.params['sigma_fwhm'].value

                    area = areas.get(peak_name, 0)
                    pct_quant = (selected_areas_quant.get(peak_name, 0) / total_quant * 100) if total_quant > 0 else 0
                    in_fit = peak_sel['fit_peaks'].get(peak_name, True)
                    in_quant = peak_sel['quant_peaks'].get(peak_name, True)

                    peak_row = {
                        'Sample': name,
                        'Peak': display_name,
                        'Internal_Name': peak_name,
                        'Center_eV': round(center, 4),
                        'Height': round(height, 6),
                        'FWHM_eV': round(fwhm, 4),
                        'Area': round(area, 6),
                        'Percentage': round(pct_quant, 2),
                        'In_Fit': in_fit,
                        'In_Quantification': in_quant
                    }
                    peak_data.append(peak_row)
                    all_peak_data.append(peak_row)

                # Export peak parameters CSV for this sample
                pd.DataFrame(peak_data).to_csv(sample_dir / f"{name}_peak_parameters.csv", index=False)

                # Export raw + fitted spectrum data
                # Safely get parameter values, excluding computed params
                excluded = {'baseline_total'}
                params_dict = {p: result.params[p].value for p in result.params if p not in excluded}

                baseline = double_arctangent(energy,
                                            params_dict['arc1_center'], params_dict['arc1_height'], params_dict['arc1_width'],
                                            params_dict['arc2_center'], params_dict['arc2_height'], params_dict['arc2_width'])

                total_fit = total_model(energy, **params_dict)

                # Get the original intensity from the result
                intensity = result.data
                residual = intensity - total_fit

                spectrum_data = {
                    'Energy_eV': energy,
                    'Raw_Intensity': intensity,
                    'Total_Fit': total_fit,
                    'Baseline': baseline,
                    'Residual': residual
                }

                # Add individual peak contributions
                for i, peak_name in enumerate(self.peak_names, 1):
                    if i <= 6:
                        fwhm = params_dict['main_fwhm']
                    elif i == 7:
                        fwhm = params_dict['carb_fwhm']
                    else:
                        fwhm = params_dict['sigma_fwhm']

                    center = params_dict[f'c{i}']
                    height = params_dict[f'h{i}']
                    peak = gaussian(energy, center, height, fwhm)
                    spectrum_data[f'Peak_{peak_name}'] = peak

                pd.DataFrame(spectrum_data).to_csv(sample_dir / f"{name}_spectrum_data.csv", index=False)

                # Export metadata JSON
                metadata = {
                    'sample_name': name,
                    'export_timestamp': datetime.now().isoformat(),
                    'fit_quality': {
                        'r_squared': round(r_squared, 6),
                        'reduced_chi_square': round(result.redchi, 6),
                        'aic': round(result.aic, 2),
                        'bic': round(result.bic, 2),
                        'n_data_points': len(energy),
                        'n_variables': result.nvarys
                    },
                    'baseline_parameters': {
                        'arc1_center': round(result.params['arc1_center'].value, 4),
                        'arc1_height': round(result.params['arc1_height'].value, 6),
                        'arc1_width': round(result.params['arc1_width'].value, 4),
                        'arc2_center': round(result.params['arc2_center'].value, 4),
                        'arc2_height': round(result.params['arc2_height'].value, 6),
                        'arc2_width': round(result.params['arc2_width'].value, 4)
                    },
                    'fwhm_parameters': {
                        'main_fwhm': round(result.params['main_fwhm'].value, 4),
                        'carbonate_fwhm': round(result.params['carb_fwhm'].value, 4),
                        'sigma_star_fwhm': round(result.params['sigma_fwhm'].value, 4)
                    },
                    'peak_summary': {
                        peak_name: {
                            'center': round(result.params[f'c{i}'].value, 4),
                            'height': round(result.params[f'h{i}'].value, 6),
                            'area': round(areas.get(peak_name, 0), 6),
                            'percentage': round((selected_areas_quant.get(peak_name, 0) / total_quant * 100) if total_quant > 0 else 0, 2)
                        }
                        for i, peak_name in enumerate(self.peak_names, 1)
                    },
                    'settings': {
                        'peaks_in_fit': [pn for pn in self.peak_names if peak_sel['fit_peaks'].get(pn, True)],
                        'peaks_in_quantification': [pn for pn in self.peak_names if peak_sel['quant_peaks'].get(pn, True)]
                    }
                }

                with open(sample_dir / f"{name}_metadata.json", 'w') as f:
                    json.dump(metadata, f, indent=2)

                # Build summary row
                summary_row = {
                    'Sample': name,
                    'R_squared': round(r_squared, 6),
                    'Reduced_Chi_Sq': round(result.redchi, 6),
                    'AIC': round(result.aic, 2),
                    'BIC': round(result.bic, 2),
                    'Main_FWHM': round(result.params['main_fwhm'].value, 4),
                    'Carb_FWHM': round(result.params['carb_fwhm'].value, 4),
                    'Sigma_FWHM': round(result.params['sigma_fwhm'].value, 4)
                }

                # Add percentages for this sample's enabled peaks
                for peak_name, display_name in zip(self.peak_names, self.peak_display_names):
                    if peak_sel['quant_peaks'].get(peak_name, False):
                        pct = (selected_areas_quant.get(peak_name, 0) / total_quant * 100) if total_quant > 0 else 0
                        summary_row[f'{display_name}_%'] = round(pct, 2)

                # Add peak centers
                for i, peak_name in enumerate(self.peak_names, 1):
                    summary_row[f'{peak_name}_center_eV'] = round(result.params[f'c{i}'].value, 4)

                summary_data.append(summary_row)

            # 3. Export combined summary CSV
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(export_dir / "analysis_summary.csv", index=False)

            # 4. Export all peak data combined
            all_peaks_df = pd.DataFrame(all_peak_data)
            all_peaks_df.to_csv(export_dir / "all_samples_peak_parameters.csv", index=False)

            # 5. Generate comprehensive report
            report_path = export_dir / "analysis_report.txt"
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write("=" * 70 + "\n")
                f.write("C1s XAS MULTI-SAMPLE ANALYSIS REPORT\n")
                f.write("=" * 70 + "\n\n")
                f.write(f"Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Number of Samples: {len(self.comparison_samples)}\n\n")

                f.write("-" * 70 + "\n")
                f.write("SAMPLES INCLUDED\n")
                f.write("-" * 70 + "\n")
                for i, (name, _, _, _) in enumerate(self.comparison_samples, 1):
                    f.write(f"  {i}. {name}\n")
                f.write("\n")

                f.write("-" * 70 + "\n")
                f.write("PEAKS USED FOR QUANTIFICATION (per sample)\n")
                f.write("-" * 70 + "\n")
                for sname, _r, _e, ps in self.comparison_samples:
                    quant_list = [dn for pn, dn in zip(self.peak_names, self.peak_display_names)
                                  if ps['quant_peaks'].get(pn, False)]
                    f.write(f"  {sname}: {', '.join(quant_list)}\n")
                f.write("\n")

                f.write("-" * 70 + "\n")
                f.write("FIT QUALITY SUMMARY\n")
                f.write("-" * 70 + "\n")
                f.write(f"{'Sample':<30} {'R²':>10} {'Red.χ²':>12} {'AIC':>10}\n")
                f.write("-" * 70 + "\n")
                for row in summary_data:
                    f.write(f"{row['Sample']:<30} {row['R_squared']:>10.6f} {row['Reduced_Chi_Sq']:>12.2e} {row['AIC']:>10.1f}\n")
                f.write("\n")

                # Average statistics
                f.write("-" * 70 + "\n")
                f.write("AVERAGE COMPOSITION (%)\n")
                f.write("-" * 70 + "\n")

                for peak_name, display_name in zip(self.peak_names, self.peak_display_names):
                    col_name = f'{display_name}_%'
                    if col_name in summary_df.columns:
                        avg = summary_df[col_name].mean()
                        std = summary_df[col_name].std() if len(summary_df) > 1 else 0.0
                        if pd.isna(std):
                            std = 0.0
                        f.write(f"  {display_name:<20}: {avg:>6.2f} +/- {std:>5.2f} %\n")
                f.write("\n")

                f.write("-" * 70 + "\n")
                f.write("EXPORTED FILES\n")
                f.write("-" * 70 + "\n")
                f.write("plots/\n")
                f.write("  - statistics_summary.pdf/png\n")
                f.write("  - comparison_charts.pdf/png\n")
                f.write("  - spectral_overlay.pdf/png\n")
                f.write("sample_data/\n")
                f.write("  - [sample_name]/\n")
                f.write("    - [sample_name]_peak_parameters.csv\n")
                f.write("    - [sample_name]_spectrum_data.csv\n")
                f.write("    - [sample_name]_metadata.json\n")
                f.write("analysis_summary.csv\n")
                f.write("all_samples_peak_parameters.csv\n")
                f.write("analysis_report.txt\n")

            # Show success message
            n_samples = len(self.comparison_samples)
            messagebox.showinfo("Export Complete",
                              f"Complete analysis exported to:\n{export_dir}\n\n"
                              f"Included:\n"
                              f"• 3 plot types (PDF & PNG)\n"
                              f"• {n_samples} sample data packages\n"
                              f"• Analysis summary CSV\n"
                              f"• Combined peak parameters CSV\n"
                              f"• Analysis report")

        except Exception as e:
            messagebox.showerror("Export Error", f"Export failed:\n{str(e)}")
            import traceback
            traceback.print_exc()

    def generate_spectral_overlay(self):
        """Generate spectral overlay plot with y-axis offset."""
        if len(self.comparison_samples) < 1:
            messagebox.showwarning("Warning", "Add at least 1 sample to generate overlay")
            return

        self.comparison_fig.clear()

        # Publication quality settings
        plt.rcParams.update({
            'font.family': 'Arial',
            'font.size': 11,
            'axes.linewidth': 1.5,
            'axes.labelweight': 'bold',
            'xtick.major.width': 1.5,
            'ytick.major.width': 1.5,
        })

        ax = self.comparison_fig.add_subplot(111)

        y_offset_factor = self.y_offset_var.get()
        show_peaks = self.show_overlay_peaks.get()
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.comparison_samples)))

        # Store intensity data for each sample
        sample_data = []
        for i, (name, result, energy, _peak_sel) in enumerate(self.comparison_samples):
            # Find the corresponding file to get original intensity
            sample_path = None
            for f in self.spectra_files:
                if f.stem == name:
                    sample_path = f
                    break

            if sample_path:
                try:
                    orig_energy, orig_intensity = load_spectrum(sample_path)
                    sample_data.append((name, result, orig_energy, orig_intensity))
                except Exception as e:
                    print(f"Error loading {name}: {e}")
                    # Use model data as fallback
                    from c1s_fitter_optimized import total_model
                    excluded = {'baseline_total'}
                    params = {p: result.params[p].value for p in result.params if p not in excluded}
                    model_intensity = total_model(energy, **params)
                    sample_data.append((name, result, energy, model_intensity))
            else:
                # Use model data
                from c1s_fitter_optimized import total_model
                excluded = {'baseline_total'}
                params = {p: result.params[p].value for p in result.params if p not in excluded}
                model_intensity = total_model(energy, **params)
                sample_data.append((name, result, energy, model_intensity))

        # Calculate data range for proper offset scaling
        all_intensities = [intensity for _, _, _, intensity in sample_data]
        if all_intensities:
            max_intensity = max(np.max(i) for i in all_intensities)
            min_intensity = min(np.min(i) for i in all_intensities)
            data_range = max_intensity - min_intensity
        else:
            data_range = 1.0

        legend_handles = []
        legend_labels = []

        for i, (name, result, energy, intensity) in enumerate(sample_data):
            # Scale offset by data range so spectra don't overlap
            offset = i * y_offset_factor * data_range

            # Plot original data
            line, = ax.plot(energy, intensity + offset, 'o', markersize=3,
                           color=colors[i], alpha=0.5, markeredgewidth=0)

            # Calculate total fit
            from c1s_fitter_optimized import total_model, double_arctangent, gaussian
            # Filter out computed parameters like baseline_total that aren't function arguments
            excluded_params = {'baseline_total'}
            params_dict = {p: result.params[p].value for p in result.params if p not in excluded_params}
            total_fit = total_model(energy, **params_dict)

            # Plot total fit
            fit_line, = ax.plot(energy, total_fit + offset, '-', linewidth=2,
                               color=colors[i], label=name)
            legend_handles.append(fit_line)
            legend_labels.append(name)

            # Plot individual peaks if enabled
            if show_peaks:
                # Baseline
                baseline = double_arctangent(energy,
                                            params_dict['arc1_center'], params_dict['arc1_height'], params_dict['arc1_width'],
                                            params_dict['arc2_center'], params_dict['arc2_height'], params_dict['arc2_width'])

                # Individual peaks
                for j, peak_name in enumerate(self.peak_names):
                    if not self.peak_in_fit[peak_name].get():
                        continue

                    if peak_name == 'Carbonate':
                        fwhm = params_dict.get('carb_fwhm', 0.8)
                    elif peak_name == 'Sigma_star':
                        fwhm = params_dict.get('sigma_fwhm', 0.6)
                    else:
                        fwhm = params_dict.get('main_fwhm', 1.5)

                    center = params_dict.get(f'c{j+1}', self.default_peak_centers[peak_name])
                    height = params_dict.get(f'h{j+1}', 0)

                    if height > 0.01:  # Only plot significant peaks
                        peak = gaussian(energy, center, height, fwhm)
                        ax.fill_between(energy, baseline + offset, baseline + peak + offset,
                                        color=self.peak_colors[j], alpha=0.3, linewidth=0)

            # Add sample label at the right side
            ax.text(energy.max() + 0.2, offset + intensity.max() * 0.5,
                   name, fontsize=9, va='center', color=colors[i], fontweight='bold')

        # Axis labels and styling
        ax.set_xlabel('Energy (eV)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Normalized Intensity (offset)', fontsize=13, fontweight='bold')
        ax.set_title('Spectral Overlay Comparison', fontsize=14, fontweight='bold', pad=15)
        ax.tick_params(axis='both', labelsize=11)

        # Clean up axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Legend
        if legend_handles:
            ax.legend(legend_handles, legend_labels, loc='upper left', fontsize=9)

        # Set x-axis limits
        all_energies = [e for _, _, e, _ in sample_data]
        if all_energies:
            e_min = min(e.min() for e in all_energies)
            e_max = max(e.max() for e in all_energies)
            ax.set_xlim(e_min - 0.5, e_max + 2)  # Extra space for labels

        self.comparison_fig.tight_layout(pad=2.0)
        self.comparison_canvas.draw()

    def export_all_data(self):
        """Export complete data package for current sample (plot, data, metadata)."""
        if self.current_result is None:
            messagebox.showwarning("Warning", "No fit data to export")
            return

        # Ask for output directory
        output_dir = filedialog.askdirectory(
            title="Select Output Directory for Export"
        )

        if not output_dir:
            return

        output_dir = Path(output_dir)
        sample_name = self.spectra_files[self.current_index].stem

        try:
            result = self.current_result
            energy = self.current_energy
            intensity = self.current_intensity

            # 1. Export plot (PDF and PNG)
            plot_pdf = output_dir / f"{sample_name}_fit_plot.pdf"
            plot_png = output_dir / f"{sample_name}_fit_plot.png"
            self.fig.savefig(plot_pdf, dpi=300, bbox_inches='tight')
            self.fig.savefig(plot_png, dpi=300, bbox_inches='tight')

            # 2. Calculate areas and percentages
            areas = calculate_peak_areas(energy, result)
            selected_areas_quant = {k: v for k, v in areas.items() if self.peak_in_quant[k].get()}
            total_quant = sum(selected_areas_quant.values())

            # 3. Export peak parameters CSV
            peak_data = []
            for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
                center = result.params[f'c{i}'].value
                height = result.params[f'h{i}'].value

                if i <= 6:
                    fwhm = result.params['main_fwhm'].value
                elif i == 7:
                    fwhm = result.params['carb_fwhm'].value
                else:
                    fwhm = result.params['sigma_fwhm'].value

                area = areas.get(name, 0)
                pct_quant = (selected_areas_quant.get(name, 0) / total_quant * 100) if total_quant > 0 else 0
                in_fit = self.peak_in_fit[name].get()
                in_quant = self.peak_in_quant[name].get()

                peak_data.append({
                    'Peak': display_name,
                    'Internal_Name': name,
                    'Center_eV': round(center, 4),
                    'Height': round(height, 6),
                    'FWHM_eV': round(fwhm, 4),
                    'Area': round(area, 6),
                    'Percentage': round(pct_quant, 2),
                    'In_Fit': in_fit,
                    'In_Quantification': in_quant
                })

            peaks_csv = output_dir / f"{sample_name}_peak_parameters.csv"
            pd.DataFrame(peak_data).to_csv(peaks_csv, index=False)

            # 4. Export raw + fitted data CSV
            from c1s_fitter_optimized import total_model, double_arctangent, gaussian

            # Calculate model components
            excluded = {'baseline_total'}
            params_dict = {p: result.params[p].value for p in result.params if p not in excluded}

            baseline = double_arctangent(energy,
                                        params_dict['arc1_center'], params_dict['arc1_height'], params_dict['arc1_width'],
                                        params_dict['arc2_center'], params_dict['arc2_height'], params_dict['arc2_width'])

            total_fit = total_model(energy, **params_dict)
            residual = intensity - total_fit

            spectrum_data = {
                'Energy_eV': energy,
                'Raw_Intensity': intensity,
                'Total_Fit': total_fit,
                'Baseline': baseline,
                'Residual': residual
            }

            # Add individual peak contributions
            for i, name in enumerate(self.peak_names, 1):
                if i <= 6:
                    fwhm = params_dict['main_fwhm']
                elif i == 7:
                    fwhm = params_dict['carb_fwhm']
                else:
                    fwhm = params_dict['sigma_fwhm']

                center = params_dict[f'c{i}']
                height = params_dict[f'h{i}']
                peak = gaussian(energy, center, height, fwhm)
                spectrum_data[f'Peak_{name}'] = peak

            spectrum_csv = output_dir / f"{sample_name}_spectrum_data.csv"
            pd.DataFrame(spectrum_data).to_csv(spectrum_csv, index=False)

            # 5. Export metadata JSON
            import json
            from datetime import datetime

            r_squared = 1 - result.residual.var() / np.var(result.data)

            metadata = {
                'sample_name': sample_name,
                'export_timestamp': datetime.now().isoformat(),
                'source_file': str(self.spectra_files[self.current_index]),
                'fit_quality': {
                    'r_squared': round(r_squared, 6),
                    'reduced_chi_square': round(result.redchi, 6),
                    'aic': round(result.aic, 2),
                    'bic': round(result.bic, 2),
                    'n_data_points': len(energy),
                    'n_variables': result.nvarys
                },
                'baseline_parameters': {
                    'arc1_center': round(result.params['arc1_center'].value, 4),
                    'arc1_height': round(result.params['arc1_height'].value, 6),
                    'arc1_width': round(result.params['arc1_width'].value, 4),
                    'arc2_center': round(result.params['arc2_center'].value, 4),
                    'arc2_height': round(result.params['arc2_height'].value, 6),
                    'arc2_width': round(result.params['arc2_width'].value, 4)
                },
                'fwhm_parameters': {
                    'main_fwhm': round(result.params['main_fwhm'].value, 4),
                    'carbonate_fwhm': round(result.params['carb_fwhm'].value, 4),
                    'sigma_star_fwhm': round(result.params['sigma_fwhm'].value, 4)
                },
                'peak_summary': {
                    name: {
                        'center': round(result.params[f'c{i}'].value, 4),
                        'height': round(result.params[f'h{i}'].value, 6),
                        'area': round(areas.get(name, 0), 6),
                        'percentage': round((selected_areas_quant.get(name, 0) / total_quant * 100) if total_quant > 0 else 0, 2)
                    }
                    for i, name in enumerate(self.peak_names, 1)
                },
                'settings': {
                    'peaks_in_fit': [name for name in self.peak_names if self.peak_in_fit[name].get()],
                    'peaks_in_quantification': [name for name in self.peak_names if self.peak_in_quant[name].get()],
                    'manual_baseline': self.manual_baseline.get(),
                    'custom_peak_centers': self.custom_peak_centers.get()
                }
            }

            metadata_json = output_dir / f"{sample_name}_metadata.json"
            with open(metadata_json, 'w') as f:
                json.dump(metadata, f, indent=2)

            # 6. Export fit report (text)
            report_txt = output_dir / f"{sample_name}_fit_report.txt"
            with open(report_txt, 'w') as f:
                f.write(f"{'='*60}\n")
                f.write(f"C1s XAS Peak Fitting Report\n")
                f.write(f"{'='*60}\n\n")
                f.write(f"Sample: {sample_name}\n")
                f.write(f"Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Source File: {self.spectra_files[self.current_index]}\n\n")

                f.write(f"{'-'*60}\n")
                f.write(f"FIT QUALITY METRICS\n")
                f.write(f"{'-'*60}\n")
                f.write(f"R-squared:           {r_squared:.6f}\n")
                f.write(f"Reduced Chi-square:  {result.redchi:.6e}\n")
                f.write(f"AIC:                 {result.aic:.2f}\n")
                f.write(f"BIC:                 {result.bic:.2f}\n")
                f.write(f"Data Points:         {len(energy)}\n")
                f.write(f"Variables:           {result.nvarys}\n\n")

                f.write(f"{'-'*60}\n")
                f.write(f"BASELINE PARAMETERS\n")
                f.write(f"{'-'*60}\n")
                f.write(f"Step 1: Center={result.params['arc1_center'].value:.3f} eV, ")
                f.write(f"Height={result.params['arc1_height'].value:.4f}, ")
                f.write(f"Width={result.params['arc1_width'].value:.3f}\n")
                f.write(f"Step 2: Center={result.params['arc2_center'].value:.3f} eV, ")
                f.write(f"Height={result.params['arc2_height'].value:.4f}, ")
                f.write(f"Width={result.params['arc2_width'].value:.3f}\n\n")

                f.write(f"{'-'*60}\n")
                f.write(f"PEAK PARAMETERS\n")
                f.write(f"{'-'*60}\n")
                f.write(f"{'Peak':<12} {'Center':>8} {'Height':>10} {'FWHM':>8} {'Area':>10} {'%':>8}\n")
                f.write(f"{'-'*60}\n")

                for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
                    center = result.params[f'c{i}'].value
                    height = result.params[f'h{i}'].value
                    if i <= 6:
                        fwhm = result.params['main_fwhm'].value
                    elif i == 7:
                        fwhm = result.params['carb_fwhm'].value
                    else:
                        fwhm = result.params['sigma_fwhm'].value
                    area = areas.get(name, 0)
                    pct = (selected_areas_quant.get(name, 0) / total_quant * 100) if total_quant > 0 else 0

                    fit_mark = '*' if self.peak_in_fit[name].get() else ' '
                    quant_mark = '+' if self.peak_in_quant[name].get() else ' '

                    f.write(f"{fit_mark}{quant_mark}{display_name:<10} {center:>8.3f} {height:>10.5f} {fwhm:>8.3f} {area:>10.5f} {pct:>7.1f}%\n")

                f.write(f"\n* = included in fit, + = included in quantification\n\n")

                f.write(f"{'-'*60}\n")
                f.write(f"EXPORTED FILES\n")
                f.write(f"{'-'*60}\n")
                f.write(f"- {sample_name}_fit_plot.pdf\n")
                f.write(f"- {sample_name}_fit_plot.png\n")
                f.write(f"- {sample_name}_peak_parameters.csv\n")
                f.write(f"- {sample_name}_spectrum_data.csv\n")
                f.write(f"- {sample_name}_metadata.json\n")
                f.write(f"- {sample_name}_fit_report.txt\n")

            messagebox.showinfo("Export Complete",
                              f"Exported to {output_dir}:\n\n"
                              f"• {sample_name}_fit_plot.pdf/png\n"
                              f"• {sample_name}_peak_parameters.csv\n"
                              f"• {sample_name}_spectrum_data.csv\n"
                              f"• {sample_name}_metadata.json\n"
                              f"• {sample_name}_fit_report.txt")

        except Exception as e:
            messagebox.showerror("Export Error", f"Export failed:\n{str(e)}")
            import traceback
            traceback.print_exc()

    def batch_export_all(self):
        """Batch export all loaded samples to a directory."""
        if not self.spectra_files:
            messagebox.showwarning("Warning", "No samples loaded")
            return

        # Ask for output directory
        output_dir = filedialog.askdirectory(
            title="Select Output Directory for Batch Export"
        )

        if not output_dir:
            return

        output_dir = Path(output_dir)

        # Create summary data structure
        summary_data = []

        # Progress tracking
        total = len(self.spectra_files)
        errors = []

        # Save current index to restore later
        original_index = self.current_index

        try:
            for idx, file_path in enumerate(self.spectra_files):
                sample_name = file_path.stem

                try:
                    # Load and fit the spectrum
                    energy, intensity = load_spectrum(file_path)

                    if self.custom_peak_centers.get():
                        result = self.fit_spectrum_with_custom_centers(energy, intensity)
                    else:
                        result = fit_spectrum(energy, intensity)

                    # Calculate areas
                    areas = calculate_peak_areas(energy, result)
                    selected_areas = {k: v for k, v in areas.items() if self.peak_in_quant[k].get()}
                    total_area = sum(selected_areas.values())

                    r_squared = 1 - result.residual.var() / np.var(result.data)

                    # Build summary row
                    row = {
                        'Sample': sample_name,
                        'R_squared': round(r_squared, 6),
                        'Chi_squared': round(result.redchi, 6),
                        'AIC': round(result.aic, 2)
                    }

                    # Add percentages for each peak
                    for name, display_name in zip(self.peak_names, self.peak_display_names):
                        if self.peak_in_quant[name].get():
                            pct = (selected_areas.get(name, 0) / total_area * 100) if total_area > 0 else 0
                            row[f'{display_name}_%'] = round(pct, 2)

                    # Add peak centers
                    for i, name in enumerate(self.peak_names, 1):
                        row[f'{name}_center'] = round(result.params[f'c{i}'].value, 3)

                    summary_data.append(row)

                    # Export individual sample data
                    sample_dir = output_dir / sample_name
                    sample_dir.mkdir(exist_ok=True)

                    # Quick export: peak parameters CSV
                    peak_data = []
                    for i, (name, display_name) in enumerate(zip(self.peak_names, self.peak_display_names), 1):
                        center = result.params[f'c{i}'].value
                        height = result.params[f'h{i}'].value
                        if i <= 6:
                            fwhm = result.params['main_fwhm'].value
                        elif i == 7:
                            fwhm = result.params['carb_fwhm'].value
                        else:
                            fwhm = result.params['sigma_fwhm'].value
                        area = areas.get(name, 0)
                        pct = (selected_areas.get(name, 0) / total_area * 100) if total_area > 0 else 0

                        peak_data.append({
                            'Peak': display_name,
                            'Center_eV': round(center, 4),
                            'Height': round(height, 6),
                            'FWHM_eV': round(fwhm, 4),
                            'Area': round(area, 6),
                            'Percentage': round(pct, 2)
                        })

                    pd.DataFrame(peak_data).to_csv(sample_dir / f"{sample_name}_peaks.csv", index=False)

                except Exception as e:
                    errors.append(f"{sample_name}: {str(e)}")
                    continue

            # Export summary CSV
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_csv = output_dir / "batch_summary.csv"
                summary_df.to_csv(summary_csv, index=False)

            # Restore original sample
            self.current_index = original_index
            self.load_current_spectrum()

            # Show completion message
            msg = f"Batch export complete!\n\n"
            msg += f"Processed: {len(summary_data)}/{total} samples\n"
            msg += f"Output: {output_dir}\n\n"
            msg += f"Files created:\n"
            msg += f"• batch_summary.csv (all samples)\n"
            msg += f"• Individual folders with peak data"

            if errors:
                msg += f"\n\nErrors ({len(errors)}):\n"
                for err in errors[:5]:  # Show first 5 errors
                    msg += f"• {err}\n"
                if len(errors) > 5:
                    msg += f"• ... and {len(errors) - 5} more"

            messagebox.showinfo("Batch Export Complete", msg)

        except Exception as e:
            # Restore original sample
            self.current_index = original_index
            self.load_current_spectrum()
            messagebox.showerror("Batch Export Error", f"Export failed:\n{str(e)}")

    def _set_app_icon(self):
        """Set application icon using the best available method."""
        # Handle both development and PyInstaller compiled mode
        if getattr(sys, 'frozen', False):
            # Running as compiled executable
            base_path = Path(sys._MEIPASS)
        else:
            # Running in development
            base_path = Path(__file__).parent

        # Try PNG with wm_iconphoto first (often renders more crisply)
        png_sizes = [256, 128, 64, 48, 32]  # Try larger sizes first
        for size in png_sizes:
            png_path = base_path / f"carbonpeaks_{size}.png"
            if png_path.exists():
                try:
                    icon_image = tk.PhotoImage(file=str(png_path))
                    self.root.wm_iconphoto(True, icon_image)
                    # Keep reference to prevent garbage collection
                    self._icon_image = icon_image
                    return
                except tk.TclError:
                    continue

        # Fallback to ICO file with iconbitmap
        ico_path = base_path / "carbonpeaks.ico"
        if ico_path.exists():
            try:
                self.root.iconbitmap(str(ico_path))
            except tk.TclError:
                pass


def _get_base_path():
    """Get base path for assets (handles both dev and PyInstaller)."""
    if getattr(sys, 'frozen', False):
        return Path(sys._MEIPASS)
    return Path(__file__).parent


def _create_splash(root):
    """Create a splash screen showing the app icon during loading."""
    splash = tk.Toplevel(root)
    splash.overrideredirect(True)  # No title bar or borders

    base_path = _get_base_path()

    # Load the largest available icon for the splash
    icon_image = None
    for size in [256, 128, 64]:
        png_path = base_path / f"carbonpeaks_{size}.png"
        if png_path.exists():
            try:
                icon_image = tk.PhotoImage(file=str(png_path))
                break
            except tk.TclError:
                continue

    if icon_image is None:
        splash.destroy()
        return None, None

    # Build splash content
    splash.configure(bg='#2b2b2b')
    img_label = tk.Label(splash, image=icon_image, bg='#2b2b2b')
    img_label.pack(padx=30, pady=(25, 10))
    text_label = tk.Label(splash, text="CarbonKPeaks", font=("Segoe UI", 16, "bold"),
                          fg='white', bg='#2b2b2b')
    text_label.pack(pady=(0, 5))
    sub_label = tk.Label(splash, text="Loading...", font=("Segoe UI", 10),
                         fg='#aaaaaa', bg='#2b2b2b')
    sub_label.pack(pady=(0, 20))

    # Center the splash on screen
    # Force geometry calculation, then place at center
    splash.update()
    w = splash.winfo_width()
    h = splash.winfo_height()
    screen_w = splash.winfo_screenwidth()
    screen_h = splash.winfo_screenheight()
    x = (screen_w - w) // 2
    y = (screen_h - h) // 2
    splash.geometry(f"{w}x{h}+{x}+{y}")
    splash.lift()
    splash.update()

    return splash, icon_image


def main():
    """Main function."""
    import argparse
    parser = argparse.ArgumentParser(description='C1s XAS Peak Viewer')
    parser.add_argument('--spectra_dir', type=str, default=None,
                       help='Directory containing normalized spectra')
    args = parser.parse_args()

    # Close PyInstaller bootloader splash if present (shows during extraction)
    try:
        import pyi_splash  # noqa: F401 - only exists in PyInstaller builds
        pyi_splash.close()
    except ImportError:
        pass

    root = tk.Tk()
    root.withdraw()  # Hide main window while loading

    # Show splash screen with icon
    splash, _splash_icon = _create_splash(root)

    # Build the application (this is the slow part)
    app = C1sPeakViewerFinal(root, spectra_dir=args.spectra_dir)

    # Close splash and show main window
    if splash is not None:
        splash.destroy()
    root.deiconify()

    root.mainloop()


if __name__ == '__main__':
    main()

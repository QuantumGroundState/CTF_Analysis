# -*- coding: utf-8 -*-
"""
2D Image CTF Analyzer (Spatial + Temporal Envelopes)
=============================================================
"""

import sys
import os
import numpy as np
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons
from matplotlib.patches import Circle
from PIL import Image

# -- Global Font Settings -----------------------------------------------------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 8  

# -- Electron wavelengths (relativistic) --------------------------------------
ELECTRON_WAVELENGTHS = {
    "300 kV": 1.9687e-3,   # nm
    "200 kV": 2.5079e-3,
    "120 kV": 3.3492e-3,
    "80 kV":  4.1757e-3,
    "60 kV":  4.8661e-3,
    "30 kV":  6.9791e-3,
}

# -- Analysis functions (SPEED OPTIMIZED) -------------------------------------
def perform_binning(image, bin_factor):
    if bin_factor == 1:
        return image.copy()
    h, w = image.shape
    new_h, new_w = h // bin_factor, w // bin_factor
    trimmed = image[:new_h * bin_factor, :new_w * bin_factor]
    return trimmed.reshape(new_h, bin_factor, new_w, bin_factor).mean(axis=(1, 3))

def compute_radial_profile_fast(image, cx, cy, max_radius):
    h, w = image.shape
    y, x = np.ogrid[0:h, 0:w]
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    r_int = np.round(r).astype(int)
    
    mask = r_int <= max_radius
    r_int_masked = r_int[mask]
    vals_masked = image[mask]
    
    profile = np.bincount(r_int_masked, weights=vals_masked)
    counts = np.bincount(r_int_masked)
    
    valid = counts > 0
    profile[valid] /= counts[valid]
    
    radii = np.arange(len(profile))
    return radii, profile

# -- Interactive GUI ----------------------------------------------------------
class CircularFeatureAnalyzer:
    def __init__(self, image_path):
        img = Image.open(image_path).convert("L")
        self.base_image = np.array(img, dtype=float) / 255.0
        
        self.bin_factor = 2 
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape
        
        self.max_valid_radius = int(np.hypot(self.w / 2.0, self.h / 2.0))
        self.fit_start_radius = min(35.0, max(1, self.max_valid_radius - 10))
        self.fit_end_radius = min(347.0, self.max_valid_radius)
        
        self.base_cal_nm_per_px = 0.114
        self.wavelength_nm = ELECTRON_WAVELENGTHS["300 kV"]
        self.current_kv_ev = 300000.0  
        
        # Fixing States
        self.fix_cs = True
        self.fix_cc = True
        self.current_cs = 1.2
        self.current_cc = 1.4  

        self.rad_avg_fft = None
        self.radii_px = None
        self.k_axis = None
        self.fft_img_log = None
        self.N_dimension = max(self.h, self.w)

        self._build_ui()
        self._compute_fft_and_profile()
        self._update_slider_text()
        self._run_fit()
        
        plt.show()

    def _build_ui(self):
        self.fig = plt.figure(figsize=(18, 10), facecolor="#1a1a2e")
        self.fig.canvas.manager.set_window_title("CTF Analyzer")

        gs = self.fig.add_gridspec(
            2, 4,
            left=0.02, right=0.98, top=0.88, bottom=0.04, 
            hspace=0.12, wspace=0.10, 
            height_ratios=[1.2, 1.1],
            width_ratios=[1.2, 1.2, 0.35, 0.65]
        )

        # -- Top Controls --
        ax_slider_start = self.fig.add_axes([0.08, 0.95, 0.18, 0.02], facecolor="#2a2a4a")
        self.slider_start = Slider(ax_slider_start, "Fit Start", 1, self.max_valid_radius, valinit=self.fit_start_radius, valstep=1, color="#3b82f6")
        self.slider_start.label.set_color("w")
        
        ax_slider_end = self.fig.add_axes([0.08, 0.91, 0.18, 0.02], facecolor="#2a2a4a")
        self.slider_end = Slider(ax_slider_end, "Fit End", 10, self.max_valid_radius, valinit=self.fit_end_radius, valstep=1, color="#00ff88")
        self.slider_end.label.set_color("w")

        self.slider_start.on_changed(self._on_slider_change)
        self.slider_end.on_changed(self._on_slider_change)

        ax_btn = self.fig.add_axes([0.35, 0.92, 0.10, 0.04])
        self.btn_analyze = Button(ax_btn, "UPDATE FIT", color="#ff3366", hovercolor="#ff6699")
        self.btn_analyze.label.set_color("white")
        self.btn_analyze.label.set_fontweight("bold")
        self.btn_analyze.on_clicked(self._run_fit)

        ax_btn_load = self.fig.add_axes([0.47, 0.92, 0.10, 0.04])
        self.btn_load = Button(ax_btn_load, "LOAD IMAGE", color="#3b82f6", hovercolor="#60a5fa")
        self.btn_load.label.set_color("white")
        self.btn_load.label.set_fontweight("bold")
        self.btn_load.on_clicked(self._on_load_image)

        ax_px = self.fig.add_axes([0.65, 0.925, 0.08, 0.03], facecolor="#2a2a4a")
        self.tb_px = TextBox(ax_px, "Base Px (nm): ", initial=str(self.base_cal_nm_per_px), color="0.9", hovercolor="1.0")
        self.tb_px.label.set_color("white")
        self.tb_px.on_submit(self._on_px_change)

        # 1. Real Image
        self.ax_img = self.fig.add_subplot(gs[0, 0])
        self.ax_img.set_title("Real Space Image", color="w", pad=5, fontsize=9)
        self.img_plot = self.ax_img.imshow(self.image, cmap="gray", origin="lower")
        self.ax_img.axis("off")

        # 2. 2D FFT Image
        self.ax_fft = self.fig.add_subplot(gs[0, 1])
        self.ax_fft.set_facecolor("#12122a")
        self.ax_fft.set_title(f"2D Spectrum (True Bin {self.bin_factor}x)", color="w", pad=5, fontsize=9)
        self.ax_fft.axis("off")
        
        self.fft_start_circle = Circle((0, 0), self.fit_start_radius, fill=False, edgecolor="#3b82f6", linewidth=1.5, linestyle="--")
        self.fft_end_circle = Circle((0, 0), self.fit_end_radius, fill=False, edgecolor="#00ff88", linewidth=1.5, linestyle="--")
        self.scherzer_circle = Circle((0, 0), 10.0, fill=False, edgecolor="red", linewidth=1.5, linestyle=":")
        
        self.ax_fft.add_patch(self.fft_start_circle)
        self.ax_fft.add_patch(self.fft_end_circle)
        self.ax_fft.add_patch(self.scherzer_circle)
        self.scherzer_circle.set_visible(False)

        self.text_fft_start = self.ax_fft.text(0.02, 0.98, "", color="#3b82f6", transform=self.ax_fft.transAxes, va="top", fontsize=9, fontweight="bold")
        self.text_fft_end = self.ax_fft.text(0.02, 0.92, "", color="#00ff88", transform=self.ax_fft.transAxes, va="top", fontsize=9, fontweight="bold")
        
        # 3. Voltage & Binning Stack
        gs_radio = gs[0, 2].subgridspec(2, 1, hspace=0.15)
        
        self.ax_kv = self.fig.add_subplot(gs_radio[0], facecolor="#1a1a2e")
        self.ax_kv.set_title("Voltage (kV)", color="w", pad=5, fontsize=9)
        self.ax_kv.set_aspect('equal', adjustable='box')
        self.ax_kv.axis("off") 
        
        self.radio_kv = RadioButtons(self.ax_kv, list(ELECTRON_WAVELENGTHS.keys()), active=0, activecolor='blue')
        for label in self.radio_kv.labels: label.set_color("white")
        for circle in self.radio_kv.circles:
            circle.set_edgecolor('#00ff88')
            circle.set_linewidth(1.5)
        self.radio_kv.on_clicked(self._on_kv_select)

        self.ax_bin = self.fig.add_subplot(gs_radio[1], facecolor="#1a1a2e")
        self.ax_bin.set_title("Binning", color="w", pad=5, fontsize=9)
        self.ax_bin.set_aspect('equal', adjustable='box')
        self.ax_bin.axis("off")
        
        initial_bin_index = {1:0, 2:1, 4:2, 8:3}.get(self.bin_factor, 1)
        self.radio_bin = RadioButtons(self.ax_bin, ['1', '2', '4', '8'], active=initial_bin_index, activecolor='blue')
        for label in self.radio_bin.labels: label.set_color("white")
        for circle in self.radio_bin.circles:
            circle.set_edgecolor('#3b82f6')
            circle.set_linewidth(1.5)
        self.radio_bin.on_clicked(self._on_bin_select)

        # 4. Inputs & Guesses
        self.ax_input_bg = self.fig.add_subplot(gs[0, 3], facecolor="#1a1a2e")
        self.ax_input_bg.set_title("Inputs & Guesses", color="w", pad=2, fontsize=9)
        self.ax_input_bg.axis("off")
        
        # Grid updated for 8 inputs + 1 toggle row
        gs_in = gs[0, 3].subgridspec(9, 2, wspace=0.05, hspace=0.10, width_ratios=[0.55, 0.45])
        
        # Check toggles in the top row
        ax_check = self.fig.add_subplot(gs_in[0, :], facecolor="#1a1a2e")
        ax_check.axis("off")
        self.check_fix = CheckButtons(ax_check, ['Fix Cs', 'Fix Cc'], [self.fix_cs, self.fix_cc])
        for label in self.check_fix.labels:
            label.set_color("white")
            label.set_fontsize(9)
        for rect in self.check_fix.rectangles:
            rect.set_edgecolor("white")
            rect.set_facecolor("#2a2a4a")
        for line_list in self.check_fix.lines:
            for line in line_list:
                line.set_color("#00ff88")
        self.check_fix.on_clicked(self._on_check_click)

        labels = ["Defocus (nm)", "Cs (mm)", "Cc (mm)", "dE (eV)", "Alpha (mrad)", "Phase (rad)", "CTF Amp", "Bg Offset"]
        initials = ["100.0", str(self.current_cs), str(self.current_cc), "0.6", "0.5", "0.07", "auto", "auto"]
        self.text_boxes = []

        for i in range(8):
            ax_lbl = self.fig.add_subplot(gs_in[i+1, 0])
            ax_lbl.axis("off")
            ax_lbl.text(1.0, 0.5, labels[i], color="w", va="center", ha="right", fontsize=9)

            ax_tb = self.fig.add_subplot(gs_in[i+1, 1])
            tb = TextBox(ax_tb, "", initial=initials[i], color="0.9", hovercolor="1.0")
            self.text_boxes.append(tb)

        self.tb_def, self.tb_cs, self.tb_cc, self.tb_dE, self.tb_alpha, self.tb_phase, self.tb_ctfamp, self.tb_bgoff = self.text_boxes

        # 5. 1D CTF Fit Plot & Residuals
        gs_ctf = gs[1, 0:3].subgridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        
        self.ax_ctf = self.fig.add_subplot(gs_ctf[0])
        self.ax_ctf.set_facecolor("#12122a")
        self.ax_ctf.set_title("1D Amplitude Spectrum & Absolute CTF Fit", color="w", pad=5, fontsize=9)
        self.ax_ctf.tick_params(colors="#666", labelsize=8, labelbottom=False)
        
        self.line_rad_fft, = self.ax_ctf.plot([], [], color="gray", alpha=0.7, label="Data")
        self.line_bg, = self.ax_ctf.plot([], [], color="yellow", linestyle="--", linewidth=1.5, label="Background")
        self.line_env, = self.ax_ctf.plot([], [], color="cyan", linestyle="--", linewidth=1.5, label="Total Env")
        self.line_ctf_fit, = self.ax_ctf.plot([], [], color="#ff3366", label="CTF Fit", linewidth=2.0)
        self.scherzer_vline = self.ax_ctf.axvline(x=0, color='red', linestyle=':', linewidth=1.5, label='Scherzer Freq')
        self.ax_ctf.legend(loc="upper right", facecolor="#1a1a2e", edgecolor="#333", labelcolor="#cccccc", borderpad=0.2, fontsize=8)

        self.ax_res = self.fig.add_subplot(gs_ctf[1], sharex=self.ax_ctf, facecolor="#12122a")
        self.ax_res.set_xlabel(r"Spatial Frequency k (1/nm)", color="#aaa", labelpad=2, fontsize=9)
        self.ax_res.tick_params(colors="#666", labelsize=8)
        self.ax_res.set_ylabel("Resid", color="#aaa", fontsize=8)
        self.line_res, = self.ax_res.plot([], [], color="#00ff88", linewidth=1.0)
        self.ax_res.axhline(0, color="gray", linestyle="--", linewidth=1.0)

        # 6. Scherzer & Fitted Parameters
        gs_results = gs[1, 3].subgridspec(2, 1, height_ratios=[0.45, 1.55], hspace=0.1)
        
        self.ax_sch = self.fig.add_subplot(gs_results[0], facecolor="#1a1a2e")
        self.ax_sch.set_title("Optimal Scherzer Params", color="w", pad=5, fontsize=9)
        self.ax_sch.axis("off")
        self.sch_text = self.ax_sch.text(0.0, 0.5, "Calculating...", color="#fb923c", va="center", ha="left", fontfamily="monospace", fontsize=9)

        self.ax_table = self.fig.add_subplot(gs_results[1], facecolor="#1a1a2e")
        self.ax_table.axis("off")
        self.ax_table.set_title("Fitted Parameters", color="w", pad=5, fontsize=9)
        self.table_text = self.ax_table.text(
            0.0, 0.5, "Computing...", 
            va='center', ha='left', fontfamily='monospace', color='#e2e8f0', fontsize=8
        )

    def _on_check_click(self, label):
        if label == 'Fix Cs':
            self.fix_cs = not self.fix_cs
        elif label == 'Fix Cc':
            self.fix_cc = not self.fix_cc
        self._run_fit()

    def _update_slider_text(self):
        try:
            self.base_cal_nm_per_px = float(self.tb_px.text)
        except ValueError: pass
            
        effective_px_size = self.base_cal_nm_per_px * self.bin_factor
        k_start = self.fit_start_radius / (self.N_dimension * effective_px_size)
        k_end = self.fit_end_radius / (self.N_dimension * effective_px_size)
        
        self.slider_start.valtext.set_text(f"{int(self.fit_start_radius)} px")
        self.slider_end.valtext.set_text(f"{int(self.fit_end_radius)} px")

        self.text_fft_start.set_text(f"Start: {int(self.fit_start_radius)} px [{k_start:.4f} 1/nm]")
        self.text_fft_end.set_text(f"End: {int(self.fit_end_radius)} px [{k_end:.4f} 1/nm]")

    def _on_px_change(self, text):
        self._update_slider_text()
        self._run_fit()

    def _on_load_image(self, event):
        try:
            from tkinter import Tk, filedialog
            root = Tk()
            root.withdraw()
            path = filedialog.askopenfilename(title="Select an image", filetypes=[("Images", "*.png *.jpg *.tif *.tiff"), ("All", "*.*")])
            root.destroy()
            if path and os.path.exists(path):
                img = Image.open(path).convert("L")
                self.base_image = np.array(img, dtype=float) / 255.0
                self.bin_factor = 2 
                self.radio_bin.set_active(1) 
                self._apply_binning(update_fit=True)
        except Exception as e:
            print(f"Error loading image: {e}")

    def _on_bin_select(self, label):
        self.bin_factor = int(label)
        self._apply_binning(update_fit=True)

    def _apply_binning(self, update_fit=False):
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape
        self.img_plot.set_data(self.image)
        self.img_plot.set_extent([-0.5, self.w - 0.5, -0.5, self.h - 0.5])
        self.ax_img.set_xlim(0, self.w)
        self.ax_img.set_ylim(0, self.h)
        
        old_max_radius = self.max_valid_radius
        self.max_valid_radius = int(np.hypot(self.w / 2.0, self.h / 2.0))
        boundary_ratio = self.max_valid_radius / max(1, old_max_radius)
        
        self.slider_start.valmax = self.max_valid_radius
        self.slider_start.ax.set_xlim(self.slider_start.valmin, self.slider_start.valmax)
        self.slider_end.valmax = self.max_valid_radius
        self.slider_end.ax.set_xlim(self.slider_end.valmin, self.slider_end.valmax)
        
        new_start = max(1, int(self.fit_start_radius * boundary_ratio))
        new_end = max(2, min(self.max_valid_radius, int(self.fit_end_radius * boundary_ratio)))
        
        self.slider_start.eventson = False
        self.slider_end.eventson = False
        self.slider_start.set_val(new_start)
        self.slider_end.set_val(new_end)
        self.slider_start.eventson = True
        self.slider_end.eventson = True
        
        self._update_slider_text()
        if update_fit:
            self._compute_fft_and_profile()
            self._run_fit()

    def _compute_fft_and_profile(self):
        fft_img = np.fft.fft2(self.image)
        fft_shifted = np.fft.fftshift(fft_img)
        
        amplitude_spectrum = np.abs(fft_shifted)
        self.fft_img_log = np.log(amplitude_spectrum + 1e-10)
        
        h, w = amplitude_spectrum.shape
        fft_cx, fft_cy = w // 2, h // 2
        self.N_dimension = max(h, w)
        
        self.ax_fft.clear()
        self.ax_fft.set_title(f"2D Spectrum (True Bin {self.bin_factor}x)", color="w", pad=5, fontsize=9)
        self.ax_fft.imshow(self.fft_img_log, cmap="gray")
        self.ax_fft.axis("off")
        
        self.ax_fft.add_patch(self.fft_start_circle)
        self.ax_fft.add_patch(self.fft_end_circle)
        self.ax_fft.add_patch(self.scherzer_circle)
        
        self.fft_start_circle.set_center((fft_cx, fft_cy))
        self.fft_end_circle.set_center((fft_cx, fft_cy))
        self.scherzer_circle.set_center((fft_cx, fft_cy))
        
        self.fft_start_circle.set_visible(True)
        self.fft_end_circle.set_visible(True)

        self.text_fft_start = self.ax_fft.text(0.02, 0.98, "", color="#3b82f6", transform=self.ax_fft.transAxes, va="top", fontsize=9, fontweight="bold")
        self.text_fft_end = self.ax_fft.text(0.02, 0.92, "", color="#00ff88", transform=self.ax_fft.transAxes, va="top", fontsize=9, fontweight="bold")
        self._update_slider_text()

        max_rad = int(np.hypot(w / 2.0, h / 2.0))
        self.radii_px, self.rad_avg_fft = compute_radial_profile_fast(amplitude_spectrum, fft_cx, fft_cy, max_rad)

    def _on_slider_change(self, val):
        self.fit_start_radius = self.slider_start.val
        self.fit_end_radius = self.slider_end.val
        self.fft_start_circle.set_radius(self.fit_start_radius)
        self.fft_end_circle.set_radius(self.fit_end_radius)
        self._update_slider_text()
        self._run_fit()

    def _on_kv_select(self, label):
        self.wavelength_nm = ELECTRON_WAVELENGTHS[label]
        self.current_kv_ev = float(label.split()[0]) * 1000.0

    def _ctf_model_base(self, k, defocus, cs_mm, cc_mm, dE_ev, alpha_mrad, phase_rad, ctf_amp, bg_offset):
        cs_nm = cs_mm * 1e6
        cc_nm = cc_mm * 1e6
        alpha_rad = alpha_mrad * 1e-3
        
        # Phase Argument
        chi = np.pi * self.wavelength_nm * defocus * (k**2) - 0.5 * np.pi * cs_nm * (self.wavelength_nm**3) * (k**4) + phase_rad
        
        # Temporal Envelope (Chromatic & Energy Spread)
        delta_nm = cc_nm * (dE_ev / self.current_kv_ev)
        temporal_env = np.exp(-0.5 * (np.pi**2) * (self.wavelength_nm**2) * (delta_nm**2) * (k**4))
        
        # Spatial Envelope (Beam Convergence)
        spatial_env = np.exp(-(np.pi**2) * (alpha_rad**2) * (defocus * k - cs_nm * (self.wavelength_nm**2) * (k**3))**2)
        
        total_envelope = spatial_env * temporal_env
        ctf_abs = np.abs(np.sin(chi)) * total_envelope
        
        return bg_offset + ctf_amp * ctf_abs

    def _get_dynamic_model(self):
        def dynamic_model(k, *free_params):
            pidx = 0
            defocus = free_params[pidx]; pidx += 1
            
            if self.fix_cs: cs_val = self.current_cs
            else: cs_val = free_params[pidx]; pidx += 1
            
            if self.fix_cc: cc_val = self.current_cc
            else: cc_val = free_params[pidx]; pidx += 1
            
            dE = free_params[pidx]; pidx += 1
            alpha = free_params[pidx]; pidx += 1
            phase = free_params[pidx]; pidx += 1
            ctf_amp = free_params[pidx]; pidx += 1
            bg_offset = free_params[pidx]
            
            return self._ctf_model_base(k, defocus, cs_val, cc_val, dE, alpha, phase, ctf_amp, bg_offset)
        return dynamic_model

    def _run_fit(self, event=None):
        if self.rad_avg_fft is None or self.radii_px is None:
            return
            
        try: self.base_cal_nm_per_px = float(self.tb_px.text)
        except ValueError: pass
        
        effective_px_size = self.base_cal_nm_per_px * self.bin_factor
        self.k_axis = self.radii_px / (self.N_dimension * effective_px_size)
        self._update_1d_fit()

    def _update_1d_fit(self):
        start_idx = int(self.fit_start_radius)
        end_idx = int(self.fit_end_radius)
        
        if start_idx >= end_idx - 5:
            self.table_text.set_text("Invalid Fit Range.\nEnsure End > Start.")
            self.table_text.set_color("#ef4444")
            self.fig.canvas.draw_idle()
            return

        # Fetch inputs
        try: defocus_guess = float(self.tb_def.text)
        except ValueError: defocus_guess = 100.0
        try: self.current_cs = float(self.tb_cs.text)
        except ValueError: self.current_cs = 1.2
        try: self.current_cc = float(self.tb_cc.text)
        except ValueError: self.current_cc = 1.4
        try: dE_guess = float(self.tb_dE.text)
        except ValueError: dE_guess = 0.6
        try: alpha_guess = float(self.tb_alpha.text)
        except ValueError: alpha_guess = 0.5
        try: phase_guess = float(self.tb_phase.text)
        except ValueError: phase_guess = 0.07

        k_fit = self.k_axis[start_idx:end_idx]
        data_fit = self.rad_avg_fft[start_idx:end_idx]

        if self.tb_ctfamp.text.lower() == 'auto': ctf_amp_guess = np.max(data_fit) - np.min(data_fit)
        else:
            try: ctf_amp_guess = float(self.tb_ctfamp.text)
            except ValueError: ctf_amp_guess = np.max(data_fit) - np.min(data_fit)

        if self.tb_bgoff.text.lower() == 'auto': bg_offset_guess = max(0.0, np.min(data_fit))
        else:
            try: bg_offset_guess = max(0.0, float(self.tb_bgoff.text))
            except ValueError: bg_offset_guess = max(0.0, np.min(data_fit))

        p0 = [defocus_guess]
        lower_bounds = [0.0]
        upper_bounds = [np.inf]

        if not self.fix_cs:
            p0.append(self.current_cs)
            lower_bounds.append(0.001)
            upper_bounds.append(20.0)
            
        if not self.fix_cc:
            p0.append(self.current_cc)
            lower_bounds.append(0.0)
            upper_bounds.append(20.0)

        p0.extend([dE_guess, alpha_guess, phase_guess, ctf_amp_guess, bg_offset_guess])
        lower_bounds.extend([0.0, 0.0, -np.pi, 0.0, 0.0])
        upper_bounds.extend([50.0, 10.0, np.pi, np.inf, np.inf])

        model_func = self._get_dynamic_model()
        fit_success = False

        try:
            popt, pcov = curve_fit(
                model_func, k_fit, data_fit, p0=p0, bounds=(lower_bounds, upper_bounds), 
                maxfev=10000, absolute_sigma=True
            )
            
            fit_curve_full = model_func(self.k_axis, *popt)
            residuals_fit = data_fit - fit_curve_full[start_idx:end_idx]
            
            dof = max(1, len(data_fit) - len(popt))
            reduced_chi_sq = np.sum(residuals_fit**2) / dof
            pcov_scaled = pcov * reduced_chi_sq
            perr = np.sqrt(np.diag(pcov_scaled))
            
            fit_success = True

            popt_full, perr_full = {}, {}
            pidx = 0
            
            popt_full['defocus'] = popt[pidx]; perr_full['defocus'] = perr[pidx]; pidx += 1
            if self.fix_cs: popt_full['cs'] = self.current_cs; perr_full['cs'] = None
            else: popt_full['cs'] = popt[pidx]; perr_full['cs'] = perr[pidx]; pidx += 1
            if self.fix_cc: popt_full['cc'] = self.current_cc; perr_full['cc'] = None
            else: popt_full['cc'] = popt[pidx]; perr_full['cc'] = perr[pidx]; pidx += 1
                
            popt_full['dE'] = popt[pidx]; perr_full['dE'] = perr[pidx]; pidx += 1
            popt_full['alpha'] = popt[pidx]; perr_full['alpha'] = perr[pidx]; pidx += 1
            popt_full['phase'] = popt[pidx]; perr_full['phase'] = perr[pidx]; pidx += 1
            popt_full['amp'] = popt[pidx]; perr_full['amp'] = perr[pidx]; pidx += 1
            popt_full['bg'] = popt[pidx]; perr_full['bg'] = perr[pidx]

            cs_err_str = f"{perr_full['cs']:<9.3f}" if perr_full['cs'] is not None else "[Fixed]  "
            cc_err_str = f"{perr_full['cc']:<9.3f}" if perr_full['cc'] is not None else "[Fixed]  "

            table_str = (
                f"{'Parameter':<14} | {'Fitted':<9} | {'Error':<9}\n"
                f"{'-'*41}\n"
                f"{'Defocus (nm)':<14} | {popt_full['defocus']:<9.2f} | {perr_full['defocus']:<9.2f}\n"
                f"{'Cs (mm)':<14} | {popt_full['cs']:<9.3f} | {cs_err_str}\n"
                f"{'Cc (mm)':<14} | {popt_full['cc']:<9.3f} | {cc_err_str}\n"
                f"{'dE (eV)':<14} | {popt_full['dE']:<9.3f} | {perr_full['dE']:<9.3f}\n"
                f"{'Alpha (mrad)':<14} | {popt_full['alpha']:<9.3f} | {perr_full['alpha']:<9.3f}\n"
                f"{'Phase (rad)':<14} | {popt_full['phase']:<9.3f} | {perr_full['phase']:<9.3f}\n"
                f"{'CTF Amp':<14} | {popt_full['amp']:<9.2e} | {perr_full['amp']:<9.2e}\n"
                f"{'Bg Offset':<14} | {popt_full['bg']:<9.2e} | {perr_full['bg']:<9.2e}\n"
            )

        except (RuntimeError, ValueError) as e:
            table_str = f"Curve fit failed.\nAdjust guesses or ROI.\n{'-'*41}\n"
            try: fit_curve_full = model_func(self.k_axis, *p0)
            except: fit_curve_full = np.full_like(self.k_axis, np.nan)
            popt_full = {
                'defocus': defocus_guess, 'cs': self.current_cs, 'cc': self.current_cc, 
                'dE': dE_guess, 'alpha': alpha_guess, 'phase': phase_guess, 'amp': ctf_amp_guess, 'bg': bg_offset_guess
            }

        def compute_scherzer(cs_mm):
            if cs_mm > 0:
                cs_nm = cs_mm * 1e6
                sch_def = 1.2 * np.sqrt(cs_nm * self.wavelength_nm)
                k_sch = 1.515 * (cs_nm**-0.25) * (self.wavelength_nm**-0.75)
                sch_res = 1.0 / k_sch
                return sch_def, k_sch, sch_res
            return 0.0, 0.0, 0.0

        g_def, g_k, g_res = compute_scherzer(self.current_cs)
        f_def, f_k, f_res = compute_scherzer(popt_full['cs'])

        sch_str = (
            f"{'':<14} | {'Guess':<8} | {'Fit':<8}\n"
            f"{'-'*35}\n"
            f"{'Opt Def (nm)':<14} | {g_def:<8.1f} | {f_def:<8.1f}\n"
            f"{'Sp. Freq(1/nm)':<14} | {g_k:<8.3f} | {f_k:<8.3f}\n"
            f"{'Res Limit(nm)':<14} | {g_res:<8.3f} | {f_res:<8.3f}"
        )
        self.sch_text.set_text(sch_str)
        self.sch_text.set_color("#fb923c" if fit_success else "#ef4444")
        
        if f_k > 0:
            effective_px_size = self.base_cal_nm_per_px * self.bin_factor
            r_px_sch = f_k * self.N_dimension * effective_px_size
            self.scherzer_circle.set_radius(r_px_sch)
            self.scherzer_circle.set_visible(True)
            self.scherzer_vline.set_xdata([f_k, f_k])
            self.scherzer_vline.set_visible(True)
        else:
            self.scherzer_circle.set_visible(False)
            self.scherzer_vline.set_visible(False)

        self.table_text.set_text(table_str)
        self.table_text.set_color("#fb923c" if fit_success else "#ef4444") 

        bg_curve_full = np.full_like(self.k_axis, popt_full['bg'])
        cc_nm = popt_full['cc'] * 1e6
        cs_nm = popt_full['cs'] * 1e6
        alpha_rad = popt_full['alpha'] * 1e-3
        
        delta_nm = cc_nm * (popt_full['dE'] / self.current_kv_ev)
        temporal_env = np.exp(-0.5 * (np.pi**2) * (self.wavelength_nm**2) * (delta_nm**2) * (self.k_axis**4))
        spatial_env = np.exp(-(np.pi**2) * (alpha_rad**2) * (popt_full['defocus'] * self.k_axis - cs_nm * (self.wavelength_nm**2) * (self.k_axis**3))**2)
        
        env_curve_full = bg_curve_full + popt_full['amp'] * (temporal_env * spatial_env)

        self.line_rad_fft.set_data(self.k_axis, self.rad_avg_fft)
        
        fit_plot_curve = np.full_like(self.k_axis, np.nan)
        fit_plot_curve[start_idx:end_idx] = fit_curve_full[start_idx:end_idx]
        self.line_ctf_fit.set_data(self.k_axis, fit_plot_curve)
        
        fit_plot_bg = np.full_like(self.k_axis, np.nan)
        fit_plot_bg[start_idx:end_idx] = bg_curve_full[start_idx:end_idx]
        self.line_bg.set_data(self.k_axis, fit_plot_bg)
        
        fit_plot_env = np.full_like(self.k_axis, np.nan)
        fit_plot_env[start_idx:end_idx] = env_curve_full[start_idx:end_idx]
        self.line_env.set_data(self.k_axis, fit_plot_env)

        if fit_success:
            self.line_res.set_data(k_fit, residuals_fit)
            res_min, res_max = np.min(residuals_fit), np.max(residuals_fit)
            res_pad = max(abs(res_max - res_min) * 0.1, 1e-5)
            self.ax_res.set_ylim(res_min - res_pad, res_max + res_pad)
        else:
            self.line_res.set_data([], [])

        x_max = max(self.k_axis[0:end_idx+20]) if len(self.k_axis) > end_idx+20 else max(self.k_axis)
        self.ax_ctf.set_xlim(0, x_max) 
        self.ax_ctf.set_ylim(0, np.percentile(data_fit, 99) * 1.5)
        self.fig.canvas.draw_idle()

if __name__ == "__main__":
    def generate_test_image(size=2048):  
        y, x = np.ogrid[-size // 2 : size // 2, -size // 2 : size // 2]
        r = np.sqrt(x ** 2 + y ** 2)
        k = r / (size * 0.114) 
        chi = np.pi * 1.9687e-3 * 100.0 * (k**2) 
        img = np.abs(np.sin(chi)) 
        img *= np.exp(-(r ** 2) / (2 * (size * 0.25) ** 2))
        img = np.clip(img, 0, 1)
        return (img * 255).astype(np.uint8)

    path = None
    if len(sys.argv) > 1:
        path = sys.argv[1]

    if not path or not os.path.exists(path):
        test_img_array = generate_test_image()
        path = "temp_ctf_test_image.png"
        Image.fromarray(test_img_array).save(path)

    CircularFeatureAnalyzer(path)
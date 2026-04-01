# -*- coding: utf-8 -*-
"""
2D Image CTF Analyzer (Spatial + Temporal Envelopes)
=============================================================

Version Changelog
-----------------
v1.96 (current)
  - Fixed moving-average edge artifact: replaced np.convolve(..., mode='same')
    (which zero-pads beyond the array boundary, pulling the MA below the true
    signal at the ROI edge and creating a spurious peak as it recovers) with
    explicit edge-value replication via np.pad(..., mode='edge') followed by
    np.convolve(..., mode='valid'); the MA now starts at the correct level and
    no artificial gradient or false maximum is introduced at the ROI boundary

v1.95
  - Fixed Reimer n-assignment: all detected maxima and minima are now merged
    into one list sorted jointly by q; consecutive integers n_start, n_start+1,
    … are assigned in q-order rather than odd/even independently per list,
    eliminating the interleaving error that arose when ring counts were unequal
    or rings were missing
  - Auto n-start mode: iterates n_start = 1…19, keeps only physically valid
    fits (Cs > 0, Δz > 0), and selects the starting index that maximises R²;
    Manual mode lets the user specify n_start directly via a Tkinter Spinbox;
    a radio-button pair (Auto / Manual) in the side panel switches modes
  - Layout redesign: Reimer n/q² vs q² scatter+fit plot moved to the top row
    as a third equal-width panel to the right of the FFT (top row is now
    Image | FFT | Reimer plot); wspace between image and FFT reduced for a
    tighter pairing; table columns in the bottom row widened to prevent overlap;
    Reimer text table stays in bottom-row col 5 (text only, no nested gridspec)

v1.9
  - Added "Reimer Est." panel (separate table, bottom row col 5) that applies
    Reimer eq. 6.50 (Cs λ³ q² − 2Δz λ = n/q²) to the detected RA-ring maxima
    and minima: maxima assigned odd n (1,3,5,…), minima even n (2,4,6,…);
    a linear polyfit of n/q² vs q² yields slope = Cs λ³ and intercept = −2Δzλ,
    giving direct estimates of Cs (mm) and defocus (nm) with R² quality metric
    and per-ring residual table; panel color-codes green (R²>0.99), orange
    (R²>0.95), or red (poor fit); updates automatically whenever rings change

v1.8
  - Added version number display in upper-left corner of window
  - Added this version changelog comment block at the top of the file
  - Added red (max) and blue (min) dot markers along the x-axis direction of
    the 2D FFT image, marking local maxima and minima of the horizontal center
    row profile for both positive and negative x directions from center
  - Improved extrema detection: replaced plateau-prone >= window comparison with
    first-difference sign-change algorithm (rising→falling = max, falling→rising
    = min); refines peak position with argmax/argmin in the neighbourhood; also
    enforces minimum inter-peak spacing equal to the adaptive order parameter
  - Added peaks table text box overlaid on the CTF fit plot showing k (1/nm)
    and d-spacing (nm) for every detected maximum and minimum ring; table
    updates whenever the fit or calibration changes
  - Markers now restricted to the fit ROI (between blue fit-start and green
    fit-end circles) instead of the full FFT extent
  - Switched extrema source from raw FFT row pixels to the rotational average
    profile (smoother, rotationally representative); markers placed at
    fft_cx ± r on the x-axis for every extremum radius r found in the ROI
  - UPDATE FIT now fully refreshes marker positions and peaks table (ROI and
    calibration changes both reflected immediately)
  - Exponential background model now the default on startup
  - Peaks table moved out of the CTF plot into a dedicated narrow column
    (ax_peaks) immediately to the right of the CTF/FFT area, not overlaid
  - Red downward-triangle / blue upward-triangle markers added to the bottom
    of the 1D CTF plot at each detected maximum / minimum spatial frequency
  - Moving-average pre-smoothing applied to the radial profile ROI before
    extrema detection (window ≈ order/3) to suppress high-frequency noise
  - "Lin + Exp" background model added: offset + amp*exp(-decay*k) + lin_slope*k
    (combines exponential and linear terms; adds Fix Bg Lin Slope checkbox and
    Bg Lin Slope initial-guess field, both visible only when this model is active)
  - "GaussExp" background model added: offset + A*exp[-(k/(2sigma))^2]/(k^2+Xi^2)
    (Gaussian envelope divided by Lorentzian-type denominator; adds Fix Bg Sigma
    checkbox and Bg Sigma (nm) field; Bg Slope relabeled Bg Xi when selected)
  - Layout changed to nested gridspec: top row = two equal-width image panels,
    bottom row = CTF plot + Ring Peaks panel + Fit Results panel (panels now
    occupy only the lower portion of the figure)
  - Moving-average smoothed ROI profile ("Ring Est (MA)") plotted as a dashed
    green line on the CTF plot (shows the wider minima-finding MA)
  - Separate MA windows for max vs min detection: narrow (order/3) for maxima,
    wider (2*order/3) for minima — minima ride on a sloping background and
    benefit from extra smoothing
  - New "RA Rings" table (col 3 of bottom row) shows extrema directly from the
    rotational-average profile without FFT ± mirroring; "FFT Rings" table (col 2)
    retains the deduped ± mirror result used for CTF marker placement
  - Bottom row expanded to 5 columns to accommodate both ring tables side-by-side
    without overlap; top row remains two equal-width image panels

v1.7
  - Blitting system for FFT overlay circles/text (faster circle drag feedback)
  - DPI selector (Low/Medium/High) added to toolbar
  - Fit sliders moved to toolbar; live px/k readout labels added
  - Exponential background model option added
  - Fix Bg Amp, Fix Bg Offset, Fix Bg Slope checkboxes added
  - Scherzer circle and vertical line on CTF plot
  - Residuals sub-panel below CTF plot

v1.6
  - Initial public version with 2D FFT panel, radial profile, CTF fit
  - Spatial and temporal envelope curves
  - Tkinter side panel with fix/free parameter checkboxes and initial guess entries
"""

import sys
import os
import threading
import numpy as np
from scipy.optimize import curve_fit
try:
    from scipy.fft import fft2, fftshift, next_fast_len
    _HAS_SCIPY_FFT = True
except ImportError:
    _HAS_SCIPY_FFT = False

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.transforms import blended_transform_factory
import tkinter as tk
from tkinter import ttk
from PIL import Image
try:
    import tifffile
    HAS_TIFFFILE = True
except ImportError:
    HAS_TIFFFILE = False

# -- Global Font Settings -----------------------------------------------------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

# Default radio/check button size (doubled from matplotlib default)
_RADIO_RADIUS = 0.10


def _style_radio(radio, edge_color='#00ff88', label_color='white',
                 active_color='#00ff88', radius=_RADIO_RADIUS):
    """Style a RadioButtons widget for visibility across matplotlib versions."""
    for label in radio.labels:
        label.set_color(label_color)
        label.set_fontsize(12)
        label.set_fontfamily('Arial')

    # Newer matplotlib (>=3.7): buttons stored as AxesImage / _buttons
    # Older matplotlib: .circles list
    if hasattr(radio, 'circles'):
        for circle in radio.circles:
            circle.set_edgecolor(edge_color)
            circle.set_linewidth(2)
            circle.set_radius(radius)
    # For newer matplotlib, try to scale via the axes artists
    for artist in radio.ax.patches:
        if hasattr(artist, 'set_radius'):
            artist.set_radius(radius)
            artist.set_edgecolor(edge_color)
            artist.set_linewidth(2)


def _style_check(check, edge_color='white', face_color='#3a3a3a',
                 check_color='#00ff88'):
    """Style a CheckButtons widget for visibility across matplotlib versions."""
    for label in check.labels:
        label.set_color('white')
        label.set_fontsize(12)
        label.set_fontfamily('Arial')

    if hasattr(check, 'rectangles'):
        for rect in check.rectangles:
            rect.set_edgecolor(edge_color)
            rect.set_facecolor(face_color)
            rect.set_linewidth(2)
            # Double the checkbox size
            bbox = rect.get_bbox()
            cx, cy = bbox.x0 + bbox.width / 2, bbox.y0 + bbox.height / 2
            new_w, new_h = bbox.width * 2, bbox.height * 2
            rect.set_bounds(cx - new_w / 2, cy - new_h / 2, new_w, new_h)
    if hasattr(check, 'lines'):
        for line_list in check.lines:
            for line in line_list:
                line.set_color(check_color)
                line.set_linewidth(2)
    # Newer matplotlib: style frame patches
    for artist in check.ax.patches:
        if hasattr(artist, 'set_edgecolor'):
            artist.set_edgecolor(edge_color)
            artist.set_facecolor(face_color)
            artist.set_linewidth(2)

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

_radial_cache = {}   # keyed by (h, w, cx, cy) → (r_int, counts_for_max)

def compute_radial_profile_fast(image, cx, cy, max_radius):
    h, w = image.shape
    key = (h, w, cx, cy)
    if key in _radial_cache:
        r_int = _radial_cache[key]
    else:
        y, x = np.ogrid[0:h, 0:w]
        r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
        r_int = np.round(r).astype(np.int32)
        _radial_cache.clear()      # keep only one entry to avoid memory growth
        _radial_cache[key] = r_int

    mask = r_int <= max_radius
    r_int_masked = r_int[mask]
    vals_masked = image[mask]

    profile = np.bincount(r_int_masked, weights=vals_masked)
    counts = np.bincount(r_int_masked)

    valid = counts > 0
    profile[valid] /= counts[valid]

    radii = np.arange(len(profile))
    return radii, profile

def load_image_as_float(path):
    """Load an image file and return a 2D float64 array normalized to [0, 1].

    Handles 8/16/32-bit TIF/TIFF (including multi-page stacks where
    the first plane is used), as well as standard PNG/JPG via PIL.
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in ('.tif', '.tiff'):
        # --- Try tifffile first (best for scientific TIFFs) ---
        data = None
        if HAS_TIFFFILE:
            try:
                data = tifffile.imread(path)
            except Exception:
                pass  # fall through to PIL
        if data is None:
            # PIL fallback (handles LZW and other common compressions)
            pil_img = Image.open(path)
            data = np.array(pil_img)

        # If multi-page / 3-D stack, take the first plane
        if data.ndim == 3:
            data = data[0]
        elif data.ndim > 3:
            data = data.reshape(-1, *data.shape[-2:])[0]

        # Convert to float64
        data = data.astype(np.float64)

        # Normalise to [0, 1] based on actual data range
        dmin, dmax = data.min(), data.max()
        if dmax - dmin > 0:
            data = (data - dmin) / (dmax - dmin)
        else:
            data = np.zeros_like(data)
        return data

    else:
        # PNG, JPG, etc. – 8-bit pathway via PIL
        img = Image.open(path).convert("L")
        return np.array(img, dtype=np.float64) / 255.0


def _moving_average(arr, w):
    """Box-car moving average with edge-value replication boundary handling.

    Returns an array the same length as ``arr``.  If ``w < 2`` or the array
    is shorter than the window the original array is returned unchanged.

    The array is padded on both sides by replicating the boundary value
    (np.pad mode='edge') before convolving, so the MA at the start and end
    of the array reflects the actual data level rather than being pulled
    toward zero by implicit zero-padding.  This prevents the false local
    maximum that arises when a high-valued signal at the ROI boundary causes
    the zero-padded MA to dip below the data and then recover, creating a
    spurious peak that would be misidentified as a CTF ring.
    """
    a = np.asarray(arr, dtype=np.float64)
    if w < 2 or len(a) < w:
        return a.copy()
    kernel = np.ones(w, dtype=np.float64) / w
    pad = w // 2
    a_padded = np.pad(a, (pad, pad), mode='edge')
    result = np.convolve(a_padded, kernel, mode='valid')
    return result[:len(a)]


def _local_extrema_1d(arr, order=5):
    """Find local maxima and minima via first-difference sign changes.

    A local maximum is where the finite difference changes from positive to
    negative (rising then falling).  A local minimum is the reverse.  Once a
    candidate transition is found, the precise peak/trough is located with
    argmax/argmin inside a window of ``order`` samples on each side so that
    the result is not biased by the exact sign-change index.  A minimum
    inter-peak spacing of ``order`` samples is enforced to suppress
    closely-clustered detections on broad, flat extrema.
    """
    a = np.asarray(arr, dtype=np.float64)
    n = len(a)
    if n < 3:
        return np.array([], dtype=int), np.array([], dtype=int)

    d = np.diff(a)          # length n-1
    maxima, minima = [], []
    i = 1
    while i < n - 1:
        if d[i - 1] > 0 and d[i] < 0:          # rising → falling: maximum
            lo = max(0, i - order)
            hi = min(n, i + order + 1)
            peak = int(lo + np.argmax(a[lo:hi]))
            if not maxima or peak - maxima[-1] >= order:
                maxima.append(peak)
            i += max(1, order)
        elif d[i - 1] < 0 and d[i] > 0:        # falling → rising: minimum
            lo = max(0, i - order)
            hi = min(n, i + order + 1)
            trough = int(lo + np.argmin(a[lo:hi]))
            if not minima or trough - minima[-1] >= order:
                minima.append(trough)
            i += max(1, order)
        else:
            i += 1

    return np.array(maxima, dtype=int), np.array(minima, dtype=int)


# -- Interactive GUI ----------------------------------------------------------
class CircularFeatureAnalyzer:
    def __init__(self, image_path):
        self.base_image = load_image_as_float(image_path)

        self.bin_factor = 2
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape

        self.max_valid_radius = int(np.hypot(self.w / 2.0, self.h / 2.0))
        self.fit_start_radius = min(90.0, max(1, self.max_valid_radius - 10))
        self.fit_end_radius = self.max_valid_radius

        self.base_cal_nm_per_px = 0.114
        self.wavelength_nm = ELECTRON_WAVELENGTHS["300 kV"]
        self.current_kv_ev = 300000.0

        # Fixing States
        self.fix_defocus = False
        self.fix_cs = True
        self.fix_cc = True
        self.fix_phase = False
        self.fix_dE = False
        self.fix_alpha = False
        self.fix_bg_amp = False
        self.fix_bg_offset = False
        self.fix_bg_slope = False
        self.fix_bg_linslope = False
        self.fix_bg_sigma = False
        self.use_spatial = True
        self.full_yrange = False
        self.current_cs = 1.2
        self.current_cc = 1.4
        self.bg_model = "Exponential"
        self.reimer_auto = True
        self.reimer_n_start_manual = 1
        self.bg_amp_fixed = 0.0
        self.bg_offset_fixed = 0.0
        self.bg_slope_fixed = 0.0
        self.bg_linslope_fixed = 0.0
        self.bg_sigma_fixed = 1.0

        self.rad_avg_fft = None
        self.radii_px = None
        self.k_axis = None
        self.fft_img_log = None
        self._fft_imshow = None
        self.N_dimension = max(self.h, self.w)
        self._cached_r_int = None
        self._cached_r_shape = None


        # --- FFT row-extrema storage (pixel columns, for table updates) ---
        self._fft_row_max_px = []
        self._fft_row_min_px = []
        self._fft_cx_last = 0
        self._fft_cy_last = 0
        # --- RA-direct extrema (k values in 1/nm, for the RA Extrema table) ---
        self._ra_max_ks = []
        self._ra_min_ks = []

        # --- Blitting state ---
        self._bg_fft = None
        self._fft_overlays = []

        self._build_ui()

        # Mark FFT overlays as animated (excluded from full draw → enables blitting)
        self._fft_overlays = [
            self.fft_start_circle, self.fft_end_circle, self.scherzer_circle,
            self.text_fft_start, self.text_fft_end,
            self.fft_row_max, self.fft_row_min,
        ]
        for a in self._fft_overlays:
            a.set_animated(True)

        # Re-cache FFT background + blit overlays after every full draw
        self.fig.canvas.mpl_connect('draw_event', self._on_draw_event)

        self._compute_fft_and_profile()
        self._update_slider_text()
        self._run_fit()

        plt.show()

    def _build_ui(self):
        self.fig = plt.figure(figsize=(18, 10), facecolor="#2d2d2d", dpi=80)
        self.fig.canvas.manager.set_window_title("CTF Analyzer  v1.96")
        self.fig.text(0.01, 0.99, "v1.96", color="#00ff88", fontsize=12,
                      fontweight="bold", ha="left", va="top")

        # Outer 2-row grid: top row = images + Reimer plot, bottom row = CTF + tables
        gs = self.fig.add_gridspec(
            2, 1,
            left=0.02, right=0.98, top=0.96, bottom=0.04,
            hspace=0.10,
            height_ratios=[1.35, 1.0]
        )
        # Top row: Image | FFT (tight) | Reimer plot
        gs_top = gs[0].subgridspec(1, 3, wspace=0.04, width_ratios=[1, 1, 1])
        # Bottom row: CTF(0:2) | FFT Rings | RA Rings | Fit Results | Reimer Table
        gs_bot = gs[1].subgridspec(
            1, 6,
            width_ratios=[1.2, 1.2, 0.50, 0.50, 0.82, 0.76],
            wspace=0.06
        )

        # 1. Real Image
        self.ax_img = self.fig.add_subplot(gs_top[0, 0])
        self.ax_img.set_title("Real Space Image", color="w", pad=5, fontsize=12)
        self.img_plot = self.ax_img.imshow(self.image, cmap="gray", origin="lower")
        self.ax_img.axis("off")

        # 2. 2D FFT Image
        self.ax_fft = self.fig.add_subplot(gs_top[0, 1])
        self.ax_fft.set_facecolor("#383838")
        self.ax_fft.set_title(f"2D Spectrum (True Bin {self.bin_factor}x)", color="w", pad=5, fontsize=12)
        self.ax_fft.axis("off")
        
        self.fft_start_circle = Circle((0, 0), self.fit_start_radius, fill=False, edgecolor="#3b82f6", linewidth=1.5, linestyle="--")
        self.fft_end_circle = Circle((0, 0), self.fit_end_radius, fill=False, edgecolor="#00ff88", linewidth=1.5, linestyle="--")
        self.scherzer_circle = Circle((0, 0), 10.0, fill=False, edgecolor="red", linewidth=1.5, linestyle=":")
        
        self.ax_fft.add_patch(self.fft_start_circle)
        self.ax_fft.add_patch(self.fft_end_circle)
        self.ax_fft.add_patch(self.scherzer_circle)
        self.scherzer_circle.set_visible(False)

        self.text_fft_start = self.ax_fft.text(0.02, 0.98, "", color="#3b82f6", transform=self.ax_fft.transAxes, va="top", fontsize=12, fontweight="bold")
        self.text_fft_end = self.ax_fft.text(0.02, 0.92, "", color="#00ff88", transform=self.ax_fft.transAxes, va="top", fontsize=12, fontweight="bold")

        # Horizontal x-axis max/min markers on the FFT image
        self.fft_row_max, = self.ax_fft.plot(
            [], [], 'o', color='red', markersize=5, zorder=6,
            markeredgewidth=0, alpha=0.9, linestyle='none'
        )
        self.fft_row_min, = self.ax_fft.plot(
            [], [], 'o', color='#3399ff', markersize=5, zorder=6,
            markeredgewidth=0, alpha=0.9, linestyle='none'
        )

        self._dpi_map = {'Low': 80, 'Medium': 120, 'High': 150}
        self._current_dpi = self._dpi_map['Low']

        # 5. 1D CTF Fit Plot & Residuals  (spans bottom cols 0-1)
        gs_ctf = gs_bot[0, 0:2].subgridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
        
        self.ax_ctf = self.fig.add_subplot(gs_ctf[0])
        self.ax_ctf.set_facecolor("#383838")
        self.ax_ctf.set_title("1D Amplitude Spectrum & Absolute CTF Fit", color="w", pad=5, fontsize=12)
        self.ax_ctf.tick_params(colors="#666", labelsize=12, labelbottom=False)
        
        self.line_rad_fft, = self.ax_ctf.plot([], [], color="gray", alpha=0.7, label="Data")
        self.line_ma, = self.ax_ctf.plot([], [], color="#44dd99", linestyle="--", linewidth=1.1,
                                          alpha=0.85, label="Ring Est (MA)", zorder=4)
        self.line_bg, = self.ax_ctf.plot([], [], color="yellow", linestyle="--", linewidth=1.5, label="Background")
        self.line_env_spatial, = self.ax_ctf.plot([], [], color="#ff9f1c", linestyle="--", linewidth=1.2, alpha=0.8, label="Spatial Env")
        self.line_env_temporal, = self.ax_ctf.plot([], [], color="#a78bfa", linestyle="--", linewidth=1.2, alpha=0.8, label="Temporal Env")
        self.line_env, = self.ax_ctf.plot([], [], color="cyan", linestyle="-", linewidth=1.5, label="Total Env")
        self.line_ctf_fit, = self.ax_ctf.plot([], [], color="#ff3366", label="CTF Fit", linewidth=2.0)
        self.scherzer_vline = self.ax_ctf.axvline(x=0, color='red', linestyle=':', linewidth=1.5, label='Scherzer Freq')
        self.ax_ctf.legend(loc="upper right", facecolor="#2d2d2d", edgecolor="#333", labelcolor="#cccccc", borderpad=0.2, fontsize=10)

        # Triangle tick-marks along the bottom of the CTF plot for ring peaks.
        # Blended transform: x in data (spatial-frequency) coords, y in axes
        # fraction so the markers sit at a fixed distance above the x-axis
        # regardless of y-scale changes.
        _btrans = blended_transform_factory(self.ax_ctf.transData, self.ax_ctf.transAxes)
        self.ctf_max_marks, = self.ax_ctf.plot(
            [], [], 'v', color='red', markersize=7, zorder=8,
            transform=_btrans, clip_on=False, linestyle='none',
            markeredgewidth=0.6, markeredgecolor='#550000'
        )
        self.ctf_min_marks, = self.ax_ctf.plot(
            [], [], '^', color='#3399ff', markersize=7, zorder=8,
            transform=_btrans, clip_on=False, linestyle='none',
            markeredgewidth=0.6, markeredgecolor='#002255'
        )

        self.ax_res = self.fig.add_subplot(gs_ctf[1], sharex=self.ax_ctf, facecolor="#383838")
        self.ax_res.set_xlabel(r"Spatial Frequency k (1/nm)", color="#aaa", labelpad=2, fontsize=12)
        self.ax_res.tick_params(colors="#666", labelsize=12)
        self.ax_res.set_ylabel("Resid", color="#aaa", fontsize=12)
        self.line_res, = self.ax_res.plot([], [], color="#00ff88", linewidth=1.0)
        self.ax_res.axhline(0, color="gray", linestyle="--", linewidth=1.0)

        # 6a. FFT-based Ring Peaks table (col 2)
        self.ax_peaks = self.fig.add_subplot(gs_bot[0, 2], facecolor="#2d2d2d")
        self.ax_peaks.set_title("FFT Rings", color="w", pad=4, fontsize=10)
        self.ax_peaks.axis("off")
        self.peaks_text = self.ax_peaks.text(
            0.04, 0.97, "",
            color='white', va='top', ha='left',
            fontfamily='monospace', fontsize=8
        )

        # 6b. RA-direct Extrema table (col 3)
        self.ax_ra_peaks = self.fig.add_subplot(gs_bot[0, 3], facecolor="#2d2d2d")
        self.ax_ra_peaks.set_title("RA Rings", color="#44dd99", pad=4, fontsize=10)
        self.ax_ra_peaks.axis("off")
        self.ra_peaks_text = self.ax_ra_peaks.text(
            0.04, 0.97, "",
            color='white', va='top', ha='left',
            fontfamily='monospace', fontsize=8
        )

        # 7. Combined Results Panel (bottom row only, col 4)
        self.ax_results = self.fig.add_subplot(gs_bot[0, 4], facecolor="#2d2d2d")
        self.ax_results.set_title("Fit Results", color="w", pad=5, fontsize=12)
        self.ax_results.axis("off")

        # Scherzer text (orange, top portion)
        self.sch_text = self.ax_results.text(
            0.02, 0.97, "Calculating...",
            color="#fb923c", va="top", ha="left",
            fontfamily="monospace", fontsize=10, fontweight="bold"
        )
        # Fitted params text (white, below Scherzer)
        self.table_text = self.ax_results.text(
            0.02, 0.72, "Computing...",
            va='top', ha='left', fontfamily='monospace',
            color='#e2e8f0', fontsize=10
        )

        # 3. Reimer scatter+fit plot — top row col 2
        self.ax_reimer_plot = self.fig.add_subplot(gs_top[0, 2], facecolor="#383838")
        self.ax_reimer_plot.set_facecolor("#383838")
        self.ax_reimer_plot.set_title("Reimer  n/q²  vs  q²", color="#44dd99",
                                      pad=4, fontsize=12)
        self.ax_reimer_plot.tick_params(colors="#888", labelsize=9,
                                        length=3, pad=2, direction='in')
        self.ax_reimer_plot.set_xlabel("q²  (nm⁻²)", color="#aaa",
                                       fontsize=10, labelpad=2)
        self.ax_reimer_plot.set_ylabel("n / q²  (nm²)", color="#aaa",
                                       fontsize=10, labelpad=2)
        for sp in self.ax_reimer_plot.spines.values():
            sp.set_edgecolor("#666")
        self.reimer_scat_max, = self.ax_reimer_plot.plot(
            [], [], 'v', color='red', markersize=7, zorder=5,
            markeredgewidth=0.6, markeredgecolor='#550000', label='Max (odd n)'
        )
        self.reimer_scat_min, = self.ax_reimer_plot.plot(
            [], [], '^', color='#3399ff', markersize=7, zorder=5,
            markeredgewidth=0.6, markeredgecolor='#002255', label='Min (even n)'
        )
        self.reimer_fit_line, = self.ax_reimer_plot.plot(
            [], [], '--', color='#00ff88', linewidth=1.5, zorder=4, label='Linear fit'
        )
        self.ax_reimer_plot.legend(
            fontsize=9, facecolor='#2d2d2d', edgecolor='#555',
            labelcolor='white', borderpad=0.4, loc='upper right',
            handlelength=1.2, handletextpad=0.4, labelspacing=0.3
        )
        self._reimer_annots = []   # cleared and rebuilt on each ring update

        # 8. Reimer text summary table — bottom row col 5
        self.ax_reimer = self.fig.add_subplot(gs_bot[0, 5], facecolor="#2d2d2d")
        self.ax_reimer.set_title("Reimer Est.", color="#44dd99", pad=4, fontsize=10)
        self.ax_reimer.axis("off")
        self.reimer_text = self.ax_reimer.text(
            0.03, 0.99, "No rings\ndetected.",
            color='#888888', va='top', ha='left',
            fontfamily='monospace', fontsize=8
        )

        self._build_tk_controls()

    def _build_tk_controls(self):
        canvas_widget = self.fig.canvas.get_tk_widget()
        container = canvas_widget.master

        for child in container.pack_slaves():
            child.pack_forget()

        self.tk_toolbar = tk.Frame(container, bg="#202020", padx=8, pady=6)
        self.tk_toolbar.pack(fill="x", side="top")

        self.tk_side = tk.Frame(container, bg="#1f1f1f", width=360, padx=10, pady=10)
        self.tk_side.pack(side="right", fill="y")
        self.tk_side.pack_propagate(False)

        canvas_widget.pack(side="left", fill="both", expand=True)

        self.var_kv = tk.StringVar(value="300 kV")
        self.var_bin = tk.StringVar(value=str(self.bin_factor))
        self.var_dpi = tk.StringVar(value="Low")
        self.var_px = tk.StringVar(value=str(self.base_cal_nm_per_px))

        self.var_fit_start = tk.IntVar(value=int(self.fit_start_radius))
        self.var_fit_end = tk.IntVar(value=int(self.fit_end_radius))

        self.var_fix_cs = tk.BooleanVar(value=self.fix_cs)
        self.var_fix_cc = tk.BooleanVar(value=self.fix_cc)
        self.var_fix_phase = tk.BooleanVar(value=self.fix_phase)
        self.var_fix_dE = tk.BooleanVar(value=self.fix_dE)
        self.var_fix_defocus = tk.BooleanVar(value=self.fix_defocus)
        self.var_fix_alpha = tk.BooleanVar(value=self.fix_alpha)
        self.var_fix_bg_amp = tk.BooleanVar(value=self.fix_bg_amp)
        self.var_fix_bg_offset = tk.BooleanVar(value=self.fix_bg_offset)
        self.var_fix_bg_slope = tk.BooleanVar(value=self.fix_bg_slope)
        self.var_fix_bg_linslope = tk.BooleanVar(value=self.fix_bg_linslope)
        self.var_fix_bg_sigma = tk.BooleanVar(value=self.fix_bg_sigma)
        self.var_use_spatial = tk.BooleanVar(value=self.use_spatial)
        self.var_full_yrange = tk.BooleanVar(value=self.full_yrange)
        self.var_bg_model = tk.StringVar(value=self.bg_model)

        self.var_defocus = tk.StringVar(value="500.0")
        self.var_cs = tk.StringVar(value=str(self.current_cs))
        self.var_cc = tk.StringVar(value=str(self.current_cc))
        self.var_dE = tk.StringVar(value="0.6")
        self.var_alpha = tk.StringVar(value="0.5")
        self.var_phase = tk.StringVar(value="0.07")
        self.var_ctfamp = tk.StringVar(value="auto")
        self.var_bgamp = tk.StringVar(value="auto")
        self.var_bgoff = tk.StringVar(value="auto")
        self.var_bgslope = tk.StringVar(value="auto")
        self.var_bglinslope = tk.StringVar(value="0.0")
        self.var_bgsigma = tk.StringVar(value="auto")

        tk.Button(
            self.tk_toolbar, text="LOAD IMAGE", command=self._on_load_image,
            bg="#3b82f6", fg="white", activebackground="#60a5fa", relief="flat", padx=8
        ).pack(side="left", padx=(0, 8))
        tk.Button(
            self.tk_toolbar, text="UPDATE FIT", command=self._run_fit,
            bg="#ff3366", fg="white", activebackground="#ff6699", relief="flat", padx=8
        ).pack(side="left", padx=(0, 12))

        tk.Label(self.tk_toolbar, text="Voltage", bg="#202020", fg="#dddddd").pack(side="left")
        ttk.Combobox(
            self.tk_toolbar, textvariable=self.var_kv, state="readonly", width=8,
            values=list(ELECTRON_WAVELENGTHS.keys())
        ).pack(side="left", padx=(4, 8))

        tk.Label(self.tk_toolbar, text="Binning", bg="#202020", fg="#dddddd").pack(side="left")
        ttk.Combobox(
            self.tk_toolbar, textvariable=self.var_bin, state="readonly", width=4,
            values=["1", "2", "4", "8"]
        ).pack(side="left", padx=(4, 8))

        tk.Label(self.tk_toolbar, text="DPI", bg="#202020", fg="#dddddd").pack(side="left")
        ttk.Combobox(
            self.tk_toolbar, textvariable=self.var_dpi, state="readonly", width=8,
            values=list(self._dpi_map.keys())
        ).pack(side="left", padx=(4, 8))

        tk.Label(self.tk_toolbar, text="Base Px (nm)", bg="#202020", fg="#dddddd").pack(side="left")
        ttk.Entry(self.tk_toolbar, textvariable=self.var_px, width=8).pack(side="left", padx=(4, 16))

        tk.Label(self.tk_toolbar, text="Fit Start", bg="#202020", fg="#3b82f6").pack(side="left")
        self.scale_fit_start = tk.Scale(
            self.tk_toolbar, variable=self.var_fit_start, from_=1, to=self.max_valid_radius,
            orient="horizontal", command=self._on_radius_scale_change,
            length=160, resolution=1, showvalue=False, bg="#202020", fg="#dddddd",
            troughcolor="#3a3a3a", highlightthickness=0
        )
        self.scale_fit_start.pack(side="left", padx=(4, 4))
        self.lbl_fit_start_info = tk.Label(self.tk_toolbar, text="", bg="#202020", fg="#3b82f6", width=21, anchor="w")
        self.lbl_fit_start_info.pack(side="left", padx=(0, 8))

        tk.Label(self.tk_toolbar, text="Fit End", bg="#202020", fg="#00ff88").pack(side="left")
        self.scale_fit_end = tk.Scale(
            self.tk_toolbar, variable=self.var_fit_end, from_=10, to=self.max_valid_radius,
            orient="horizontal", command=self._on_radius_scale_change,
            length=160, resolution=1, showvalue=False, bg="#202020", fg="#dddddd",
            troughcolor="#3a3a3a", highlightthickness=0
        )
        self.scale_fit_end.pack(side="left", padx=(4, 4))
        self.lbl_fit_end_info = tk.Label(self.tk_toolbar, text="", bg="#202020", fg="#00ff88", width=21, anchor="w")
        self.lbl_fit_end_info.pack(side="left", padx=(0, 8))

        tk.Label(self.tk_side, text="Fit Options", bg="#1f1f1f", fg="#ffffff", font=("Arial", 12, "bold")).pack(anchor="w", pady=(0, 6))

        bg_model_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        bg_model_row.pack(fill="x", pady=(0, 6))
        tk.Label(bg_model_row, text="Background Model", bg="#1f1f1f", fg="#e5e7eb").pack(side="left")
        self.cmb_bg_model = ttk.Combobox(
            bg_model_row, textvariable=self.var_bg_model, state="readonly", width=12,
            values=["Linear", "Exponential", "Lin + Exp", "GaussExp"]
        )
        self.cmb_bg_model.pack(side="right")
        self.cmb_bg_model.bind("<<ComboboxSelected>>", self._on_bg_model_change)

        for txt, var in [
            ("Fix Defocus", self.var_fix_defocus),
            ("Fix Cs", self.var_fix_cs),
            ("Fix Cc", self.var_fix_cc),
            ("Fix Phase", self.var_fix_phase),
            ("Fix dE", self.var_fix_dE),
            ("Fix Alpha", self.var_fix_alpha),
            ("Use Spatial Envelope", self.var_use_spatial),
            ("Full Y-Range", self.var_full_yrange),
        ]:
            tk.Checkbutton(
                self.tk_side, text=txt, variable=var,
                bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
                activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
            ).pack(fill="x", pady=1)

        self.chk_fix_bg_offset = tk.Checkbutton(
            self.tk_side, text="Fix Bg Offset", variable=self.var_fix_bg_offset,
            bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
            activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
        )
        self.chk_fix_bg_offset.pack(fill="x", pady=1)

        self.chk_fix_bg_amp = tk.Checkbutton(
            self.tk_side, text="Fix Bg Amp", variable=self.var_fix_bg_amp,
            bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
            activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
        )
        self.chk_fix_bg_amp.pack(fill="x", pady=1)

        self.chk_fix_bg_slope = tk.Checkbutton(
            self.tk_side, text="Fix Bg Slope", variable=self.var_fix_bg_slope,
            bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
            activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
        )
        self.chk_fix_bg_slope.pack(fill="x", pady=1)

        self.chk_fix_bg_linslope = tk.Checkbutton(
            self.tk_side, text="Fix Bg Lin Slope", variable=self.var_fix_bg_linslope,
            bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
            activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
        )
        self.chk_fix_bg_linslope.pack(fill="x", pady=1)
        self.chk_fix_bg_linslope.pack_forget()  # hidden unless Lin + Exp model

        self.chk_fix_bg_sigma = tk.Checkbutton(
            self.tk_side, text="Fix Bg Sigma", variable=self.var_fix_bg_sigma,
            bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
            activebackground="#1f1f1f", activeforeground="#ffffff", anchor="w"
        )
        self.chk_fix_bg_sigma.pack(fill="x", pady=1)
        self.chk_fix_bg_sigma.pack_forget()  # hidden unless GaussExp model

        tk.Label(self.tk_side, text="Initial Guesses", bg="#1f1f1f", fg="#ffffff", font=("Arial", 12, "bold")).pack(anchor="w", pady=(10, 4))

        self._make_labeled_entry(self.tk_side, "Defocus (nm)", self.var_defocus)
        self._make_labeled_entry(self.tk_side, "Cs (mm)", self.var_cs)
        self._make_labeled_entry(self.tk_side, "Cc (mm)", self.var_cc)
        self._make_labeled_entry(self.tk_side, "dE (eV)", self.var_dE)
        self._make_labeled_entry(self.tk_side, "Alpha (mrad)", self.var_alpha)
        self._make_labeled_entry(self.tk_side, "Phase (rad)", self.var_phase)
        self._make_labeled_entry(self.tk_side, "CTF Amp", self.var_ctfamp)
        self.lbl_bg_param1, _ = self._make_labeled_entry(self.tk_side, "Bg Offset", self.var_bgoff)
        self.lbl_bg_param2, _ = self._make_labeled_entry(self.tk_side, "Bg Amp", self.var_bgamp)
        self.lbl_bg_param3, _ = self._make_labeled_entry(self.tk_side, "Bg Slope", self.var_bgslope)
        self.lbl_bg_param4, _ = self._make_labeled_entry(self.tk_side, "Bg Lin Slope", self.var_bglinslope)
        self.lbl_bg_param4.master.pack_forget()  # hidden unless Lin + Exp model
        self.lbl_bg_param5, _ = self._make_labeled_entry(self.tk_side, "Bg Sigma (nm)", self.var_bgsigma)
        self.lbl_bg_param5.master.pack_forget()  # hidden unless GaussExp model

        self._on_bg_model_change()

        # ── Reimer n-start controls ──────────────────────────────────────────
        tk.Label(self.tk_side, text="Reimer n-Start", bg="#1f1f1f", fg="#ffffff",
                 font=("Arial", 11, "bold")).pack(anchor="w", pady=(10, 3))

        self.var_reimer_mode = tk.StringVar(value="Auto")
        self.var_reimer_n_start = tk.IntVar(value=1)

        reimer_radio_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        reimer_radio_row.pack(fill="x", pady=2)
        for mode_label in ("Auto", "Manual"):
            tk.Radiobutton(
                reimer_radio_row, text=mode_label,
                variable=self.var_reimer_mode, value=mode_label,
                command=self._on_reimer_mode_change,
                bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
                activebackground="#1f1f1f", activeforeground="#ffffff"
            ).pack(side="left", padx=(0, 6))

        reimer_n_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        reimer_n_row.pack(fill="x", pady=2)
        tk.Label(reimer_n_row, text="n-start value:", width=13, anchor="w",
                 bg="#1f1f1f", fg="#d1d5db").pack(side="left")
        self.spin_reimer_n = tk.Spinbox(
            reimer_n_row, from_=1, to=30,
            textvariable=self.var_reimer_n_start,
            width=5, bg="#2b2b2b", fg="white",
            buttonbackground="#3a3a3a", disabledbackground="#1a1a1a",
            disabledforeground="#555555", state="disabled"
        )
        self.spin_reimer_n.pack(side="left", padx=4)

        self.fig.canvas.draw_idle()

    def _on_reimer_mode_change(self):
        """Enable/disable the n-start spinbox based on Auto/Manual selection."""
        if self.var_reimer_mode.get() == "Manual":
            self.spin_reimer_n.config(state="normal")
        else:
            self.spin_reimer_n.config(state="disabled")

    def _make_labeled_entry(self, parent, label, variable):
        row = tk.Frame(parent, bg="#1f1f1f")
        row.pack(fill="x", pady=2)
        lbl = tk.Label(row, text=label, width=14, anchor="w", bg="#1f1f1f", fg="#d1d5db")
        lbl.pack(side="left")
        ent = ttk.Entry(row, textvariable=variable, width=14)
        ent.pack(side="right", fill="x", expand=True)
        return lbl, ent

    def _on_bg_model_change(self, _event=None):
        model = self.var_bg_model.get()
        prev_model = getattr(self, '_last_bg_model', model)
        if model == "Exponential":
            self.chk_fix_bg_offset.config(text="Fix Bg Offset")
            self.chk_fix_bg_amp.config(state="normal")
            self.chk_fix_bg_slope.config(text="Fix Bg Decay")
            self.lbl_bg_param1.config(text="Bg Offset")
            self.lbl_bg_param2.config(text="Bg Amp")
            self.lbl_bg_param3.config(text="Bg Decay (1/nm)")
            self.chk_fix_bg_linslope.pack_forget()
            self.lbl_bg_param4.master.pack_forget()
            self.chk_fix_bg_sigma.pack_forget()
            self.lbl_bg_param5.master.pack_forget()
            if self.var_bgoff.get().strip() == "":
                self.var_bgoff.set("auto")
            if self.var_bgamp.get().strip() == "":
                self.var_bgamp.set("auto")
            decay_text = self.var_bgslope.get().strip().lower()
            if decay_text == "":
                self.var_bgslope.set("1.0")
            elif prev_model in ("Linear", "Lin + Exp") and decay_text in ("auto", "0", "0.0", "-0", "-0.0"):
                self.var_bgslope.set("1.0")
        elif model == "Lin + Exp":
            self.chk_fix_bg_offset.config(text="Fix Bg Offset")
            self.chk_fix_bg_amp.config(state="normal")
            self.chk_fix_bg_slope.config(text="Fix Bg Decay")
            self.lbl_bg_param1.config(text="Bg Offset")
            self.lbl_bg_param2.config(text="Bg Amp")
            self.lbl_bg_param3.config(text="Bg Decay (1/nm)")
            self.chk_fix_bg_linslope.pack(fill="x", pady=1, after=self.chk_fix_bg_slope)
            self.lbl_bg_param4.master.pack(fill="x", pady=2, after=self.lbl_bg_param3.master)
            self.chk_fix_bg_sigma.pack_forget()
            self.lbl_bg_param5.master.pack_forget()
            if self.var_bgoff.get().strip() == "":
                self.var_bgoff.set("auto")
            if self.var_bgamp.get().strip() == "":
                self.var_bgamp.set("auto")
            decay_text = self.var_bgslope.get().strip().lower()
            if decay_text == "":
                self.var_bgslope.set("1.0")
            elif prev_model in ("Linear",) and decay_text in ("auto", "0", "0.0", "-0", "-0.0"):
                self.var_bgslope.set("1.0")
            linslope_text = self.var_bglinslope.get().strip().lower()
            if linslope_text == "":
                self.var_bglinslope.set("0.0")
        elif model == "GaussExp":
            self.chk_fix_bg_offset.config(text="Fix Bg Offset")
            self.chk_fix_bg_amp.config(state="normal")
            self.chk_fix_bg_slope.config(text="Fix Bg Xi")
            self.lbl_bg_param1.config(text="Bg Offset")
            self.lbl_bg_param2.config(text="Bg Amp")
            self.lbl_bg_param3.config(text="Bg Xi (1/nm)")
            self.chk_fix_bg_linslope.pack_forget()
            self.lbl_bg_param4.master.pack_forget()
            self.chk_fix_bg_sigma.pack(fill="x", pady=1, after=self.chk_fix_bg_slope)
            self.lbl_bg_param5.master.pack(fill="x", pady=2, after=self.lbl_bg_param3.master)
            if self.var_bgoff.get().strip() == "":
                self.var_bgoff.set("auto")
            if self.var_bgamp.get().strip() == "":
                self.var_bgamp.set("auto")
            xi_text = self.var_bgslope.get().strip().lower()
            if xi_text in ("", "0", "0.0", "-0", "-0.0"):
                self.var_bgslope.set("1.0")
            sigma_text = self.var_bgsigma.get().strip().lower()
            if sigma_text == "":
                self.var_bgsigma.set("auto")
        else:  # Linear
            self.chk_fix_bg_offset.config(text="Fix Bg Offset")
            self.chk_fix_bg_amp.config(state="disabled")
            self.var_fix_bg_amp.set(False)
            self.chk_fix_bg_slope.config(text="Fix Bg Slope")
            self.lbl_bg_param1.config(text="Bg Offset")
            self.lbl_bg_param2.config(text="Bg Amp (exp only)")
            self.lbl_bg_param3.config(text="Bg Slope")
            self.chk_fix_bg_linslope.pack_forget()
            self.lbl_bg_param4.master.pack_forget()
            self.chk_fix_bg_sigma.pack_forget()
            self.lbl_bg_param5.master.pack_forget()
            if self.var_bgoff.get().strip() == "":
                self.var_bgoff.set("auto")
            slope_text = self.var_bgslope.get().strip().lower()
            if slope_text == "":
                self.var_bgslope.set("0.0")
            elif prev_model in ("Exponential", "Lin + Exp", "GaussExp") and slope_text in ("auto", "1", "1.0"):
                self.var_bgslope.set("0.0")

        self._last_bg_model = model

    # -- Blitting helpers --------------------------------------------------------
    def _on_draw_event(self, event):
        """Cache clean FFT background after every full draw, then re-blit overlays."""
        self._bg_fft = self.fig.canvas.copy_from_bbox(self.ax_fft.bbox)
        # Schedule overlay blit after draw finishes (can't blit inside draw)
        try:
            self.fig.canvas.get_tk_widget().after(1, self._blit_fft)
        except Exception:
            pass

    def _blit_fft(self):
        """Blit only the FFT overlay artists (circles + text) — very fast."""
        if self._bg_fft is None:
            return
        self.fig.canvas.restore_region(self._bg_fft)
        for a in self._fft_overlays:
            if a.get_visible():
                self.ax_fft.draw_artist(a)
        self.fig.canvas.blit(self.ax_fft.bbox)
        self.fig.canvas.flush_events()

    def _update_slider_text(self):
        try:
            self.base_cal_nm_per_px = float(self.var_px.get())
        except ValueError: pass

        effective_px_size = self.base_cal_nm_per_px * self.bin_factor
        k_start = self.fit_start_radius / (self.N_dimension * effective_px_size)
        k_end = self.fit_end_radius / (self.N_dimension * effective_px_size)

        if hasattr(self, 'lbl_fit_start_info'):
            self.lbl_fit_start_info.config(text=f"{int(self.fit_start_radius)} px | {k_start:.4f} 1/nm")
        if hasattr(self, 'lbl_fit_end_info'):
            self.lbl_fit_end_info.config(text=f"{int(self.fit_end_radius)} px | {k_end:.4f} 1/nm")

        self.text_fft_start.set_text(f"Start: {int(self.fit_start_radius)} px  [{k_start:.4f} 1/nm]")
        self.text_fft_end.set_text(f"End: {int(self.fit_end_radius)} px  [{k_end:.4f} 1/nm]")

    def _on_load_image(self, event=None):
        try:
            from tkinter import filedialog
            # Use the existing Tk root from matplotlib's TkAgg backend
            tk_root = self.fig.canvas.get_tk_widget().winfo_toplevel()
            path = filedialog.askopenfilename(
                parent=tk_root,
                title="Select an image",
                filetypes=[("Images", "*.png *.jpg *.tif *.tiff"), ("All", "*.*")]
            )
            if path and os.path.exists(path):
                self.base_image = load_image_as_float(path)
                self.bin_factor = 2
                self.var_bin.set("2")

                # Invalidate caches
                self._bg_fft = None
                _radial_cache.clear()

                # Recompute everything with the new image
                self._apply_binning()
                self._compute_fft_and_profile()

                # Force a synchronous redraw so the new image/FFT appear immediately
                self.fig.canvas.draw()
                self.fig.canvas.flush_events()

                self.fig.canvas.manager.set_window_title(f"CTF Analyzer  v1.96  \u2014  {os.path.basename(path)}")
        except Exception as e:
            import traceback
            print(f"Error loading image: {e}")
            traceback.print_exc()

    def _apply_binning(self):
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape
        self.img_plot.set_data(self.image)
        self.img_plot.set_clim(self.image.min(), self.image.max())
        self.img_plot.set_extent([-0.5, self.w - 0.5, -0.5, self.h - 0.5])
        self.ax_img.set_xlim(0, self.w)
        self.ax_img.set_ylim(0, self.h)
        
        old_max_radius = self.max_valid_radius
        self.max_valid_radius = int(np.hypot(self.w / 2.0, self.h / 2.0))
        boundary_ratio = self.max_valid_radius / max(1, old_max_radius)

        self.scale_fit_start.config(to=self.max_valid_radius)
        self.scale_fit_end.config(to=self.max_valid_radius)
        
        new_start = max(1, int(self.fit_start_radius * boundary_ratio))
        new_end = max(2, min(self.max_valid_radius, int(self.fit_end_radius * boundary_ratio)))

        self.var_fit_start.set(new_start)
        self.var_fit_end.set(new_end)
        self._on_radius_scale_change()
        
        self._update_slider_text()

    def _compute_fft_and_profile(self):
        h_img, w_img = self.image.shape

        # Use scipy.fft (multithreaded) with fast-length padding when available
        if _HAS_SCIPY_FFT:
            fast_h = next_fast_len(h_img)
            fast_w = next_fast_len(w_img)
            fft_img = fft2(self.image, s=(fast_h, fast_w), workers=-1)
            fft_shifted = fftshift(fft_img)
        else:
            fft_img = np.fft.fft2(self.image)
            fft_shifted = np.fft.fftshift(fft_img)

        amplitude_spectrum = np.abs(fft_shifted)
        self.fft_img_log = np.log(amplitude_spectrum + 1e-10).astype(np.float32)

        h, w = amplitude_spectrum.shape
        fft_cx, fft_cy = w // 2, h // 2
        self.N_dimension = max(h_img, w_img)   # use original, not padded

        # Update existing imshow artist instead of clearing the axes
        if hasattr(self, '_fft_imshow') and self._fft_imshow is not None:
            self._fft_imshow.set_data(self.fft_img_log)
            self._fft_imshow.set_clim(self.fft_img_log.min(), self.fft_img_log.max())
            self._fft_imshow.set_extent([-0.5, w - 0.5, h - 0.5, -0.5])
            self.ax_fft.set_xlim(-0.5, w - 0.5)
            self.ax_fft.set_ylim(h - 0.5, -0.5)
        else:
            self._fft_imshow = self.ax_fft.imshow(self.fft_img_log, cmap="gray")

        self.ax_fft.set_title(f"2D Spectrum (True Bin {self.bin_factor}x)", color="w", pad=5, fontsize=12)
        self.ax_fft.axis("off")

        self.fft_start_circle.set_center((fft_cx, fft_cy))
        self.fft_end_circle.set_center((fft_cx, fft_cy))
        self.scherzer_circle.set_center((fft_cx, fft_cy))

        self.fft_start_circle.set_visible(True)
        self.fft_end_circle.set_visible(True)

        self.text_fft_start.set_text("")
        self.text_fft_end.set_text("")
        self._update_slider_text()

        max_rad = int(np.hypot(w / 2.0, h / 2.0))
        self.radii_px, self.rad_avg_fft = compute_radial_profile_fast(amplitude_spectrum, fft_cx, fft_cy, max_rad)

        # Persist centre coords so _apply_fit_results can call the marker updater
        self._fft_cx_last = fft_cx
        self._fft_cy_last = fft_cy
        self._update_fft_row_markers(fft_cx, fft_cy)

    def _update_fft_row_markers(self, fft_cx, fft_cy):
        """Find max/min rings in the radial average within the current fit ROI
        (between blue fit-start and green fit-end circles) and place red/blue
        dot markers at fft_cx ± r on the FFT x-axis for each ring radius r.

        Using the rotational average instead of a raw pixel row gives a much
        smoother, angularly-representative profile so the extrema correspond to
        genuine CTF rings rather than noise spikes in one direction.
        """
        if self.rad_avg_fft is None or self.fft_img_log is None:
            return
        h, w = self.fft_img_log.shape
        if fft_cy < 0 or fft_cy >= h:
            return

        # --- ROI in radial pixels (matching the blue/green circles) ---
        r_lo = int(self.fit_start_radius)
        r_hi = int(self.fit_end_radius)
        max_r = len(self.rad_avg_fft) - 1
        r_lo = max(1, min(r_lo, max_r - 2))
        r_hi = max(r_lo + 2, min(r_hi, max_r))

        # Slice the rotational average to the ROI only
        roi_raw = self.rad_avg_fft[r_lo : r_hi + 1].astype(np.float64)
        roi_len = len(roi_raw)

        # Adaptive minimum ring separation: ~1/20 of the ROI width
        order = max(3, roi_len // 20)

        # Two smoothing levels: narrower for maxima, wider for minima.
        # Minima in the CTF (zeros) sit on a steeply falling background and
        # are harder to isolate, so more smoothing stabilises their detection.
        ma_w_max = max(3, order // 3)          # e.g. window ≈ ROI/60
        ma_w_min = max(5, order * 2 // 3)      # e.g. window ≈ ROI/30 — 2× wider
        roi_for_max = _moving_average(roi_raw, ma_w_max)
        roi_for_min = _moving_average(roi_raw, ma_w_min)

        max_rel = np.array([], dtype=int)
        min_rel = np.array([], dtype=int)
        if roi_len > order * 2 + 2:
            max_rel, _ = _local_extrema_1d(roi_for_max, order)
            _, min_rel = _local_extrema_1d(roi_for_min, order)

        # Show the minima-finding MA on the CTF plot (wider window, so it
        # represents the more aggressive smoothing used for the harder minima).
        if hasattr(self, 'line_ma') and self.k_axis is not None and len(self.k_axis) > r_hi:
            ma_full = np.full(len(self.k_axis), np.nan)
            ma_full[r_lo : r_hi + 1] = roi_for_min
            self.line_ma.set_data(self.k_axis, ma_full)
        elif hasattr(self, 'line_ma'):
            self.line_ma.set_data([], [])

        # Absolute radii in pixels (offset back from the ROI start)
        max_radii = (r_lo + max_rel).tolist()
        min_radii = (r_lo + min_rel).tolist()

        # Direct RA k-values (no ±symmetry step) for the RA Rings table
        eff_px = self.base_cal_nm_per_px * self.bin_factor
        N_dim  = self.N_dimension
        self._ra_max_ks = [r / (N_dim * eff_px) for r in max_radii if r > 0]
        self._ra_min_ks = [r / (N_dim * eff_px) for r in min_radii if r > 0]

        # Place markers at fft_cx ± r (positive AND negative x side)
        max_cols, min_cols = [], []
        for r in max_radii:
            for col in (fft_cx + r, fft_cx - r):
                if 0 <= col < w:
                    max_cols.append(col)
        for r in min_radii:
            for col in (fft_cx + r, fft_cx - r):
                if 0 <= col < w:
                    min_cols.append(col)

        # Persist for _update_peaks_table (k-value conversion)
        self._fft_row_max_px = max_cols
        self._fft_row_min_px = min_cols
        self._fft_cx_last    = fft_cx

        if max_cols:
            mx = np.array(max_cols, dtype=float)
            self.fft_row_max.set_data(mx, np.full_like(mx, fft_cy))
        else:
            self.fft_row_max.set_data([], [])

        if min_cols:
            mn = np.array(min_cols, dtype=float)
            self.fft_row_min.set_data(mn, np.full_like(mn, fft_cy))
        else:
            self.fft_row_min.set_data([], [])

        self._update_peaks_table()

    def _update_peaks_table(self):
        """Refresh both ring-peak tables and the CTF triangle markers.

        FFT Rings table: pixel columns from FFT x-axis dots → |k|, then
        de-duplicated across the ± mirror pair.

        RA Rings table: k-values computed directly from radial-average
        extrema radii, no mirroring or deduplication needed.
        """
        if not hasattr(self, 'peaks_text'):
            return

        eff_px = self.base_cal_nm_per_px * self.bin_factor
        N = self.N_dimension
        cx = self._fft_cx_last

        def px_to_k(cols):
            if not cols:
                return np.array([])
            return np.sort(np.abs(np.array(cols, dtype=float) - cx) / (N * eff_px))

        def dedup(ks, tol=0.08):
            if len(ks) == 0:
                return []
            out = [float(ks[0])]
            for k in ks[1:]:
                if k <= 1e-9:
                    continue
                if k - out[-1] < tol * out[-1]:
                    out[-1] = (out[-1] + k) / 2.0
                else:
                    out.append(float(k))
            return out

        max_ks = dedup(px_to_k(self._fft_row_max_px))
        min_ks = dedup(px_to_k(self._fft_row_min_px))

        # ── CTF triangle markers (use FFT-deduped k values) ──────────────
        if hasattr(self, 'ctf_max_marks') and hasattr(self, 'ctf_min_marks'):
            y_tick = 0.015
            mk_arr = np.array(max_ks, dtype=float)
            mn_arr = np.array(min_ks, dtype=float)
            self.ctf_max_marks.set_data(mk_arr, np.full_like(mk_arr, y_tick))
            self.ctf_min_marks.set_data(mn_arr, np.full_like(mn_arr, y_tick))

        # ── FFT Rings table (deduped ± mirror pairs) ─────────────────────
        def _fmt_row(mk_list, nk_list, i):
            mk = f"{mk_list[i]:6.4f} {1.0/mk_list[i]:5.3f}" if i < len(mk_list) and mk_list[i] > 1e-9 else f"{'—':>6} {'—':>5}"
            nk = f"{nk_list[i]:6.4f} {1.0/nk_list[i]:5.3f}" if i < len(nk_list) and nk_list[i] > 1e-9 else f"{'—':>6} {'—':>5}"
            return mk, nk

        hdr  = ["# Mx-k  d(nm) Mn-k  d(nm)", "─" * 27]
        n_fft = max(len(max_ks), len(min_ks))
        for i in range(n_fft):
            mk, nk = _fmt_row(max_ks, min_ks, i)
            hdr.append(f"{i+1:<2} {mk} {nk}")
        self.peaks_text.set_text("\n".join(hdr))

        # ── RA Rings table (direct from rotational average extrema) ───────
        ra_max = list(self._ra_max_ks)
        ra_min = list(self._ra_min_ks)

        def _ra_section(label, ks, color_tag=""):
            rows = [label]
            rows.append("  # k(1/nm)  d(nm)")
            rows.append("  " + "─" * 18)
            if ks:
                for i, k in enumerate(ks):
                    if k > 1e-9:
                        rows.append(f"  {i+1:<2} {k:6.4f}  {1.0/k:6.3f}")
            else:
                rows.append("  (none)")
            return rows

        ra_lines = _ra_section("Maxima (MA narrow):", ra_max)
        ra_lines.append("")
        ra_lines += _ra_section("Minima (MA wide):",  ra_min)
        if hasattr(self, 'ra_peaks_text'):
            self.ra_peaks_text.set_text("\n".join(ra_lines))

        self._estimate_from_rings()

    def _estimate_from_rings(self):
        """Estimate Cs and defocus from detected ring positions via Reimer eq. 6.50.

        All detected maxima and minima are merged into one list sorted jointly
        by q.  Consecutive integers  n_start, n_start+1, n_start+2, …  are
        assigned in ascending-q order (no separate odd/even lists).

        In the code's CTF convention  chi = π λ Δz q² − ½ π Cs λ³ q⁴  with
        Δz > 0 for underfocus, setting chi = n π/2 and dividing by q² gives:

            n / q²  =  2 λ Δz  −  Cs λ³ · q²

        so  slope m = −Cs λ³  (negative)  and  intercept b = +2 λ Δz.

        Auto mode:  tries n_start = 1 … 19, keeps only physically valid results
        (Cs > 0, Δz > 0), picks the n_start giving the highest R².
        Manual mode: uses the user-supplied n_start value directly.
        """
        if not hasattr(self, 'reimer_text'):
            return

        def _clear_plot():
            if hasattr(self, 'reimer_scat_max'):
                self.reimer_scat_max.set_data([], [])
                self.reimer_scat_min.set_data([], [])
                self.reimer_fit_line.set_data([], [])
                for ann in self._reimer_annots:
                    ann.remove()
                self._reimer_annots.clear()

        ra_max = sorted(k for k in self._ra_max_ks if k > 1e-9)
        ra_min = sorted(k for k in self._ra_min_ks if k > 1e-9)
        n_rings = len(ra_max) + len(ra_min)

        if n_rings < 2:
            self.reimer_text.set_text(
                "── Reimer Eq. 6.50 ──────\n"
                "  Need ≥2 detected rings.\n"
                "  Adjust ROI or smoothing."
            )
            self.reimer_text.set_color('#888888')
            _clear_plot()
            return

        wl = self.wavelength_nm   # nm

        # ── Merge ALL extrema into one list sorted by q ───────────────────────
        all_pts = [(q, 'max') for q in ra_max] + [(q, 'min') for q in ra_min]
        all_pts.sort(key=lambda p: p[0])
        q_all  = np.array([p[0] for p in all_pts])   # 1/nm
        types  = [p[1] for p in all_pts]
        q2_all = q_all ** 2

        # ── Search over n_start values ────────────────────────────────────────
        n_start_range = (range(1, 20) if self.reimer_auto
                         else [max(1, int(self.reimer_n_start_manual))])

        best_r2 = -np.inf
        best    = None

        for n_start in n_start_range:
            n_arr  = np.arange(n_start, n_start + len(all_pts), dtype=float)
            y      = n_arr / q2_all

            coeffs = np.polyfit(q2_all, y, 1)
            m_try, b_try = float(coeffs[0]), float(coeffs[1])

            Cs_try = -m_try / (wl ** 3)      # nm  (negative slope → positive Cs)
            dz_try =  b_try / (2.0 * wl)     # nm

            # Reject physically impossible results
            if Cs_try <= 0 or dz_try <= 0:
                continue

            y_pred = m_try * q2_all + b_try
            ss_res = float(np.sum((y - y_pred) ** 2))
            ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
            r2_try = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0

            if r2_try > best_r2:
                best_r2 = r2_try
                best = dict(
                    n_start=n_start, m=m_try, b=b_try,
                    Cs_nm=Cs_try, dz_nm=dz_try, r2=r2_try,
                    n_arr=n_arr, y=y, y_pred=y_pred
                )

        if best is None:
            self.reimer_text.set_text(
                "── Reimer Eq. 6.50 ──────\n"
                "  No physically valid fit.\n"
                "  (Cs≤0 or Δz≤0 for all\n"
                "   tried n_start values)"
            )
            self.reimer_text.set_color('#ef4444')
            _clear_plot()
            return

        # ── Unpack best result ────────────────────────────────────────────────
        n_start = best['n_start']
        m, b    = best['m'], best['b']
        Cs_nm   = best['Cs_nm'];  Cs_mm = Cs_nm / 1e6
        dz_nm   = best['dz_nm']
        r2      = best['r2']
        n_arr   = best['n_arr']
        y       = best['y']
        y_pred  = best['y_pred']

        mode_tag = "(auto)" if self.reimer_auto else "(manual)"

        # Colour-code by R²
        if r2 > 0.99:
            color = '#00ff88'
        elif r2 > 0.95:
            color = '#fb923c'
        else:
            color = '#ef4444'

        # ── Text table ───────────────────────────────────────────────────────
        lines = [
            "── Reimer Eq. 6.50 ──────",
            f"  n₀={n_start} {mode_tag}",
            "",
            f"  {'Cs (mm)':<13} {Cs_mm:>8.4f}",
            f"  {'Defocus (nm)':<13} {dz_nm:>8.1f}",
            f"  {'R²':<13} {r2:>8.5f}",
            f"  {'Rings':<13} {len(q_all):>8d}",
            f"  {'λ (nm)':<13} {wl:>8.5f}",
            "",
            f"  {'n':>3}  {'type':3} {'q(1/nm)':>7} {'resid':>7}",
            "  " + "─" * 27,
        ]
        for ni, qi, ti, yi, ypi in zip(n_arr, q_all, types, y, y_pred):
            lbl = "Mx" if ti == 'max' else "Mn"
            lines.append(f"  {int(ni):>3}  {lbl:<3} {qi:>7.4f} {yi-ypi:>+7.4f}")

        self.reimer_text.set_text("\n".join(lines))
        self.reimer_text.set_color(color)

        # ── Scatter + fit plot ────────────────────────────────────────────────
        if not hasattr(self, 'reimer_scat_max'):
            return

        max_x = [qi ** 2 for qi, ti in zip(q_all, types) if ti == 'max']
        max_y = [ni / qi ** 2 for ni, qi, ti in zip(n_arr, q_all, types) if ti == 'max']
        min_x = [qi ** 2 for qi, ti in zip(q_all, types) if ti == 'min']
        min_y = [ni / qi ** 2 for ni, qi, ti in zip(n_arr, q_all, types) if ti == 'min']

        self.reimer_scat_max.set_data(max_x, max_y)
        self.reimer_scat_min.set_data(min_x, min_y)

        if len(q2_all) >= 2:
            x_lo   = max(0.0, q2_all.min() * 0.85)
            x_hi   = q2_all.max() * 1.15
            x_line = np.linspace(x_lo, x_hi, 150)
            self.reimer_fit_line.set_data(x_line, m * x_line + b)
        else:
            self.reimer_fit_line.set_data([], [])

        # Refresh n-index annotations
        for ann in self._reimer_annots:
            ann.remove()
        self._reimer_annots.clear()
        for ni, qi in zip(n_arr, q_all):
            ann = self.ax_reimer_plot.annotate(
                str(int(ni)),
                xy=(qi ** 2, ni / qi ** 2),
                xytext=(3, 2), textcoords='offset points',
                color='white', fontsize=8, zorder=6,
                annotation_clip=True
            )
            self._reimer_annots.append(ann)

        self.ax_reimer_plot.relim()
        self.ax_reimer_plot.autoscale_view()

    def _on_radius_scale_change(self, _val=None):
        self.fit_start_radius = self.var_fit_start.get()
        self.fit_end_radius = self.var_fit_end.get()
        self.fft_start_circle.set_radius(self.fit_start_radius)
        self.fft_end_circle.set_radius(self.fit_end_radius)
        self._update_slider_text()
        self._blit_fft()                   # instant circle feedback (no full redraw!)

    def _get_dynamic_model(self):
        """Build a fully self-contained CTF model closure (thread-safe, no self refs)."""
        k_all = self.k_axis
        k2_all = k_all * k_all
        k3_all = k2_all * k_all
        k4_all = k2_all * k2_all

        fix_defocus = self.fix_defocus
        fix_cs = self.fix_cs
        fix_cc = self.fix_cc
        fix_phase = self.fix_phase
        fix_dE = self.fix_dE
        fix_alpha = self.fix_alpha
        fix_bg_amp = self.fix_bg_amp
        fix_bg_offset = self.fix_bg_offset
        fix_bg_slope = self.fix_bg_slope
        bg_model = self.bg_model
        use_spatial = self.use_spatial
        try: defocus_fixed = float(self.var_defocus.get())
        except ValueError: defocus_fixed = 100.0
        cs_fixed = self.current_cs * 1e6
        cc_fixed = self.current_cc * 1e6
        try: phase_fixed = float(self.var_phase.get())
        except ValueError: phase_fixed = 0.07
        try: dE_fixed = float(self.var_dE.get())
        except ValueError: dE_fixed = 0.6
        try: alpha_fixed = float(self.var_alpha.get()) * 1e-3
        except ValueError: alpha_fixed = 0.5e-3
        bg_amp_fixed = self.bg_amp_fixed
        bg_offset_fixed = self.bg_offset_fixed
        bg_slope_fixed = self.bg_slope_fixed
        bg_linslope_fixed = self.bg_linslope_fixed
        fix_bg_linslope = self.fix_bg_linslope
        bg_sigma_fixed = self.bg_sigma_fixed
        fix_bg_sigma = self.fix_bg_sigma
        wl = self.wavelength_nm
        kv_ev = self.current_kv_ev
        pi = np.pi

        def dynamic_model(k, *free_params):
            n = len(k)
            start = 0
            if n < len(k_all):
                start = int(round((k[0] - k_all[0]) / max(k_all[1] - k_all[0], 1e-30)))
                start = max(0, min(start, len(k_all) - n))
            sl = slice(start, start + n)
            k2 = k2_all[sl]; k3 = k3_all[sl]; k4 = k4_all[sl]

            pidx = 0
            if fix_defocus: defocus = defocus_fixed
            else: defocus = free_params[pidx]; pidx += 1
            if fix_cs: cs_nm = cs_fixed
            else: cs_nm = free_params[pidx] * 1e6; pidx += 1
            if fix_cc: cc_nm = cc_fixed
            else: cc_nm = free_params[pidx] * 1e6; pidx += 1
            if fix_dE: dE = dE_fixed
            else: dE = free_params[pidx]; pidx += 1
            if fix_alpha: alpha_rad = alpha_fixed
            else: alpha_rad = free_params[pidx] * 1e-3; pidx += 1
            if fix_phase: phase = phase_fixed
            else: phase = free_params[pidx]; pidx += 1
            ctf_amp = free_params[pidx]; pidx += 1
            if fix_bg_offset: bg_offset = bg_offset_fixed
            else: bg_offset = free_params[pidx]; pidx += 1
            if bg_model in ("Exponential", "Lin + Exp", "GaussExp"):
                if fix_bg_amp: bg_amp = bg_amp_fixed
                else: bg_amp = free_params[pidx]; pidx += 1
                if fix_bg_slope: bg_decay = bg_slope_fixed
                else: bg_decay = free_params[pidx]; pidx += 1
                if bg_model == "Lin + Exp":
                    if fix_bg_linslope: bg_linslope = bg_linslope_fixed
                    else: bg_linslope = free_params[pidx]; pidx += 1
                elif bg_model == "GaussExp":
                    if fix_bg_sigma: bg_sigma = bg_sigma_fixed
                    else: bg_sigma = free_params[pidx]; pidx += 1
                    bg_linslope = 0.0
                else:
                    bg_linslope = 0.0
            else:
                bg_amp = 0.0
                bg_linslope = 0.0
                bg_sigma = 1.0
                if fix_bg_slope: bg_slope = bg_slope_fixed
                else: bg_slope = free_params[pidx]; pidx += 1

            chi = pi * wl * defocus * k2 - 0.5 * pi * cs_nm * (wl**3) * k4 + phase
            delta_nm = cc_nm * (dE / kv_ev)
            temporal = np.exp(-0.5 * (pi**2) * (wl**2) * (delta_nm**2) * k4)
            t = defocus * k - cs_nm * (wl**2) * k3
            spatial = np.exp(-(pi**2) * (alpha_rad**2) * (t * t)) if use_spatial else np.ones_like(k)
            if bg_model == "Exponential":
                bg_term = bg_offset + bg_amp * np.exp(-bg_decay * k)
            elif bg_model == "Lin + Exp":
                bg_term = bg_offset + bg_amp * np.exp(-bg_decay * k) + bg_linslope * k
            elif bg_model == "GaussExp":
                sig2 = max(bg_sigma, 1e-12) ** 2
                xi2  = max(bg_decay, 1e-12) ** 2
                bg_term = bg_offset + bg_amp * np.exp(-(k * k) / (4.0 * sig2)) / (k * k + xi2)
            else:
                bg_term = bg_offset + bg_slope * k
            return bg_term + ctf_amp * np.sin(chi) * (spatial * temporal)

        return dynamic_model

    def _run_fit(self, event=None):
        """Read UI, prepare parameters, and run curve_fit synchronously."""
        self._sync_controls_from_ui()

        if self.rad_avg_fft is None or self.radii_px is None:
            return

        try: self.base_cal_nm_per_px = float(self.var_px.get())
        except ValueError: pass

        effective_px_size = self.base_cal_nm_per_px * self.bin_factor
        self.k_axis = self.radii_px / (self.N_dimension * effective_px_size)

        start_idx = int(self.fit_start_radius)
        end_idx = int(self.fit_end_radius)

        if start_idx >= end_idx - 5:
            self.table_text.set_text("Invalid Fit Range.\nEnsure End > Start.")
            self.table_text.set_color("#ef4444")
            self.fig.canvas.draw_idle()
            return

        # Read all text box values on main thread
        try: defocus_guess = float(self.var_defocus.get())
        except ValueError: defocus_guess = 100.0
        try: self.current_cs = float(self.var_cs.get())
        except ValueError: self.current_cs = 1.2
        try: self.current_cc = float(self.var_cc.get())
        except ValueError: self.current_cc = 1.4
        try: dE_guess = float(self.var_dE.get())
        except ValueError: dE_guess = 0.6
        try: alpha_guess = float(self.var_alpha.get())
        except ValueError: alpha_guess = 0.5
        try: phase_guess = float(self.var_phase.get())
        except ValueError: phase_guess = 0.07

        k_fit = self.k_axis[start_idx:end_idx]
        data_fit = self.rad_avg_fft[start_idx:end_idx]

        ctfamp_text = self.var_ctfamp.get().strip().lower()
        if ctfamp_text == 'auto':
            ctf_amp_guess = self.rad_avg_fft[start_idx]
        else:
            try: ctf_amp_guess = float(self.var_ctfamp.get())
            except ValueError: ctf_amp_guess = self.rad_avg_fft[start_idx]

        bgamp_text = self.var_bgamp.get().strip().lower()
        bgoff_text = self.var_bgoff.get().strip().lower()
        bgslope_text = self.var_bgslope.get().strip().lower()
        bglinslope_text = self.var_bglinslope.get().strip().lower()
        bgsigma_text = self.var_bgsigma.get().strip().lower()
        if self.bg_model in ("Exponential", "Lin + Exp", "GaussExp"):
            if bgoff_text == 'auto':
                bg_offset_guess = max(0.0, np.min(data_fit) * 0.5)
            else:
                try: bg_offset_guess = max(0.0, float(self.var_bgoff.get()))
                except ValueError: bg_offset_guess = max(0.0, np.min(data_fit) * 0.5)

            if bgamp_text == 'auto':
                bg_amp_guess = max(0.0, np.max(data_fit) - bg_offset_guess)
            else:
                try: bg_amp_guess = max(0.0, float(self.var_bgamp.get()))
                except ValueError: bg_amp_guess = max(0.0, np.max(data_fit) - bg_offset_guess)

            if bgslope_text == 'auto':
                bg_slope_guess = 1.0
            else:
                try: bg_slope_guess = max(1e-6, float(self.var_bgslope.get()))
                except ValueError: bg_slope_guess = 1.0

            if self.bg_model == "Lin + Exp":
                if bglinslope_text in ('auto', ''):
                    bg_linslope_guess = 0.0
                else:
                    try: bg_linslope_guess = float(self.var_bglinslope.get())
                    except ValueError: bg_linslope_guess = 0.0
                bg_sigma_guess = 1.0
            elif self.bg_model == "GaussExp":
                # Auto-guess sigma as half the fitted k-range
                if bgsigma_text in ('auto', ''):
                    k_span = float(k_fit[-1] - k_fit[0]) if len(k_fit) > 1 else 1.0
                    bg_sigma_guess = max(1e-4, k_span * 0.5)
                else:
                    try: bg_sigma_guess = max(1e-4, float(self.var_bgsigma.get()))
                    except ValueError: bg_sigma_guess = 1.0
                # bg(0) = offset + A/Xi^2 → A ≈ amplitude * Xi^2
                bg_amp_guess = max(0.0, (np.max(data_fit) - bg_offset_guess) * bg_slope_guess ** 2)
                bg_linslope_guess = 0.0
            else:
                bg_linslope_guess = 0.0
                bg_sigma_guess = 1.0
        else:  # Linear
            if bgoff_text == 'auto':
                bg_offset_guess = max(0.0, np.min(data_fit))
            else:
                try: bg_offset_guess = max(0.0, float(self.var_bgoff.get()))
                except ValueError: bg_offset_guess = max(0.0, np.min(data_fit))

            if bgslope_text == 'auto':
                bg_slope_guess = 0.0
            else:
                try: bg_slope_guess = min(0.0, float(self.var_bgslope.get()))
                except ValueError: bg_slope_guess = 0.0
            bg_amp_guess = 0.0
            bg_linslope_guess = 0.0
            bg_sigma_guess = 1.0

        self.bg_amp_fixed = bg_amp_guess
        self.bg_offset_fixed = bg_offset_guess
        self.bg_slope_fixed = bg_slope_guess
        self.bg_linslope_fixed = bg_linslope_guess
        self.bg_sigma_fixed = bg_sigma_guess

        p0 = []
        lower_bounds = []
        upper_bounds = []
        if not self.fix_defocus:
            p0.append(defocus_guess); lower_bounds.append(0.0); upper_bounds.append(np.inf)
        if not self.fix_cs:
            p0.append(self.current_cs); lower_bounds.append(0.0); upper_bounds.append(np.inf)
        if not self.fix_cc:
            p0.append(self.current_cc); lower_bounds.append(0.0); upper_bounds.append(np.inf)
        if not self.fix_dE:
            p0.append(dE_guess); lower_bounds.append(0.0); upper_bounds.append(np.inf)
        if not self.fix_alpha:
            p0.append(alpha_guess); lower_bounds.append(0.0); upper_bounds.append(np.inf)
        if not self.fix_phase:
            p0.append(phase_guess); lower_bounds.append(-np.inf); upper_bounds.append(np.inf)
        p0.append(ctf_amp_guess)
        lower_bounds.append(0.0)
        upper_bounds.append(np.inf)
        if not self.fix_bg_offset:
            p0.append(bg_offset_guess)
            lower_bounds.append(0.0)
            upper_bounds.append(np.inf)
        if self.bg_model in ("Exponential", "Lin + Exp", "GaussExp") and not self.fix_bg_amp:
            p0.append(bg_amp_guess)
            lower_bounds.append(0.0)
            upper_bounds.append(np.inf)
        if not self.fix_bg_slope:
            p0.append(bg_slope_guess)
            if self.bg_model in ("Exponential", "Lin + Exp", "GaussExp"):
                lower_bounds.append(1e-6); upper_bounds.append(np.inf)
            else:
                lower_bounds.append(-np.inf); upper_bounds.append(0.0)
        if self.bg_model == "Lin + Exp" and not self.fix_bg_linslope:
            p0.append(bg_linslope_guess)
            lower_bounds.append(-np.inf); upper_bounds.append(np.inf)
        if self.bg_model == "GaussExp" and not self.fix_bg_sigma:
            p0.append(bg_sigma_guess)
            lower_bounds.append(1e-4); upper_bounds.append(np.inf)

        model_func = self._get_dynamic_model()

        # Snapshot everything the background thread + callback need
        ctx = dict(
            start_idx=start_idx, end_idx=end_idx,
            k_fit=k_fit, data_fit=data_fit,
            k_axis=self.k_axis.copy(), rad_avg_fft=self.rad_avg_fft.copy(),
            p0=p0, lower_bounds=lower_bounds, upper_bounds=upper_bounds,
            model_func=model_func,
            defocus_guess=defocus_guess,
            cs=self.current_cs, cc=self.current_cc,
            dE_guess=dE_guess, alpha_guess=alpha_guess,
            phase_guess=phase_guess, ctf_amp_guess=ctf_amp_guess,
            bg_offset_guess=bg_offset_guess, bg_amp_guess=bg_amp_guess, bg_slope_guess=bg_slope_guess,
            bg_linslope_guess=bg_linslope_guess,
            bg_sigma_guess=bg_sigma_guess,
            fix_defocus=self.fix_defocus, fix_cs=self.fix_cs, fix_cc=self.fix_cc,
            fix_phase=self.fix_phase, fix_dE=self.fix_dE, fix_alpha=self.fix_alpha,
            fix_bg_amp=self.fix_bg_amp, fix_bg_offset=self.fix_bg_offset, fix_bg_slope=self.fix_bg_slope,
            fix_bg_linslope=self.fix_bg_linslope, fix_bg_sigma=self.fix_bg_sigma,
            bg_model=self.bg_model,
            use_spatial=self.use_spatial,
            wavelength_nm=self.wavelength_nm, current_kv_ev=self.current_kv_ev,
            base_cal_nm_per_px=self.base_cal_nm_per_px,
            bin_factor=self.bin_factor, N_dimension=self.N_dimension,
        )

        self.table_text.set_text("Fitting...")
        self.table_text.set_color("#888888")

        mf = ctx['model_func']
        try:
            popt, pcov = curve_fit(
                mf, ctx['k_fit'], ctx['data_fit'],
                p0=ctx['p0'], bounds=(ctx['lower_bounds'], ctx['upper_bounds']),
                maxfev=3000, absolute_sigma=True, ftol=1e-6, xtol=1e-6,
                diff_step=1e-5,
            )
            ctx['fit_curve_full'] = mf(ctx['k_axis'], *popt)
            ctx['residuals_fit'] = ctx['data_fit'] - ctx['fit_curve_full'][ctx['start_idx']:ctx['end_idx']]
            dof = max(1, len(ctx['data_fit']) - len(popt))
            rchi2 = np.sum(ctx['residuals_fit']**2) / dof
            ctx['perr'] = np.sqrt(np.diag(pcov * rchi2))
            ctx['popt'] = popt
            ctx['fit_success'] = True
        except (RuntimeError, ValueError):
            ctx['fit_success'] = False
            try: ctx['fit_curve_full'] = mf(ctx['k_axis'], *ctx['p0'])
            except: ctx['fit_curve_full'] = np.full_like(ctx['k_axis'], np.nan)

        self._apply_fit_results(ctx)

    def _sync_controls_from_ui(self):
        """Apply UI control state only when UPDATE FIT is pressed."""
        # Text boxes
        try:
            self.base_cal_nm_per_px = float(self.var_px.get())
        except ValueError:
            pass

        kv_label = self.var_kv.get()
        if kv_label in ELECTRON_WAVELENGTHS:
            self.wavelength_nm = ELECTRON_WAVELENGTHS[kv_label]
            self.current_kv_ev = float(kv_label.split()[0]) * 1000.0

        self.fix_defocus = self.var_fix_defocus.get()
        self.fix_cs = self.var_fix_cs.get()
        self.fix_cc = self.var_fix_cc.get()
        self.fix_phase = self.var_fix_phase.get()
        self.fix_dE = self.var_fix_dE.get()
        self.fix_alpha = self.var_fix_alpha.get()
        self.fix_bg_amp = self.var_fix_bg_amp.get()
        self.fix_bg_offset = self.var_fix_bg_offset.get()
        self.fix_bg_slope = self.var_fix_bg_slope.get()
        self.fix_bg_linslope = self.var_fix_bg_linslope.get()
        self.fix_bg_sigma = self.var_fix_bg_sigma.get()
        self.use_spatial = self.var_use_spatial.get()
        self.full_yrange = self.var_full_yrange.get()
        self.bg_model = self.var_bg_model.get() if self.var_bg_model.get() in ("Linear", "Exponential", "Lin + Exp", "GaussExp") else "Linear"
        # Reimer n-start
        if hasattr(self, 'var_reimer_mode'):
            self.reimer_auto = (self.var_reimer_mode.get() == "Auto")
            try:
                self.reimer_n_start_manual = max(1, int(self.var_reimer_n_start.get()))
            except (ValueError, tk.TclError):
                self.reimer_n_start_manual = 1

        selected_bin = self.var_bin.get()
        try:
            selected_bin_int = int(selected_bin)
        except (TypeError, ValueError):
            selected_bin_int = self.bin_factor

        if selected_bin_int != self.bin_factor:
            self.bin_factor = selected_bin_int
            self._apply_binning()
            self._compute_fft_and_profile()

        selected_dpi = self.var_dpi.get()
        new_dpi = self._dpi_map.get(selected_dpi, self._current_dpi)
        if new_dpi != self._current_dpi:
            self._current_dpi = new_dpi
            self.fig.set_dpi(new_dpi)
            self._bg_fft = None

        self.fit_start_radius = self.var_fit_start.get()
        self.fit_end_radius = self.var_fit_end.get()

        self._update_slider_text()

    def _apply_fit_results(self, ctx):
        si, ei = ctx['start_idx'], ctx['end_idx']
        k_axis = ctx['k_axis']
        data_fit = ctx['data_fit']
        k_fit = ctx['k_fit']
        fit_curve_full = ctx['fit_curve_full']
        fit_success = ctx['fit_success']

        # --- Build popt_full dict ---
        if fit_success:
            popt, perr = ctx['popt'], ctx['perr']
            popt_full, perr_full = {}, {}
            pidx = 0
            if ctx['fix_defocus']:
                popt_full['defocus'] = ctx['defocus_guess']; perr_full['defocus'] = None
            else:
                popt_full['defocus'] = popt[pidx]; perr_full['defocus'] = perr[pidx]; pidx += 1
            if ctx['fix_cs']: popt_full['cs'] = ctx['cs']; perr_full['cs'] = None
            else: popt_full['cs'] = popt[pidx]; perr_full['cs'] = perr[pidx]; pidx += 1
            if ctx['fix_cc']: popt_full['cc'] = ctx['cc']; perr_full['cc'] = None
            else: popt_full['cc'] = popt[pidx]; perr_full['cc'] = perr[pidx]; pidx += 1
            if ctx['fix_dE']:
                popt_full['dE'] = ctx['dE_guess']; perr_full['dE'] = None
            else:
                popt_full['dE'] = popt[pidx]; perr_full['dE'] = perr[pidx]; pidx += 1
            if ctx['fix_alpha']:
                popt_full['alpha'] = ctx['alpha_guess']; perr_full['alpha'] = None
            else:
                popt_full['alpha'] = popt[pidx]; perr_full['alpha'] = perr[pidx]; pidx += 1
            if ctx['fix_phase']: popt_full['phase'] = ctx['phase_guess']; perr_full['phase'] = None
            else: popt_full['phase'] = popt[pidx]; perr_full['phase'] = perr[pidx]; pidx += 1
            popt_full['amp'] = popt[pidx]; perr_full['amp'] = perr[pidx]; pidx += 1
            if ctx['fix_bg_offset']:
                popt_full['bg'] = ctx['bg_offset_guess']; perr_full['bg'] = None
            else:
                popt_full['bg'] = popt[pidx]; perr_full['bg'] = perr[pidx]; pidx += 1
            if ctx.get('bg_model') in ("Exponential", "Lin + Exp", "GaussExp"):
                if ctx['fix_bg_amp']:
                    popt_full['bg_amp'] = ctx['bg_amp_guess']; perr_full['bg_amp'] = None
                else:
                    popt_full['bg_amp'] = popt[pidx]; perr_full['bg_amp'] = perr[pidx]; pidx += 1
            else:
                popt_full['bg_amp'] = 0.0; perr_full['bg_amp'] = None
            if ctx['fix_bg_slope']:
                popt_full['bg_slope'] = ctx['bg_slope_guess']; perr_full['bg_slope'] = None
            else:
                popt_full['bg_slope'] = popt[pidx]; perr_full['bg_slope'] = perr[pidx]; pidx += 1
            if ctx.get('bg_model') == "Lin + Exp":
                if ctx.get('fix_bg_linslope', False):
                    popt_full['bg_linslope'] = ctx.get('bg_linslope_guess', 0.0); perr_full['bg_linslope'] = None
                else:
                    popt_full['bg_linslope'] = popt[pidx]; perr_full['bg_linslope'] = perr[pidx]; pidx += 1
                popt_full['bg_sigma'] = 1.0; perr_full['bg_sigma'] = None
            elif ctx.get('bg_model') == "GaussExp":
                popt_full['bg_linslope'] = 0.0; perr_full['bg_linslope'] = None
                if ctx.get('fix_bg_sigma', False):
                    popt_full['bg_sigma'] = ctx.get('bg_sigma_guess', 1.0); perr_full['bg_sigma'] = None
                else:
                    popt_full['bg_sigma'] = popt[pidx]; perr_full['bg_sigma'] = perr[pidx]; pidx += 1
            else:
                popt_full['bg_linslope'] = 0.0; perr_full['bg_linslope'] = None
                popt_full['bg_sigma'] = 1.0; perr_full['bg_sigma'] = None

            def_e = f"{perr_full['defocus']:>9.2f}" if perr_full['defocus'] is not None else "  [Fixed]"
            cs_e = f"{perr_full['cs']:>9.3f}" if perr_full['cs'] is not None else "  [Fixed]"
            cc_e = f"{perr_full['cc']:>9.3f}" if perr_full['cc'] is not None else "  [Fixed]"
            de_e = f"{perr_full['dE']:>9.3f}" if perr_full['dE'] is not None else "  [Fixed]"
            alpha_e = f"{perr_full['alpha']:>9.3f}" if perr_full['alpha'] is not None else "  [Fixed]"
            ph_e = f"{perr_full['phase']*1000:>9.1f}" if perr_full['phase'] is not None else "  [Fixed]"
            ph_deg_e = f"{np.degrees(perr_full['phase']):>9.2f}" if perr_full['phase'] is not None else "  [Fixed]"
            bga_e = f"{perr_full['bg_amp']:>9.2e}" if perr_full['bg_amp'] is not None else "  [Fixed]"
            bg_e = f"{perr_full['bg']:>9.2e}" if perr_full['bg'] is not None else "  [Fixed]"
            bgs_e = f"{perr_full['bg_slope']:>9.2e}" if perr_full['bg_slope'] is not None else "  [Fixed]"
            table_lines = [
                f"── Fitted Parameters ─────────",
                f"  {'':.<12} {'Fitted':>9}  {'Error':>9}",
                f"  {'Defocus':<12} {popt_full['defocus']:>9.2f}  {def_e}",
                f"  {'Cs (mm)':<12} {popt_full['cs']:>9.3f}  {cs_e}",
                f"  {'Cc (mm)':<12} {popt_full['cc']:>9.3f}  {cc_e}",
                f"  {'dE (eV)':<12} {popt_full['dE']:>9.3f}  {de_e}",
                f"  {'Alpha (mrad)':<12} {popt_full['alpha']:>9.3f}  {alpha_e}",
                f"  {'Phase (mrad)':<12} {popt_full['phase']*1000:>9.1f}  {ph_e}",
                f"  {'Phase (deg)':<12} {np.degrees(popt_full['phase']):>9.2f}  {ph_deg_e}",
                f"  {'CTF Amp':<12} {popt_full['amp']:>9.2e}  {perr_full['amp']:>9.2e}",
                f"  {'Bg Offset':<12} {popt_full['bg']:>9.2e}  {bg_e}",
            ]
            if ctx.get('bg_model') in ("Exponential", "Lin + Exp", "GaussExp"):
                table_lines.append(f"  {'Bg Amp':<12} {popt_full['bg_amp']:>9.2e}  {bga_e}")
                slope_lbl = "Bg Xi (1/nm)" if ctx.get('bg_model') == "GaussExp" else "Bg Decay"
                table_lines.append(f"  {slope_lbl:<12} {popt_full['bg_slope']:>9.2e}  {bgs_e}")
            else:
                table_lines.append(f"  {'Bg Slope':<12} {popt_full['bg_slope']:>9.2e}  {bgs_e}")
            if ctx.get('bg_model') == "Lin + Exp":
                bgls_e = f"{perr_full['bg_linslope']:>9.2e}" if perr_full.get('bg_linslope') is not None else "  [Fixed]"
                table_lines.append(f"  {'Bg Lin Slp':<12} {popt_full['bg_linslope']:>9.2e}  {bgls_e}")
            elif ctx.get('bg_model') == "GaussExp":
                bgsig_e = f"{perr_full['bg_sigma']:>9.2e}" if perr_full.get('bg_sigma') is not None else "  [Fixed]"
                table_lines.append(f"  {'Bg Sigma':<12} {popt_full['bg_sigma']:>9.2e}  {bgsig_e}")
            table_str = "\n".join(table_lines)
        else:
            table_str = "── Fitted Parameters ─────────\n  Curve fit failed.\n  Adjust guesses or ROI."
            popt_full = dict(defocus=ctx['defocus_guess'], cs=ctx['cs'], cc=ctx['cc'],
                             dE=ctx['dE_guess'], alpha=ctx['alpha_guess'],
                             phase=ctx['phase_guess'], amp=ctx['ctf_amp_guess'],
                             bg_amp=ctx.get('bg_amp_guess', 0.0), bg=ctx['bg_offset_guess'],
                             bg_slope=ctx['bg_slope_guess'], bg_linslope=ctx.get('bg_linslope_guess', 0.0),
                             bg_sigma=ctx.get('bg_sigma_guess', 1.0))

        # --- Scherzer ---
        wl = ctx['wavelength_nm']
        def scherzer(cs_mm):
            if cs_mm > 0:
                cs_nm = cs_mm * 1e6
                return (1.2 * np.sqrt(cs_nm * wl),
                        1.515 * cs_nm**-0.25 * wl**-0.75,
                        1.0 / (1.515 * cs_nm**-0.25 * wl**-0.75))
            return 0., 0., 0.

        # Scherzer column: theoretical optimum from fitted Cs
        s_d, s_k, s_r = scherzer(popt_full['cs'])

        # Fit column: actual fitted defocus + its corresponding resolution
        fit_def = popt_full['defocus']
        # Point resolution at fitted defocus: k where chi crosses pi/4
        fit_cs_nm = popt_full['cs'] * 1e6
        # Information limit from total envelope at 1/e²
        fit_k = s_k   # Scherzer freq as reference
        fit_r = s_r

        self.sch_text.set_text(
            f"── Scherzer ──────────────────\n"
            f"  {'Defocus':<12} {s_d:>9.1f} nm\n"
            f"  {'Sp. Freq':<12} {s_k:>9.3f} 1/nm\n"
            f"  {'Res Limit':<12} {s_r:>9.3f} nm")
        self.sch_text.set_color("#fb923c" if fit_success else "#ef4444")

        if s_k > 0:
            r_px = s_k * ctx['N_dimension'] * ctx['base_cal_nm_per_px'] * ctx['bin_factor']
            self.scherzer_circle.set_radius(r_px)
            self.scherzer_circle.set_visible(True)
            self.scherzer_vline.set_xdata([s_k, s_k])
            self.scherzer_vline.set_visible(True)
        else:
            self.scherzer_circle.set_visible(False)
            self.scherzer_vline.set_visible(False)

        self.table_text.set_text(table_str)
        self.table_text.set_color("#fb923c" if fit_success else "#ef4444")

        # --- Envelope curves ---
        k = k_axis
        k2 = k * k; k3 = k2 * k; k4 = k2 * k2
        cc_nm = popt_full['cc'] * 1e6; cs_nm = popt_full['cs'] * 1e6
        alpha_rad = popt_full['alpha'] * 1e-3
        dn = cc_nm * (popt_full['dE'] / ctx['current_kv_ev'])
        temporal_env = np.exp(-0.5 * (np.pi**2) * (wl**2) * (dn**2) * k4)
        tt = popt_full['defocus'] * k - cs_nm * (wl**2) * k3
        spatial_env = np.exp(-(np.pi**2) * (alpha_rad**2) * (tt * tt)) if ctx['use_spatial'] else np.ones_like(k)
        _bg_model = ctx.get('bg_model', 'Linear')
        def _eval_bg(kv, pf):
            """Evaluate background curve at spatial frequencies kv."""
            if _bg_model == "Exponential":
                return pf['bg'] + pf['bg_amp'] * np.exp(-pf['bg_slope'] * kv)
            elif _bg_model == "Lin + Exp":
                return (pf['bg'] + pf['bg_amp'] * np.exp(-pf['bg_slope'] * kv)
                        + pf['bg_linslope'] * kv)
            elif _bg_model == "GaussExp":
                sig2 = max(pf['bg_sigma'], 1e-12) ** 2
                xi2  = max(pf['bg_slope'],  1e-12) ** 2
                return pf['bg'] + pf['bg_amp'] * np.exp(-(kv * kv) / (4.0 * sig2)) / (kv * kv + xi2)
            else:
                return pf['bg'] + pf['bg_slope'] * kv

        bg_line = _eval_bg(k, popt_full)
        total_env = bg_line + popt_full['amp'] * (temporal_env * spatial_env)
        spatial_curve = bg_line + popt_full['amp'] * spatial_env
        temporal_curve = bg_line + popt_full['amp'] * temporal_env

        # --- Update line data ---
        self.line_rad_fft.set_data(k_axis, ctx['rad_avg_fft'])

        fc = np.full_like(k_axis, np.nan); fc[si:ei] = fit_curve_full[si:ei]
        self.line_ctf_fit.set_data(k_axis, fc)

        bg = np.full_like(k_axis, np.nan)
        bg[si:ei] = _eval_bg(k_axis[si:ei], popt_full)
        self.line_bg.set_data(k_axis, bg)

        se_plot = np.full_like(k_axis, np.nan); se_plot[si:ei] = spatial_curve[si:ei]
        self.line_env_spatial.set_data(k_axis, se_plot)

        te_plot = np.full_like(k_axis, np.nan); te_plot[si:ei] = temporal_curve[si:ei]
        self.line_env_temporal.set_data(k_axis, te_plot)

        ev = np.full_like(k_axis, np.nan); ev[si:ei] = total_env[si:ei]
        self.line_env.set_data(k_axis, ev)

        if fit_success:
            res = ctx['residuals_fit']
            self.line_res.set_data(k_fit, res)
            rmin, rmax = res.min(), res.max()
            self.ax_res.set_ylim(rmin - max(abs(rmax - rmin) * 0.1, 1e-5),
                                 rmax + max(abs(rmax - rmin) * 0.1, 1e-5))
        else:
            self.line_res.set_data([], [])

        x_max = k_axis[min(ei + 20, len(k_axis) - 1)]
        self.ax_ctf.set_xlim(0, x_max)
        if self.full_yrange:
            self.ax_ctf.set_ylim(0, np.max(ctx['rad_avg_fft']) * 1.1)
        else:
            y_max = ctx['rad_avg_fft'][si] * 1.3
            self.ax_ctf.set_ylim(0, y_max)

        # Re-run the full marker update: ROI may have moved, calibration may
        # have changed, so recompute extrema within the new fit_start/end range
        # and refresh both the FFT dot markers and the peaks table together.
        self._update_fft_row_markers(self._fft_cx_last, self._fft_cy_last)

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
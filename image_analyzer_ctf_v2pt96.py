# -*- coding: utf-8 -*-
"""
2D Image CTF Analyzer (Spatial + Temporal Envelopes)
=============================================================

Version Changelog
-----------------
v2.97 (current)
  - Information limit formula updated to d_info = (πλΔ · 1/(2√ln2))^½
    (Reimer formula from user spec) where Δ = Cc·(ΔE/E₀); focal spread Δ
    now shown explicitly in the Info Limit panel
  - Information limit drawn as a purple dashed ring on the 2D spectrum
  - Version bumped from v2.95 → v2.97

v2.95
  - Information limit corrected to ρ_c = √(πλΔ/2) (Reimer / globalsino
    formula) where Δ = Cc·(ΔE/E₀); voltage and current fluctuations excluded
  - Spacing-consistency filter added after alternation enforcement: ring
    spacings are fit to a linear trend (monotonically shrinking with q); any
    ring producing a gap less than 0.35× the predicted spacing is removed
    (less prominent of the pair); iterates to convergence
  - Phase χ vs q plot: CTF-minimum (χ=−π/2, sin=−1) and first zero-crossing
    (χ=−π) markers drawn as explicit bounded line segments inside the plot,
    each with a dot at the intersection and a q-value readout; markers are
    correctly removed when switching back to Reimer n/q² view
  - Information limit from chromatic envelope added to the Scherzer panel:
    q_info = (2 / (π² λ² δ²))^(1/4) where δ = Cc·(ΔE/E); shown as
    spatial frequency (1/nm) and real-space distance (nm)
  - "EXPORT RINGS" toolbar button saves detected ring data (n, type, q,
    d-spacing, χ) as a two-column xlsx workbook via openpyxl

v2.94
  - C5 (5th-order spherical aberration) term removed from the CTF model,
    phase-fit model, spatial-envelope gradient, and all UI elements; the
    CTF is now χ = −π λ Δz q² + ½ π Cs λ³ q⁴  (standard 2-term Scherzer)
  - Fit-parameter errors rounded to 2 significant figures throughout;
    fitted values displayed to the same decimal precision as their error
  - Phase plot result text now states the CTF formula and defocus sign
    convention (Δz > 0 underfocus, Δz < 0 overfocus)
  - Phase χ vs q plot: when underfocused (A < 0), CTF minimum and first
    zero crossing annotated with q read-offs

v2.93
  - Boundary guard: ring detections within 5% of the fit-start or fit-end
    radius are rejected (top-hat filter is unreliable near ROI edges)
  - Reduced chi-square (χ²_red = RSS/(N-2)) added to the Reimer table
  - Parameter errors on Cs and Defocus via the χ²_red=1 trick: np.polyfit
    cov=True already scales by RSS/(N-2), so √diag(cov) gives the correct
    standard error; σ_Cs and σ_Δz are shown as ± values in the table
  - Uniform data-point error bars (σ_y = √χ²_red) drawn as gray vertical
    stems on the n/q² vs q² scatter plot
  - Color-code legend appended to the Reimer text table explaining green /
    orange / red R² thresholds

v2.93
  - Reimer n-Start index range clamped to ±15 (was ±50)
  - Phase-χ plot x-axis now always starts at spatial frequency zero; fit curve
    also drawn from q=0 so the intercept and low-q behaviour are visible

v2.91
  - Fixed Reimer n-Start spinbox in Manual mode: added command= callback
    (_on_reimer_n_change) so arrow-button clicks immediately sync
    reimer_n_start_manual and re-run _estimate_from_rings(); added <Return>
    and <FocusOut> bindings so typed values are also applied without needing
    to press UPDATE FIT

v2.90
  - DM3 calibration now reads from the same ImageList that owns the largest
    data array (two-pass: preferred ImageList first, then any calibrated entry);
    uncalibrated placeholder entries (Scale=1.0, no units) are skipped, so
    the actual 0.11431 nm/px value is found instead of the thumbnail's 1.0
  - Radial-average data line (gray) is updated directly in
    _compute_fft_and_profile so it always refreshes when a new slice is loaded,
    even if _run_fit exits early due to an invalid fit range or failed curve_fit

v2.35
  - Fit is re-run automatically whenever the DM3 slice changes (◀/▶ or slider)
  - Fit is re-run automatically when the "Base Px (nm)" entry is committed
    (press Return or click away), so k-axis and CTF curve always reflect the
    current pixel size without needing to press "UPDATE FIT"
  - Default pixel size is now 1.0 nm/px before any image is loaded
  - DM3 Scale tag now also parsed as float64 (type 7) in addition to float32;
    sorted iteration over cal_scales ensures Dimension[0] (X-axis) is always
    preferred; zero-valued calibrations are ignored so a meaningful value is
    never overwritten by an uncalibrated dimension

v2.30
  - "Fit Range" Auto/Manual radio in toolbar: Auto disables the sliders and
    auto-computes fit-start/end from the first and last prominent peak in the
    radial-average profile (80 %–115 % margins); Manual re-enables sliders for
    full user control; auto-range re-runs whenever a new image/slice is loaded
  - Fixed "Reimer n-Start Manual": switching the radio now immediately updates
    self.reimer_auto, syncs the spinbox value into self.reimer_n_start_manual,
    and re-runs _estimate_from_rings() so the table reflects the chosen n₀

v2.25
  - Arrow ◀/▶ buttons flanking the DM3 slice slider for one-click step navigation
  - Radio button "Reimer Auto" vs "Manual" initial guesses: in Auto mode the
    Defocus and Cs entries are auto-populated from the Reimer linear fit; in
    Manual mode the user's entered values are never overwritten by the fit
  - Fit-start (blue dashed) and fit-end (green dashed) vertical lines overlaid
    on the CTF plot, updated live as the fit-range sliders move

v2.20
  - Accelerating voltage auto-read from DM3: ImageList[N].ImageTags.Microscope
    Info.Voltage (float32, stored in volts); rounded to nearest standard value
    (30/60/80/120/200/300 kV) and pushed into the Voltage combobox so the
    correct electron wavelength is applied immediately on load
  - Cs auto-read from DM3: ImageList[N].ImageTags.Microscope Info.Cs(mm)
    (float32 or float64); pushed into the Cs entry and current_cs state
  - Reimer n/q² vs q² scatter plot axis limits now set explicitly from the
    actual data range (±15 % padding) instead of relying on relim/autoscale
    which failed to rescale Line2D artists; scatter points now always visible

v2.15
  - Reimer panel no longer goes blank when no physically valid fit exists:
    the n/q² scatter points (maxima red ▼, minima blue ▲) are always plotted
    using n_start=1 as a display fallback; the text table also lists every
    detected ring with its q (1/nm) and n/q² value so they can be read off
    even without a converged fit; "No physically valid fit" header remains in
    red to make the status clear
  - Pixel size (nm/px) extraction from DM3 calibration tags is attempted on
    every load; if the DM3 stores Scale=1.0 / blank units (uncalibrated file,
    as with the focal-series stack) the field retains its previous value so
    the user's manual entry is not overwritten; default startup value remains
    0.11431 nm/px (correct for the 300 kV 80k× focal series)

v2.10
  - FFT always taken on the largest possible centre-square crop of the (binned)
    image: min(H, W) × min(H, W) pixels extracted from the image centre;
    a yellow dashed rectangle overlaid on the Real Space Image panel shows
    exactly which region is used, updating whenever the image or binning changes
  - N_dimension now equals the square side rather than max(H, W), so all
    k-axis and spatial-frequency calculations are correct for non-square images
  - max_valid_radius capped at half the square side (inscribed circle) instead
    of the full diagonal, preventing the fit-ROI from reaching outside the FFT
  - Left/right arrow keys step the DM3 stack slice slider ±1 frame (keyboard
    focus must be on the matplotlib figure window)
  - Pixel size (nm/px) is now read automatically from the DM3 calibration tags
    (Calibrations.Dimension[0].Scale / .Units) and pushed into the "Base Px (nm)"
    entry on load; units in nm, µm, Å, pm, m, mm are all handled

v2.00
  - DM3 file support: LOAD IMAGE now accepts Gatan Digital Micrograph 3 (.dm3)
    files in addition to PNG/JPG/TIFF; a pure-Python binary walker (_dm3_read)
    parses the DM3 tag tree, locates the largest float32 Data array, and reads
    the Dimensions group to determine the image shape without any external library
  - Stack / focal-series support: when a DM3 file contains a 3-D stack
    (e.g. shape n_slices × ny × nx) a "Slice" slider appears in the toolbar;
    dragging it extracts the selected frame, normalises it to [0, 1], and
    immediately recomputes the binning, FFT, and radial profile so every frame
    in a focal series can be analysed interactively
  - Slice slider is hidden automatically when a 2-D image (DM3 or otherwise) is
    loaded, keeping the toolbar uncluttered for normal use

v1.99
  - Reimer Cs and defocus estimates are automatically pushed into the CTF fit
    parameter input fields (Defocus (nm) and Cs (mm)) each time the Reimer fit
    updates, so they become the starting values for the next UPDATE FIT run;
    internal current_cs state is kept in sync so fixed-Cs mode also picks up
    the new value immediately

v1.98
  - Interactive Reimer plot point selection: left-click any point to toggle it
    in/out of the linear fit; included points are ringed with a white circle,
    excluded points lose their circle and the fit updates immediately; a
    "Clear" button in the Reimer controls resets all exclusions; the excluded
    set resets automatically whenever the detected ring count changes
  - Excluded points are still plotted (as dim markers) so their position
    relative to the fit line remains visible

v1.97
  - Replaced moving-average ring detection entirely with morphological top-hat:
    white_tophat (signal − morphological opening) isolates local peaks;
    black_tophat (morphological closing − signal) isolates local troughs;
    both applied with a flat structuring element of ~2 ring-period width so
    the slowly-varying background is fully removed before extrema are sought
  - Top-hat background estimate (morphological opening) plotted as a solid
    dark-orange line on the CTF plot for background-removal quality inspection
  - Detected ring positions shown with red ▼ (maxima) and blue ▲ (minima)
    tick-marks along the bottom of the CTF plot
  - Dashed vertical lines drawn from each detected ring marker up to the data
    curve: red dashes for maxima, blue dashes for minima — making it visually
    unambiguous which data points correspond to each ring
  - RA Rings table simplified to show only top-hat maxima and minima sections
  - Reimer estimation uses top-hat extrema exclusively

v1.96
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

try:
    from scipy.ndimage import white_tophat, black_tophat, grey_opening
    _HAS_MORPHOLOGY = True
except ImportError:
    _HAS_MORPHOLOGY = False

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

    if ext == '.dm3':
        data, _, _cal, _kv, _cs = _dm3_read(path)
        if data.ndim == 3:
            data = data[0]
        data = data.astype(np.float64)
        dmin, dmax = data.min(), data.max()
        if dmax - dmin > 0:
            data = (data - dmin) / (dmax - dmin)
        else:
            data = np.zeros_like(data)
        return data

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


def _dm3_read(path):
    """Read a DM3 file and return the image data plus spatial calibration.

    Handles both 2-D images and 3-D stacks (e.g. focal series).

    Returns
    -------
    data       : ndarray  float32.
                 Shape (ny, nx) for a single image or (n_slices, ny, nx) for a stack.
    n_slices   : int   1 for a 2-D image; number of frames for a 3-D stack.
    cal_nm_per_px : float  pixel size in nm (0.0 if not found in file).
    kv            : float  accelerating voltage in kV (0.0 if not found).
    cs_mm         : float  Cs in mm (−1.0 if not found).
    """
    import struct as _st

    with open(path, 'rb') as fh:
        blob = fh.read()

    p = [12]   # mutable position pointer; skip 12-byte file header
    ESIZES = {2: 2, 3: 4, 4: 2, 5: 4, 6: 4, 7: 8, 8: 1, 9: 1, 10: 1, 11: 8, 12: 8}

    # Maps DM3 type code → (numpy dtype string, bytes-per-element)
    # Covers every simple numeric type GMS may store image data as
    _NP_DTYPE = {
        2:  ('<i2', 2),   # int16
        3:  ('<i4', 4),   # int32
        4:  ('<u2', 2),   # uint16  ← most common for TEM images
        5:  ('<u4', 4),   # uint32
        6:  ('<f4', 4),   # float32
        7:  ('<f8', 8),   # float64
        8:  ('u1',  1),   # bool8 / uint8
        9:  ('i1',  1),   # int8
        10: ('u1',  1),   # uint8
    }

    best = [None]          # (offset, n_elem, etype) of the largest numeric Data tag
    dims_by_parent = {}    # parent_path -> [d0, d1, ...] for each Dimensions group
    cal_scales = {}        # "ImageList[N].ImageData.Calibrations.Dimension[D]" -> scale_nm
    microscope = {}        # 'Voltage' -> volts (float),  'Cs' -> mm (float)

    def rb(n):    v = blob[p[0]:p[0]+n]; p[0] += n; return v
    def ru8():    v = blob[p[0]];        p[0] += 1; return v
    def ri16be(): v, = _st.unpack_from('>h', blob, p[0]); p[0] += 2; return v
    def ri32be(): v, = _st.unpack_from('>i', blob, p[0]); p[0] += 4; return v
    def ru32be(): v, = _st.unpack_from('>I', blob, p[0]); p[0] += 4; return v

    def skip_val(descs):
        enc = descs[0]
        if enc == 20:                         # typed array
            etype = descs[1]; alen = descs[-1]
            if etype == 15:
                esz = sum(ESIZES.get(descs[4 + i*2 + 1], 0) for i in range(descs[3]))
            else:
                esz = ESIZES.get(etype, 0)
            p[0] += esz * alen
        elif enc == 15:                       # struct
            p[0] += sum(ESIZES.get(descs[3 + i*2 + 1], 0) for i in range(descs[2]))
        elif enc == 18:                       # variable-length string
            n = ru32be(); p[0] += n * 2
        else:
            p[0] += ESIZES.get(enc, 0)

    def walk_group(path):
        ru8(); ru8()                          # sorted / open flag bytes
        n = ru32be()
        anon = [0]
        for _ in range(n):
            walk_entry(path, anon)

    def walk_entry(parent, anon):
        ttype = ru8()
        nl    = ri16be()
        if nl > 0:
            raw = rb(nl).decode('latin-1', errors='replace')
            path = f'{parent}.{raw}' if parent else raw
        else:
            raw  = ''
            path = f'{parent}[{anon[0]}]' if parent else f'[{anon[0]}]'
            anon[0] += 1

        if ttype == 20:
            walk_group(path)
        elif ttype == 21:
            magic = rb(4)
            if magic != b'%%%%':
                raise ValueError(f'DM3 bad magic {magic!r} at offset {p[0]-4}')
            nd    = ri32be()
            descs = [ri32be() for _ in range(nd)]
            if not descs:
                return
            if (raw == 'Data' and len(descs) == 3
                    and descs[0] == 20 and descs[1] in _NP_DTYPE):
                # Any simple numeric array (int16, uint16, float32, …)
                etype = descs[1]
                esz   = _NP_DTYPE[etype][1]
                n_el  = descs[2]
                nbytes = n_el * esz
                if best[0] is None or n_el > best[0][1]:
                    best[0] = (p[0], n_el, etype, parent)  # include parent path
                p[0] += nbytes
            elif (raw == '' and len(descs) == 1
                  and descs[0] in (2, 3, 4, 5)   # int16 / int32 / uint16 / uint32
                  and (parent.endswith('.Dimensions') or parent == 'Dimensions')):
                # Anonymous integer inside an ImageData.Dimensions group
                esz_d = ESIZES[descs[0]]
                fmt   = {2: '<h', 3: '<i', 4: '<H', 5: '<I'}[descs[0]]
                v,    = _st.unpack_from(fmt, blob, p[0])
                p[0] += esz_d
                dims_by_parent.setdefault(parent, []).append(int(v))
            elif (raw == 'Scale' and len(descs) == 1 and descs[0] in (6, 7)
                  and 'Calibrations.Dimension' in parent):
                # float32 (type 6) or float64 (type 7) Scale tag
                if descs[0] == 6:
                    v, = _st.unpack_from('<f', blob, p[0]); p[0] += 4
                else:
                    v, = _st.unpack_from('<d', blob, p[0]); p[0] += 8
                cal_scales[parent] = float(v)
            elif (raw == 'Units' and 'Calibrations.Dimension' in parent):
                # Units string — several encoding variants exist in DM3 files.
                # Variant A: array-of-uint16 (descs=[20, 4, n_chars])
                if len(descs) == 3 and descs[0] == 20 and descs[1] == 4:
                    n_chars = descs[2]
                    raw_str = blob[p[0]:p[0] + n_chars * 2].decode('utf-16-le', errors='replace').rstrip('\x00')
                    p[0] += n_chars * 2
                    cal_scales[parent + '.__units__'] = raw_str
                # Variant B: simple string type 18 (length already consumed — skip)
                else:
                    skip_val(descs)
            elif (raw == 'Voltage' and len(descs) == 1 and descs[0] == 6
                  and 'Microscope Info' in parent):
                # Accelerating voltage — float32, stored in volts
                v, = _st.unpack_from('<f', blob, p[0]); p[0] += 4
                microscope['Voltage'] = float(v)
            elif (raw == 'Cs(mm)' and len(descs) == 1 and descs[0] in (6, 7)
                  and 'Microscope Info' in parent):
                # Cs — float32 or float64, stored in mm
                if descs[0] == 6:
                    v, = _st.unpack_from('<f', blob, p[0]); p[0] += 4
                else:
                    v, = _st.unpack_from('<d', blob, p[0]); p[0] += 8
                microscope['Cs'] = float(v)
            else:
                skip_val(descs)

    walk_group('')

    if best[0] is None:
        raise ValueError('No numeric image data found in DM3 file.')

    offset, n_elem, etype, best_parent = best[0]
    np_dtype, esz = _NP_DTYPE[etype]
    data = np.frombuffer(blob[offset:offset + n_elem * esz], dtype=np_dtype).copy().astype(np.float32)

    # Find the Dimensions set whose element product matches the data size.
    # Prefer longer dim lists (3-D over 2-D) when multiple sets match.
    best_dims = None
    for dims in dims_by_parent.values():
        prod = 1
        for d in dims:
            prod *= d
        if prod == n_elem:
            if best_dims is None or len(dims) > len(best_dims):
                best_dims = dims

    if best_dims is not None:
        # DM3 stores dimensions fastest-first: [nx, ny] or [nx, ny, nz]
        # C-order on-disk layout is therefore (ny, nx) or (nz, ny, nx)
        if len(best_dims) == 2:
            nx, ny = best_dims
            shaped = data.reshape(ny, nx)
        elif len(best_dims) >= 3:
            nx, ny, nz = best_dims[0], best_dims[1], best_dims[2]
            shaped = data.reshape(nz, ny, nx)
        else:
            shaped = None
    else:
        shaped = None

    if shaped is None:
        # Fallback: try square 2-D
        side = int(np.sqrt(n_elem))
        if side * side == n_elem:
            shaped = data.reshape(side, side)
        else:
            raise ValueError(
                f'_dm3_read: could not determine image shape for {n_elem:,} elements. '
                f'Dimension sets found: {list(dims_by_parent.values())}'
            )

    n_out = 1 if shaped.ndim == 2 else shaped.shape[0]

    # ── Extract kV and Cs ─────────────────────────────────────────────────────
    kv_out   = microscope.get('Voltage', 0.0) / 1000.0   # V → kV
    cs_mm_out = microscope.get('Cs', -1.0)

    # ── Extract pixel size in nm from calibration tags ────────────────────────
    # DM3 Dimension[0] = X (fastest axis = columns); Dimension[1] = Y (rows).
    # We want the X spatial scale converted to nm.
    cal_nm_per_px = 0.0
    unit_factors = {'nm': 1.0, 'um': 1e3, 'µm': 1e3, 'Å': 0.1, 'A': 0.1,
                    'pm': 0.001, 'm': 1e9, 'mm': 1e6}

    # Determine which ImageList owns the best (largest) data array so we can
    # prefer its calibration over a thumbnail / preview ImageList[0].
    import re as _re
    _m = _re.match(r'(ImageList\[\d+\])', best_parent)
    best_il_prefix = _m.group(1) if _m else ''

    def _try_cal(keys_iterable):
        """Return calibrated nm/px from the first valid Dimension[0] key."""
        for key in keys_iterable:
            if '__units__' in key or 'Dimension[0]' not in key:
                continue
            scale = cal_scales[key]
            units_key = key + '.__units__'
            units = cal_scales.get(units_key, '').strip()
            # Skip uncalibrated DM3 placeholder (scale==1.0, no units)
            if abs(scale - 1.0) < 1e-9 and not units:
                continue
            factor = unit_factors.get(units) or unit_factors.get(units, None)
            if factor is None:
                factor = 1.0  # assume nm if units unrecognised
            val = abs(scale) * factor
            if val > 0.0:
                return val
        return 0.0

    # Pass 1: look only within the same ImageList as the data block
    if best_il_prefix:
        cal_nm_per_px = _try_cal(
            k for k in sorted(cal_scales.keys()) if k.startswith(best_il_prefix)
        )
    # Pass 2: fall back to any calibrated entry in the file
    if cal_nm_per_px == 0.0:
        cal_nm_per_px = _try_cal(sorted(cal_scales.keys()))

    return shaped, n_out, cal_nm_per_px, kv_out, cs_mm_out


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

        self.bin_factor = 1
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape

        self.max_valid_radius = int(np.hypot(self.w / 2.0, self.h / 2.0))
        self.fit_start_radius = min(90.0, max(1, self.max_valid_radius - 10))
        self.fit_end_radius = self.max_valid_radius

        self.base_cal_nm_per_px = 1.0
        self.wavelength_nm = ELECTRON_WAVELENGTHS["80 kV"]
        self.current_kv_ev = 80000.0

        # Fixing States
        self.fix_defocus = False
        self.fix_cs = True
        self.fix_cc = True
        self.fix_phase = False
        self.fix_dE = True
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
        self.bg_model = "Lin + Exp"
        self.reimer_auto = True
        self.reimer_n_start_manual = 1
        self.reimer_push_guesses = True   # True = Reimer auto-fills Defocus/Cs
        self._reimer_excluded  = set()

        # DM3 stack state (populated when a 3-D DM3 file is loaded)
        self._dm3_stack = None           # ndarray (n_slices, ny, nx) float32, or None
        self._dm3_n_slices = 1
        self._dm3_current_slice = 0   # indices into current q_all to skip
        self._reimer_pt_count  = 0       # last ring count; reset excluded on change
        self._reimer_all_q2    = np.array([])
        self._reimer_all_y     = np.array([])
        # Last ring data cached for the phase-χ view
        self._last_ring_q_all  = np.array([])
        self._last_ring_n_arr  = np.array([])
        self._last_ring_types  = []
        self._last_ring_incl   = np.array([], dtype=bool)
        self._last_ring_wl     = 0.0
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
        # --- Top-hat extrema (k values in 1/nm) ---
        self._th_max_ks = []
        self._th_min_ks = []

        # --- Blitting state ---
        self._bg_fft = None
        self._fft_overlays = []

        self._build_ui()

        # Mark FFT overlays as animated (excluded from full draw → enables blitting)
        self._fft_overlays = [
            self.fft_start_circle, self.fft_end_circle, self.scherzer_circle,
            self.info_limit_circle,
            self.text_fft_start, self.text_fft_end,
            self.fft_row_max, self.fft_row_min,
        ]
        for a in self._fft_overlays:
            a.set_animated(True)

        # Re-cache FFT background + blit overlays after every full draw
        self.fig.canvas.mpl_connect('draw_event', self._on_draw_event)

        # Arrow keys step through DM3 stack slices
        self.fig.canvas.mpl_connect('key_press_event', self._on_key_press)

        self._compute_fft_and_profile()
        self._update_slider_text()
        self._run_fit()

        plt.show()

    def _build_ui(self):
        self.fig = plt.figure(figsize=(18, 10), facecolor="#2d2d2d", dpi=80)
        self.fig.canvas.manager.set_window_title("CTF Analyzer  v2.97")
        self.fig.text(0.01, 0.99, "v2.97", color="#00ff88", fontsize=12,
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

        # Yellow dashed rectangle showing the square crop used for the FFT
        from matplotlib.patches import Rectangle as _Rect
        self._fft_crop_rect = _Rect((0, 0), 1, 1, linewidth=1.5,
                                    edgecolor='#ffdd00', facecolor='none',
                                    linestyle='--', visible=False)
        self.ax_img.add_patch(self._fft_crop_rect)

        # 2. 2D FFT Image
        self.ax_fft = self.fig.add_subplot(gs_top[0, 1])
        self.ax_fft.set_facecolor("#383838")
        self.ax_fft.set_title(f"2D Spectrum (True Bin {self.bin_factor}x)", color="w", pad=5, fontsize=12)
        self.ax_fft.axis("off")
        
        self.fft_start_circle = Circle((0, 0), self.fit_start_radius, fill=False, edgecolor="#3b82f6", linewidth=1.5, linestyle="--")
        self.fft_end_circle = Circle((0, 0), self.fit_end_radius, fill=False, edgecolor="#00ff88", linewidth=1.5, linestyle="--")
        self.scherzer_circle = Circle((0, 0), 10.0, fill=False, edgecolor="red", linewidth=1.5, linestyle=":")
        self.info_limit_circle = Circle((0, 0), 10.0, fill=False, edgecolor="#a855f7", linewidth=1.5, linestyle="--")

        self.ax_fft.add_patch(self.fft_start_circle)
        self.ax_fft.add_patch(self.fft_end_circle)
        self.ax_fft.add_patch(self.scherzer_circle)
        self.ax_fft.add_patch(self.info_limit_circle)
        self.scherzer_circle.set_visible(False)
        self.info_limit_circle.set_visible(False)

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
        # Light moving average fed into the top-hat filter
        self.line_th_ma, = self.ax_ctf.plot([], [], color="#44dd99", linestyle="--",
                                             linewidth=1.1, alpha=0.85,
                                             label="TH Input (MA)", zorder=4)
        # Top-hat background estimate (morphological opening)
        self.line_th_open, = self.ax_ctf.plot([], [], color="#cc6600", linestyle="-",
                                               linewidth=1.3, alpha=0.85,
                                               label="TH Background", zorder=4)
        self.line_bg, = self.ax_ctf.plot([], [], color="yellow", linestyle="--", linewidth=1.5, label="Background")
        self.line_env_spatial, = self.ax_ctf.plot([], [], color="#ff9f1c", linestyle="--", linewidth=1.2, alpha=0.8, label="Spatial Env")
        self.line_env_temporal, = self.ax_ctf.plot([], [], color="#a78bfa", linestyle="--", linewidth=1.2, alpha=0.8, label="Temporal Env")
        self.line_env, = self.ax_ctf.plot([], [], color="cyan", linestyle="-", linewidth=1.5, label="Total Env")
        self.line_ctf_fit, = self.ax_ctf.plot([], [], color="#ff3366", label="CTF Fit", linewidth=2.0)
        self.scherzer_vline = self.ax_ctf.axvline(x=0, color='red', linestyle=':', linewidth=1.5, label='Scherzer Freq')
        self.fit_start_vline = self.ax_ctf.axvline(x=0, color='#3b82f6', linestyle='--', linewidth=1.2, alpha=0.8, label='Fit Start')
        self.fit_end_vline   = self.ax_ctf.axvline(x=0, color='#00ff88', linestyle='--', linewidth=1.2, alpha=0.8, label='Fit End')
        self.ax_ctf.legend(loc="upper right", facecolor="#2d2d2d", edgecolor="#333", labelcolor="#cccccc", borderpad=0.2, fontsize=10)

        # Triangle tick-marks along the bottom of the CTF plot for ring peaks.
        # Blended transform: x in data (spatial-frequency) coords, y in axes
        # fraction so the markers sit at a fixed distance above the x-axis
        # regardless of y-scale changes.
        _btrans = blended_transform_factory(self.ax_ctf.transData, self.ax_ctf.transAxes)
        # Top-hat detections — red ▼ for maxima, blue ▲ for minima
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
        # Container for vertical dashed lines from markers to data curve
        self._th_vlines = []

        self.ax_res = self.fig.add_subplot(gs_ctf[1], sharex=self.ax_ctf, facecolor="#383838")
        self.ax_res.set_xlabel(r"Spatial Frequency k (1/nm)", color="#aaa", labelpad=2, fontsize=12)
        self.ax_res.tick_params(colors="#666", labelsize=12)
        self.ax_res.set_ylabel("Resid", color="#aaa", fontsize=12)
        # Raw CTF data shown in the residuals panel (scaled to residual range)
        self.line_res_data, = self.ax_res.plot([], [], color="gray", linewidth=0.8,
                                                alpha=0.5, zorder=1, label="Data")
        self.line_res, = self.ax_res.plot([], [], color="#00ff88", linewidth=1.0,
                                           zorder=2, label="Residual")
        self.ax_res.axhline(0, color="gray", linestyle="--", linewidth=1.0, zorder=3)

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

        # Hollow white circles drawn around every point that is IN the fit
        self._reimer_incl_circles, = self.ax_reimer_plot.plot(
            [], [], 'o', color='white', markersize=14, zorder=7,
            markerfacecolor='none', markeredgewidth=1.5, linestyle='none'
        )
        # Error bar stems (vertical lines, one per scatter point, uniform σ_y)
        self._reimer_errorbars, = self.ax_reimer_plot.plot(
            [], [], '-', color='#aaaaaa', linewidth=1.2, zorder=4, alpha=0.75
        )
        self._reimer_hyperbola_lines = []  # light reference curves, rebuilt each update

        # ── Phase-χ view artists (initially hidden) ──────────────────────────
        self.phase_scat_max, = self.ax_reimer_plot.plot(
            [], [], 'v', color='red', markersize=7, zorder=5,
            markeredgewidth=0.6, markeredgecolor='#550000',
            label='Max (odd n)', visible=False
        )
        self.phase_scat_min, = self.ax_reimer_plot.plot(
            [], [], '^', color='#3399ff', markersize=7, zorder=5,
            markeredgewidth=0.6, markeredgecolor='#002255',
            label='Min (even n)', visible=False
        )
        self.phase_fit_curve, = self.ax_reimer_plot.plot(
            [], [], '-', color='#ffbb44', linewidth=1.8, zorder=4,
            label='Quadratic fit', visible=False
        )
        self._phase_incl_circles, = self.ax_reimer_plot.plot(
            [], [], 'o', color='white', markersize=14, zorder=7,
            markerfacecolor='none', markeredgewidth=1.5,
            linestyle='none', visible=False
        )
        self._phase_result_text = self.ax_reimer_plot.text(
            0.03, 0.97, '', transform=self.ax_reimer_plot.transAxes,
            color='#ffcc66', va='top', ha='left', fontsize=9,
            fontfamily='monospace', visible=False, zorder=10
        )

        # Connect click handler for point inclusion/exclusion toggling
        self.fig.canvas.mpl_connect('button_press_event', self._on_reimer_click)

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

        self.var_kv = tk.StringVar(value="80 kV")
        self.var_bin = tk.StringVar(value=str(self.bin_factor))
        self.var_dpi = tk.StringVar(value="Low")
        self.var_px = tk.StringVar(value="1.0")

        self.var_fit_start = tk.IntVar(value=int(self.fit_start_radius))
        self.var_fit_end = tk.IntVar(value=int(self.fit_end_radius))

        self.var_fix_cs = tk.BooleanVar(value=self.fix_cs)
        self.var_fix_cc = tk.BooleanVar(value=self.fix_cc)
        self.var_fix_phase = tk.BooleanVar(value=self.fix_phase)
        self.var_fix_dE = tk.BooleanVar(value=True)
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
        self.var_dE = tk.StringVar(value="1.2")
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
        ).pack(side="left", padx=(0, 8))
        tk.Button(
            self.tk_toolbar, text="EXPORT RINGS", command=self._on_export_rings,
            bg="#22aa66", fg="white", activebackground="#44cc88", relief="flat", padx=8
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
        self.entry_px = ttk.Entry(self.tk_toolbar, textvariable=self.var_px, width=8)
        self.entry_px.pack(side="left", padx=(4, 16))
        self.entry_px.bind("<Return>",   self._on_px_commit)
        self.entry_px.bind("<FocusOut>", self._on_px_commit)

        # Fit Range: Auto (best-guess from profile) vs Manual (user sliders)
        tk.Label(self.tk_toolbar, text="Fit Range:", bg="#202020", fg="#dddddd").pack(side="left", padx=(8, 2))
        self.var_fit_range_mode = tk.StringVar(value="Auto")
        for _frl in ("Auto", "Manual"):
            tk.Radiobutton(
                self.tk_toolbar, text=_frl,
                variable=self.var_fit_range_mode, value=_frl,
                command=self._on_fit_range_mode_change,
                bg="#202020", fg="#e5e7eb", selectcolor="#3a3a3a",
                activebackground="#202020", activeforeground="#ffffff",
                indicatoron=1
            ).pack(side="left", padx=(0, 3))

        tk.Label(self.tk_toolbar, text="Fit Start", bg="#202020", fg="#3b82f6").pack(side="left", padx=(6, 0))
        self.scale_fit_start = tk.Scale(
            self.tk_toolbar, variable=self.var_fit_start, from_=1, to=self.max_valid_radius,
            orient="horizontal", command=self._on_radius_scale_change,
            length=160, resolution=1, showvalue=False, bg="#202020", fg="#dddddd",
            troughcolor="#3a3a3a", highlightthickness=0, state="disabled"
        )
        self.scale_fit_start.pack(side="left", padx=(4, 4))
        self.lbl_fit_start_info = tk.Label(self.tk_toolbar, text="", bg="#202020", fg="#3b82f6", width=21, anchor="w")
        self.lbl_fit_start_info.pack(side="left", padx=(0, 8))

        tk.Label(self.tk_toolbar, text="Fit End", bg="#202020", fg="#00ff88").pack(side="left")
        self.scale_fit_end = tk.Scale(
            self.tk_toolbar, variable=self.var_fit_end, from_=10, to=self.max_valid_radius,
            orient="horizontal", command=self._on_radius_scale_change,
            length=160, resolution=1, showvalue=False, bg="#202020", fg="#dddddd",
            troughcolor="#3a3a3a", highlightthickness=0, state="disabled"
        )
        self.scale_fit_end.pack(side="left", padx=(4, 4))
        self.lbl_fit_end_info = tk.Label(self.tk_toolbar, text="", bg="#202020", fg="#00ff88", width=21, anchor="w")
        self.lbl_fit_end_info.pack(side="left", padx=(0, 8))

        # --- DM3 stack slice selector (hidden until a multi-slice DM3 is loaded) ---
        self.dm3_slice_frame = tk.Frame(self.tk_toolbar, bg="#202020")
        self.dm3_slice_frame.pack(side="left", padx=(8, 0))
        tk.Label(self.dm3_slice_frame, text="Slice", bg="#202020", fg="#ffaa44").pack(side="left")
        tk.Button(
            self.dm3_slice_frame, text="◀", command=self._on_dm3_prev,
            bg="#3a3a3a", fg="#ffaa44", activebackground="#555", relief="flat",
            padx=3, pady=0, font=("Arial", 10, "bold")
        ).pack(side="left", padx=(4, 1))
        self.var_dm3_slice = tk.IntVar(value=0)
        self.var_dm3_slice_str = tk.StringVar(value="0")
        self.cmb_dm3_slice = ttk.Combobox(
            self.dm3_slice_frame, textvariable=self.var_dm3_slice_str,
            state="readonly", width=5, values=["0"]
        )
        self.cmb_dm3_slice.pack(side="left", padx=(1, 1))
        self.cmb_dm3_slice.bind("<<ComboboxSelected>>", self._on_dm3_slice_combo)
        tk.Button(
            self.dm3_slice_frame, text="▶", command=self._on_dm3_next,
            bg="#3a3a3a", fg="#ffaa44", activebackground="#555", relief="flat",
            padx=3, pady=0, font=("Arial", 10, "bold")
        ).pack(side="left", padx=(1, 4))
        self.dm3_slice_frame.pack_forget()   # hidden until a DM3 stack is loaded

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

        tk.Label(self.tk_side, text="Initial Guesses", bg="#1f1f1f", fg="#ffffff", font=("Arial", 12, "bold")).pack(anchor="w", pady=(10, 2))

        # Radio: Reimer Auto-fill vs Manual for Defocus / Cs guesses
        self.var_guess_mode = tk.StringVar(value="Reimer Auto")
        guess_radio_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        guess_radio_row.pack(fill="x", pady=(0, 4))
        for lbl in ("Reimer Auto", "Manual"):
            tk.Radiobutton(
                guess_radio_row, text=lbl,
                variable=self.var_guess_mode, value=lbl,
                command=self._on_guess_mode_change,
                bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
                activebackground="#1f1f1f", activeforeground="#ffffff",
                indicatoron=1
            ).pack(side="left", padx=(0, 8))

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
            reimer_n_row, from_=-15, to=15,
            textvariable=self.var_reimer_n_start,
            command=self._on_reimer_n_change,
            width=5, bg="#2b2b2b", fg="white",
            buttonbackground="#3a3a3a", disabledbackground="#1a1a1a",
            disabledforeground="#555555", state="disabled"
        )
        self.spin_reimer_n.pack(side="left", padx=4)
        self.spin_reimer_n.bind("<Return>",   lambda e: self._on_reimer_n_change())
        self.spin_reimer_n.bind("<FocusOut>", lambda e: self._on_reimer_n_change())

        # "Clear" button resets Reimer point exclusions
        reimer_clear_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        reimer_clear_row.pack(fill="x", pady=2)
        tk.Button(
            reimer_clear_row, text="Clear Reimer Exclusions",
            bg="#3a3a3a", fg="white", activebackground="#555",
            relief="flat", padx=6, pady=2,
            command=self._clear_reimer_exclusions
        ).pack(side="left", padx=4)

        # Line source: Reimer estimate vs CTF fit parameters
        self.var_reimer_curve_src = tk.StringVar(value="fit")
        reimer_curve_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        reimer_curve_row.pack(fill="x", pady=2)
        tk.Label(reimer_curve_row, text="Line src:", width=10, anchor="w",
                 bg="#1f1f1f", fg="#d1d5db").pack(side="left")
        for lbl, val in (("Reimer", "reimer"), ("Fit", "fit")):
            tk.Radiobutton(
                reimer_curve_row, text=lbl,
                variable=self.var_reimer_curve_src, value=val,
                command=self._estimate_from_rings,
                bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
                activebackground="#1f1f1f", activeforeground="#ffffff"
            ).pack(side="left", padx=(0, 4))

        # Plot view: Reimer (n/q²) vs Phase (χ vs q, quadratic fit)
        self.var_reimer_view = tk.StringVar(value="Reimer")
        reimer_view_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        reimer_view_row.pack(fill="x", pady=(4, 2))
        tk.Label(reimer_view_row, text="Plot view:", width=10, anchor="w",
                 bg="#1f1f1f", fg="#d1d5db").pack(side="left")
        for lbl, val in (("n/q²", "Reimer"), ("Phase fit", "Phase")):
            tk.Radiobutton(
                reimer_view_row, text=lbl,
                variable=self.var_reimer_view, value=val,
                command=self._estimate_from_rings,
                bg="#1f1f1f", fg="#e5e7eb", selectcolor="#2b2b2b",
                activebackground="#1f1f1f", activeforeground="#ffffff"
            ).pack(side="left", padx=(0, 4))

        # Phase fit initial guesses
        tk.Label(self.tk_side, text="Phase fit p₀", bg="#1f1f1f", fg="#aaaaaa",
                 font=("Arial", 9, "italic")).pack(anchor="w", pady=(6, 1))
        phase_p0_dz_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        phase_p0_dz_row.pack(fill="x", pady=1)
        tk.Label(phase_p0_dz_row, text="Δz guess (nm):", width=14, anchor="w",
                 bg="#1f1f1f", fg="#d1d5db").pack(side="left")
        self.var_phase_dz_guess = tk.StringVar(value="500")
        tk.Entry(phase_p0_dz_row, textvariable=self.var_phase_dz_guess,
                 width=8, bg="#2b2b2b", fg="white",
                 insertbackground="white").pack(side="left", padx=4)

        phase_p0_cs_row = tk.Frame(self.tk_side, bg="#1f1f1f")
        phase_p0_cs_row.pack(fill="x", pady=1)
        tk.Label(phase_p0_cs_row, text="Cs guess (mm):", width=14, anchor="w",
                 bg="#1f1f1f", fg="#d1d5db").pack(side="left")
        self.var_phase_cs_guess = tk.StringVar(value="1.2")
        tk.Entry(phase_p0_cs_row, textvariable=self.var_phase_cs_guess,
                 width=8, bg="#2b2b2b", fg="white",
                 insertbackground="white").pack(side="left", padx=4)

        self.fig.canvas.draw_idle()

    def _on_reimer_mode_change(self):
        """Enable/disable the n-start spinbox based on Auto/Manual selection; re-run Reimer."""
        is_manual = (self.var_reimer_mode.get() == "Manual")
        self.spin_reimer_n.config(state="normal" if is_manual else "disabled")
        self.reimer_auto = not is_manual
        if is_manual:
            try:
                self.reimer_n_start_manual = int(self.var_reimer_n_start.get())
            except (ValueError, tk.TclError):
                self.reimer_n_start_manual = 1
        self._estimate_from_rings()

    def _on_reimer_n_change(self):
        """Called by Spinbox arrow clicks, Return, or FocusOut in Manual mode.
        Syncs var_reimer_n_start → reimer_n_start_manual and re-runs Reimer."""
        if self.reimer_auto:
            return   # Auto mode — spinbox is disabled, nothing to do
        try:
            self.reimer_n_start_manual = int(self.var_reimer_n_start.get())
        except (ValueError, tk.TclError):
            return
        self._estimate_from_rings()
        self.fig.canvas.draw_idle()

    def _on_guess_mode_change(self):
        """Toggle whether Reimer fit results auto-fill Defocus/Cs guesses."""
        self.reimer_push_guesses = (self.var_guess_mode.get() == "Reimer Auto")

    def _on_px_commit(self, event=None):
        """Re-run the full fit when the pixel size entry is committed (Return or FocusOut)."""
        try:
            val = float(self.var_px.get())
            if val <= 0:
                return
            self.base_cal_nm_per_px = val
        except ValueError:
            return
        self._update_slider_text()
        if hasattr(self, 'var_fit_range_mode') and self.var_fit_range_mode.get() == "Auto":
            self._auto_fit_range()
        self._run_fit()

    def _on_fit_range_mode_change(self):
        """Enable/disable fit-range sliders; trigger auto-range when Auto is chosen."""
        auto = (self.var_fit_range_mode.get() == "Auto")
        state = "disabled" if auto else "normal"
        if hasattr(self, 'scale_fit_start'):
            self.scale_fit_start.config(state=state)
        if hasattr(self, 'scale_fit_end'):
            self.scale_fit_end.config(state=state)
        if auto:
            self._auto_fit_range()
            self._blit_fft()

    def _auto_fit_range(self):
        """Auto-compute optimal fit-start / fit-end radii from the radial FFT profile.

        Strategy:
          1. Skip the DC region (~4 % of max radius, min 5 px).
          2. Smooth the radial-average profile with a running-mean window.
          3. Find local maxima above the median of the search region.
          4. fit_start = 80 % of first-peak radius (at least dc_skip).
          5. fit_end   = 115 % of last-peak radius (capped at max_valid_radius).
          6. Fallback when no peaks detected: 12 % – 75 % of max radius.
        """
        if not hasattr(self, 'rad_avg_fft') or self.rad_avg_fft is None:
            return
        profile = np.asarray(self.rad_avg_fft, dtype=float)
        n = len(profile)
        if n < 10:
            return

        dc_skip = max(5, int(0.04 * n))
        win = max(3, n // 40)
        kernel = np.ones(win) / win
        smoothed = np.convolve(profile[dc_skip:], kernel, mode='same')

        threshold = np.median(smoothed)
        peaks_rel = [i for i in range(1, len(smoothed) - 1)
                     if smoothed[i] >= smoothed[i - 1] and smoothed[i] >= smoothed[i + 1]
                     and smoothed[i] > threshold]
        peaks = [p + dc_skip for p in peaks_rel]

        if len(peaks) >= 2:
            fit_start = max(dc_skip, int(peaks[0]  * 0.80))
            fit_end   = min(self.max_valid_radius, int(peaks[-1] * 1.15))
        elif len(peaks) == 1:
            fit_start = max(dc_skip, int(peaks[0]  * 0.80))
            fit_end   = min(self.max_valid_radius, int(peaks[0]  * 1.50))
        else:
            fit_start = max(1, int(0.12 * self.max_valid_radius))
            fit_end   = max(fit_start + 20, int(0.75 * self.max_valid_radius))

        fit_end   = min(fit_end,  self.max_valid_radius)
        fit_start = max(1, min(fit_start, fit_end - 5))

        self.fit_start_radius = float(fit_start)
        self.fit_end_radius   = float(fit_end)
        self.var_fit_start.set(fit_start)
        self.var_fit_end.set(fit_end)
        if hasattr(self, 'fft_start_circle'):
            self.fft_start_circle.set_radius(self.fit_start_radius)
        if hasattr(self, 'fft_end_circle'):
            self.fft_end_circle.set_radius(self.fit_end_radius)
        self._update_slider_text()

    def _draw_reimer_hyperbolas(self, n_arr, q2_all):
        """Draw light reference hyperbolas y = n/x for each ring's n value."""
        for line in self._reimer_hyperbola_lines:
            line.remove()
        self._reimer_hyperbola_lines.clear()
        if len(q2_all) < 1:
            return
        x_min = max(1e-4, q2_all.min() * 0.4)
        x_max = q2_all.max() * 1.5
        x_plot = np.linspace(x_min, x_max, 400)
        for n_i in n_arr:
            if abs(n_i) < 1e-9:
                continue
            y_hyp = float(n_i) / x_plot
            line, = self.ax_reimer_plot.plot(
                x_plot, y_hyp, '-', color='white', alpha=0.15,
                linewidth=0.5, zorder=1
            )
            self._reimer_hyperbola_lines.append(line)

    def _clear_reimer_exclusions(self):
        """Reset all Reimer point exclusions and refit."""
        self._reimer_excluded.clear()
        self._estimate_from_rings()
        self.fig.canvas.draw_idle()

    def _on_reimer_click(self, event):
        """Toggle a Reimer point in/out of the fit when clicked."""
        if event.inaxes is not self.ax_reimer_plot:
            return
        if len(self._reimer_all_q2) == 0:
            return
        xc, yc = event.xdata, event.ydata
        if xc is None or yc is None:
            return

        q2 = self._reimer_all_q2
        y  = self._reimer_all_y

        # Normalise distances so x and y axes are treated equally
        xr = float(q2.max() - q2.min()) or 1.0
        yr = float(y.max()  - y.min())  or 1.0
        dists = np.hypot((q2 - xc) / xr, (y - yc) / yr)
        nearest = int(np.argmin(dists))

        # Ignore clicks that are far from any point (>15 % of axis range)
        if dists[nearest] > 0.15:
            return

        if nearest in self._reimer_excluded:
            self._reimer_excluded.discard(nearest)
        else:
            self._reimer_excluded.add(nearest)

        self._estimate_from_rings()
        self.fig.canvas.draw_idle()

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

        # Update CTF plot vertical lines for fit-start (blue) and fit-end (green)
        if hasattr(self, 'fit_start_vline'):
            self.fit_start_vline.set_xdata([k_start, k_start])
        if hasattr(self, 'fit_end_vline'):
            self.fit_end_vline.set_xdata([k_end, k_end])

    def _on_load_image(self, event=None):
        try:
            from tkinter import filedialog
            tk_root = self.fig.canvas.get_tk_widget().winfo_toplevel()
            path = filedialog.askopenfilename(
                parent=tk_root,
                title="Select an image",
                filetypes=[
                    ("Images & DM3", "*.png *.jpg *.tif *.tiff *.dm3"),
                    ("DM3 files",    "*.dm3"),
                    ("All",          "*.*"),
                ]
            )
            if not path or not os.path.exists(path):
                return

            ext = os.path.splitext(path)[1].lower()
            if ext == '.dm3':
                stack, n_slices, cal_nm, kv_dm3, cs_dm3 = _dm3_read(path)
                # Push calibration into the pixel-size entry if found
                if cal_nm > 0.0:
                    self.base_cal_nm_per_px = cal_nm
                    self.var_px.set(f"{cal_nm:.6g}")
                # Push kV into the voltage combobox if found
                if kv_dm3 > 0.0:
                    _kv_map = {300: "300 kV", 200: "200 kV", 120: "120 kV",
                               80: "80 kV", 60: "60 kV", 30: "30 kV"}
                    kv_round = min(_kv_map.keys(), key=lambda k: abs(k - kv_dm3))
                    kv_str = _kv_map[kv_round]
                    self.var_kv.set(kv_str)
                    self.wavelength_nm = ELECTRON_WAVELENGTHS[kv_str]
                    self.current_kv_ev = kv_round * 1000.0
                # Push Cs into the Cs entry if found
                if cs_dm3 >= 0.0:
                    self.current_cs = cs_dm3
                    if hasattr(self, 'var_cs'):
                        self.var_cs.set(f"{cs_dm3:.4f}")
                if stack.ndim == 3 and n_slices > 1:
                    # Store the full 3-D stack and configure the slice slider
                    self._dm3_stack = stack
                    self._dm3_n_slices = n_slices
                    self._dm3_current_slice = 0
                    slice_vals = [str(i) for i in range(n_slices)]
                    self.cmb_dm3_slice.config(values=slice_vals)
                    self.var_dm3_slice_str.set("0")
                    self.var_dm3_slice.set(0)
                    self.dm3_slice_frame.pack(side="left", padx=(8, 0))
                    # Extract first slice as the working image
                    sl = stack[0].astype(np.float64)
                else:
                    # Single 2-D DM3
                    self._dm3_stack = None
                    self.dm3_slice_frame.pack_forget()
                    sl = (stack if stack.ndim == 2 else stack[0]).astype(np.float64)
                dmin, dmax = sl.min(), sl.max()
                self.base_image = (sl - dmin) / (dmax - dmin) if dmax > dmin else np.zeros_like(sl)
            else:
                # Non-DM3 file — hide any previously shown slice slider
                self._dm3_stack = None
                self.dm3_slice_frame.pack_forget()
                self.base_image = load_image_as_float(path)

            self.bin_factor = 2
            self.var_bin.set("2")
            self._bg_fft = None
            _radial_cache.clear()
            self._apply_binning()
            self._compute_fft_and_profile()
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
            self.fig.canvas.manager.set_window_title(
                f"CTF Analyzer  v2.97  \u2014  {os.path.basename(path)}")
        except Exception as e:
            import traceback
            print(f"Error loading image: {e}")
            traceback.print_exc()

    def _on_export_rings(self, event=None):
        """Export detected ring positions (q, d, type, χ) to an xlsx file."""
        try:
            import openpyxl
        except ImportError:
            from tkinter import messagebox as _mb
            tk_root = self.fig.canvas.get_tk_widget().winfo_toplevel()
            _mb.showerror("Missing library",
                          "openpyxl is required for xlsx export.\n"
                          "Install with:  pip install openpyxl",
                          parent=tk_root)
            return

        q_all  = self._last_ring_q_all
        n_arr  = self._last_ring_n_arr
        types  = self._last_ring_types

        if len(q_all) == 0:
            from tkinter import messagebox as _mb
            tk_root = self.fig.canvas.get_tk_widget().winfo_toplevel()
            _mb.showinfo("No data", "No ring data available. Load an image and run the fit first.",
                         parent=tk_root)
            return

        from tkinter import filedialog as _fd
        tk_root = self.fig.canvas.get_tk_widget().winfo_toplevel()
        path = _fd.asksaveasfilename(
            parent=tk_root,
            title="Save ring data as xlsx",
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", "*.xlsx"), ("All files", "*.*")],
        )
        if not path:
            return

        wl = self._last_ring_wl
        chi_n = n_arr * (np.pi / 2.0)

        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Ring data"
        ws.append(["n", "Type", "q (1/nm)", "d (nm)", "chi (rad)"])
        for n, t, q, chi in zip(n_arr, types, q_all, chi_n):
            d = 1.0 / q if q > 0 else 0.0
            ws.append([int(n), t, float(q), float(d), float(chi)])

        # Auto-width columns
        for col in ws.columns:
            max_len = max(len(str(cell.value)) for cell in col if cell.value is not None)
            ws.column_dimensions[col[0].column_letter].width = max_len + 4

        wb.save(path)

        from tkinter import messagebox as _mb
        _mb.showinfo("Exported", f"Ring data saved to:\n{path}", parent=tk_root)

    def _on_key_press(self, event):
        """Left/right arrow keys step the DM3 stack slice slider by one frame."""
        if self._dm3_stack is None:
            return
        if event.key not in ('left', 'right'):
            return
        idx = int(self.var_dm3_slice.get())
        idx += 1 if event.key == 'right' else -1
        idx = max(0, min(idx, self._dm3_n_slices - 1))
        self.var_dm3_slice.set(idx)
        self.var_dm3_slice_str.set(str(idx))
        self._on_dm3_slice_change()

    def _on_dm3_prev(self):
        """◀ button — step one slice back."""
        if self._dm3_stack is None:
            return
        idx = max(0, int(self.var_dm3_slice.get()) - 1)
        self.var_dm3_slice.set(idx)
        self.var_dm3_slice_str.set(str(idx))
        self._on_dm3_slice_change()

    def _on_dm3_next(self):
        """▶ button — step one slice forward."""
        if self._dm3_stack is None:
            return
        idx = min(self._dm3_n_slices - 1, int(self.var_dm3_slice.get()) + 1)
        self.var_dm3_slice.set(idx)
        self.var_dm3_slice_str.set(str(idx))
        self._on_dm3_slice_change()

    def _on_dm3_slice_combo(self, event=None):
        """Called when the slice combobox selection changes."""
        if self._dm3_stack is None:
            return
        try:
            idx = int(self.var_dm3_slice_str.get())
        except ValueError:
            return
        idx = max(0, min(idx, self._dm3_n_slices - 1))
        self.var_dm3_slice.set(idx)
        self._on_dm3_slice_change()

    def _on_dm3_slice_change(self, val=None):
        """Called when the DM3 slice selection changes; re-extracts the selected frame."""
        if self._dm3_stack is None:
            return
        idx = int(self.var_dm3_slice.get())
        idx = max(0, min(idx, self._dm3_n_slices - 1))
        self._dm3_current_slice = idx
        self.var_dm3_slice_str.set(str(idx))

        sl = self._dm3_stack[idx].astype(np.float64)
        dmin, dmax = sl.min(), sl.max()
        self.base_image = (sl - dmin) / (dmax - dmin) if dmax > dmin else np.zeros_like(sl)

        self._bg_fft = None
        _radial_cache.clear()
        self._apply_binning()
        self._compute_fft_and_profile()
        self._run_fit()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def _apply_binning(self):
        self.image = perform_binning(self.base_image, self.bin_factor)
        self.h, self.w = self.image.shape
        self.img_plot.set_data(self.image)
        self.img_plot.set_clim(self.image.min(), self.image.max())
        self.img_plot.set_extent([-0.5, self.w - 0.5, -0.5, self.h - 0.5])
        self.ax_img.set_xlim(0, self.w)
        self.ax_img.set_ylim(0, self.h)
        
        old_max_radius = self.max_valid_radius
        # Largest circle that fits inside the square crop (half-side)
        sq_side = min(self.h, self.w)
        self.max_valid_radius = sq_side // 2
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

        # ── Crop to largest center square ────────────────────────────────────
        sq = min(h_img, w_img)
        y0 = (h_img - sq) // 2
        x0 = (w_img - sq) // 2
        crop = self.image[y0:y0 + sq, x0:x0 + sq]

        # Update the yellow crop-square overlay on the real-space image panel.
        # imshow with origin='lower' maps array row 0 to y=0 at the bottom,
        # so the rectangle corner in axes coords is (x0, y0) as expected.
        if hasattr(self, '_fft_crop_rect'):
            self._fft_crop_rect.set_bounds(x0 - 0.5, y0 - 0.5, sq, sq)
            self._fft_crop_rect.set_visible(True)

        # Use scipy.fft (multithreaded) with fast-length padding when available
        if _HAS_SCIPY_FFT:
            fast_sq = next_fast_len(sq)
            fft_img = fft2(crop, s=(fast_sq, fast_sq), workers=-1)
            fft_shifted = fftshift(fft_img)
        else:
            fft_img = np.fft.fft2(crop)
            fft_shifted = np.fft.fftshift(fft_img)

        amplitude_spectrum = np.abs(fft_shifted)
        self.fft_img_log = np.log(amplitude_spectrum + 1e-10).astype(np.float32)

        h, w = amplitude_spectrum.shape
        fft_cx, fft_cy = w // 2, h // 2
        self.N_dimension = sq   # square side = the spatial scale for k-axis

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
        self.info_limit_circle.set_center((fft_cx, fft_cy))

        self.fft_start_circle.set_visible(True)
        self.fft_end_circle.set_visible(True)

        self.text_fft_start.set_text("")
        self.text_fft_end.set_text("")
        self._update_slider_text()

        max_rad = int(np.hypot(w / 2.0, h / 2.0))
        self.radii_px, self.rad_avg_fft = compute_radial_profile_fast(amplitude_spectrum, fft_cx, fft_cy, max_rad)

        # Auto-update fit range from profile peaks (if Auto mode is active)
        if hasattr(self, 'var_fit_range_mode') and self.var_fit_range_mode.get() == "Auto":
            self._auto_fit_range()

        # Always refresh the gray data line immediately so it reflects the new
        # slice even if _run_fit exits early (invalid range, failed curve_fit…)
        eff_px_now = self.base_cal_nm_per_px * self.bin_factor
        if eff_px_now > 0:
            k_now = self.radii_px / (self.N_dimension * eff_px_now)
            self.k_axis = k_now
            if hasattr(self, 'line_rad_fft'):
                self.line_rad_fft.set_data(k_now, self.rad_avg_fft)

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

        eff_px = self.base_cal_nm_per_px * self.bin_factor
        N_dim  = self.N_dimension

        # ── Top-hat ring detection ────────────────────────────────────────────
        # A light moving average (window ≈ ROI/40, ≥3) is applied first to
        # suppress sample noise without blurring the CTF ring structure; the
        # top-hat filter then removes the slowly-varying background.
        # Structuring element ~2 ring-periods wide so the morphological opening
        # tracks the slowly-varying background and the top-hat signals contain
        # only the CTF oscillations, independent of background slope.
        th_max_radii, th_min_radii = [], []
        if _HAS_MORPHOLOGY and roi_len > 10:
            se_size  = max(7, roi_len // 7)
            th_order = max(3, roi_len // 35)

            # Light pre-smooth: small window to kill noise, preserve ring shape
            ma_w_light = max(3, roi_len // 40)
            roi_smooth = _moving_average(roi_raw, ma_w_light)

            # Plot the smoothed signal used as top-hat input
            if hasattr(self, 'line_th_ma') and self.k_axis is not None and len(self.k_axis) > r_hi:
                th_ma_full = np.full(len(self.k_axis), np.nan)
                th_ma_full[r_lo : r_hi + 1] = roi_smooth
                self.line_th_ma.set_data(self.k_axis, th_ma_full)

            opening = grey_opening(roi_smooth, size=se_size)
            wth = roi_smooth - opening                     # white top-hat: peaks
            bth = grey_opening(-roi_smooth, size=se_size)  # black top-hat: troughs
            bth = -bth - roi_smooth                        # = closing − signal > 0 at troughs
            # Clip negatives that arise from floating-point rounding
            wth = np.clip(wth, 0, None)
            bth = np.clip(bth, 0, None)

            # Plot background estimate on the CTF axes
            if hasattr(self, 'line_th_open') and self.k_axis is not None and len(self.k_axis) > r_hi:
                th_open_full = np.full(len(self.k_axis), np.nan)
                th_open_full[r_lo : r_hi + 1] = opening
                self.line_th_open.set_data(self.k_axis, th_open_full)

            # Threshold: ignore peaks below 5 % of the respective top-hat max
            wth_thresh = float(np.max(wth)) * 0.05 if wth.max() > 0 else 0.0
            bth_thresh = float(np.max(bth)) * 0.05 if bth.max() > 0 else 0.0

            if roi_len > th_order * 2 + 2:
                th_max_rel_raw, _ = _local_extrema_1d(wth, th_order)
                th_min_rel_raw, _ = _local_extrema_1d(bth, th_order)
                th_max_rel = [i for i in th_max_rel_raw if wth[i] >= wth_thresh]
                th_min_rel = [i for i in th_min_rel_raw if bth[i] >= bth_thresh]
                th_max_radii = (r_lo + np.array(th_max_rel, dtype=int)).tolist()
                th_min_radii = (r_lo + np.array(th_min_rel, dtype=int)).tolist()
                # Reject detections within 5 % of either ROI boundary; the
                # rolling average and top-hat filter produce unreliable extrema
                # near the edges of the ROI window.
                _guard = max(3, int((r_hi - r_lo) * 0.05))
                th_max_radii = [r for r in th_max_radii
                                if r_lo + _guard <= r <= r_hi - _guard]
                th_min_radii = [r for r in th_min_radii
                                if r_lo + _guard <= r <= r_hi - _guard]

                # ── Enforce max–min–max alternation ──────────────────────────
                # Merge all detections in q-order, then walk the sequence and
                # keep only those that form a strictly alternating pattern.
                # The first detection is always kept; each subsequent one must
                # be the opposite type to the previously accepted one.
                # Among runs of the same type, keep the most prominent (highest
                # top-hat value for maxima, highest bth value for minima).
                _tagged = (
                    [(r, 'max') for r in th_max_radii] +
                    [(r, 'min') for r in th_min_radii]
                )
                _tagged.sort(key=lambda x: x[0])

                _alt_max, _alt_min = [], []
                _last_type = None
                _pending = None   # (radius, type, prominence)

                def _prom(r, t):
                    rel = r - r_lo
                    if 0 <= rel < len(wth):
                        return float(wth[rel]) if t == 'max' else float(bth[rel])
                    return 0.0

                for _r, _t in _tagged:
                    if _last_type is None:
                        _pending = (_r, _t, _prom(_r, _t))
                        _last_type = _t
                    elif _t == _last_type:
                        # Same type as pending: keep the more prominent one
                        _p = _prom(_r, _t)
                        if _p > _pending[2]:
                            _pending = (_r, _t, _p)
                    else:
                        # Different type: commit pending, start new pending
                        if _pending[1] == 'max':
                            _alt_max.append(_pending[0])
                        else:
                            _alt_min.append(_pending[0])
                        _pending = (_r, _t, _prom(_r, _t))
                        _last_type = _t
                if _pending is not None:
                    if _pending[1] == 'max':
                        _alt_max.append(_pending[0])
                    else:
                        _alt_min.append(_pending[0])

                th_max_radii = sorted(_alt_max)
                th_min_radii = sorted(_alt_min)

                # ── Spacing-consistency filter ────────────────────────────────
                # Build the merged ring sequence (all maxima + minima in q-order).
                # CTF ring spacings shrink monotonically with q, so we fit a
                # linear trend to spacing vs mid-point index and reject any ring
                # whose gap to its predecessor is less than 0.35× the predicted
                # spacing at that position (too close → spurious split) OR whose
                # gap is more than 2.8× the predicted spacing (large gap → likely
                # missed ring, but we don't fabricate rings, so just leave it).
                # Rejection always removes the *less prominent* of the two
                # neighbours producing the bad gap.
                if len(th_max_radii) + len(th_min_radii) >= 4:
                    _all_rings = sorted(
                        [(r, 'max') for r in th_max_radii] +
                        [(r, 'min') for r in th_min_radii]
                    )
                    _radii_seq = [r for r, _ in _all_rings]
                    _types_seq = [t for _, t in _all_rings]

                    # Iteratively remove spacing outliers until stable
                    for _iter in range(10):
                        if len(_radii_seq) < 3:
                            break
                        _gaps = [_radii_seq[i+1] - _radii_seq[i]
                                 for i in range(len(_radii_seq) - 1)]
                        # Fit linear trend: spacing = a + b*i
                        _xs = np.arange(len(_gaps), dtype=float)
                        if len(_gaps) >= 2:
                            _pb = np.polyfit(_xs, _gaps, 1)
                            _pred = np.polyval(_pb, _xs)
                        else:
                            _pred = np.full(len(_gaps), float(np.median(_gaps)))
                        # Clamp predictions to positive
                        _pred = np.clip(_pred, max(1, np.min(_gaps) * 0.3), None)
                        _removed = False
                        for _gi in range(len(_gaps)):
                            ratio = _gaps[_gi] / _pred[_gi]
                            if ratio < 0.35:
                                # Too close: remove less prominent neighbour
                                _ia, _ib = _gi, _gi + 1
                                _pa = _prom(_radii_seq[_ia], _types_seq[_ia])
                                _pb2 = _prom(_radii_seq[_ib], _types_seq[_ib])
                                _drop = _ib if _pa >= _pb2 else _ia
                                _radii_seq.pop(_drop)
                                _types_seq.pop(_drop)
                                _removed = True
                                break  # restart with updated sequence
                        if not _removed:
                            break

                    th_max_radii = sorted(r for r, t in zip(_radii_seq, _types_seq) if t == 'max')
                    th_min_radii = sorted(r for r, t in zip(_radii_seq, _types_seq) if t == 'min')
        else:
            if hasattr(self, 'line_th_ma'):
                self.line_th_ma.set_data([], [])
            if hasattr(self, 'line_th_open'):
                self.line_th_open.set_data([], [])

        self._th_max_ks = [r / (N_dim * eff_px) for r in th_max_radii if r > 0]
        self._th_min_ks = [r / (N_dim * eff_px) for r in th_min_radii if r > 0]

        # Drive the shared RA ring lists from top-hat results so Reimer
        # estimation and the FFT image dots use the top-hat detections.
        self._ra_max_ks = list(self._th_max_ks)
        self._ra_min_ks = list(self._th_min_ks)

        # Place markers at fft_cx ± r (positive AND negative x side)
        max_cols, min_cols = [], []
        for r in th_max_radii:
            for col in (fft_cx + r, fft_cx - r):
                if 0 <= col < w:
                    max_cols.append(col)
        for r in th_min_radii:
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

        # ── CTF triangle markers (top-hat detections) ────────────────────────
        if hasattr(self, 'ctf_max_marks') and hasattr(self, 'ctf_min_marks'):
            y_row = 0.015
            th_mk = np.array(self._th_max_ks, dtype=float)
            th_mn = np.array(self._th_min_ks, dtype=float)
            self.ctf_max_marks.set_data(th_mk, np.full_like(th_mk, y_row))
            self.ctf_min_marks.set_data(th_mn, np.full_like(th_mn, y_row))

        # ── Vertical dashed lines from each marker up to the data curve ───────
        if hasattr(self, '_th_vlines'):
            for vl in self._th_vlines:
                try:
                    vl.remove()
                except Exception:
                    pass
            self._th_vlines.clear()
            if self.k_axis is not None and self.rad_avg_fft is not None:
                for k_val in self._th_max_ks:
                    idx = int(np.argmin(np.abs(self.k_axis - k_val)))
                    y_val = float(self.rad_avg_fft[idx])
                    ln, = self.ax_ctf.plot(
                        [k_val, k_val], [0, y_val],
                        color='red', linestyle='--', linewidth=0.9,
                        alpha=0.6, zorder=3
                    )
                    self._th_vlines.append(ln)
                for k_val in self._th_min_ks:
                    idx = int(np.argmin(np.abs(self.k_axis - k_val)))
                    y_val = float(self.rad_avg_fft[idx])
                    ln, = self.ax_ctf.plot(
                        [k_val, k_val], [0, y_val],
                        color='#3399ff', linestyle='--', linewidth=0.9,
                        alpha=0.6, zorder=3
                    )
                    self._th_vlines.append(ln)

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

        # ── RA Rings table (top-hat detections only) ──────────────────────
        def _ra_section(label, ks):
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

        ra_lines = _ra_section("Maxima (Top-hat):", list(self._th_max_ks))
        ra_lines.append("")
        ra_lines += _ra_section("Minima (Top-hat):", list(self._th_min_ks))
        if hasattr(self, 'ra_peaks_text'):
            self.ra_peaks_text.set_text("\n".join(ra_lines))

        self._estimate_from_rings()

    def _estimate_from_rings(self):
        """Estimate Cs and defocus from detected ring positions via Reimer eq. 6.50.

        All detected maxima and minima are merged into one list sorted jointly
        by q.  Consecutive integers  n_start, n_start+1, n_start+2, …  are
        assigned in ascending-q order (no separate odd/even lists).

        In the code's CTF convention  chi = -π λ Δz q² + ½ π Cs λ³ q⁴
        (Reimer/textbook convention; Δz > 0 underfocus, < 0 overfocus),
        setting chi = n π/2 and dividing by q² gives:

            n / q²  =  Cs λ³ · q²  -  2 λ Δz

        so  slope m = +Cs λ³  (positive)  and  intercept b = -2 λ Δz.

        Auto mode:  tries n_start = 1 … 19, keeps only physically valid results
        (Cs > 0), picks the n_start giving the highest R².
        Manual mode: uses the user-supplied n_start value directly.
        """
        if not hasattr(self, 'reimer_text'):
            return

        def _clear_plot():
            if hasattr(self, 'reimer_scat_max'):
                self.reimer_scat_max.set_data([], [])
                self.reimer_scat_min.set_data([], [])
                self.reimer_fit_line.set_data([], [])
            if hasattr(self, '_reimer_errorbars'):
                self._reimer_errorbars.set_data([], [])
                self._reimer_errorbars.set_visible(False)
                for ann in self._reimer_annots:
                    ann.remove()
                self._reimer_annots.clear()
            if hasattr(self, '_reimer_incl_circles'):
                self._reimer_incl_circles.set_data([], [])

        # Use top-hat extrema (background-independent ring detection)
        ra_max = sorted(k for k in self._th_max_ks if k > 1e-9)
        ra_min = sorted(k for k in self._th_min_ks if k > 1e-9)
        src_label = "TH"
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
        n_total = len(q_all)

        # Reset exclusion set when the ring population changes
        if n_total != self._reimer_pt_count:
            self._reimer_excluded.clear()
            self._reimer_pt_count = n_total

        # Inclusion mask — True for points used in the fit
        incl = np.array([i not in self._reimer_excluded
                         for i in range(n_total)], dtype=bool)

        # Cache geometry for phase-χ view (n_arr cached per-branch below)
        self._last_ring_q_all = q_all.copy()
        self._last_ring_types = types[:]
        self._last_ring_incl  = incl.copy()
        self._last_ring_wl    = wl

        # ── Search over n_start values ────────────────────────────────────────
        n_start_range = (range(-15, 16) if self.reimer_auto
                         else [int(self.reimer_n_start_manual)])

        best_r2 = -np.inf
        best    = None

        for n_start in n_start_range:
            n_arr  = np.arange(n_start, n_start + n_total, dtype=float)
            y      = n_arr / q2_all

            # Fit only on included points (need ≥2)
            q2_fit = q2_all[incl];  y_fit = y[incl]
            if len(q2_fit) < 2:
                continue

            try:
                coeffs, cov_try = np.polyfit(q2_fit, y_fit, 1, cov=True)
            except (np.linalg.LinAlgError, ValueError):
                coeffs = np.polyfit(q2_fit, y_fit, 1)
                cov_try = None
            m_try, b_try = float(coeffs[0]), float(coeffs[1])

            Cs_try = m_try / (wl ** 3)       # nm  (positive slope → positive Cs)
            dz_try = -b_try / (2.0 * wl)     # nm  (b = −2λΔz → Δz = −b/2λ)

            # Reject physically impossible results
            if Cs_try <= 0:
                continue

            y_pred_fit = m_try * q2_fit + b_try
            ss_res = float(np.sum((y_fit - y_pred_fit) ** 2))
            ss_tot = float(np.sum((y_fit - float(np.mean(y_fit))) ** 2))
            r2_try = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0
            # Reduced chi-square: set σ=1 per point → χ²_red = RSS/(N-2)
            chi2r_try = ss_res / max(1, len(q2_fit) - 2)

            if r2_try > best_r2:
                best_r2 = r2_try
                best = dict(
                    n_start=n_start, m=m_try, b=b_try,
                    Cs_nm=Cs_try, dz_nm=dz_try, r2=r2_try,
                    chi2r=chi2r_try, cov=cov_try,
                    n_arr=n_arr, y=y, incl=incl
                )

        if best is None:
            # ── No valid fit, but still show the ring points ──────────────────
            # Use n_start=1 as a display fallback so the scatter plot and table
            # remain informative even though no physically valid line was found.
            n_start_fb = 1
            n_arr_fb   = np.arange(n_start_fb, n_start_fb + n_total, dtype=float)
            y_fb       = n_arr_fb / q2_all
            self._last_ring_n_arr = n_arr_fb.copy()   # cache for phase view

            # Store for click-handler
            self._reimer_all_q2 = q2_all.copy()
            self._reimer_all_y  = y_fb.copy()

            # Scatter plot (no fit line)
            if hasattr(self, 'reimer_scat_max'):
                max_x = [qi**2 for qi, ti in zip(q_all, types) if ti == 'max']
                max_y = [ni/qi**2 for ni, qi, ti in zip(n_arr_fb, q_all, types) if ti == 'max']
                min_x = [qi**2 for qi, ti in zip(q_all, types) if ti == 'min']
                min_y = [ni/qi**2 for ni, qi, ti in zip(n_arr_fb, q_all, types) if ti == 'min']
                self.reimer_scat_max.set_data(max_x, max_y)
                self.reimer_scat_min.set_data(min_x, min_y)
                # Draw fit-parameter line if "Fit" source selected, else blank
                if (getattr(self, 'var_reimer_curve_src', None) and
                        self.var_reimer_curve_src.get() == 'fit' and
                        len(q2_all) >= 2):
                    x_lo = max(0.0, q2_all.min() * 0.85)
                    x_hi = q2_all.max() * 1.15
                    x_line = np.linspace(x_lo, x_hi, 150)
                    try:
                        dz_fit = float(self.var_defocus.get())
                        cs_fit_nm = float(self.var_cs.get()) * 1e6
                        m_fit = cs_fit_nm * (wl ** 3)
                        b_fit = -2.0 * wl * dz_fit
                        self.reimer_fit_line.set_data(x_line, m_fit * x_line + b_fit)
                    except (ValueError, AttributeError):
                        self.reimer_fit_line.set_data([], [])
                else:
                    self.reimer_fit_line.set_data([], [])
                # Inclusion circles for all included points
                incl_q2 = q2_all[incl]; incl_y = y_fb[incl]
                self._reimer_incl_circles.set_data(incl_q2, incl_y)
                for ann in self._reimer_annots: ann.remove()
                self._reimer_annots.clear()
                # Auto-scale axes
                all_x = (max_x or [0]) + (min_x or [0])
                all_y_pts = (max_y or [0]) + (min_y or [0])
                if all_x and all_y_pts:
                    px = np.array(all_x); py = np.array(all_y_pts)
                    self.ax_reimer_plot.set_xlim(max(0, px.min()*0.8), px.max()*1.2)
                    self.ax_reimer_plot.set_ylim(py.min()*0.8, py.max()*1.2)
                self._draw_reimer_hyperbolas(n_arr_fb, q2_all)

            # Text table — show rings but flag no valid fit
            lines = [
                "── Reimer Eq. 6.50 ──────",
                "  No physically valid fit.",
                "  (Cs≤0 for all tried",
                f"   n_start values)",
                "",
                f"  n₀=1 (display only)",
                "",
                f"  {'n':>3}  {'type':3} {'q(1/nm)':>7} {'n/q²':>7}",
                "  " + "─" * 27,
            ]
            for i, (ni, qi, ti, yi) in enumerate(zip(n_arr_fb, q_all, types, y_fb)):
                lbl = "Mx" if ti == 'max' else "Mn"
                excl_tag = " *" if i in self._reimer_excluded else ""
                lines.append(
                    f"  {int(ni):>3}  {lbl:<3} {qi:>7.4f} {yi:>7.4f}{excl_tag}"
                )
            lines += ["", "  Δz: + underfocus  − overfocus"]
            self.reimer_text.set_text("\n".join(lines))
            self.reimer_text.set_color('#ef4444')
            self._apply_reimer_view(q_all, n_arr_fb, types, incl, wl)
            return

        # ── Unpack best result ────────────────────────────────────────────────
        n_start = best['n_start']
        m, b    = best['m'], best['b']
        Cs_nm   = best['Cs_nm'];  Cs_mm = Cs_nm / 1e6
        dz_nm   = best['dz_nm']
        r2      = best['r2']
        chi2r   = best.get('chi2r', float('nan'))
        cov_mat = best.get('cov', None)
        n_arr   = best['n_arr']
        y       = best['y']
        incl    = best['incl']
        self._last_ring_n_arr = n_arr.copy()   # cache for phase view
        # Residuals evaluated for all points against the fit line
        y_pred  = m * q2_all + b

        # ── Parameter errors via the χ²_red=1 trick ───────────────────────────
        # numpy polyfit cov=True already scales by RSS/(N-2), so √diag(cov) is
        # the standard error on each coefficient directly.
        if cov_mat is not None:
            sigma_Cs_mm = float(np.sqrt(max(0.0, cov_mat[0, 0]))) / (wl ** 3) / 1e6
            sigma_dz_nm = float(np.sqrt(max(0.0, cov_mat[1, 1]))) / (2.0 * wl)
        else:
            sigma_Cs_mm = sigma_dz_nm = float('nan')
        # Uniform data-point error bar: σ_y = √(χ²_red)  (i.e. setting χ²_red=1)
        sigma_y = float(np.sqrt(chi2r)) if np.isfinite(chi2r) else float('nan')

        def _fmt_err(sigma, nsig=2):
            """Round sigma to nsig significant figures; return (value_decimals, sigma_str)."""
            import math
            if not np.isfinite(sigma) or sigma <= 0:
                return 4, "±  n/a"
            mag = math.floor(math.log10(abs(sigma)))
            dec = max(0, nsig - 1 - mag)
            sigma_r = round(sigma, dec)
            return dec, f"±{sigma_r:.{dec}f}"

        _cs_dec,  _cs_err_str  = _fmt_err(sigma_Cs_mm)
        _dz_dec,  _dz_err_str  = _fmt_err(sigma_dz_nm)

        # Store for click-handler proximity detection
        self._reimer_all_q2 = q2_all.copy()
        self._reimer_all_y  = y.copy()

        # ── Push Reimer results into the fit parameter input fields ──────────
        if self.reimer_push_guesses and hasattr(self, 'var_defocus') and hasattr(self, 'var_cs'):
            self.var_defocus.set(f"{dz_nm:.1f}")
            self.var_cs.set(f"{Cs_mm:.4f}")
            self.current_cs = Cs_mm   # keep internal state in sync

        n_used = int(incl.sum())
        mode_tag = f"({'auto' if self.reimer_auto else 'manual'},{src_label})"

        # Colour-code by R²
        if r2 > 0.99:
            color = '#00ff88'
        elif r2 > 0.95:
            color = '#fb923c'
        else:
            color = '#ef4444'

        # ── Text table ───────────────────────────────────────────────────────
        n_excl = n_total - n_used
        lines = [
            "── Reimer Eq. 6.50 ──────",
            f"  n₀={n_start} {mode_tag}",
            "",
            f"  {'Cs (mm)':<13} {format(Cs_mm, f'>8.{_cs_dec}f')} {_cs_err_str}",
            f"  {'Defocus (nm)':<13} {format(dz_nm, f'>8.{_dz_dec}f')} {_dz_err_str}",
            f"  {'R²':<13} {r2:>8.5f}",
            f"  {'χ²_red':<13} {chi2r:>8.4f}",
            f"  {'Used/Total':<13} {n_used:>4d}/{n_total:<4d}",
            f"  {'λ (nm)':<13} {wl:>8.5f}",
            "",
            f"  {'n':>3}  {'type':3} {'q(1/nm)':>7} {'resid':>7}",
            "  " + "─" * 27,
        ]
        for i, (ni, qi, ti, yi, ypi) in enumerate(
                zip(n_arr, q_all, types, y, y_pred)):
            lbl = "Mx" if ti == 'max' else "Mn"
            excl_tag = " *" if i in self._reimer_excluded else ""
            lines.append(
                f"  {int(ni):>3}  {lbl:<3} {qi:>7.4f} {yi-ypi:>+7.4f}{excl_tag}"
            )

        lines += [
            "",
            "  Δz: + underfocus  − overfocus",
            "  ── R² color key ──────────────",
            "  ≥ 0.99 → green (excellent)",
            "  ≥ 0.95 → orange (acceptable)",
            "  < 0.95 → red (poor fit)",
        ]
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

        # Draw uniform error bars (σ_y = √χ²_red) as vertical stems
        if hasattr(self, '_reimer_errorbars') and np.isfinite(sigma_y) and sigma_y > 0:
            eb_x, eb_y = [], []
            for xi, yi in zip(max_x + min_x, max_y + min_y):
                eb_x += [xi, xi, None]
                eb_y += [yi - sigma_y, yi + sigma_y, None]
            self._reimer_errorbars.set_data(eb_x, eb_y)
            self._reimer_errorbars.set_visible(True)
        elif hasattr(self, '_reimer_errorbars'):
            self._reimer_errorbars.set_data([], [])
            self._reimer_errorbars.set_visible(False)

        if len(q2_all) >= 2:
            x_lo   = max(0.0, q2_all.min() * 0.85)
            x_hi   = q2_all.max() * 1.15
            x_line = np.linspace(x_lo, x_hi, 150)
            if (getattr(self, 'var_reimer_curve_src', None) and
                    self.var_reimer_curve_src.get() == 'fit'):
                try:
                    dz_fit = float(self.var_defocus.get())
                    cs_fit_nm = float(self.var_cs.get()) * 1e6
                    m_line = cs_fit_nm * (wl ** 3)
                    b_line = -2.0 * wl * dz_fit
                except (ValueError, AttributeError):
                    m_line, b_line = m, b
            else:
                m_line, b_line = m, b
            self.reimer_fit_line.set_data(x_line, m_line * x_line + b_line)
        else:
            self.reimer_fit_line.set_data([], [])

        # ── Inclusion circles (white hollow ring around every included point) ──
        if hasattr(self, '_reimer_incl_circles'):
            circ_x = q2_all[incl]
            circ_y = y[incl]
            self._reimer_incl_circles.set_data(circ_x, circ_y)

        # Refresh n-index annotations (dim excluded labels)
        for ann in self._reimer_annots:
            ann.remove()
        self._reimer_annots.clear()
        for i, (ni, qi) in enumerate(zip(n_arr, q_all)):
            ann_color = '#888888' if i in self._reimer_excluded else 'white'
            ann = self.ax_reimer_plot.annotate(
                str(int(ni)),
                xy=(qi ** 2, ni / qi ** 2),
                xytext=(3, 2), textcoords='offset points',
                color=ann_color, fontsize=8, zorder=6,
                annotation_clip=True
            )
            self._reimer_annots.append(ann)

        # ── Explicit axis limits (relim/autoscale unreliable for Line2D) ───────
        all_x = np.concatenate([q2_all])
        all_y = y   # n/q² for all points
        if len(all_x) >= 2:
            xpad = (all_x.max() - all_x.min()) * 0.15 or all_x.max() * 0.15
            ypad = (all_y.max() - all_y.min()) * 0.15 or all_y.max() * 0.15
            self.ax_reimer_plot.set_xlim(max(0.0, all_x.min() - xpad),
                                         all_x.max() + xpad)
            self.ax_reimer_plot.set_ylim(all_y.min() - ypad,
                                         all_y.max() + ypad)
        self._draw_reimer_hyperbolas(n_arr, q2_all)
        self._apply_reimer_view(q_all, n_arr, types, incl, wl)

    # ── Phase-χ view helpers ──────────────────────────────────────────────────

    def _hide_phase_artists(self):
        """Make all phase-χ view artists invisible and remove CTF markers."""
        for a in (self.phase_scat_max, self.phase_scat_min,
                  self.phase_fit_curve, self._phase_incl_circles,
                  self._phase_result_text):
            a.set_visible(False)
        for _a in getattr(self, '_phase_ctf_markers', []):
            try: _a.remove()
            except Exception: pass
        self._phase_ctf_markers = []

    def _show_reimer_artists(self):
        """Restore visibility of all Reimer n/q² artists."""
        for a in (self.reimer_scat_max, self.reimer_scat_min,
                  self.reimer_fit_line, self._reimer_incl_circles):
            a.set_visible(True)
        for ann in self._reimer_annots:
            ann.set_visible(True)
        for ln in self._reimer_hyperbola_lines:
            ln.set_visible(True)

    def _apply_reimer_view(self, q_all, n_arr, types, incl, wl):
        """Route to the correct plot view based on var_reimer_view."""
        view = getattr(self, 'var_reimer_view', None)
        if view and view.get() == 'Phase':
            # Hide Reimer artists, update labels, draw phase fit
            for a in (self.reimer_scat_max, self.reimer_scat_min,
                      self.reimer_fit_line, self._reimer_incl_circles):
                a.set_visible(False)
            for ann in self._reimer_annots:
                ann.set_visible(False)
            for ln in self._reimer_hyperbola_lines:
                ln.set_visible(False)
            self._draw_phase_view(q_all, n_arr, types, incl, wl)
        else:
            # Reimer view: ensure Reimer artists visible, phase artists hidden
            self._hide_phase_artists()
            self._show_reimer_artists()
            self.ax_reimer_plot.set_title("Reimer  n/q²  vs  q²",
                                          color="#44dd99", pad=4, fontsize=12)
            self.ax_reimer_plot.set_xlabel("q²  (nm⁻²)", color="#aaa",
                                           fontsize=10, labelpad=2)
            self.ax_reimer_plot.set_ylabel("n / q²  (nm²)", color="#aaa",
                                           fontsize=10, labelpad=2)

    def _draw_phase_view(self, q_all, n_arr, types, incl, wl):
        """Draw χ = nπ/2 vs q and fit to χ = A·q² + B·q⁴.

        CTF convention (Reimer/Williams & Carter):
            CTF = sin(χ),  χ = −π λ Δz q² + ½ π Cs λ³ q⁴
            Δz > 0 : underfocus  |  Δz < 0 : overfocus

        Extracts defocus and Cs from the quadratic coefficients:
            A = −π λ Δz  →  Δz = −A / (π λ)
            B = ½ π Cs λ³  →  Cs = 2B / (π λ³)
        """
        ax = self.ax_reimer_plot
        ax.set_title("Phase  χ  vs  q", color="#44dd99", pad=4, fontsize=12)
        ax.set_xlabel("q  (nm⁻¹)", color="#aaa", fontsize=10, labelpad=2)
        ax.set_ylabel("χ = nπ/2  (rad)", color="#aaa", fontsize=10, labelpad=2)

        # Measured phase at each ring position
        chi_n = n_arr * (np.pi / 2.0)

        # Scatter: max (odd n) and min (even n)
        max_q   = [q for q, t in zip(q_all, types) if t == 'max']
        max_chi = [c for c, t in zip(chi_n, types) if t == 'max']
        min_q   = [q for q, t in zip(q_all, types) if t == 'min']
        min_chi = [c for c, t in zip(chi_n, types) if t == 'min']
        self.phase_scat_max.set_data(max_q, max_chi)
        self.phase_scat_min.set_data(min_q, min_chi)
        self.phase_scat_max.set_visible(True)
        self.phase_scat_min.set_visible(True)

        # Inclusion circles
        self._phase_incl_circles.set_data(q_all[incl], chi_n[incl])
        self._phase_incl_circles.set_visible(True)

        # ── Quadratic fit: χ = A·q² + B·q⁴  (no constant term) ──────────────
        from scipy.optimize import curve_fit as _curve_fit

        def _chi_model(q, A, B):
            return A * q ** 2 + B * q ** 4

        # Clear any CTF-annotation artists left from a previous draw
        for _a in getattr(self, '_phase_ctf_markers', []):
            try: _a.remove()
            except Exception: pass
        self._phase_ctf_markers = []
        _ctf_annot = {}   # filled when underfocus fit is available

        result_str = ""
        try:
            q_fit   = q_all[incl]
            chi_fit = chi_n[incl]
            if len(q_fit) >= 2:
                # Build p0 from user-supplied guesses (Dz in nm, Cs in mm)
                try:
                    dz_guess_nm = float(self.var_phase_dz_guess.get())
                except (ValueError, AttributeError):
                    dz_guess_nm = 500.0
                try:
                    cs_guess_mm = float(self.var_phase_cs_guess.get())
                except (ValueError, AttributeError):
                    cs_guess_mm = 1.2
                cs_guess_nm = cs_guess_mm * 1e6
                A0 = -np.pi * wl * dz_guess_nm
                B0 =  0.5 * np.pi * cs_guess_nm * wl ** 3
                popt, pcov = _curve_fit(_chi_model, q_fit, chi_fit, p0=[A0, B0])
                A, B = float(popt[0]), float(popt[1])
                var = np.clip(np.diag(pcov), 0.0, None)
                perr = np.sqrt(var)
                Dz_nm = -A / (np.pi * wl)              # nm
                Cs_nm = 2.0 * B / (np.pi * wl ** 3)    # nm
                Cs_mm = Cs_nm / 1e6
                dDz   = perr[0] / (np.pi * wl)
                dCs   = 2.0 * perr[1] / (np.pi * wl ** 3) / 1e6

                import math as _math2
                def _ve2(val, err):
                    """Round error to 2 sig figs; match value precision."""
                    if not np.isfinite(err) or err <= 0 or abs(err) > 1e6:
                        return f"{val:>+9.4f}", "  [n/a]  "
                    mag = _math2.floor(_math2.log10(abs(err)))
                    dec = max(0, 2 - 1 - mag)
                    return f"{round(val, dec):>+9.{dec}f}", f"{round(err, dec):>9.{dec}f}"

                Av, Ae = _ve2(A,     perr[0])
                Bv, Be = _ve2(B,     perr[1])
                Dv, De = _ve2(Dz_nm, dDz)
                Cv, Ce = _ve2(Cs_mm, dCs)

                # Fit curve from q=0 to beyond last ring
                q_lo    = 0.0
                q_hi    = q_all.max() * 1.08
                q_curve = np.linspace(q_lo, q_hi, 400)
                chi_curve = _chi_model(q_curve, A, B)
                self.phase_fit_curve.set_data(q_curve, chi_curve)
                self.phase_fit_curve.set_visible(True)

                # Solve A·u + B·u² = χ_target for u = q² (only meaningful for underfocus, A<0, B>0)
                if A < 0 and B > 0:
                    def _q_for_chi(chi_t):
                        disc = A ** 2 + 4.0 * B * chi_t
                        if disc < 0:
                            return None
                        u1 = (-A - np.sqrt(disc)) / (2.0 * B)
                        u2 = (-A + np.sqrt(disc)) / (2.0 * B)
                        cands = [u for u in (u1, u2) if u > 1e-20]
                        return float(np.sqrt(min(cands))) if cands else None
                    _ctf_annot['q_min']  = _q_for_chi(-np.pi / 2.0)   # sin(χ) = −1
                    _ctf_annot['q_zero'] = _q_for_chi(-np.pi)           # sin(χ) = 0

                result_str = (
                    f"  χ = A·q² + B·q⁴\n"
                    f"  CTF = sin(χ)\n"
                    f"  Δz>0 underfocus, Δz<0 overfocus\n"
                    f"  {'':.<10} {'Fitted':>9}  {'±Err':>9}\n"
                    f"  {'A (r·nm²)':<10} {Av}  {Ae}\n"
                    f"  {'B (r·nm⁴)':<10} {Bv}  {Be}\n"
                    f"  {'─'*32}\n"
                    f"  {'Dz (nm)':<10} {Dv}  {De}\n"
                    f"  {'Cs (mm)':<10} {Cv}  {Ce}"
                )
            else:
                self.phase_fit_curve.set_data([], [])
                self.phase_fit_curve.set_visible(False)
                result_str = "  Need ≥2 included\n  points to fit."
        except Exception as exc:
            self.phase_fit_curve.set_data([], [])
            self.phase_fit_curve.set_visible(False)
            result_str = f"  Fit error:\n  {exc}"

        self._phase_result_text.set_text(result_str)
        self._phase_result_text.set_visible(bool(result_str))

        # Auto-scale axes to data
        all_q   = q_all
        all_chi = chi_n
        if len(all_q) >= 2:
            xpad = all_q.max() * 0.08
            chi_range = all_chi.max() - all_chi.min()
            ypad = chi_range * 0.20 if chi_range > 1e-9 else abs(all_chi.max()) * 0.20
            # Extend y-limit to accommodate −π level if underfocused
            y_lo = all_chi.min() - ypad
            if _ctf_annot:
                y_lo = min(y_lo, -np.pi - ypad * 0.5)
            ax.set_xlim(0.0, all_q.max() + xpad)
            ax.set_ylim(y_lo, all_chi.max() + ypad)

        # Draw CTF minimum and first zero-crossing markers for underfocus
        # Use explicit line segments (not axhline/axvline) so they are
        # bounded to the data range and scoped only to this axes.
        if _ctf_annot:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x0, x1 = xlim
            y0, y1 = ylim
            yspan = y1 - y0
            for chi_level, q_val, clr, lbl in (
                (-np.pi / 2.0, _ctf_annot.get('q_min'),  '#ff5555', 'CTF min  χ=−π/2'),
                (-np.pi,        _ctf_annot.get('q_zero'), '#5599ff', 'zero  χ=−π'),
            ):
                if q_val is None or q_val <= x0 or q_val >= x1:
                    continue
                # Horizontal dashed line from x=0 to x=q_val at height chi_level
                hl, = ax.plot([x0, q_val], [chi_level, chi_level],
                              color=clr, linewidth=1.1, linestyle='--',
                              alpha=0.85, zorder=4)
                self._phase_ctf_markers.append(hl)
                # Vertical dotted drop-line from y=y0 up to chi_level
                vl, = ax.plot([q_val, q_val], [y0, chi_level],
                              color=clr, linewidth=1.1, linestyle=':',
                              alpha=0.85, zorder=4)
                self._phase_ctf_markers.append(vl)
                # Dot at the intersection on the fit curve
                dot, = ax.plot([q_val], [chi_level], 'o',
                               color=clr, markersize=5, zorder=5)
                self._phase_ctf_markers.append(dot)
                # q label just to the right of the vertical line, near bottom
                qt = ax.text(q_val + (x1 - x0) * 0.015, y0 + yspan * 0.03,
                             f'q = {q_val:.3f} nm⁻¹', color=clr, fontsize=8,
                             va='bottom', ha='left', rotation=90,
                             fontfamily='monospace', zorder=6)
                self._phase_ctf_markers.append(qt)
                # χ-level label at the left end of horizontal line
                lt = ax.text(x0 + (x1 - x0) * 0.01, chi_level + yspan * 0.012,
                             lbl, color=clr, fontsize=8,
                             va='bottom', ha='left',
                             fontfamily='monospace', zorder=6)
                self._phase_ctf_markers.append(lt)

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

            chi = (-pi * wl * defocus * k2
                   + 0.5 * pi * cs_nm * (wl**3) * k4
                   + phase)
            delta_nm = cc_nm * (dE / kv_ev)
            temporal = np.exp(-0.5 * (pi**2) * (wl**2) * (delta_nm**2) * k4)
            t = -defocus * k + cs_nm * (wl**2) * k3
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
            p0.append(defocus_guess); lower_bounds.append(-np.inf); upper_bounds.append(np.inf)
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
        # Reimer n-start and guess mode
        if hasattr(self, 'var_guess_mode'):
            self.reimer_push_guesses = (self.var_guess_mode.get() == "Reimer Auto")
        if hasattr(self, 'var_reimer_mode'):
            self.reimer_auto = (self.var_reimer_mode.get() == "Auto")
            try:
                self.reimer_n_start_manual = int(self.var_reimer_n_start.get())
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

            import math as _math

            def _ve(val, err, use_exp=False):
                """Format value and error: error rounded to 2 sig figs,
                value shown to the same decimal precision as the error."""
                if err is None:
                    if use_exp:
                        return f"{val:>9.2e}", "  [Fixed]"
                    mag = _math.floor(_math.log10(abs(val))) if val != 0 else 0
                    dec = max(0, 2 - 1 - mag)
                    return f"{val:>9.{dec}f}", "  [Fixed]"
                if not np.isfinite(err) or err <= 0:
                    if use_exp:
                        return f"{val:>9.2e}", "  [n/a]  "
                    return f"{val:>9.4f}", "  [n/a]  "
                if use_exp:
                    # For scientific notation parameters (bg terms) keep 2e format
                    return f"{val:>9.2e}", f"{err:>9.2e}"
                mag = _math.floor(_math.log10(abs(err)))
                dec = max(0, 2 - 1 - mag)  # 2 sig figs in error
                val_r = round(val, dec)
                err_r = round(err, dec)
                return f"{val_r:>9.{dec}f}", f"{err_r:>9.{dec}f}"

            dz_v,  dz_e  = _ve(popt_full['defocus'],       perr_full['defocus'])
            cs_v,  cs_e  = _ve(popt_full['cs'],             perr_full['cs'])
            cc_v,  cc_e  = _ve(popt_full['cc'],             perr_full['cc'])
            de_v,  de_e  = _ve(popt_full['dE'],             perr_full['dE'])
            al_v,  al_e  = _ve(popt_full['alpha'],          perr_full['alpha'])
            ph_v,  ph_e  = _ve(popt_full['phase'] * 1000,
                                perr_full['phase'] * 1000 if perr_full['phase'] is not None else None)
            phd_v, phd_e = _ve(np.degrees(popt_full['phase']),
                                np.degrees(perr_full['phase']) if perr_full['phase'] is not None else None)
            amp_v, amp_e = _ve(popt_full['amp'],            perr_full['amp'],  use_exp=True)
            bg_v,  bg_e  = _ve(popt_full['bg'],             perr_full['bg'],   use_exp=True)
            bga_v, bga_e = _ve(popt_full['bg_amp'],         perr_full['bg_amp'], use_exp=True)
            bgs_v, bgs_e = _ve(popt_full['bg_slope'],       perr_full['bg_slope'], use_exp=True)
            table_lines = [
                f"── Fitted Parameters ─────────",
                f"  {'':.<12} {'Fitted':>9}  {'Error':>9}",
                f"  {'Defocus':<12} {dz_v}  {dz_e}",
                f"  {'Cs (mm)':<12} {cs_v}  {cs_e}",
                f"  {'Cc (mm)':<12} {cc_v}  {cc_e}",
                f"  {'dE (eV)':<12} {de_v}  {de_e}",
                f"  {'Alpha (mrad)':<12} {al_v}  {al_e}",
                f"  {'Phase (mrad)':<12} {ph_v}  {ph_e}",
                f"  {'Phase (deg)':<12} {phd_v}  {phd_e}",
                f"  {'CTF Amp':<12} {amp_v}  {amp_e}",
                f"  {'Bg Offset':<12} {bg_v}  {bg_e}",
            ]
            if ctx.get('bg_model') in ("Exponential", "Lin + Exp", "GaussExp"):
                table_lines.append(f"  {'Bg Amp':<12} {bga_v}  {bga_e}")
                slope_lbl = "Bg Xi (1/nm)" if ctx.get('bg_model') == "GaussExp" else "Bg Decay"
                table_lines.append(f"  {slope_lbl:<12} {bgs_v}  {bgs_e}")
            else:
                table_lines.append(f"  {'Bg Slope':<12} {bgs_v}  {bgs_e}")
            if ctx.get('bg_model') == "Lin + Exp":
                bgls_v, bgls_e = _ve(popt_full['bg_linslope'], perr_full.get('bg_linslope'), use_exp=True)
                table_lines.append(f"  {'Bg Lin Slp':<12} {bgls_v}  {bgls_e}")
            elif ctx.get('bg_model') == "GaussExp":
                bgsig_v, bgsig_e = _ve(popt_full['bg_sigma'], perr_full.get('bg_sigma'), use_exp=True)
                table_lines.append(f"  {'Bg Sigma':<12} {bgsig_v}  {bgsig_e}")
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

        # Information limit (chromatic, energy-spread only):
        #   d_info = (π λ Δ · 1/(2√ln2))^(1/2)  [Reimer formula]
        #   where Δ = Cc · (ΔE / E₀)  (focal spread; voltage & current fluctuations excluded)
        #   → q_info = 1/d_info = (2√ln2 / (π λ Δ))^(1/2)
        cc_nm_sch  = popt_full['cc'] * 1e6
        dE_sch     = popt_full['dE']
        kv_ev_sch  = ctx['current_kv_ev']
        delta_nm   = cc_nm_sch * (dE_sch / kv_ev_sch)   # focal spread Δ (nm)
        if delta_nm > 0 and wl > 0:
            d_info  = np.sqrt(np.pi * wl * delta_nm / (2.0 * np.sqrt(np.log(2))))
            q_info  = 1.0 / d_info if d_info > 0 else 0.0
        else:
            q_info = 0.0; d_info = 0.0; delta_nm = 0.0

        sch_lines = [
            f"── Scherzer ──────────────────",
            f"  {'Defocus':<12} {s_d:>9.1f} nm",
            f"  {'Sp. Freq':<12} {s_k:>9.3f} 1/nm",
            f"  {'Res Limit':<12} {s_r:>9.3f} nm",
        ]
        if q_info > 0:
            _lbl_delta = "Focal Spr \u0394"
            _lbl_dinf  = "d\u1d62\u207f\u2099\u2092 (nm)"
            sch_lines += [
                "\u2500\u2500 Info Limit (\u03c0\u03bb\u0394/2\u221aln2)\u00bd \u2500\u2500",
                f"  {_lbl_delta:<12} {delta_nm:>9.3f} nm",
                f"  {'Sp. Freq':<12} {q_info:>9.3f} 1/nm",
                f"  {_lbl_dinf:<12} {d_info:>9.3f} nm",
            ]
        self.sch_text.set_text("\n".join(sch_lines))
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

        if q_info > 0:
            r_info_px = q_info * ctx['N_dimension'] * ctx['base_cal_nm_per_px'] * ctx['bin_factor']
            self.info_limit_circle.set_radius(r_info_px)
            self.info_limit_circle.set_visible(True)
        else:
            self.info_limit_circle.set_visible(False)

        self.table_text.set_text(table_str)
        self.table_text.set_color("#fb923c" if fit_success else "#ef4444")

        # --- Envelope curves ---
        k = k_axis
        k2 = k * k; k3 = k2 * k; k4 = k2 * k2
        cc_nm = popt_full['cc'] * 1e6; cs_nm = popt_full['cs'] * 1e6
        alpha_rad = popt_full['alpha'] * 1e-3
        dn = cc_nm * (popt_full['dE'] / ctx['current_kv_ev'])
        temporal_env = np.exp(-0.5 * (np.pi**2) * (wl**2) * (dn**2) * k4)
        tt = -popt_full['defocus'] * k + cs_nm * (wl**2) * k3
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
            rmin, rmax = float(res.min()), float(res.max())
            pad = max(abs(rmax - rmin) * 0.1, 1e-5)
            self.ax_res.set_ylim(rmin - pad, rmax + pad)

            # Overlay the raw CTF data scaled to the residual y-range so both
            # curves share the same panel without axis distortion.
            data_fit = ctx['data_fit']
            if len(data_fit) > 0:
                d_min, d_max = float(data_fit.min()), float(data_fit.max())
                d_range = d_max - d_min if d_max > d_min else 1.0
                r_range = (rmax - rmin) if rmax > rmin else 1.0
                data_scaled = (data_fit - d_min) / d_range * r_range + rmin
                self.line_res_data.set_data(k_fit, data_scaled)
            else:
                self.line_res_data.set_data([], [])
        else:
            self.line_res.set_data([], [])
            self.line_res_data.set_data([], [])

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
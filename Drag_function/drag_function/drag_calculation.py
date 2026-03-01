import matplotlib.pyplot as plt
import numpy as np

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
from config_utils_local import config as C
from drag_function import io


dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Actual Run Tests
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

from dataclasses import dataclass
from scipy.signal import savgol_filter

# ---------------------------
# Unit conversion constants
# ---------------------------
AMU_TO_KG = 1.66053906660e-27
A_TO_M = 1.0e-10
PS_TO_S = 1.0e-12
EV_TO_J = 1.602176634e-19

# Convert (amu * Å/ps^2) to (eV/Å)
# Start: amu * Å/ps^2
# -> kg * m/s^2 = N
# -> eV/Å
AMU_A_PS2_TO_EV_PER_A = (
    AMU_TO_KG * (A_TO_M / (PS_TO_S**2))   # kg * m/s^2 per (amu * Å/ps^2) = N
    * (1.0 / EV_TO_J)                     # eV/m per N
    * A_TO_M                              # eV/Å per (eV/m)
)

# Coulomb constant in eV·Å for elementary charges
KE_EV_A = 14.3996454784255


@dataclass
class DragExtractionParams:
    window_length: int = 1001  # must be odd
    polyorder: int = 3
    q1: int = 1
    q2: int = 1
    mu_eff_amu: float = None   # effective reduced mass in amu (must be set)
    t_min_ps: float = 0.0      # optional transient exclusion
    r_min_A: float = 0.0       # optional: exclude very small R
    r_max_A: float = np.inf    # optional: exclude very large R
    show_plot: bool = True    # if True, display control plot of raw vs smoothed R


def _check_uniform_dt(t: np.ndarray, tol: float = 1e-3) -> float:
    dt = np.diff(t)
    rel = np.std(dt) / np.mean(dt)
    if rel > tol:
        raise ValueError(f"Non-uniform dt detected: std/mean = {rel:.2e}. "
                         "SG derivative assumes uniform sampling.")
    return float(np.mean(dt))


def _ensure_odd(n: int) -> int:
    return n if (n % 2 == 1) else n + 1


def extract_drag_from_R(dict_data: dict, params: DragExtractionParams):
    """
    Expects dict_data to contain:
      - 't' : time in ps
      - 'R' : separation in Å
    Optional:
      - could also include v1/v2 but not needed here

    Returns:
      arrays for:
        t_sel, R_smooth, vR, aR, F_drive, F_drag
    """
    if params.mu_eff_amu is None:
        raise ValueError("params.mu_eff_amu must be set (effective reduced mass in amu).")

    t = np.asarray(dict_data["t"], dtype=float)
    R = np.asarray(dict_data["R"], dtype=float)

    dt = _check_uniform_dt(t)

    wl = _ensure_odd(params.window_length)
    if wl >= len(t):
        raise ValueError(f"window_length={wl} must be smaller than data length={len(t)}.")
    if params.polyorder >= wl:
        raise ValueError("polyorder must be < window_length.")

    # Smooth R(t)
    R_s = savgol_filter(R, window_length=wl, polyorder=params.polyorder, deriv=0, delta=dt, mode="interp")

    # Optional control plot comparing raw R and smoothed R (non-fatal)
    if getattr(params, "show_plot", False):
        try:
            plt.figure(figsize=(8, 4))
            plt.plot(t, R, 'r-', label="R (raw)")
            plt.plot(t, R_s, 'b--', lw=1.5, label="R (smoothed)")
            plt.xlabel("t (ps)")
            plt.ylabel("R (Å)")
            plt.title("Control: raw R vs smoothed R")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()
        except Exception:
            # Don't let plotting errors break the extraction
            pass

    # First derivative: vR = dR/dt
    vR = savgol_filter(R, window_length=wl, polyorder=params.polyorder, deriv=1, delta=dt, mode="interp")

    # Compare vR to sum of velocities from the input dict (v1 + v2)
    if ("v1" in dict_data) and ("v2" in dict_data):
        v1_arr = np.asarray(dict_data["v1"], dtype=float)
        v2_arr = np.asarray(dict_data["v2"], dtype=float)
        # Align lengths: if lengths match t, sum directly; otherwise interpolate across the same time span
        if len(v1_arr) == len(t) and len(v2_arr) == len(t):
            v_sum = v1_arr + v2_arr
        else:
            # Map v1/v2 linearly onto t assuming they span the same time interval
            try:
                v1_on_t = np.interp(t, np.linspace(t[0], t[-1], len(v1_arr)), v1_arr)
                v2_on_t = np.interp(t, np.linspace(t[0], t[-1], len(v2_arr)), v2_arr)
                v_sum = v1_on_t + v2_on_t
            except Exception:
                v_sum = None

        if (v_sum is not None) and getattr(params, "show_plot", False):
            try:
                plt.figure(figsize=(8, 4))
                plt.plot(t, vR, color="#ff7f0e", lw=1.5, label="vR (smoothed)")
                plt.plot(t, v_sum, color="#2ca02c", alpha=0.9, label="v1+v2 (dict)")
                plt.xlabel("t (ps)")
                plt.ylabel("velocity (Å/ps)")
                plt.title("Control: vR vs v1+v2")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.show()
            except Exception:
                pass

    # Second derivative: aR = d2R/dt2
    aR = savgol_filter(R, window_length=wl, polyorder=params.polyorder, deriv=2, delta=dt, mode="interp")

    # Selection mask (transient and R-range)
    mask = (t >= params.t_min_ps) & (R_s >= params.r_min_A) & (R_s <= params.r_max_A)

    t_sel = t[mask]
    R_sel = R_s[mask]
    vR_sel = vR[mask]
    aR_sel = aR[mask]

    # Coulomb drive along R (magnitude), in eV/Å
    F_drive = (KE_EV_A * params.q1 * params.q2) / (R_sel**2)

    # Inertial term mu * aR converted to eV/Å
    F_inert = (params.mu_eff_amu * aR_sel) * AMU_A_PS2_TO_EV_PER_A

    # Drag (magnitude along the separation coordinate)
    F_drag = F_drive - F_inert

    return {
        "t": t_sel,
        "R": R_sel,
        "vR": vR_sel,
        "aR": aR_sel,
        "F_drive": F_drive,
        "F_inert": F_inert,
        "F_drag": F_drag,
        "dt": dt,
        "window_length": wl,
        "polyorder": params.polyorder,
    }

my_settings = DragExtractionParams(
    mu_eff_amu=105,
    t_min_ps=2,
)

#test = extract_drag_from_R(dict18, my_settings)


def sg_smooth_v(t, v, window_length, polyorder=3):
    t = np.asarray(t, float)
    v = np.asarray(v, float)

    # Ensure window_length is a plain int for savgol_filter and for parity checks
    window_length = int(window_length)

    dt = np.mean(np.diff(t))
    if not np.allclose(np.diff(t), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("t is not uniformly spaced enough for SG derivative.")

    if window_length % 2 == 0:
        window_length += 1
    if window_length >= len(v):
        raise ValueError("window_length must be smaller than data length.")
    if polyorder >= window_length:
        raise ValueError("polyorder must be < window_length.")

    v_s = savgol_filter(v, window_length=window_length, polyorder=polyorder, deriv=0, delta=dt, mode="interp")
    resid = v - v_s
    return v_s, resid, dt



def plot_test_smoothing(dictionary, atom, z_component = False,  polyorder = 3):
    v_list = []
    resid_list = []
    t_list = []
    wl = []
    if z_component:
        y = dictionary[f"v{atom}_z"]
    else:
        y = dictionary[f"v{atom}"]

    for i in range(500, 8500, 2000):
        v_s, resid, dt = sg_smooth_v(dictionary["t"], y, window_length=i, polyorder = polyorder)
        v_list.append(v_s)
        wl.append(i)
        resid_list.append(resid)
        t_list.append(dt)
        t_full = np.asarray(dictionary["t"], float)
        v_full = np.asarray(y, float)
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))

        # Left: overlay raw v1 with smoothed signals
        axs[0].plot(t_full, v_full, color='k', lw=0.8, label=f'v{atom} (raw)')
        axs[0].plot(t_full, v_s, 'b--', lw=1.2, alpha=0.9, label=f'v{atom}_s wl={i}')
        axs[0].set_xlabel('t (ps)')
        axs[0].set_ylabel('velocity (Å/ps)')
        axs[0].set_title(f'v{atom}: raw vs smoothed')
        axs[0].legend(fontsize='small')
        axs[0].grid(True)

        # Right: residuals (v1 - v1_s) for each smoothing; include RMS in legend
        rms = np.sqrt(np.mean(resid ** 2))
        axs[1].plot(t_full, resid, lw=1.0, alpha=0.9, label=f'wl={i}  RMS={rms:.3e}')
        axs[1].axhline(0.0, color='k', lw=0.8, alpha=0.6)
        axs[1].set_xlabel('t (ps)')
        axs[1].set_ylabel('residual (Å/ps)')
        axs[1].set_title(f'Residuals: v{atom} - v{atom}_s')
        axs[1].legend(fontsize='small')
        axs[1].grid(True)

        plt.tight_layout()
        plt.show()


    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    t_full = np.asarray(dictionary["t"], float)
    axs[0].plot(t_full, v_full, color='k', lw=0.7, alpha=0.95, label='v1 (raw)')
    for i in range(len(v_list)):
        axs[0].plot(t_full, v_list[i], ls = '--', lw=1.2, label=f'v{atom}_s wl={wl[i]}')
        axs[0].set_xlabel('t (ps)')
        axs[0].set_ylabel('velocity (Å/ps)')
        axs[0].set_title('v1: smoothed')
        axs[0].legend(fontsize='small')
        axs[0].grid(True)

        rms = np.sqrt(np.mean(resid_list[i] ** 2))
        axs[1].plot(t_full, resid_list[i], lw=1.0, alpha=0.9, label=f'wl={wl[i]}  RMS={rms:.3e}')
        axs[1].axhline(0.0, color='k', lw=0.8, alpha=0.6)
        axs[1].set_xlabel('t (ps)')
        axs[1].set_ylabel('residual (Å/ps)')
        axs[1].set_title('Residuals: v1 - v1_s')
        axs[1].legend(fontsize='small')
        axs[1].grid(True)
    plt.tight_layout()
    plt.show()

#plot_test_smoothing(dict18, 2, z_component = False, polyorder = 3)


def dominant_periods_fft(t, y, t_min=None, t_max=None, detrend='linear', n_peaks=3, fmin=0.05, fmax=None):
    """
    Returns dominant periods (ps) from FFT power spectrum of y(t).

    Parameters
    ----------
    t : array, ps
    y : array
    t_min, t_max : float or None
        restrict time window (helps focus on quasi-stationary region)
    detrend : {'none','mean','linear'}
    n_peaks : int
        number of peaks to report
    fmin : float
        ignore frequencies below this (1/ps) to avoid drift dominating
    fmax : float or None
        optional max frequency (1/ps), default Nyquist
    show : bool
        plot spectrum and mark peaks

    Returns
    -------
    peaks : list of (freq[1/ps], period[ps], power)
    freqs, power : arrays for further analysis
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)

    # Time window
    mask = np.ones_like(t, dtype=bool)
    if t_min is not None:
        mask &= (t >= t_min)
    if t_max is not None:
        mask &= (t <= t_max)

    tt = t[mask]
    yy = y[mask]

    # Uniform dt check
    dt = np.mean(np.diff(tt))
    if not np.allclose(np.diff(tt), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("t is not uniformly spaced enough for FFT peak detection.")

    # Detrend
    if detrend == 'mean':
        yy = yy - np.mean(yy)
    elif detrend == 'linear':
        A = np.vstack([tt, np.ones_like(tt)]).T
        m, c = np.linalg.lstsq(A, yy, rcond=None)[0]
        yy = yy - (m*tt + c)
    elif detrend == 'none':
        pass
    else:
        raise ValueError("detrend must be 'none', 'mean', or 'linear'.")

    # Apply a mild window to reduce leakage
    w = np.hanning(len(yy))
    yyw = yy * w

    # One-sided FFT
    Y = np.fft.rfft(yyw)
    freqs = np.fft.rfftfreq(len(yyw), d=dt)  # 1/ps
    power = (np.abs(Y)**2)

    # Frequency limits
    nyq = 0.5/dt
    if fmax is None:
        fmax = nyq
    band = (freqs >= fmin) & (freqs <= fmax)
    freqs_b = freqs[band]
    power_b = power[band]

    # Find peak candidates (simple: sort by power)
    idx_sorted = np.argsort(power_b)[::-1]
    peaks = []
    used = []
    for idx in idx_sorted:
        f = freqs_b[idx]
        p = power_b[idx]
        # avoid near-duplicates
        if any(abs(f - fu) < 0.02 for fu in used):
            continue
        used.append(f)
        peaks.append((float(f), float(1.0/f), float(p)))
        if len(peaks) >= n_peaks:
            break

    return peaks, freqs_b, power_b

import numpy as np
from scipy.signal import find_peaks

def dominant_periods_fft(
    t, y,
    t_min=None, t_max=None,
    detrend='linear',
    n_peaks=3,
    fmin=0.05,
    fmax=None,
    # new peak-finding controls:
    peak_min_prominence_ratio=0.05,  # relative to max(power_b)
    peak_min_distance_hz=0.05,       # min separation in 1/ps (frequency units)
):
    """
    Returns dominant periods (ps) from FFT power spectrum of y(t), using local maxima peak finding.

    Returns:
      peaks : list of (freq[1/ps], period[ps], power)
      freqs_b, power_b : arrays (already band-limited)
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)

    # Time window
    mask = np.ones_like(t, dtype=bool)
    if t_min is not None:
        mask &= (t >= t_min)
    if t_max is not None:
        mask &= (t <= t_max)

    tt = t[mask]
    yy = y[mask]

    # Uniform dt check
    dt = np.mean(np.diff(tt))
    if not np.allclose(np.diff(tt), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("t is not uniformly spaced enough for FFT peak detection.")

    # Detrend
    if detrend == 'mean':
        yy = yy - np.mean(yy)
    elif detrend == 'linear':
        A = np.vstack([tt, np.ones_like(tt)]).T
        m, c = np.linalg.lstsq(A, yy, rcond=None)[0]
        yy = yy - (m*tt + c)
    elif detrend == 'none':
        pass
    else:
        raise ValueError("detrend must be 'none', 'mean', or 'linear'.")

    # Window to reduce leakage
    w = np.hanning(len(yy))
    yyw = yy * w

    # One-sided FFT
    Y = np.fft.rfft(yyw)
    freqs = np.fft.rfftfreq(len(yyw), d=dt)  # 1/ps
    power = (np.abs(Y)**2)

    # Frequency limits
    nyq = 0.5 / dt
    if fmax is None:
        fmax = nyq
    band = (freqs >= fmin) & (freqs <= fmax)
    freqs_b = freqs[band]
    power_b = power[band]

    if len(freqs_b) < 5:
        return [], freqs_b, power_b

    # --- Local maxima peak finding ---
    # Convert desired minimum separation in frequency to bins
    df = freqs_b[1] - freqs_b[0]
    min_dist_bins = max(1, int(np.ceil(peak_min_distance_hz / df)))

    # Prominence threshold
    prom_abs = peak_min_prominence_ratio * np.max(power_b)

    peak_idx, props = find_peaks(power_b, prominence=prom_abs, distance=min_dist_bins)

    if peak_idx.size == 0:
        # fallback: return the single global max if no peak passes threshold
        i = int(np.argmax(power_b))
        f = float(freqs_b[i])
        return [(f, float(1.0/f), float(power_b[i]))], freqs_b, power_b

    # Rank peaks by prominence primarily (more robust than raw height), break ties by height
    prominences = props.get("prominences", np.zeros_like(peak_idx, dtype=float))
    heights = power_b[peak_idx]
    order = np.lexsort((heights, prominences))[::-1]  # descending
    peak_idx_sorted = peak_idx[order]

    # Build output list
    peaks = []
    for i in peak_idx_sorted[:n_peaks]:
        f = float(freqs_b[i])
        p = float(power_b[i])
        peaks.append((f, float(1.0/f), p))

    return peaks, freqs_b, power_b



def plot_spectrum(freqs, power, fmax=5.0, logy=False, t_min = None, t_max = None, title="FFT Power Spectrum"):
    title += f" (fmin={freqs[0]:.2f} 1/ps, fmax={fmax:.2f} 1/ps)"
    if t_max != None and t_min != None:
        title += f"\n (t={t_min:.1f}-{t_max:.1f} ps)"
    m = freqs <= fmax
    plt.figure(figsize=(8,4))
    plt.plot(freqs[m], power[m], lw=1.0)
    if logy:
        plt.yscale("log")
    plt.xlabel("frequency (1/ps)")
    plt.ylabel("power (arb.)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_spectrum_stack(t, y, t_min_list = (2.0, 3.0, 3.5, 4.0), t_max = 8.5, detrend="linear",
    n_peaks=3, fmin=0.05, fmax_plot=5.0, logy=False, figsize=(9, 10)):
    """
    Make a 4x1 plot: power spectra for different t_min values stacked vertically.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)

    n = len(t_min_list)
    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True)

    # If only one axis, wrap to list for uniform handling
    if n == 1:
        axes = [axes]

    for ax, t_min in zip(axes, t_min_list):
        peaks, freqs_b, power_b = dominant_periods_fft(
            t, y,
            t_min=t_min,
            t_max=t_max,
            detrend=detrend,
            n_peaks=n_peaks,
            fmin=fmin,
            fmax=None
        )

        m = freqs_b <= fmax_plot
        ax.plot(freqs_b[m], power_b[m], lw=1.0)

        if logy:
            ax.set_yscale("log")

        # Mark peaks
        for f, T, _p in peaks:
            if f <= fmax_plot:
                ax.axvline(f, ls="--", lw=1.0, alpha=0.8)
                ax.text(
                    f, 0.95, f"T≈{T:.2f} ps",
                    transform=ax.get_xaxis_transform(),
                    rotation=90, va="top", ha="right", fontsize=9
                )

        ax.set_title(f"Power spectrum (t_min={t_min} ps) | peaks: "
                     + ", ".join([f"{T:.2f} ps" for _, T, _ in peaks]))
        ax.grid(True)

    axes[-1].set_xlabel("frequency (1/ps)")
    fig.supylabel("power (arb.)")

    plt.tight_layout()
    plt.show()

# t_min = 2
# t_max = 9
# f_min = 0.05
# peaks, freqs_b, power_b = dominant_periods_fft(dict9["t"], dict9["v2"], t_min= t_min, t_max = t_max, detrend='linear', n_peaks=3, fmin=f_min)
# print(peaks)  # list of (freq, period, power)
# plot_spectrum(freqs_b, power_b, t_min= t_min, t_max = t_max, fmax=5,  logy=False)
# plot_spectrum_stack(dict9["t"], dict9["v2"], t_max = None, logy = False)

def oscillation_summary(data, t_star, t_out=9.0, fmin=0.3, fmax_plot=5.0, n_peaks=5, title=""):
    t = np.asarray(data["t"], float)
    v = np.asarray(data["v1"], float)

    peaks, freqs_b, power_b = dominant_periods_fft(
        t, v, t_min=t_star, t_max=t_out,
        detrend="linear", n_peaks=n_peaks, fmin=fmin
    )

    # plot (zoomed)
    m = freqs_b <= fmax_plot
    plt.figure(figsize=(8,4))
    plt.plot(freqs_b[m], power_b[m], lw=1.0)
    #plt.yscale("log")
    for f, T, _ in peaks:
        if f <= fmax_plot:
            plt.axvline(f, ls="--", lw=1.0)
            plt.text(f, np.max(power_b[m])*0.8, f"T≈{T:.2f} ps", rotation=90, va="top", fontsize=9)
    plt.xlabel("frequency (1/ps)")
    plt.ylabel("power (arb.)")
    plt.title(title + f" | window [{t_star},{t_out}] ps")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return peaks

# peaks_9  = oscillation_summary(dict9,  t_star=2.67, t_out=9.0, title="9Å: v1 spectrum")
# peaks_18 = oscillation_summary(dict18, t_star=4.90, t_out=8.5, title="18Å: v1 spectrum")
# print("9Å peaks:", peaks_9[:3])
# print("18Å peaks:", peaks_18[:3])


def smooth_range_only(t, v, t_start, t_end, window_length, polyorder=3, enable_plot=False):
    """
    Smooths a specific time range and optionally plots the result with boundaries.
    """
    t = np.asarray(t)
    v = np.asarray(v)

    # 1. Identify indices
    mask = (t >= t_start) & (t <= t_end)
    indices = np.where(mask)[0]

    if len(indices) == 0:
        raise ValueError(f"No data found in range {t_start} to {t_end}")

    idx_start, idx_end = indices[0], indices[-1]

    # 2. Extract and Smooth
    t_segment = t[idx_start: idx_end + 1]
    v_segment = v[idx_start: idx_end + 1]
    v_s_segment, resid, dt = sg_smooth_v(t_segment, v_segment, window_length, polyorder)

    # 3. Stitch
    v_full_smoothed = np.concatenate([
        v[:idx_start],
        v_s_segment,
        v[idx_end + 1:]
    ])

    # 4. Optional Diagnostic Plot
    if enable_plot:
        plt.figure(figsize=(10, 4))
        plt.plot(t, v, color='gray', alpha=0.4, label='Original (Full)')
        plt.plot(t, v_full_smoothed, color='#0072BD', lw=1.5, label='Stitched Result')

        # Mark boundaries
        plt.axvline(t_start, color='red', linestyle='--', alpha=0.7, label='Start Smoothing')
        plt.axvline(t_end, color='red', linestyle='--', alpha=0.7, label='End Smoothing')

        plt.title(f'Smoothing Verification: {t_start}s to {t_end}s')
        plt.xlabel('t / ps')
        plt.ylabel('v / Å/ps')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    return v_full_smoothed

v_part_smoothed = smooth_range_only(dict9["t"], dict9["v2"], t_start=2.67, t_end=8.2, window_length=4501, polyorder=3, enable_plot=True)
v_part_smoothed = smooth_range_only(dict18["t"], dict18["v2"], t_start=4.50, t_end=8.0, window_length=2501, polyorder=3, enable_plot=True)


import matplotlib.pyplot as plt
import numpy as np

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
from config_utils_local import config as C
from drag_function import io

home = True

if home:
    dict9 = io.load_data(C.PATH9A)
    dict18 = io.load_data(C.PATH18A)


else:
    dict9 = io.load_data(C.PATH9A_OFFICE)
    dict18 = io.load_data(C.PATH18A_OFFICE)

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Actual Run Tests
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
import numpy as np
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

    # First derivative: vR = dR/dt
    vR = savgol_filter(R, window_length=wl, polyorder=params.polyorder, deriv=1, delta=dt, mode="interp")

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

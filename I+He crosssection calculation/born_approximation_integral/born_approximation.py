"""
I+–He velocity-dependent differential & total cross section
MATLAB -> Python translation (SciPy/NumPy)

Notes:
- Distances R, r are in Angström unless explicitly converted.
- Energies from the Buchachenko et al. parametrization are in cm^-1 and
  converted to eV via eV_per_wavenumber.
- For the scattering integral, potential is converted from eV -> J.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import curve_fit
import physical_constants as pc

# -----------------------------
# Physical constants (SI)
# -----------------------------
hbar_SI = 1.05457182e-34  # J*s
e_charge = 1.602176634e-19  # C
epsilon_0 = 8.8541878128e-12  # F/m
eV_J = 1.602176634e-19  # J per eV
u = 1.66053906660e-27  # atomic mass unit (kg)

# -----------------------------
# Masses (amu)
# -----------------------------
m1 = 127.0  # iodine (amu) in your code
m2 = 4.0    # helium (amu)
mu_amu = m1 * m2 / (m1 + m2)  # reduced mass (amu)
mu = mu_amu * u               # reduced mass in kg


# -----------------------------
# Potential parameters (Buchachenko et al.)
# -----------------------------
alpha = 0.5132894e+1
beta = -0.8349857e+1

g = np.array([
    0.8833500e+7, -0.1624297e+8,  0.1156310e+8,
   -0.3647464e+7,  0.2590770e+6,  0.1324182e+6,
   -0.2782933e+5, -0.1278957e+2,  0.2773608e+3
], dtype=float)

delta = 0.7159647e+3
D4 = 0.1185215e+5
D6 = 0.1117964e+6
D8 = 0.7710727e+6

eV_per_wavenumber = 1.0 / 8065.54429  # eV per cm^-1


# -----------------------------
# Potential building blocks (R in Å)
# -----------------------------
l = np.arange(0, 9, dtype=int)

def VSR(R_ang: np.ndarray | float) -> np.ndarray | float:
    """
    Short-range part (in cm^-1, then you later convert).
    Matches MATLAB:
      sum( g(l+1) .* (R').^l .* exp(-alpha*R' - beta), 2)'
    """
    R = np.atleast_1d(R_ang).astype(float)
    # shape: (len(R), len(l))
    terms = (R[:, None] ** l[None, :]) * g[None, :] * np.exp(-alpha * R[:, None] - beta)
    out = np.sum(terms, axis=1)
    return out[0] if np.isscalar(R_ang) else out

def VLR(R_ang: np.ndarray | float) -> np.ndarray | float:
    """
    Long-range part in cm^-1.
    """
    R = np.asarray(R_ang, dtype=float)
    return -D4 * R**(-4) - D6 * R**(-6) - D8 * R**(-8)

def sw(R_ang: np.ndarray | float) -> np.ndarray | float:
    """
    Switch function (kept for completeness; you noted it's effectively irrelevant here).
    """
    R = np.asarray(R_ang, dtype=float)
    return 0.5 * (1.0 + np.tanh(1.0 + delta * R))

def V_original(R_ang: np.ndarray | float) -> np.ndarray | float:
    """
    Full potential in eV (MATLAB: (VLR*1 + VSR) * eV_per_wavenumber).
    """
    return (VLR(R_ang) + VSR(R_ang)) * eV_per_wavenumber

def softstep(R_ang: np.ndarray | float, R0: float, width: float = 0.01) -> np.ndarray | float:
    """
    Soft step: ~0 below R0, ~1 above R0.
    MATLAB: (erf((R-R0)/0.01)+1)/2
    """
    R = np.asarray(R_ang, dtype=float)
    return 0.5 * (erf((R - R0) / width) + 1.0)

def V_clipped(R_ang: np.ndarray | float, R_switch: float = 2.0) -> np.ndarray | float:
    """
    Modified potential well-behaved below R_switch.
    MATLAB:
      V = V_original(R).*softstep(R,R_switch) + V_original(R_switch)*(1-softstep(...))
    """
    s = softstep(R_ang, R_switch)
    return V_original(R_ang) * s + V_original(R_switch) * (1.0 - s)


# -----------------------------
# Lennard-Jones regularizations (r in Å)
# -----------------------------
eps = 0.01784
sig = 3.25 / (2.0 ** (1.0 / 6.0))
regpar = 0.05 * sig
rc = 0.05 * sig
n_damp = 8

import numpy as np

R_EPS = 1e-12  # Å, purely numerical guard

def _asarray_r(r_ang):
    r = np.asarray(r_ang, dtype=float)
    return np.maximum(r, R_EPS)

def V_LJ(eps_, sig_, r_ang):
    r = _asarray_r(r_ang)
    return 4.0 * eps_ * ((sig_ / r)**12 - (sig_ / r)**6)

def V_exp(eps_, sig_, rc_, n_, r_ang):
    """
    Exponential cutoff / damping:
      V_LJ * (1 - exp(-(r/rc)^n))
    Matches your MATLAB intent.
    """
    r = _asarray_r(r_ang)
    return V_LJ(eps_, sig_, r) * (1.0 - np.exp(- (r / rc_)**n_))


def V_soft(eps_: float, sig_: float, a: float, r_ang: np.ndarray | float) -> np.ndarray | float:
    """
    Soft-core regularization (kept from your MATLAB).
    """
    r = np.asarray(r_ang, dtype=float)
    num = (sig_**2 + a**2)
    den = (r**2 + a**2)
    return 4.0 * eps_ * ((num / den)**6 - (num / den)**3)

# -----------------------------
# Scattering integral setup
# -----------------------------
def scattering_amplitude_f(theta, v, potential_func, mu, hbar,
                           rmin_ang=0.0, rmax_ang=200.0, abstol=1e-8):
    k0 = (mu * v) / hbar
    K = 2.0 * k0 * np.sin(theta / 2.0)

    prefactor = -2.0 * mu / (hbar**2 * K)

    def integrand(z_ang):
        z_m = z_ang * pc.angstrom
        return z_m * np.sin(K * z_m) * potential_func(z_ang) * eV_J

    val, _ = quad(integrand, rmin_ang, rmax_ang, epsabs=abstol, limit=500)

    return prefactor * val



def compute_sigma(v_array: np.ndarray,
                  theta_array: np.ndarray,
                  potential_func,
                  mu: float,
                  hbar: float,
                  rmax_ang: float = 300.0) -> np.ndarray:
    """
    sigma(v,theta) = |f|^2 on a grid.
    """
    sigma = np.zeros((len(v_array), len(theta_array)), dtype=float)

    total = sigma.size
    counter = 0

    for i, v in enumerate(v_array):
        for j, theta in enumerate(theta_array):
            f = scattering_amplitude_f(theta, v, potential_func,
                                       mu=mu, hbar=hbar,
                                       rmax_ang=rmax_ang)
            sigma[i, j] = np.abs(f)**2

            counter += 1
            if counter % int(np.ceil(total / 20)) == 0:
                pct = 100.0 * counter / total
                print(f"calculation progress: {pct:.0f} percent")

    sigma[~np.isfinite(sigma)] = 0.0
    return sigma



# -----------------------------
# Main run (mirrors your MATLAB choices)
# -----------------------------
def main() -> None:
    # Plot potentials (modified and original)
    x = np.arange(0.01, 25.0 + 1e-12, 0.01)

    plt.figure()
    plt.plot(x, V_clipped(x, R_switch=2.0) * 1000.0, label="clipped")
    plt.plot(x, V_original(x) * 1000.0, label="original")
    plt.xlabel("R / Å")
    plt.ylabel("E / meV")
    plt.legend()
    plt.tight_layout()
    plt.close()

    # LJ regularization plots (matching your intention)
    plt.figure()
    plt.plot(x, V_LJ(eps, sig, x) * 1000.0, label="Lennard-Jones")
    plt.plot(x, V_soft(eps, sig, regpar, x) * 1000.0, label="soft-core")
    plt.plot(x, V_exp(eps, sig, rc, n_damp, x) * 1000.0, label="exp cutoff")
    plt.xlim(0.0, 5.0)
    plt.ylim(-20.0, 100.0)
    plt.xlabel("r / Å")
    plt.ylabel("E / meV")
    plt.legend()
    plt.tight_layout()
    plt.close()

    # Cross section grid
    v_array = np.arange(10.0, 2200.0 + 1e-12, 100.0)
    theta_array = np.linspace(-np.pi, np.pi, 100)
    theta_array = theta_array[np.abs(theta_array) > 1e-3]

    # Choose potential for scattering integral (your MATLAB uses V_exp in the integral)
    potential = lambda r_ang: V_exp(eps, sig, rc, n_damp, r_ang)

    sigma = compute_sigma(v_array, theta_array, potential_func=potential, mu=mu, hbar=hbar_SI, rmax_ang=1e12)

    # Total cross section by integrating over angles
    sigma_total = np.trapezoid(sigma, theta_array, axis=1)

    # Fit power law in log space: log(sigma) = a + b log(v)
    # (same structure as your MATLAB nlinfit block)
    def log_fit_func(logv, a, b):
        return a + b * logv

    # Avoid zeros in log
    mask = np.isfinite(sigma_total) & (sigma_total > 0)
    print("sigma_total finite:", np.isfinite(sigma_total).sum(), "/", sigma_total.size)
    print("sigma_total > 0:", (sigma_total > 0).sum(), "/", sigma_total.size)
    print("sigma_total min/max:", np.nanmin(sigma_total), np.nanmax(sigma_total))
    print("fit points:", mask.sum())

    if mask.sum() < 3:
        print("Not enough valid points for a 2-parameter fit. Skipping fit.")
    else:
        popt, _ = curve_fit(log_fit_func, np.log(v_array[mask]), np.log(sigma_total[mask]), p0=(np.log(np.max(sigma_total[mask])), -1.0))
        a_fit, b_fit = popt

    plt.figure()
    plt.plot(v_array, sigma_total, label="sigma_total")
    if mask.sum() > 3:
        plt.plot(v_array, np.exp(log_fit_func(np.log(v_array), a_fit, b_fit)),
                 label=f"fit: σ ≈ {np.exp(a_fit):.2e} · v^{b_fit:.2f}")
    plt.yscale("log")
    plt.xlabel("v / m/s")
    plt.ylabel("σ / m$^2$")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()

    # # Surface plot (v, theta, sigma)
    # # Matplotlib surface needs meshgrid
    # from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    #
    # Vv, Tt = np.meshgrid(v_array, theta_array, indexing="ij")
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection="3d")
    # ax.plot_surface(Vv, Tt, sigma, rstride=1, cstride=1, linewidth=0, antialiased=True)
    # ax.set_xlabel("v / m/s")
    # ax.set_ylabel("θ / rad")
    # ax.set_zlabel("σ / m$^2$")
    # ax.view_init(elev=45, azim=45)
    # plt.tight_layout()
    #
    # plt.show()


if __name__ == "__main__":
    main()

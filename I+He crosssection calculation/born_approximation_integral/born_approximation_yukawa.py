"""
I+–He velocity-dependent differential & total cross section
MATLAB -> Python (NumPy/SciPy)

This is a direct conversion of your Yukawa-test MATLAB script.

Units:
- r, z, x are in Angström (Å) unless converted
- alpha is in 1/Å in the definition of V(r)
- V0 is in eV
- In the Born integral we convert Å -> m and eV -> J to keep SI consistency
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import curve_fit


# -----------------------------
# Physical constants (SI)
# -----------------------------
hbar_SI = 1.05457182e-34  # J*s
fs = 1e-15                # s (unused here, kept to mirror MATLAB)
u = 1.66053906660e-27     # atomic mass unit (kg)
eV_J = 1.602176634e-19    # J per eV
ANGSTROM_TO_M = 1e-10     # m per Å


# -----------------------------
# Masses (amu) and reduced mass
# -----------------------------
m1 = 127.0
m2 = 4.0
mu_amu = m1 * m2 / (m1 + m2)   # reduced mass in amu
mu = mu_amu * u               # reduced mass in kg


# -----------------------------
# Yukawa potential for testing
# -----------------------------
alpha = 1.0   # 1/Å
V0 = 20.0     # eV*Å

def V_yukawa_eV(r_ang: np.ndarray | float) -> np.ndarray | float:
    """Yukawa potential V(r) in eV, with r in Å."""
    r = np.asarray(r_ang, dtype=float)
    # avoid division by zero
    r = np.maximum(r, 1e-12)
    return V0 * np.exp(-alpha * r) / r


# -----------------------------
# Plot the potential
# -----------------------------
x = np.arange(1.9, 25.0 + 1e-12, 0.01)
plt.figure()
plt.plot(x, V_yukawa_eV(x) * 1000.0)
plt.xlabel("r / Å")
plt.ylabel("E / meV")
plt.tight_layout()


# -----------------------------
# Direct integral (Born approximation)
# -----------------------------
rmin_ang = 0.0
rmax_ang = 100.0

v_array = np.arange(1.0, 2200.0 + 1e-12, 100.0)
theta_array = np.linspace(-np.pi, np.pi, 100)

sigma = np.zeros((len(v_array), len(theta_array)), dtype=float)
counter = 0
total = sigma.size


def scattering_amplitude(theta: float, v: float) -> float:
    """
    Implements the MATLAB:
      k0 = mu*u*v/hbar
      K = 2*k0*sin(theta/2)
      prefactor = -2*mu*u/hbar^2 / K
      f = prefactor * integral_0^rmax ( (z*1e-10)*sin((z*1e-10)*K)*V(z)*eV ) dz

    Here we keep full SI consistency:
      - z_m = z_ang * 1e-10
      - V(z) in eV -> Joule by * eV_J
      - IMPORTANT: dz_ang is in Å, so dr = (1e-10) dz_ang must be included.
        That means the integrand must include an extra factor ANGSTROM_TO_M.
    """
    k0 = mu * v / hbar_SI                 # 1/m
    K = 2.0 * k0 * np.sin(theta / 2.0)    # 1/m

    # handle the K->0 case robustly (prevents blow-ups)
    if np.abs(K) < 1e-20:
        # small-K limit: sin(Kr) ~ Kr, so prefactor*(...) -> finite
        # f(K->0) = -(2 mu / hbar^2) * ∫ r^2 V(r) dr
        def integrand0(z_ang: float) -> float:
            r_m = z_ang * ANGSTROM_TO_M
            V_J = V_yukawa_eV(z_ang) * eV_J
            return (r_m**2) * V_J * ANGSTROM_TO_M  # dr = 1e-10 dz
        I0, _ = quad(integrand0, rmin_ang, rmax_ang, epsabs=1e-8, limit=500)
        return (-2.0 * mu / (hbar_SI**2)) * I0

    prefactor = -2.0 * mu / (hbar_SI**2 * K)

    def integrand(z_ang: float) -> float:
        r_m = z_ang * ANGSTROM_TO_M
        V_J = V_yukawa_eV(z_ang) * eV_J
        # integral is over z_ang (Å), so include dr = 1e-10 dz
        return (r_m * np.sin(K * r_m) * V_J) * ANGSTROM_TO_M

    I, _ = quad(integrand, rmin_ang, rmax_ang, epsabs=1e-8, limit=500)
    return prefactor * I


# compute sigma(v,theta)
for i, v in enumerate(v_array):
    for j, theta in enumerate(theta_array):
        f = scattering_amplitude(theta, v)
        sigma[i, j] = np.abs(f)**2

        counter += 1
        if counter % int(np.ceil(total / 20)) == 0:
            pct = 100.0 * counter / total
            print(f"calculation progress: {pct:.0f} percent")

# remove NaNs/Infs
sigma[~np.isfinite(sigma)] = 0.0

# total cross section by integrating over angles
sigma_total = np.trapezoid(sigma, theta_array, axis=1)


# -----------------------------
# Fit v-dependence: log(sigma) = a + b log(v)
# -----------------------------
plt.figure()
normalization = 1.0

plt.plot(v_array, sigma_total / normalization, label="sigma_total")

mask = (sigma_total > 0) & np.isfinite(sigma_total)

def log_fit_func(logv, a, b):
    return a + b * logv

popt, _ = curve_fit(
    log_fit_func,
    np.log(v_array[mask]),
    np.log(sigma_total[mask]),
    p0=(np.log(np.max(sigma_total[mask])), -1.0)
)
a_fit, b_fit = popt

plt.plot(
    v_array,
    np.exp(log_fit_func(np.log(v_array), a_fit, b_fit)),
    label=f"fit: σ ≈ {np.exp(a_fit):.1e} × v^{b_fit:.2f}"
)

plt.yscale("log")
plt.xlabel("v / m/s")
plt.ylabel("σ / m$^{-2}$")
plt.legend()
plt.tight_layout()


# -----------------------------
# Compare to analytical Yukawa (Born) expression in your MATLAB
# -----------------------------
plt.figure()
plt.plot(v_array, sigma_total, label="numerical")
plt.xlabel("v / m/s")
plt.ylabel("σ / m$^{-2}$")

k = mu * v_array / hbar_SI  # 1/m

# Build sigma_yukawa with shape (len(v), len(theta)) to match MATLAB
# MATLAB:
# sigma_yukawa = 4*(mu*u)^2*(V0*eV)^2/hbar^4 * 1 / ( (alpha*1e10)^2 + 4*(k^2)'*sin(theta/2)^2 )^2
alpha_SI = alpha * 1e10  # 1/m (since alpha in 1/Å)
# Actually V_0 is in [eV*A], so we need to multiply by ANGSTROM_TO_M when converting to Joules*m
pref = 4.0 * (mu**2) * ((V0 * eV_J*ANGSTROM_TO_M)**2) / (hbar_SI**4)


den = (alpha_SI**2) + 4.0 * (k[:, None]**2) * (np.sin(theta_array[None, :] / 2.0)**2)
sigma_yukawa = pref / (den**2)

sigma_total_yukawa = np.trapezoid(sigma_yukawa, theta_array, axis=1)

plt.plot(v_array, sigma_total_yukawa, ":", label="analytical")
plt.legend()
plt.tight_layout()
plt.show()

print("mu [kg] =", mu)
print("alpha [1/Å] =", alpha, "alpha_SI [1/m] =", alpha_SI)
print("V0 [eV] =", V0, "V0 [J] =", V0 * eV_J)
print("k range [1/m] =", k.min(), k.max())
print("prefactor pref =", pref)

# -----------------------------
# Surface plots (numerical, analytical, difference)
# -----------------------------
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
#
# Vv, Tt = np.meshgrid(v_array, theta_array, indexing="ij")
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, projection="3d")
# ax1.plot_surface(Vv, Tt, sigma, linewidth=0, antialiased=True)
# ax1.set_xlabel("v / m/s")
# ax1.set_ylabel("θ / rad")
# ax1.set_zlabel("σ / m$^{-2}$")
# ax1.set_title("numerical solution")
# plt.tight_layout()
#
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111, projection="3d")
# ax2.plot_surface(Vv, Tt, sigma_yukawa, linewidth=0, antialiased=True)
# ax2.set_xlabel("v / m/s")
# ax2.set_ylabel("θ / rad")
# ax2.set_zlabel("σ / m$^{-2}$")
# ax2.set_title("analytical solution")
# plt.tight_layout()
#
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111, projection="3d")
# ax3.plot_surface(Vv, Tt, (sigma_yukawa - sigma), linewidth=0, antialiased=True)
# ax3.set_xlabel("v / m/s")
# ax3.set_ylabel("θ / rad")
# ax3.set_zlabel("Δσ / m$^{-2}$")
# ax3.set_title("difference")
# plt.tight_layout()
#
# plt.show()

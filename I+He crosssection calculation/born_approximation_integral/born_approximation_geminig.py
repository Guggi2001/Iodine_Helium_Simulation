import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# -----------------------------
# 1. Physical Constants
# -----------------------------
hbar_SI = 1.05457182e-34
eV_J = 1.602176634e-19
u = 1.66053906660e-27
angstrom = 1e-10

# Masses (I+ and He)
m1, m2 = 127.0, 4.0
mu = (m1 * m2 / (m1 + m2)) * u

# Soft-Core LJ Parameters
eps = 0.01784  # eV
sig = 3.25  # Angstrom
softening_a = 0.5  # Angstroms


def V_soft_LJ_eV_part(r_ang):
    """V(r) without the singularity."""
    r_sq = r_ang ** 2
    x_sq = r_sq + softening_a ** 2
    term6 = (sig ** 2 / x_sq) ** 3
    term12 = term6 ** 2
    return 4.0 * eps * (term12 - term6)


# -----------------------------
# 2. High-Limit Born Integrator
# -----------------------------
def get_sigma_mt_converged(v, rmax):
    k0 = (mu * v) / hbar_SI

    # Angular grid
    theta_array = np.linspace(0.01, np.pi, 100)
    d_sigmas = []

    for theta in theta_array:
        K = 2.0 * k0 * np.sin(theta / 2.0)
        prefactor = -2.0 * mu / (hbar_SI ** 2 * K)

        def integrand_core(z_ang):
            z_m = z_ang * angstrom
            V_J = V_soft_LJ_eV_part(z_ang) * eV_J
            # Return r * V(r)
            return z_m * V_J * angstrom

        # --- FIX: Handling Integration Bounds ---
        wvar_val = K * angstrom

        if K * angstrom < 1e-6:
            # Forward scattering limit (non-oscillatory)
            # If rmax is infinite, this part might diverge for pure 1/r^6
            # but we are doing MT cross section so (1-cos) kills this part anyway.
            # For safety, we cap forward scattering RMAX at 5000 if infinity is requested.
            eff_rmax = 5000.0 if np.isinf(rmax) else rmax

            def integrand_lowk(z_ang):
                return (z_ang * angstrom) ** 2 * V_soft_LJ_eV_part(z_ang) * eV_J * angstrom

            val, _ = quad(integrand_lowk, 0.0, eff_rmax, limit=1000)
            f = (-2.0 * mu / hbar_SI ** 2) * val

        else:
            # Oscillatory Integration
            # CRITICAL CHANGE: limit=10000 ensures we track high-frequency waves
            val, _ = quad(integrand_core, 0.0, rmax,
                          weight='sin', wvar=wvar_val,
                          limit=10000)  # Increased from 200
            f = prefactor * val

        d_sigmas.append(np.abs(f) ** 2)

    d_sigmas = np.array(d_sigmas)
    weight_mt = (1 - np.cos(theta_array)) * np.sin(theta_array)
    sigma_mt = 2 * np.pi * np.trapz(d_sigmas * weight_mt, theta_array)

    return sigma_mt


# -----------------------------
# 3. Final Execution
# -----------------------------
if __name__ == "__main__":
    v_test = 500.0
    print(f"--- Final Stability Check (v={v_test} m/s) ---")

    # We include np.inf to show the theoretical limit
    rmax_list = [50, 100, 200, 500, 1000, 2000, 5000, np.inf]
    results = []

    for r_mx in rmax_list:
        try:
            s = get_sigma_mt_converged(v_test, r_mx)
            results.append(s)
            label = "Infinity" if np.isinf(r_mx) else f"{r_mx}"
            print(f"  rmax={label:8s} -> sigma_mt={s:.4e} m^2")
        except Exception as e:
            print(f"  rmax={r_mx} failed: {e}")
            results.append(0.0)

    # Plot (excluding Infinity for the X-axis scale, but drawing it as a line)
    valid_x = [x for x in rmax_list if not np.isinf(x)]
    valid_y = results[:-1]
    inf_val = results[-1]

    plt.figure(figsize=(8, 6))
    plt.plot(valid_x, valid_y, 'o-', linewidth=2, label='Finite Integration')
    plt.axhline(inf_val, color='r', linestyle='--', label=f'Infinite Limit ({inf_val:.2e})')

    plt.title(f"Corrected Convergence of $\sigma_{{MT}}$\n(Soft-Core LJ, limit=10,000)")
    plt.xlabel("Integration Cutoff $r_{max}$ (Å)")
    plt.ylabel("Momentum Transfer Cross Section ($m^2$)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
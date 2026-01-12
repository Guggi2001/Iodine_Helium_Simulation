import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# -----------------------------
# 1. Constants & Setup
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

# CHANGE: Weaker core for validation to stay in Born regime
# a = 3.0 makes V(0) approx equal to epsilon (weak scattering)
softening_a = 3.0


def V_soft_LJ_eV_part(r_ang):
    """V(r) without the singularity."""
    r_sq = r_ang ** 2
    x_sq = r_sq + softening_a ** 2
    term6 = (sig ** 2 / x_sq) ** 3
    term12 = term6 ** 2
    return 4.0 * eps * (term12 - term6)


# -----------------------------
# 2. Split-Domain Integrator
# -----------------------------
def get_sigma_mt_split(v, rmax):
    k0 = (mu * v) / hbar_SI

    # We define a "Core Cutoff" where the potential is non-zero
    r_core_edge = 20.0  # Angstroms

    # Angular grid
    theta_array = np.linspace(0.01, np.pi, 100)
    d_sigmas = []

    for theta in theta_array:
        K = 2.0 * k0 * np.sin(theta / 2.0)
        prefactor = -2.0 * mu / (hbar_SI ** 2 * K)

        # Integrand Function (r * V)
        def integrand_core(z_ang):
            z_m = z_ang * angstrom
            V_J = V_soft_LJ_eV_part(z_ang) * eV_J
            return z_m * V_J * angstrom

        # --- PART 1: The Core (0 to 20 A) ---
        # Standard quad is fine here, range is small
        wvar_val = K * angstrom

        # Calculate Core contribution
        # Note: We use 'weight=sin' even here to be consistent
        val_core, _ = quad(integrand_core, 0.0, min(rmax, r_core_edge),
                           weight='sin', wvar=wvar_val, limit=500)

        # --- PART 2: The Tail (20 A to rmax) ---
        val_tail = 0.0
        if rmax > r_core_edge:
            val_tail, _ = quad(integrand_core, r_core_edge, rmax,
                               weight='sin', wvar=wvar_val, limit=5000)

        f = prefactor * (val_core + val_tail)
        d_sigmas.append(np.abs(f) ** 2)

    d_sigmas = np.array(d_sigmas)
    weight_mt = (1 - np.cos(theta_array)) * np.sin(theta_array)
    sigma_mt = 2 * np.pi * np.trapz(d_sigmas * weight_mt, theta_array)

    return sigma_mt


# -----------------------------
# 3. Execution
# -----------------------------
if __name__ == "__main__":
    v_test = 500.0
    print(f"--- Final Convergence Test (Split Integration) ---")

    # Note: We include np.inf
    rmax_list = [50, 100, 200, 500, 1000, 5000, np.inf]
    results = []

    for r_mx in rmax_list:
        s = get_sigma_mt_split(v_test, r_mx)
        results.append(s)
        label = "Inf" if np.isinf(r_mx) else f"{r_mx}"
        print(f"  rmax={label:5s} -> sigma_mt={s:.4e} m^2")

    # Plot
    valid_x = [x for x in rmax_list if not np.isinf(x)]
    valid_y = results[:-1]

    plt.figure(figsize=(8, 6))
    plt.plot(valid_x, valid_y, 'o-', linewidth=2, label='Finite Integration')
    plt.axhline(results[-1], color='r', linestyle='--', label=f'Infinite Limit ({results[-1]:.2e})')

    plt.title(f"Convergence with Split Integration\n(Softening a={softening_a} $\AA$)")
    plt.xlabel("Integration Cutoff $r_{max}$ (Å)")
    plt.ylabel("$\sigma_{MT}$ ($m^2$)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
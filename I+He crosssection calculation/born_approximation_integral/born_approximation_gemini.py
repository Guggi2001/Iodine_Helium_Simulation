import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad, trapezoid

# =============================================================================
# 1. PHYSICAL CONSTANTS
# =============================================================================
HBAR = 1.0545718e-34  # J s
AMU = 1.660539e-27  # kg
ANG = 1.0e-10  # m (Scaling factor)
EV = 1.6021766e-19  # J

# Reduced Mass (I + He)
M_I = 126.90 * AMU
M_HE = 4.0026 * AMU
MU = (M_I * M_HE) / (M_I + M_HE)

# =============================================================================
# 2. LENNARD-JONES PARAMETERS (Rescaled)
# =============================================================================
# We use the parameters you provided earlier:
EPS_eV = 0.01784  # Depth ~18 meV
SIG_ANG = 3.25  # Zero crossing at 3.25 A

# Pre-calculate Epsilon in Joules
EPS_J = EPS_eV * EV

# HYBRID CUTOFF:
# Born approx is valid only when potential is "weak" (V << Kinetic Energy).
# At 10 m/s, Energy ~ 0.1 meV.
# The LJ potential V(r) drops to ~0.1 meV at roughly r = 10 Angstroms.
# We set the switch there. Inside is Hard Sphere; Outside is Born Tail.
R_SWITCH_ANG = 10


def lj_integrand_scaled(r_ang):
    """
    Returns r * V(r).
    Inputs and outputs are in ANGSTROMS and JOULES (no meters).
    """
    # Standard LJ: 4*eps * [ (sigma/r)^12 - (sigma/r)^6 ]
    # Note: We only integrate the TAIL, so r > sigma always.

    ratio = SIG_ANG / r_ang
    term6 = ratio ** 6
    term12 = ratio ** 12

    V_r = 4.0 * EPS_J * (term12 - term6)

    return r_ang * V_r


# =============================================================================
# 3. CROSS SECTION CALCULATOR
# =============================================================================
def calc_lj_sigma(v_ms):
    # 1. Wavenumber k in Angstrom^-1
    k_ang = (MU * v_ms / HBAR) * ANG

    # 2. Core Contribution (Geometric Hard Sphere)
    # Everything inside R_SWITCH is treated as a collision.
    sigma_core_ang2 = np.pi * (R_SWITCH_ANG ** 2)

    # 3. Born Tail Contribution (r > R_SWITCH)
    theta_grid = np.linspace(0.005, np.pi, 100)
    d_sigmas = []

    for theta in theta_grid:
        # Momentum transfer K in Angstrom^-1
        K_ang = 2.0 * k_ang * np.sin(theta / 2.0)

        # --- SCALED INTEGRATION ---
        # Integrate from R_SWITCH to Infinity in Angstrom-space.
        # This avoids the "nan" and "overflow" errors.

        # Split for stability
        val_near, _ = quad(lj_integrand_scaled, R_SWITCH_ANG, 100.0,
                           weight='sin', wvar=K_ang, limit=500)
        val_far, _ = quad(lj_integrand_scaled, 100.0, np.inf,
                          weight='sin', wvar=K_ang, limit=5000)

        integral_val = val_near + val_far  # Unit: J * Ang^2

        # Convert integral to SI for the prefactor math
        integral_SI = integral_val * (ANG ** 2)

        # Calculate Born Amplitude f in meters
        K_SI = K_ang / ANG
        prefactor_SI = -2.0 * MU / (HBAR ** 2 * K_SI)
        f_SI = prefactor_SI * integral_SI

        d_sigmas.append(np.abs(f_SI) ** 2)

    d_sigmas = np.array(d_sigmas)

    # Angular Integration for Momentum Transfer
    # Integral [ (1-cos)*DCS*sin ]
    integrand = (1.0 - np.cos(theta_grid)) * d_sigmas * np.sin(theta_grid)
    sigma_tail_SI = 2.0 * np.pi * trapezoid(integrand, theta_grid)

    # Final Sum in Angstroms^2
    return sigma_core_ang2 + (sigma_tail_SI / (ANG ** 2))


# =============================================================================
# 4. EXECUTION
# =============================================================================
if __name__ == "__main__":
    velocities = np.geomspace(10, 2000, 20)
    results = []

    print(f"Lennard-Jones Hybrid Cross Section (R_switch={R_SWITCH_ANG} A)")
    print("-" * 45)
    print(f"{'v (m/s)':<10} | {'Sigma (Ang^2)':<15}")

    for v in velocities:
        try:
            s = calc_lj_sigma(v)
            results.append(s)
            print(f"{v:<10.1f} | {s:<15.2f}")
        except Exception as e:
            print(f"{v:<10.1f} | Error: {e}")
            results.append(np.nan)

    # --- Plotting ---
    plt.figure(figsize=(8, 6))
    plt.plot(velocities, results, 'o-', linewidth=2, label='Hybrid Lennard-Jones')

    # Add Reference Lines
    # 1. Hard Sphere Floor
    floor = np.pi * R_SWITCH_ANG ** 2
    plt.axhline(floor, color='k', linestyle=':', label=f'Hard Sphere Limit ({floor:.0f} $\AA^2$)')

    # 2. Langevin Slope (v^-2/3 for LJ tail)
    if not np.isnan(results[0]):
        # Anchor the theoretical curve to the first calculated point
        ref_s = results[0]
        ref_v = velocities[0]
        # LJ tail (1/r^6) leads to sigma ~ v^-2/3
        plt.plot(velocities, ref_s * (velocities / ref_v) ** (-0.66), 'r--', label=r'Theoretical Tail $\sim v^{-2/3}$')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel(r'Cross Section ($\AA^2$)')
    plt.title('Corrected Lennard-Jones Cross Section')
    plt.legend()
    plt.grid(True, which='both', alpha=0.3)
    plt.show()
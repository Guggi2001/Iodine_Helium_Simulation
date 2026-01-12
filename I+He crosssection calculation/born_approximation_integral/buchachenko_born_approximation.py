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
CM_1 = 1.986e-23  # J

# System Masses
M_I = 126.90 * AMU
M_HE = 4.0026 * AMU
MU = (M_I * M_HE) / (M_I + M_HE)

# =============================================================================
# 2. BUCHACHENKO POTENTIAL (Units: Angstroms & Joules)
# =============================================================================
# Parameters adapted for Rescaled Integration (r in Angstroms)
PARAMS = {
    'A': 5.6e5 * CM_1,  # J
    'alpha': 3.6,  # Angstrom^-1 (Note: No 1/ANG factor here!)
    'C4': 1.5e4 * CM_1,  # J * Ang^4
    'C6': 8.0e5 * CM_1,  # J * Ang^6
    'delta': 1.0  # Angstrom^-1
}

# Cutoff parameters
R_SWITCH_ANG = 15.0  # We treat r < 15 A as a hard core


def buchachenko_integrand_scaled(r_ang):
    """
    Integrand for Born Approximation.
    Input r_ang is in Angstroms.
    Returns: r * V(r) * ScalingFactors
    Result Unit: [J * Angstrom]
    """
    # 1. Short Range (Born-Mayer)
    V_sr = PARAMS['A'] * np.exp(-PARAMS['alpha'] * r_ang)

    # 2. Damping
    f_damp = 0.5 * (1 + np.tanh(1 + PARAMS['delta'] * r_ang))

    # 3. Long Range (Induction + Dispersion)
    # Note: Since C4, C6 are in units of Ang^n, we just divide by r_ang^n
    V_lr = -PARAMS['C4'] / (r_ang ** 4) - PARAMS['C6'] / (r_ang ** 6)

    V_total = V_sr + f_damp * V_lr

    # Return r * V.
    # The sin(Kr) part is handled by the weight function.
    return r_ang * V_total


# =============================================================================
# 3. CROSS SECTION CALCULATOR
# =============================================================================
def calc_sigma_scaled(v_ms):
    # Wavenumber k in Angstrom^-1
    # k_SI = mu*v/hbar [m^-1] -> k_Ang = k_SI * 1e-10
    k_ang = (MU * v_ms / HBAR) * ANG

    # --- A. Core Contribution (Geometric) ---
    # Sigma = pi * R^2 [Ang^2]
    sigma_core_ang2 = np.pi * (R_SWITCH_ANG ** 2)

    # --- B. Born Tail Contribution ---
    theta_grid = np.linspace(0.005, np.pi, 100)
    d_sigmas = []

    for theta in theta_grid:
        # Momentum transfer K in Angstrom^-1
        K_ang = 2.0 * k_ang * np.sin(theta / 2.0)

        # Prefactor needs careful unit handling
        # Formula: f = - (2 mu / hbar^2 K) * Integral
        # We want f in Angstroms.
        # Let's verify units:
        # Integral = [J * Ang^2] (dr adds one Ang, integrand has one Ang)
        # Prefactor_SI = [kg / (J^2 s^2 m^-1)] ... complicated.

        # Let's compute Prefactor in strictly SI then convert to Angstroms
        # K_SI = K_ang / ANG
        # P_SI = -2 * MU / (HBAR**2 * (K_ang/ANG))
        # f_SI = P_SI * (Integral_val * ANG^2 * ANG) <--- NO.

        # SIMPLER WAY:
        # f_Born = - (2*mu / hbar^2) * (1/K) * Int(r V sin(Kr) dr)
        # Let's perform calculation in SI, but use the scaled integral result.

        # 1. Evaluate Integral (result is in J * Ang^2)
        #    Because integrand returns (J*Ang) and dr is (Ang)

        # Split integral for stability
        val_1, _ = quad(buchachenko_integrand_scaled, R_SWITCH_ANG, 100.0,
                        weight='sin', wvar=K_ang, limit=500)
        val_2, _ = quad(buchachenko_integrand_scaled, 100.0, np.inf,
                        weight='sin', wvar=K_ang, limit=5000)

        integral_scaled = val_1 + val_2  # Units: J * Ang^2

        # 2. Convert Integral to SI [J * m^2]
        integral_SI = integral_scaled * (ANG ** 2)

        # 3. Calculate K in SI
        K_SI = K_ang / ANG

        # 4. Calculate f in SI [m]
        prefactor_SI = -2.0 * MU / (HBAR ** 2 * K_SI)
        f_SI = prefactor_SI * integral_SI

        d_sigmas.append(np.abs(f_SI) ** 2)

    d_sigmas = np.array(d_sigmas)

    # Integrate over angles
    integrand = (1.0 - np.cos(theta_grid)) * d_sigmas * np.sin(theta_grid)
    sigma_tail_SI = 2.0 * np.pi * trapezoid(integrand, theta_grid)

    # Convert Total to Ang^2
    sigma_total_ang2 = sigma_core_ang2 + (sigma_tail_SI / (ANG ** 2))

    return sigma_total_ang2


# =============================================================================
# 4. EXECUTION
# =============================================================================
if __name__ == "__main__":
    velocities = np.geomspace(10, 2000, 20)
    results = []

    print(f"I+ - He Cross Section (Rescaled Integration)")
    print("-" * 40)
    print(f"{'v (m/s)':<10} | {'Sigma (Ang^2)':<15}")

    for v in velocities:
        try:
            s = calc_sigma_scaled(v)
            results.append(s)
            print(f"{v:<10.1f} | {s:<15.2f}")
        except Exception as e:
            print(f"{v:<10.1f} | Error: {e}")
            results.append(np.nan)

    # Theoretical Langevin for comparison
    # Sigma_L = pi * (4 * C4 / (mu*v^2))^(1/2) * (Scaling factors...)
    # Actually, Langevin Sigma = 2 * pi * (C4_SI / E_kin)^1/2 is for capture
    # Standard formula: Sigma(v) = pi * (4 * C4_SI / (mu * v^2))^(1/2) is not quite right dimensionally for cross section
    # Let's just plot the 1/v trend line normalized to the first point

    plt.figure(figsize=(8, 6))
    plt.plot(velocities, results, 'b-o', linewidth=2, label='Buchachenko Hybrid')

    # Reference Slope (v^-1)
    if not np.isnan(results[0]):
        ref_s = results[0]
        ref_v = velocities[0]
        plt.plot(velocities, ref_s * (ref_v / velocities), 'r--', label='Langevin Slope ($v^{-1}$)')

    # Hard Sphere Floor
    floor = np.pi * R_SWITCH_ANG ** 2
    plt.axhline(floor, color='k', linestyle=':', label=f'Core Limit ({floor:.0f} $\AA^2$)')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Cross Section ($\AA^2$)')
    plt.title('Final Converged I$^+$–He Cross Section')
    plt.legend()
    plt.grid(True, which='both', alpha=0.3)
    plt.show()

    # Save
    df = pd.DataFrame({'v_ms': velocities, 'sigma_A2': results})
    df.to_csv('I_He_Final_Rescaled.csv', index=False)
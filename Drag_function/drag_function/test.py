import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Configuration ---
base_path = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT\18A"
output_file = os.path.join(base_path, "18A_Decimated_Unified.csv")
step = 100  # Take every 100th point

def marti_loader(path):
    """Robust loader for large Marti files with scientific notation."""
    return np.genfromtxt(path, comments='#', invalid_raise=False)

# --- 1. Load Original Data ---
print("Loading large Marti files...")
v1_raw = marti_loader(os.path.join(base_path, "vimp.I_1"))
v2_raw = marti_loader(os.path.join(base_path, "vimp.I_2"))
r1_raw = marti_loader(os.path.join(base_path, "rimp.I_1"))
r2_raw = marti_loader(os.path.join(base_path, "rimp.I_2"))

# --- 2. Calculate Magnitudes and Distance ---
# Calculate R (Distance) from positions
min_len_r = min(len(r1_raw), len(r2_raw))
t_master = r1_raw[:min_len_r, 0]
R_actual = np.linalg.norm(r1_raw[:min_len_r, 1:4] - r2_raw[:min_len_r, 1:4], axis=1)

# Calculate V1 Magnitude
v1_mag = np.linalg.norm(v1_raw[:min_len_r, 1:4], axis=1)

# Handle V2 (which is shorter)
# We create an array of NaNs for V2 that matches the master time length
v2_mag_full = np.full(len(t_master), np.nan)
v2_mag_orig = np.linalg.norm(v2_raw[:, 1:4], axis=1)
v2_mag_full[:len(v2_mag_orig)] = v2_mag_orig

# --- 3. Decimate (Sub-sampling) ---
t_final  = t_master[::step]
v1_final = v1_mag[::step]
v2_final = v2_mag_full[::step]
R_final  = R_actual[::step]

# --- 4. Export to Unified CSV ---
export_df = pd.DataFrame({
    'Time_ps': t_final,
    'V1_mag': v1_final,
    'V2_mag': v2_final,
    'R_distance': R_final
})

# Save to CSV
export_df.to_csv(output_file, index=False)
print(f"Export complete: {output_file}")
print(f"Rows reduced from {len(t_master)} to {len(t_final)}")

# --- 5. Verification Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t_final, v1_final, label='Atom 1 Speed (Full)')
ax1.plot(t_final, v2_final, label='Atom 2 Speed (Shortened)')
ax1.set_ylabel('v / Å/ps')
ax1.set_title('Decimated Velocity Data (Every 100th Point)')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(t_final, R_final, color='black', label='Ground Truth Distance')
ax2.set_ylabel('R / Å')
ax2.set_xlabel('t / ps')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
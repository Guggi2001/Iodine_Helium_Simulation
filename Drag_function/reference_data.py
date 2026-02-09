import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- Configuration ---
base_path = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT"
cases = ["9A", "18A"]


def marti_loader(path):
    """
    Uses genfromtxt to handle large files with potential
    formatting errors or corrupted lines.
    """
    try:
        # invalid_raise=False skips lines that don't match the column count (like your E+ error)
        # loose=True and invalid_raise=False together make it very hard to crash
        data = np.genfromtxt(path, comments='#', invalid_raise=False)
        return pd.DataFrame(data)
    except Exception as e:
        print(f"Error loading {os.path.basename(path)}: {e}")
        return pd.DataFrame()


for case in cases:
    case_path = os.path.join(base_path, case)
    print(f"\n--- Processing {case} ---")

    try:
        if case == "9A":
            # 9A Files: Ensure paths are correct based on your previous screenshots
            df_v1 = pd.read_csv(os.path.join(case_path, "data_vabs.csv"), header=None)
            df_v2 = pd.read_csv(os.path.join(case_path, "data_vabs2.csv"), header=None)
            df_R = pd.read_csv(os.path.join(case_path, "R1-R2.csv"), header=None)
        else:
            # 18A Files using the new robust loader
            df_v1 = marti_loader(os.path.join(case_path, "vimp.I_1"))
            df_v2 = marti_loader(os.path.join(case_path, "vimp.I_2"))
            df_r1_pos = marti_loader(os.path.join(case_path, "rimp.I_1"))
            df_r2_pos = marti_loader(os.path.join(case_path, "rimp.I_2"))

            if not df_r1_pos.empty and not df_r2_pos.empty:
                min_len = min(len(df_r1_pos), len(df_r2_pos))
                p1 = df_r1_pos.iloc[:min_len, 1:4].values
                p2 = df_r2_pos.iloc[:min_len, 1:4].values
                dist = np.linalg.norm(p1 - p2, axis=1)
                df_R = pd.DataFrame({0: df_r1_pos.iloc[:min_len, 0], 1: dist})
            else:
                df_R = pd.DataFrame()

        if df_v1.empty or df_v2.empty or df_R.empty:
            print(f"Data missing for {case}. Check filenames in: {case_path}")
            continue

        # --- Visualization ---
        fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle(f'HeDFT Reference Results: {case}', fontsize=16, fontweight='bold')

        # Use iloc to ensure we grab columns by index regardless of names
        axs[0, 0].plot(df_v1.iloc[:, 0], df_v1.iloc[:, 1].abs(), color='#0072BD', lw=1.5)
        axs[0, 0].set_title('Velocity: Iodine 1')
        axs[0, 0].set_ylabel('v / Å/ps')
        axs[0, 0].grid(True, alpha=0.3)

        axs[0, 1].plot(df_v2.iloc[:, 0], df_v2.iloc[:, 1].abs(), color='#D95319', lw=1.5)
        axs[0, 1].set_title('Velocity: Iodine 2')
        axs[0, 1].grid(True, alpha=0.3)

        gs = axs[1, 0].get_gridspec()
        for ax in axs[1, :]: ax.remove()
        ax_bottom = fig.add_subplot(gs[1, :])

        ax_bottom.plot(df_R.iloc[:, 0], df_R.iloc[:, 1], color='black', lw=2)
        ax_bottom.set_title('Internuclear Distance R')
        ax_bottom.set_xlabel('t / ps')
        ax_bottom.set_ylabel('R / Å')
        ax_bottom.grid(True, alpha=0.3)

    except Exception as e:
        print(f"Fatal error in {case}: {e}")

plt.show()
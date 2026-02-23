import pandas as pd
import os
import numpy as np


def load_data(file_path):
    """
    Safely loads our pre-processed CSVs.
    Resaves panda data type to a dictionary with useful variables for physics analysis.
        - Handles missing files and parsing errors gracefully.
    Returns:
        dict: A dictionary containing time, velocities, and reconstructed distance arrays.
    """
    data = None
    if not os.path.exists(file_path):
        print(f"Warning: File not found at {file_path}")
        return data
    try:
        data = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error loading 9A data: {e}")
    data_dict = {
        't':   data['Time_ps'].values,
        'v1': data['V1_mag'].values,
        'v2': data['V2_mag'].values,
        'R':   data['R_distance'].values
    } if data is not None else None
    return data_dict


def load_data_18(folder_path):
    """
    Loads raw velocity and position data for the 18A case.
    Returns:
        dict: A dictionary containing time, velocities, and positions arrays.
    """
    def marti_loader(filename):
        path = os.path.join(folder_path, filename)
        # Skips corrupted lines and handles scientific notation
        try:
            ret = np.genfromtxt(path, comments='#', invalid_raise=False)
        except Exception as e:
            print(f"Error loading file {filename}: {e}")
            ret = None
        return ret

    data_dict = {}
    v1_data = marti_loader("vimp.I_1")
    v2_data = marti_loader("vimp.I_2")
    r1_pos = marti_loader("rimp.I_1")
    r2_pos = marti_loader("rimp.I_2")

    t1, v1 = v1_data[:, 0], np.sqrt(np.sum(v1_data[:, 1:4] ** 2, axis=1))
    t2, v2 = v2_data[:, 0], np.sqrt(np.sum(v2_data[:, 1:4] ** 2, axis=1))

    # Calculate Distance R
    min_len = min(len(r1_pos), len(r2_pos))
    r_dist = np.sqrt(np.sum((r1_pos[:min_len, 1:4] - r2_pos[:min_len, 1:4]) ** 2, axis=1))
    data_dict['t'] = t1[::100]  # Decimation for performance (every 100th point)
    data_dict['v1'] = v1[::100]
    data_dict['v2'] = v2[::100]
    data_dict['R'] = r_dist[::100]
    return data_dict

# --- Testing if data loaded correctly ---
test = True
if test:
    from config_utils_local import config as C
    import matplotlib.pyplot as plt
    # Load 9A data
    dict9 = load_data(C.PATH9A)

    fig9, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig9.suptitle('HeDFT Reference Results: 9Å', fontsize=16, fontweight='bold')

    # Velocity Iodine 1
    axs[0, 0].plot(dict9['t'], dict9['v1'], color='#0072BD', lw=1.5)
    axs[0, 0].set_title('Velocity: Iodine 1 (9Å)')
    axs[0, 0].set_ylabel('v / Å/ps')
    axs[0, 0].set_xlim(-2, 12)
    axs[0, 0].grid(True, alpha=0.3)

    # Velocity Iodine 2
    axs[0, 1].plot(dict9['t'], dict9['v2'], color='#D95319', lw=1.5)
    axs[0, 1].set_title('Velocity: Iodine 2 (9Å)')
    axs[0, 1].set_xlim(-2, 12)
    axs[0, 1].grid(True, alpha=0.3)

    # Internuclear Distance (spanning bottom)
    gs = axs[1, 0].get_gridspec()
    for ax in axs[1, :]: ax.remove()
    ax_bottom = fig9.add_subplot(gs[1, :])
    ax_bottom.plot(dict9['t'], dict9['R'], color='black', lw=2)
    ax_bottom.set_title('9Å Case: Distance Verification & Extension')
    ax_bottom.set_xlabel('t / ps')
    ax_bottom.set_ylabel('R / Å')
    ax_bottom.set_xlim(0, 12)  # Matches your MATLAB screenshot zoom
    ax_bottom.grid(True, alpha=0.3)

    plt.show()

    # Plotting 18A
    dict18 = load_data(C.PATH18A)

    fig18, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig18.suptitle('HeDFT Reference Results: 18Å', fontsize=16, fontweight='bold')


    axs[0, 0].plot(dict18['t'], dict18['v1'], color='#0072BD', lw=1.2)
    axs[0, 0].set_title('Velocity: Iodine 1 (18Å)')
    axs[0, 0].set_ylabel('v / Å/ps')

    axs[0, 1].plot(dict18['t'][:len(dict18['v2'])], dict18['v2'], color='#D95319', lw=1.2)
    axs[0, 1].set_title('Velocity: Iodine 2 (18Å)')

    gs = axs[1, 0].get_gridspec()
    for ax in axs[1, :]: ax.remove()
    ax_bottom = fig18.add_subplot(gs[1, :])
    ax_bottom.plot(dict18['t'], dict18['R'], color='black', lw=1.5)
    ax_bottom.set_title('Internuclear Distance R (18Å)')
    ax_bottom.set_xlabel('t / ps')
    ax_bottom.set_ylabel('R / Å')

    for ax in fig18.axes: ax.grid(True, alpha=0.3)

    plt.show()
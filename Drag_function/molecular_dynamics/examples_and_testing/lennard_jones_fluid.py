import matplotlib.pyplot as plt
from lammps import lammps
import numpy as np

# 1. Initialize LAMMPS and force it to write a specific log file
log_filename = "log.minimize"
lmp = lammps(cmdargs=["-log", log_filename, "-screen", "none"])

setup_cmds = [
    # 1) Initialization
    "units lj",
    "dimension 3",
    "atom_style atomic",
    "boundary p p p",
    # 2) System definition
    "region simbox block -20 20 -20 20 -20 20",
    "create_box 2 simbox",
    "create_atoms 1 random 1500 34134 simbox overlap 0.3",
    "create_atoms 2 random 100 12756 simbox overlap 0.3",
    # 3) Settings
    "mass 1 1.0",
    "mass 2 5.0",
    "pair_style lj/cut 4.0",
    "pair_coeff 1 1 1.0 1.0",
    "pair_coeff 2 2 0.5 3.0",
    # 4) Monitoring
    "thermo 10",
    "thermo_style custom step step ke pe",
    # 5) Run
    "minimize 1.0e-6 1.0e-6 1000 10000",
]

for cmd in setup_cmds:
    lmp.command(cmd)

natoms = lmp.get_natoms()


lmp.finalize()

# 3. Parse the Log File
steps, ke_data, pe_data = [], [], []
natoms = 1600  # From the 1500 + 100 atoms created

print("Parsing log file for visualization...")
with open(log_filename, "r") as f:
    in_thermo = False
    step_idx = ke_idx = pe_idx = 0

    for line in f:
        parts = line.split()
        if not parts:
            continue  # Skip empty lines

        # Detect the thermo header
        if parts[0] == "Step" and "KinEng" in parts and "PotEng" in parts:
            in_thermo = True
            step_idx = parts.index("Step")
            ke_idx = parts.index("KinEng")
            pe_idx = parts.index("PotEng")
            continue

        # Detect the end of the thermo block
        if in_thermo and (parts[0] == "Loop" or parts[0] == "Minimization"):
            in_thermo = False
            continue

        # Parse the numeric data
        if in_thermo:
            try:
                # Convert total energies to per-atom energies
                steps.append(int(parts[step_idx]))
                ke_data.append(float(parts[ke_idx]))
                pe_data.append(float(parts[pe_idx]))
            except ValueError:
                # Ignore lines that aren't strictly numeric (e.g., LAMMPS warnings)
                pass

print(f"Successfully extracted {len(steps)} data points.")

steps = np.array(steps)
ke_data = np.array(ke_data)
pe_data = np.array(pe_data)

# 4. Plot the Data
plt.figure(figsize=(12, 4))

# Kinetic Energy Plot (Will be a flat line at 0 for standard minimization)
plt.subplot(1, 2, 1)
plt.plot(steps, ke_data, color='red', linewidth=1.5)
plt.title('Minimization: Kinetic Energy (Always 0)')
plt.xlabel('Minimization Step')
plt.ylabel('KE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Potential Energy Plot (Will show the sharp descent)
plt.subplot(1, 2, 2)
plt.ylim(np.min(pe_data) * 1.1, 1)
plt.plot(steps, pe_data, color='blue', linewidth=1.5)
plt.title('Minimization: Potential Energy Descent')
plt.xlabel('Minimization Step')
plt.ylabel('PE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
# plt.savefig("minimization_evolution.png", dpi=300)
# print("Success: Plot saved as minimization_evolution.png")
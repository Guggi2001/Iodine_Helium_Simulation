from lammps import lammps
import lammps_logfile
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np

# Initialize LAMMPS instance
log_filename = "log.lj"
lmp = lammps(cmdargs=["-log", log_filename, "-screen", "none"])
# Use mpi4py to get the rank, ensuring MPI is initialized correctly
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"LAMMPS Version Date Code: {lmp.version()}")


lj_script = """
# 3D Lennard-Jones melt
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
# 2) System definition
region simbox block -20 20 -20 20 -20 20
create_box 2 simbox
create_atoms 1 random 1500 34134 simbox overlap 0.3
create_atoms 2 random 100 12756 simbox overlap 0.3
# 3) Settings
mass 1 1.0
mass 2 5.0
pair_style lj/cut 4.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
# 4) Monitoring
thermo 10
thermo_style custom step ke pe
# 5) Run

minimize 1.0e-6 1.0e-6 1000 10000

# PART B - MOLECULAR DYNAMICS
# 4) Monitoring
thermo 50
thermo_style custom step temp etotal pe ke press

# 5) Run
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 10917
timestep 0.005
run 15000

"""

# Execute the commands line-by-line
lmp.commands_string(lj_script)

# lmp.finalize() explicitly calls MPI_Finalize() before shutting down LAMMPS
lmp.finalize()


#--------------------------------------------------------------------------------------------------------------
#                                   Post processing
#--------------------------------------------------------------------------------------------------------------

# 3. Parse the Log File
steps, ke_data, pe_data = [], [], []
natoms = lmp.get_natoms()

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

# Reads the log file. If you have multiple 'run' or 'minimize' commands,
# it stores them as separate blocks.
log = lammps_logfile.File("log.lammps")

# Get data from the first run/minimize block
step = log.get("Step", run_num=0)
pe = log.get("PotEng", run_num=0)
ke = log.get("KinEng", run_num=0)


# 4. Plot the Data
plt.figure(figsize=(12, 4))

# Kinetic Energy Plot (Will be a flat line at 0 for standard minimization)
plt.subplot(1, 2, 1)
plt.plot(step, ke, color='red', linewidth=1.5)
plt.title('Minimization: Kinetic Energy (Always 0)')
plt.xlabel('Minimization Step')
plt.ylabel('KE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Potential Energy Plot (Will show the sharp descent)
plt.subplot(1, 2, 2)
plt.ylim(np.min(pe) * 1.1, 1)
plt.plot(step, pe, color='blue', linewidth=1.5)
plt.title('Minimization: Potential Energy Descent')
plt.xlabel('Minimization Step')
plt.ylabel('PE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt
from lammps import lammps

# Initialize LAMMPS silently
lmp = lammps(cmdargs=["-log", "none", "-screen", "none"])

# 1. System Setup (Matching Tutorial 1 parameters)
setup_cmds = [
    "units lj",
    "dimension 3",
    "atom_style atomic",
    "boundary p p p",
    "region simbox block -20 20 -20 20 -20 20",
    "create_box 2 simbox",
    "create_atoms 1 random 1500 34134 simbox overlap 0.3",
    "create_atoms 2 random 100 12756 simbox overlap 0.3",
    "mass 1 1.0",
    "mass 2 5.0",
    "pair_style lj/cut 4.0",
    "pair_coeff 1 1 1.0 1.0",
    "pair_coeff 2 2 0.5 3.0"
]

for cmd in setup_cmds:
    lmp.command(cmd)

# 2. Conjugate Gradient Minimization
# This resolves overlaps from the random insertion. KE remains 0 here.
print("Running energy minimization...")
lmp.command("minimize 1.0e-6 1.0e-6 1000 10000")

# 3. Equilibration Setup
lmp.command("reset_timestep 0")
lmp.command("timestep 0.005")
lmp.command("fix mynve all nve")
# The Langevin thermostat pumps KE into the system to reach T=1.0
lmp.command("fix mylgv all langevin 1.0 1.0 0.1 10917")

lmp.command("compute pe_total all pe")
lmp.command("compute ke_total all ke")

# 4. Step-by-Step Execution
steps, pe_data, ke_data = [], [], []
natoms = lmp.get_natoms()

print("Running MD equilibration to extract plot data...")
for i in range(150):  # Run 150 chunks of 10 steps (1500 total steps)
    lmp.command("run 10")

    pe = lmp.extract_compute("pe_total", 0, 0) / natoms
    ke = lmp.extract_compute("ke_total", 0, 0) / natoms

    steps.append((i + 1) * 10)
    pe_data.append(pe)
    ke_data.append(ke)

lmp.close()

# 5. Plotting Subpictures A and C
plt.figure(figsize=(12, 4))

# Subpicture A: Kinetic Energy
plt.subplot(1, 2, 1)
plt.plot(steps, ke_data, color='red', linewidth=1.5)
plt.title('Kinetic Energy vs. Timestep')
plt.xlabel('Timestep')
plt.ylabel('KE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Subpicture C: Potential Energy
plt.subplot(1, 2, 2)
plt.plot(steps, pe_data, color='blue', linewidth=1.5)
plt.title('Potential Energy vs. Timestep')
plt.xlabel('Timestep')
plt.ylabel('PE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
# plt.savefig("tutorial_energy_evolution.png", dpi=300)
# print("Success: Plot saved as tutorial_energy_evolution.png")
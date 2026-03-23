import matplotlib.pyplot as plt
from lammps import lammps

lmp = lammps()

# Define the NVT setup with tail corrections
commands = [
    "units lj",
    "atom_style atomic",
    "lattice fcc 0.8442",
    "region box block 0 10 0 10 0 10",
    "create_box 1 box",
    "create_atoms 1 box",
    "mass 1 1.0",
    "velocity all create 0.8 87287 loop geom",
    "pair_style lj/cut 2.5",
    "pair_coeff 1 1 1.0 1.0 2.5",

    # CRITICAL: Add long-range tail corrections to match NIST benchmark energies
    "pair_modify tail yes",

    "neighbor 0.3 bin",
    "neigh_modify delay 0 every 20 check no",

    # CRITICAL: Use NVT thermostat to force the system to exactly T* = 0.8
    "fix 1 all nvt temp 0.8 0.8 0.1",

    "compute t all temp",
    "compute p all pressure t",
    "compute pe all pe",
    "compute ke all ke",

    "thermo 100",

    # Equilibrate the system first before collecting our plot data
    "run 2000",
    "reset_timestep 0"
]

for cmd in commands:
    lmp.command(cmd)

# Data extraction (Production phase)
steps, temps, energies, pressures = [], [], [], []
natoms = lmp.get_natoms()

for i in range(10):
    lmp.command("run 100")

    t = lmp.extract_compute("t", 0, 0)
    pe = lmp.extract_compute("pe", 0, 0)
    ke = lmp.extract_compute("ke", 0, 0)
    p = lmp.extract_compute("p", 0, 0)

    steps.append((i + 1) * 100)
    temps.append(t)
    energies.append((pe + ke) / natoms)
    pressures.append(p)

me = lmp.extract_setting("world_rank")
if me == 0:
    # Plotting all three metrics to verify against benchmarks
    plt.figure(figsize=(15, 4))

    plt.subplot(1, 3, 1)
    plt.plot(steps, temps, marker='o', color='red')
    plt.title('Temperature (Target: ~0.8)')
    plt.ylabel('T*')

    plt.subplot(1, 3, 2)
    plt.plot(steps, energies, marker='s', color='blue')
    plt.title('Total Energy/Atom (Target: -4.4 to -4.5)')
    plt.ylabel('E*')

    plt.subplot(1, 3, 3)
    plt.plot(steps, pressures, marker='^', color='green')
    plt.title('Pressure (Target: ~5.5 to 6.0)')
    plt.ylabel('P*')

    plt.tight_layout()
    plt.savefig("benchmark_plot.png")

    print("--- Benchmark Validation ---")
    print(f"Final Temperature: {temps[-1]:.4f}")
    print(f"Final Energy/Atom: {energies[-1]:.4f}")
    print(f"Final Pressure:    {pressures[-1]:.4f}")

lmp.close()
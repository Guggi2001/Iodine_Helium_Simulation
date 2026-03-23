from lammps import lammps
from mpi4py import MPI


# Initialize LAMMPS instance
lmp = lammps()
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
thermo_style custom step etotal press
# 5) Run

minimize 1.0e-6 1.0e-6 1000 10000

# PART B - MOLECULAR DYNAMICS
# 4) Monitoring
thermo 50
thermo_style custom step temp etotal pe ke press
# 5) Run
fix mynve all nve
timestep 0.005
run 50000

"""

# Execute the commands line-by-line
lmp.commands_string(lj_script)

# lmp.finalize() explicitly calls MPI_Finalize() before shutting down LAMMPS
lmp.finalize()
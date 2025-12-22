from scipy import constants
import numpy as np


# Definitions
e_charge = -1.602e-19; # C
epsilon_0 = 8.85418781e-12; #farads  / meter

eV_per_wavenumber = 1/8065.54429

def coulomb_energy(x):
    """
    Coulomb potential energy [J]
    x : distance in Å
    """
    return e_charge**2 / (4 * np.pi * epsilon_0 * x * 1e-10)

def coulomb_velocity(x, m):
    """
    Velocity derived from Coulomb energy
    x : distance in Å
    m : mass in kg
    """
    return np.sqrt(coulomb_energy(x) / m)

u = 1.66053907e-27 # kg
eV = 1.602e-19 # joule
hc = 1240 # eV nm


# hard sphere collision parameters
#  mean free path method
bulk_density_helium = 0.0219 # Angström ^-3
density_droplet = 0.8*bulk_density_helium
# https://journals.aps.org/prb/pdf/10.1103/PhysRevB.58.3341


k_B = 1.380649E-23 # J / K

hbar = constants.hbar  # J*s
amu  = constants.physical_constants["atomic mass constant"][0]  # kg
eV   = constants.electron_volt  # J
angstrom = 1e-10  # m

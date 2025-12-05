global R0_GS
R0_GS = 2.666;
%R0_GS = 12;

global sigma_dependent_on_v
sigma_dependent_on_v = true;

global single_pulse
single_pulse = true;

global partner_interaction
partner_interaction = true; % iodine atoms interact in the ground state X potential

global additional_droplet_charges
additional_droplet_charges =0;

global highly_charged_iodine
highly_charged_iodine = false;

global num_molecules
num_molecules = 5000;

global effusive_dynamics
effusive_dynamics = false;

global hard_sphere_collision_mode
hard_sphere_collision_mode = 3; % 1: scattering rate, 2: scatter after traveling mean free path, 3: scatter based on probability sigma*dR *  rho_scatterer

% parameter for mode 1
global scattering_probability
scattering_probability = 0.0040;

% parameters for mode 2
global geometric_scattering_crosssection_I;
geometric_scattering_crosssection_I = 30; %Angström^2;

global scatter_mass_neutral
scatter_mass_neutral = 4; 
%scatter_mass_neutral = 16;

global scatter_mass_ion
scatter_mass_ion = 4;

global geometric_scattering_crosssection_Iplus;
%geometric_scattering_crosssection_Iplus = 10; %Angström^2; current best
geometric_scattering_crosssection_Iplus = 12; %Angström^2;

global binding_energy_I_ion
binding_energy_I_ion = 0.0;

global neutral_scatter_angle_std
neutral_scatter_angle_std = 40;

global ion_scatter_angle_std
ion_scatter_angle_std = 0;
%ion_scatter_angle_std = 1;

global mass_attach_probability % probability to attach a helium atom to an ion at each hard sphere collision
mass_attach_probability = 0.1;

global single_charge_ionization_allowed
single_charge_ionization_allowed = false;

global use_single_droplet_size 
global single_droplet_size
use_single_droplet_size = false;
single_droplet_size = 12800;
global p_source
global T_source

p_source = 40;
T_source = 14; 
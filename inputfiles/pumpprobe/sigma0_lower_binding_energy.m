global lambda_pump
global E_diss

lambda_pump = 634; % nm
lambda_pump = 700/2;
lambda_pump = 630/2;
E_diss = 12500*eV_per_wavenumber*eV;
E_diss = 1.556*eV;
E_diss = 2.526*eV;

global R0_GS
R0_GS = 2.666;

global E_coulomb_scale
E_coulomb_scale = 0.8;

global custom_DFT_start;
custom_DFT_start = true;

global single_initial_position
single_initial_position = false;

global T_particles
T_particles = 0.4;

global deltaR0
deltaR0 = 0.0;

global sigma_dependent_on_v
sigma_dependent_on_v = true;

global single_pulse
single_pulse = false;

global partner_interaction
partner_interaction = true; % iodine atoms interact in the ground state X potential

global additional_droplet_charges
% for uniform charge:
%additional_droplet_charges =10; % too much 
%additional_droplet_charges = 0.8;

additional_droplet_charges = 0;



global droplet_charge_model
droplet_charge_model = 2; % 1: localized charges at random positions, 2: uniformly charged sphere with droplet radius



global highly_charged_iodine
highly_charged_iodine = false;

global num_molecules
num_molecules = 1000;

global effusive_dynamics
effusive_dynamics = false;

global hard_sphere_collision_mode
hard_sphere_collision_mode = 3; % 1: scattering rate, 2: scatter after traveling mean free path

% parameter for mode 1
global scattering_probability
scattering_probability = 0.0040;

% parameters for mode 2
global geometric_scattering_crosssection_I;
geometric_scattering_crosssection_I = 60; %Angström^2;

global scatter_mass_neutral
scatter_mass_neutral = 4; 
%scatter_mass_neutral = 16;

global scatter_mass_ion
scatter_mass_ion = 4;

global binding_energy_I_ion
binding_energy_I_ion = 0.0;

global geometric_scattering_crosssection_Iplus;
%geometric_scattering_crosssection_Iplus = 10; %Angström^2; current best
geometric_scattering_crosssection_Iplus = 1900; %Angström^2; % ion single pulse
geometric_scattering_crosssection_Iplus = 800; %Angström^2; % ion pumpprobe
geometric_scattering_crosssection_Iplus = 1000; %Angström^2; % ion pumpprobe sigma = v^-1.7
geometric_scattering_crosssection_Iplus = 2000;  %Angström^2; % ion pumpprobe sigma = v^-2

%geometric_scattering_crosssection_Iplus = 2500*0.3;  %Angström^2; % ion pumpprobe sigma = v^-2

global neutral_scatter_angle_std
neutral_scatter_angle_std = 0;

global ion_scatter_angle_std
ion_scatter_angle_std = 0;

global mass_attach_probability;
%mass_attach_probability = 0.02;
mass_attach_probability = 0.00;
global single_charge_ionization_allowed
single_charge_ionization_allowed = false;

global use_single_droplet_size 
global single_droplet_size
use_single_droplet_size = true;
single_droplet_size = 12800;
global p_source
global T_source

p_source = 40;
T_source = 14; 
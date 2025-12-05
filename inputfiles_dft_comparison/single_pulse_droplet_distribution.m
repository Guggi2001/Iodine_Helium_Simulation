global R0_GS
R0_GS = 2.666;
%R0_GS = 9;

global E_coulomb_scale 
E_coulomb_scale = 0.8;

global custom_DFT_start;
custom_DFT_start = false;


global single_initial_position 
single_initial_position = false; % if true, all molecules start at the center of the droplet

global deltaR0 %width of initial interatomic distance
deltaR0 = 0;

global T_particles % this is equivalent to the translational temperature of particles moving inside the droplet (mv^2/2 = 3/2 k T)
% it is used to sample the initial positions of the particles in the
% droplet
T_particles = 0.4; % in K
%T_particles = 12.7; % this is the temperature equivalent to iodine atom moving at landau velocity (50 m/s)
% 2/3 * 127 u * (50 m/s)^2/2 / (boltzmann constant)


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
num_molecules = 8000;

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

geometric_scattering_crosssection_Iplus = 2500; %Angström^2; % for sigme = v^-2, from 9 Angström comparison

%geometric_scattering_crosssection_Iplus = 1600; %Angström^2; % for sigme = v^-2, from  18 Angström comparison
%geometric_scattering_crosssection_Iplus = 5000; %Angström^2; % for sigme = v^-2, from  18 Angström comparison

%geometric_scattering_crosssection_Iplus =100;

geometric_scattering_crosssection_Iplus = 500; % 18 Angström comparison, v^-1


global binding_energy_I_ion
%binding_energy_I_ion = 0.3; % , from 9 Angström comparison
binding_energy_I_ion = 0.05; % , from 18 Angström comparison,  v^-1
%binding_energy_I_ion = 0.00; % , from 18 Angström comparison

%binding_energy_I_ion = 0.01; % v^-1

global neutral_scatter_angle_std
neutral_scatter_angle_std = 0;

global ion_scatter_angle_std
ion_scatter_angle_std = 0;
%ion_scatter_angle_std = 1;

global mass_attach_probability % probability to attach a helium atom to an ion at each hard sphere collision
%mass_attach_probability = 0.1;
mass_attach_probability = 0.1;

mass_attach_probability = 0.09;
%mass_attach_probability = 0.4;
%mass_attach_probability = 0.005;
mass_attach_probability = 0.05; % 18 Angström, v^-1

global single_charge_ionization_allowed
single_charge_ionization_allowed = false;

global use_single_droplet_size 
global single_droplet_size
use_single_droplet_size = false;
%single_droplet_size = 2000;
%single_droplet_size =  2.3269e+04;

global p_source
global T_source

p_source = 40;
T_source = 14; 
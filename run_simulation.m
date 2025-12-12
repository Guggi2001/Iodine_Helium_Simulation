
%% define variables
close all 
clear all
clc


%setup_VMI_path_office_SSD  % establish VMI paths, needed for comparison to measurements
addpath 'T:\github synchronized\VMI_matlab'
setup_VMI_path_office_flir

physical_constants

global Xdip_active
Xdip_active = true; % include additional dip in ground state I2 potential


global DEBUG 
DEBUG =false;

cd 'T:\github synchronized\I2HeN_velocity_simulation'
addpath 'T:\github synchronized\I2HeN_velocity_simulation'
addpath additional_functions\
addpath plot_utility\ 
addpath plot_utility\colorcet
addpath 'result images'\
addpath supplementary_data\
addpath T:\github synchronized\I2HeN_velocity_simulation\inputfiles_dft_comparison


addpath 'T:\github synchronized\I2HeN_velocity_simulation\inputfiles'

set_groot_properties % graphics settings

global num_neutral_export_timestept  % number of timesteps to export from neutral dynamics
num_neutral_export_timestept = 30; % final quality

num_neutral_export_timestept =40; % trial and error

global v_limit % landau velocity
v_limit = 40; % in m/s


global sigma_ion_exponent
sigma_ion_exponent = -1;


% define these here, these variables are overwritten by the pumpprobe input
% files
global lambda_pump; lambda_pump = 630;
global E_diss; E_diss = 1.556*eV;


%% load settings from inputfiles

%% single pulse inputs (as used in thesis)

% He DFT comparison inputfiles
%run inputfiles_dft_comparison\single_pulse_N2000.m  

%run inputfiles_dft_comparison\single_pulse_droplet_distribution.m 

%run inputfiles_dft_comparison\single_pulse_gas_distribution.m

%% single pulse inputs for v^-1 tests
run inputfiles_dft_comparison\single_pulse_N2000_18Angst.m
%run inputfiles_dft_comparison\single_pulse_droplet_distribution_v_minus_one.m 

%% pumpprobe  inputs
%run inputfiles\gas_phase_pumpprobe.m


%run inputfiles\pumpprobe\HeDFT_mimic_neutral_and_ion_N2000_centerstart %
%num_neutral_export_timestept =200; % trial and error

%for comparison ofTDHeDFT and MD neutral dynamics

%run inputfiles\pumpprobe\sigma0.m  % show that the ions do not come out
%with usual single pulse parameters
 
%run inputfiles\pumpprobe\sigma0_lower_binding_energy.m %17.6.25 


%run inputfiles\pumpprobe\decreased_sigma0.m % sigme0 = 0.3*2500, e_bind_ion = 0.1

%run inputfiles\pumpprobe\sigma0_Xdip9Angst.m %22.4.25 % i cannot get this to work.. for testing currently..

% alternative, single pulse sim, to explain pumpprobe distribution at 200
% ps

%run inputfiles_dft_comparison\single_pulse_droplet_distribution_T_limit_9Angstroem.m
%run inputfiles_dft_comparison\single_pulse_droplet_distribution_T_limit_9Angstroem_zero_bind.m

%run inputfiles\pumpprobe\two_photon_excitation_Test.m % check dynamics with 1000 m/s initial velocity of neutral
%run inputfiles\pumpprobe\one_photon_excitation_Test.m

%test = 1;


%% modifiers for testing
% geometric_scattering_crosssection_Iplus = 2500*0.1;  %Angström^2; % ion pumpprobe sigma = v^-2 
% 
% binding_energy_I_ion = 0.0;
% mass_attach_probability = 0.01;
% additional_droplet_charges = 0;

%% override inputfile settings for testing
% use_single_droplet_size = true;
% single_droplet_size = 12800;
% single_initial_position = true; % if true, all molecules start at the center of the droplet





if single_pulse
cd single_pulse_simulation
end


%%
global sigma_dependent_on_v
sigma_dependent_on_v = false;
vmi_sim_3d_neutral_propa_HeDFT_mimic;


%vmi_sim_3d_neutral_propa_testing; % set v0 = 0 , same anisotropy as single pulse
%vmi_sim_post_process('neutral');





sigma_dependent_on_v = true;
vmi_sim_3d_ion_propa;



%%
global single_pulse
if single_pulse
      %  vmi_sim_post_process;

        %test = 1;
    post_process_single_pulse_paper_v3;
    %compare_dft_result;
else
    vmi_sim_post_process('pumpprobe');

end

return
figures = findall(groot,'Type','figure');

if use_single_droplet_size
    result_path = sprintf('N%.0f_sigN%.0f_sigI%.0f_N%.0f_He+_%.0f', num_molecules, geometric_scattering_crosssection_I, geometric_scattering_crosssection_Iplus, single_droplet_size,additional_droplet_charges);
else

result_path = sprintf('N%.0f_sigN%.0f_sigI%.0f_p%.0f_T%.0f_He+%.0f', num_molecules, geometric_scattering_crosssection_I, geometric_scattering_crosssection_Iplus, p_source, T_source,additional_droplet_charges);
end

mkdir("result images\"+result_path+"\");

for i=1:length(figures)
    savefig(figures(i), "result images\"+result_path+"\"+sprintf('fig%.0f', i))
end
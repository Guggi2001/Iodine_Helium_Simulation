
%% define variables
close all 
clear all
clc
setup_VMI_path_office_SSD
setup_VMI_path_office

global DEBUG
DEBUG =false;

cd 'T:\github synchronized\I2HeN_velocity_simulation'
addpath 'T:\github synchronized\I2HeN_velocity_simulation'
addpath additional_functions\
addpath plot_utility\ 
addpath plot_utility\colorcet
addpath 'result images'\
addpath supplementary_data\



set_groot_properties % graphics settings

run inputfiles/effusive_neutral

run inputfiles/droplet_neutral



%%
vmi_sim_3d_neutral_propa_HeDFT_mimic;



 vmi_sim_post_process;


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
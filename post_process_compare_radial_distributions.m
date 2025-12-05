% this script compares the radial histrogram over time to 
% the radial distribution (initially and at 60 K mean energy)
vmi_sim_post_process('pumpprobe')

%load('post_process_checkpoint.mat');
load('single_pulse_simulation/post_process_checkpoint.mat');

figure
for i=1:size(histogram_data_radius,1)
    distr = movmean(histogram_data_radius(i,:),20,2);
    plot(centers_radius, distr/max(distr)+i*0.5);

    hold on
end


recolor_lines

title('radius')

%% add p(r) r^2 for molecule r
eV = 1.602e-19; %joule
boltzmann_constant = 0.08617*1E-3*eV; % joule per kelvin

p_boltzmann = @(EE, TT, ZZ) exp(-EE./(boltzmann_constant*TT))/ZZ;
E = [0:0.001:E_max]/1000*eV; % energies in joule

normalization = trapz(E, p_boltzmann(E, T_droplet, 1)); 

p = @(E) p_boltzmann(E, T_droplet, normalization);

potential_steepness_molecule = 14.3324; % from fit of solvation potential DFT result
binding_energy_molecule = 26.9; % in meV for I2 in He from fit of solvation potential DFT result

r_max  =100;
r_step = 0.5;
r = 0:r_step:r_max;
normalization = trapz( r,  p(droplet_potential([potential_steepness_molecule, binding_energy_molecule, mean(droplet_radii)], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think

p_radius = @(r) p(droplet_potential([potential_steepness_molecule, binding_energy_molecule, mean(droplet_radii)],r)/1000*eV).*r.^2/normalization;

plot(r, p_radius(r)/max(p_radius(r)),'--')

%% 
eV = 1.602e-19; %joule
boltzmann_constant = 0.08617*1E-3*eV; % joule per kelvin

p_boltzmann = @(EE, TT, ZZ) exp(-EE./(boltzmann_constant*TT))/ZZ;
E = [0:0.001:E_max]/1000*eV; % energies in joule

normalization = trapz(E, p_boltzmann(E, 60, 1)); 

p = @(E) p_boltzmann(E,60, normalization);

potential_steepness_molecule = 14.3324; % from fit of solvation potential DFT result
binding_energy_molecule = 26.9; % in meV for I2 in He from fit of solvation potential DFT result

r_max  =100;
r_step = 0.5;
r = 0:r_step:r_max;
normalization = trapz( r,  p(droplet_potential([potential_steepness_molecule, binding_energy_molecule, mean(droplet_radii)], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think

p_radius = @(r) p(droplet_potential([potential_steepness_molecule, binding_energy_molecule, mean(droplet_radii)],r)/1000*eV).*r.^2/normalization;

plot(r, p_radius(r)/max(p_radius(r)),'-o')

%%
T_diss_atom = E_min*eV/boltzmann_constant;

normalization = trapz(E, p_boltzmann(E, T_diss_atom, 1)); 

p = @(E) p_boltzmann(E,T_diss_atom, normalization);

potential_steepness_atom = 14.2; % from fit of solvation potential DFT result


r_max  =100;
r_step = 0.5;
r = 0:r_step:r_max;
normalization = trapz( r,  p(droplet_potential([potential_steepness_atom, binding_energy_I_atom*1000, mean(droplet_radii)], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think

p_radius = @(r) p(droplet_potential([potential_steepness_atom, binding_energy_I_atom*1000, mean(droplet_radii)],r)/1000*eV).*r.^2/normalization;

plot(r, p_radius(r)/max(p_radius(r)),'-.')

vline(mean(droplet_radii));




%%
figure


for i=1:size(histogram_data_radius,1)
    distr = movmean(histogram_data_interatomic_distance(i,:),20,2);
    plot(centers_radius, distr/max(distr)+i*0.5);

    hold on
end


recolor_lines

title('interatomic distance')
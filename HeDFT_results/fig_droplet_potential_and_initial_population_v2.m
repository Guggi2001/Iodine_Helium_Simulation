addpath('T:\github synchronized\I2HeN_velocity_simulation\additional_functions')
 addpath 'T:\github synchronized\I2HeN_velocity_simulation\HeDFT_results'
addpath 'T:\github synchronized\VMI_matlab\matfile_data_scripts'
addpath('T:\github synchronized\VMI_matlab');
addpath('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures');

run a_state_paper_plot_parameters.m



% load and plot solvation energy in K from HeDFT result
load('I2_contraint_solvation.mat');
E = E - min(E);% - 574.12;
k_B = 1.381E-23; % J/K


yyaxis right
%scatter(r, E)

leg = {};

leg{end+1} =  'He-DFT';


%% random sampling of molecule initial positions and orientations
E_max = 200; % in meV
T_droplet = 0.4; % in K
%T_droplet = 4;
r_step = 0.01;
debug_plot_sample_generation = true;

R_list =  [34.14, 43.5, 55.9, 68.3]
%[N, R] = get_dropletsize(40, 17);
for radius_index = 1:length(R_list)
    R = R_list(radius_index);

rho_bulk = 0.0218;
rho_droplet = rho_bulk*0.8;
%R = (3/(4*pi*rho_droplet))^(1/3) * N.^(1/3);

N = R^3*4*pi/3*rho_droplet;

droplet_radius = R;


eV = 1.602e-19; %joule
boltzmann_constant = 0.08617*1E-3*eV; % joule per kelvin

% from ernesto dft result beta = [14.3324   26.9916   34.4431]
potential_steepness_molecule = 14.3324; % from fit of solvation potential DFT result
binding_energy_molecule = 574.12*boltzmann_constant/eV*1000; % in meV for I2 in He from fit of solvation potential DFT result


debug_plot_sample_generation = true;



% E = E*k_B; %energy in joule 
% 
% E = E/1.602E-19; % energy in eV
% E = E*1000; % energy in meV
 %16.4907   49.4918   40.2137
potential_steepness_molecule = 14.3324;
%binding_energy_molecule = 26.9916;

binding_energy_molecule =  573.3*boltzmann_constant/eV*1000;

% potential_steepness_molecule = 16.49;
% binding_energy_molecule = 49.4918 ;

steepness = potential_steepness_molecule;
binding_energy = binding_energy_molecule;


hold on
eV = 1.602e-19; %joule
boltzmann_constant = 0.08617*1E-3*eV; % joule per kelvin

p_boltzmann = @(EE, TT, ZZ) exp(-EE./(boltzmann_constant*TT))/ZZ;
E = [0:0.001:E_max]/1000*eV; % energies in joule

normalization = trapz(E, p_boltzmann(E, T_droplet, 1)); 

p = @(E) p_boltzmann(E, T_droplet, normalization);

cmap1 = colorcet('L06','N', 10);
cmap1(1:3,:) = [];


cmap2 = colorcet('L04','N', 10);
cmap2(1:3,:) = [];







    r_max = droplet_radius*2;

r = 0:r_step:r_max;
dr = 1E-6;
jacobian = @(r) (droplet_potential([steepness, binding_energy, droplet_radius], r+ dr) - droplet_potential([steepness, binding_energy, droplet_radius], r))/dr;


normalization = trapz( r,  p(droplet_potential([steepness, binding_energy, droplet_radius], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think
normalization2 = trapz( r,  jacobian(r).*p(droplet_potential([steepness, binding_energy, droplet_radius], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think


p_radius = @(r) p(droplet_potential([steepness, binding_energy, droplet_radius],r)/1000*eV).*r.^2/normalization;

%p_radius2 =  @(r) jacobian(r).*p(droplet_potential([steepness, binding_energy, droplet_radius],r)/1000*eV).*r.^2/normalization2;

y_max = max(p_radius(r));


r_accepted_total = [];






% 
% figure
% plot(E/eV*1000, p(E));
% xlim([0,30]);

%figure
yyaxis left
line = plot(r, p_radius(r), 'linestyle', ':', 'linewidth', 1.7, 'color', cmap1(radius_index,:));
line.Marker = 'none';
mean_r = trapz(r, r.*p_radius(r))./trapz(r, p_radius(r));
hold on
%vline(mean_r);

bracket_offset = -droplet_radius/700 + 0.05;

plot([mean_r, mean_r], 0.19 + 3*[-0.001, 0.001]+bracket_offset, 'HandleVisibility', 'off', 'color', cmap1(radius_index,:), 'linestyle', '-', 'Marker', 'none');
plot([mean_r, droplet_radius], 0.19 + [0, 0]+bracket_offset, 'HandleVisibility', 'off', 'color', cmap1(radius_index,:), 'linestyle', '-', 'Marker', 'none');
plot([droplet_radius, droplet_radius], 0.19 + 3*[-0.001, 0.001]+bracket_offset, 'HandleVisibility', 'off', 'color', cmap1(radius_index,:), 'linestyle', '-', 'Marker', 'none');

%hold on
%plot(r, p_radius2(r));

xlabel(['r/ ',char(197)])
ylabel('probabilitiy');

r_test = [0:r_step:100];
pot = droplet_potential([steepness, binding_energy, droplet_radius],r_test);


hold on
% bayes_hist(r_accepted_total, [0,r_max], false, 'red');
% title(sprintf('acceptance rate: %.2f\n', accept_rate));
test=1;

yyaxis right
ln = plot(r_test, pot/1000*eV/boltzmann_constant, 'linestyle', '-', 'linewidth',  1.4, 'color', cmap2(radius_index,:));
ln.Marker = 'none';
%hold on
%plot([min(r_test), max(r_test)], [0.6, 0.6]);

ylabel('T / K');
%legend('p(E(r)) r^2',  'E(r)')
%set(gca,'YScale', 'log');
set(gca,'Fontsize', 15);



end

yyaxis left
ylim([-0.00, 0.26])
%l = legend('p(E(r)) r^2',  'E(r) fit', 'He-DFT');
%l.Position = [0.5917    0.4295    0.2179    0.2114];
yyaxis right
%hline(0.6);
%hline(50.4);

a = gca;
a.XLim = [0.2513  104.0557];
a.YLim =   [-10.2503  614.8433];


global figsize
f = gcf;
f.Position(3:4) = figsize;

grid on
%exportgraphics(f, 'T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\solvation_potential\fig2_I2_solvation_potential.png', 'resolution', '300');
exportgraphics(f, 'T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\solvation_potential\fig2_I2_solvation_potential.pdf', 'ContentType','vector');


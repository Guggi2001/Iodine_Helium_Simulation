%% post processing
function vmi_sim_visualize_trajectory(input_crate)
desolve_struct(input_crate);

run physical_constants.m
% 
% global single_pulse
% if single_pulse
%   load('single_pulse_simulation\neutral_propagation_checkpoint');
%  %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);
% 
%  close all
%  tic;
% 
% load('single_pulse_simulation\ion_propagation_checkpoint.mat');
% else
% 
%   load('neutral_propagation_checkpoint');
%  %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);
% 
%  close all
%  tic;
% 
% load('ion_propagation_checkpoint.mat');
% end
% 
% 
% toc

figure

% required inputs
% mass
% mass = zeros(100,100);
% vx_components = mass;
% vy_components = mass;
% vz_components = mass;
% eV = 1;
% x_components = mass;
% y_components = mass;
% z_components = mass;
% num_molecules = 50;
% coulomb_energy = @(x)1;
% vx_total = mass;
% vy_total =mass;
% vz_total = mass;
% time = 1:length(mass);
% droplet_radii = time;
% b_ion_outside = time*0;
% % bayes hist function
% u = 1;
% E_kin = time;
% E_pot = time;
% L_droplet = mass;
% r0 = mass;
% E_initial = 1;
% E_d = 1;
% m = mass;
%effusive_dynamics = true;
% relative_energy_loss_per_ps = 1;
% binding_energy_I_atom = 1;

% vx_total = zeros(size(mass));
% vy_total = zeros(size(mass));
% vz_total = zeros(size(mass));
% b_ion_outside = zeros(size(mass));
% relative_ion_energy_loss_per_ps = 1;
% binding_energy_I_ion = 300;

E_kin_0 = sum(mass(:,end).*     (    (vx_components(:,end)*100).^2 +    (vy_components(:,end)*100).^2  + (vz_components(:,end)*100).^2  )/2)/eV;

diff_vector_x = x_components - circshift(x_components,num_molecules, 1);

diff_vector_y =  y_components - circshift( y_components,num_molecules,1);

diff_vector_z = z_components- circshift(z_components, num_molecules,1);

diff_magnitude = sqrt(diff_vector_x.^2 + diff_vector_y.^2 + diff_vector_z .^2);

diff_unit_x = diff_vector_x./diff_magnitude;
diff_unit_y = diff_vector_y./diff_magnitude;
diff_unit_z = diff_vector_z./diff_magnitude;

E_cei = sum(coulomb_energy(diff_magnitude(:,end)))/eV/2; % need to take half here to avoid double counting

% v_coulomb = coulomb_velocity(diff_magnitude, m);
% 
%  v_coulomb_x = v_coulomb/100.*diff_unit_x;
%  v_coulomb_y = v_coulomb/100.*diff_unit_y;
% 
% vx_total = vx_vectors + v_coulomb_x;
% vy_total = vy_vectors + v_coulomb_y;

E_kin_1 = sum(mass(:,end).*     (      (vx_total(:,end)*100).^2 + (vy_total(:,end)*100).^2 + (vz_total(:,end)*100).^2     )/2, 'omitnan')/eV;
defect = E_kin_1 - E_cei - E_kin_0;

disp([defect, E_kin_1, E_cei, E_kin_0, diff_magnitude(end)]);








close all

t_unique = unique(time);
edges_energy = [0:1:250]; % for energy in meV
edges_velocity = [0:0.05:22];
edges_radius = [0:1:200];

edges_mass = [126:1: 130];


histogram_data_radius = zeros(length(t_unique), length(edges_radius)-1);

histogram_data_velocity = zeros(length(t_unique), length(edges_velocity)-1);
histogram_data_total_velocity = zeros(length(t_unique), length(edges_velocity)-1);
histogram_data_kinetic_energy = zeros(length(t_unique), length(edges_energy)-1);
histogram_data_potential_energy = zeros(length(t_unique), length(edges_energy)-1);
histogram_data_total_energy =zeros(length(t_unique), length(edges_energy)-1);

histogram_data_mass = zeros(length(t_unique), length(edges_mass)-1);

if length(time)>1
t_id_interest = ceil(1/(time(2)-time(1)));
else
t_id_interest = 1;
end


v_total = sqrt(vx_total(:, t_id_interest).^2 + vy_total(:,t_id_interest).^2 + vz_total(:,t_id_interest).^2);
b_t_interest = v_total>9;


initial_depth = droplet_radii - sqrt(x_components(:,1).^2 + y_components(:, 1).^2 + z_components(:,1).^2);
depth = droplet_radii - sqrt(x_components(:,:).^2 + y_components(:, :).^2 + z_components(:,:).^2);

number_inside = sum(depth>0,1);
number_outside = sum(depth<0, 1);


%% selection of atoms to plot
t_id = length(time)


%plot_samples = sqrt(vx_vectors(:,b_t).^2 + vy_vectors(:,b_t).^2);
r_vectors = sqrt(x_components(:,t_id).^2 + y_components(:, t_id).^2 + z_components(:,t_id).^2);

b_inside = r_vectors<droplet_radii*1;


b_both_inside = b_inside & circshift(b_inside, num_molecules, 1);
b_both_outside = (~b_inside) & circshift(~b_inside, num_molecules,1);


b_one_outside = ~b_both_outside & ~b_inside;


b_recombined = sqrt((x_components(:, t_id) - circshift(x_components(:, t_id), num_molecules,1)).^2 ...
    + (y_components(:, t_id) - circshift(y_components(:, t_id), num_molecules,1)).^2 ...
    + (z_components(:, t_id) - circshift(z_components(:, t_id), num_molecules,1)).^2) < 3;


fraction_of_recombined = sum(b_recombined(1:num_molecules))/num_molecules;


v_total = sqrt(vx_total(:, t_id).^2 + vy_total(:,t_id).^2 + vz_total(:,t_id).^2);
v = sqrt(vx_components(:,t_id).^2 + vy_components(:,t_id).^2 + vz_components(:,t_id).^2);


% b_select = true(size(b_ion_outside(:, t_id)));
%b_select =
%b_select = b_ion_outside(:,t_id) & ~b_recombined & initial_depth>20;
%b_select = b_ion_outside(:,t_id) & r_vectors>droplet_radii;
%b_select = b_ion_outside(:,t_id) & r_vectors<=droplet_radii;

b_select = b_ion_outside(:,t_id) & r_vectors<droplet_radii;
b_select = b_ion_outside(:,t_id);
%b_select = true(size(b_select));

%b_select = b_ion_outside(:,t_id) & r_vectors>droplet_radii;



close all
u = 1.66053907e-27; %kg



try
b_both_ions_outside = b_ion_outside(:,end) & circshift(b_ion_outside(:,end), num_molecules)

view_both_ions_outside = false
if sum(b_both_ions_outside)>0 & view_both_ions_outside
    id_escape =  mod(find(b_both_ions_outside(:,end), num_molecules), num_molecules);
    id_escape = id_escape(1);
else
    id_escape = mod(find(b_ion_outside(:,end), num_molecules), num_molecules);
    if isempty(id_escape)
        error('no ions escaped');
    end

    %plot(id_escape, v_total(id_escape,:));

    id_escape = id_escape(4);
end

catch ME
    warning(ME.message);

    id_escape = 1;
end


%id_escape = 629;
%id_escape = id_escape(3);
%id_escape = 238;

counter = 0
while counter <5
id_escape = randi(size(x_components,1)/2, 1);
figure


plot3(x_components(id_escape,:), y_components(id_escape,:), z_components(id_escape,:), 'linewidth', 2);
hold on
plot3(x_components(id_escape+num_molecules,:), y_components(id_escape+num_molecules,:), z_components(id_escape+num_molecules,:), 'linewidth', 2);

if ~isnan(id_escape)
plot3(x_ci(id_escape,:), y_ci(id_escape,:), z_ci(id_escape,:),'LineStyle',':', 'linewidth', 2);
plot3(x_ci(id_escape+num_molecules,:), y_ci(id_escape+num_molecules,:), z_ci(id_escape+num_molecules,:), 'LineStyle',':', 'linewidth', 2);
end

legend({'neutral 1', 'neutral 2', 'ion 1', 'ion 2'});

xlim([-70,70]);
ylim([-70,70]);
zlim([-70,70]);

pbaspect([1,1,1]);

        [X,Y,Z] = sphere;
        X = X*mean(droplet_radii);
        Y = Y*mean(droplet_radii);
        Z = Z*mean(droplet_radii);

        s=  surface(X,Y,Z, 'FaceAlpha',0.1, 'FaceLighting','gouraud', 'EdgeColor',[0.1,0.3,0.9], 'EdgeAlpha',0.2, 'EdgeLighting','gouraud');
        s.CData = Z*0;
        colormap('winter')

        counter = counter + 1;
end

end

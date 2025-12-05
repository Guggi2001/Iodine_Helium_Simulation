function post_process_function(input_crate)

run physical_constants.m
desolve_struct(input_crate);


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









t_unique = unique(time);
edges_energy = [0:2:500]; % for energy in meV
edges_velocity = [0:0.05:22];
edges_radius = [0:0.2:120];
edges_interatomic_distance = edges_radius;

edges_mass = [126:1: 200];


histogram_data_radius = zeros(length(t_unique), length(edges_radius)-1);

histogram_data_velocity = zeros(length(t_unique), length(edges_velocity)-1);
histogram_data_total_velocity = zeros(length(t_unique), length(edges_velocity)-1);
histogram_data_kinetic_energy = zeros(length(t_unique), length(edges_energy)-1);
histogram_data_potential_energy = zeros(length(t_unique), length(edges_energy)-1);
histogram_data_total_energy =zeros(length(t_unique), length(edges_energy)-1);

histogram_data_mass = zeros(length(t_unique), length(edges_mass)-1);

histogram_data_interatomic_distance = zeros(length(t_unique), length(edges_interatomic_distance)-1);

%t_id_interest = ceil(1/(time(2)-time(1)));
%v_total = sqrt(vx_total(:, t_id_interest).^2 + vy_total(:,t_id_interest).^2 + vz_total(:,t_id_interest).^2);
%b_t_interest = v_total>9;


initial_depth = droplet_radii - sqrt(x_components(:,1).^2 + y_components(:, 1).^2 + z_components(:,1).^2);
depth = droplet_radii - sqrt(x_components(:,:).^2 + y_components(:, :).^2 + z_components(:,:).^2);

number_inside = sum(depth>0,1);
number_outside = sum(depth<0, 1);


% debug
    v_total = sqrt(vx_total(:, :).^2 + vy_total(:,:).^2 + vz_total(:,:).^2);

for t_id = 1:length(time)




    %plot_samples = sqrt(vx_vectors(:,b_t).^2 + vy_vectors(:,b_t).^2);
    r_vectors = sqrt(x_components(:,t_id).^2 + y_components(:, t_id).^2 + z_components(:,t_id).^2);
    if t_id ==1
        r_start = r_vectors;
    end
    
    r_interatomic = sqrt((x_components(:,t_id) - circshift(x_components(:,t_id), num_molecules)).^2 ...
        + (y_components(:,t_id) - circshift(y_components(:,t_id), num_molecules)).^2 + ...
        (z_components(:,t_id) - circshift(z_components(:,t_id), num_molecules)).^2);

    b_inside = r_vectors<droplet_radii*1;

    
    b_both_inside = b_inside & circshift(b_inside, num_molecules, 1);
    b_both_outside = (~b_inside) & circshift(~b_inside, num_molecules,1);


    b_one_outside = ~b_both_outside & ~b_inside;
    
    
    b_recombined = sqrt((x_components(:, t_id) - circshift(x_components(:, t_id), num_molecules,1)).^2 ...
        + (y_components(:, t_id) - circshift(y_components(:, t_id), num_molecules,1)).^2 ...
        + (z_components(:, t_id) - circshift(z_components(:, t_id), num_molecules,1)).^2) < 3;


    fraction_of_recombined = sum(b_recombined(1:num_molecules))/num_molecules;


    v_total = sqrt(vx_total(:, t_id).^2 + vy_total(:,t_id).^2 + vz_total(:,t_id).^2); % this is the full, unprojected velocity
    %v_total = sqrt(vx_total(:, t_id).^2 + vy_total(:,t_id).^2); % but in the experiment, we measure this quantity

    v = sqrt(vx_components(:,t_id).^2 + vy_components(:,t_id).^2 + vz_components(:,t_id).^2);
    %v = sqrt(vx_components(:,t_id).^2 + vy_components(:,t_id).^2);
    
global effusive_dynamics
    if effusive_dynamics
    b_select = true(size(b_ion_outside(:, t_id))) ;

    %b_select =
    %b_select = b_ion_outside(:,t_id) & ~b_recombined & initial_depth>20; 
    %b_select = b_ion_outside(:,t_id) & r_vectors>droplet_radii;
    %b_select = b_ion_outside(:,t_id) & r_vectors<=droplet_radii;

    %b_select = b_ion_outside(:,t_id) & r_vectors<droplet_radii*1.1;
   %b_select = b_ion_outside(:,t_id) & r_vectors<droplet_radii;
    %b_select = b_one_outside;
   %b_select = b_ion_outside(:,t_id);
      % b_select = b_ion_outside(:,t_id) & r_vectors>droplet_radii; 
    else
    b_select = b_ion_outside(:,t_id) & (~b_invalid_spawn);% & input_crate.b_v;
       % b_select = true(size(b_ion_outside(:, t_id))) & (~b_invalid_spawn) ;
   % b_select = true(size(b_select));
    
   %b_select = b_ion_outside(:,t_id) &  ( mass(:,t_id)/u>130 & mass(:,t_id)/u<132) ;

     %b_select = b_ion_outside(:,t_id);% &   mass(:,t_id)/u>315 ;
    end
    

    % velocity
    plot_samples = v;
    plot_samples = plot_samples(b_select);
    
    xinterval = [min(edges_velocity), max(edges_velocity)];
    [h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
        h = h-min(h);
    histogram_data_velocity(t_id, :) = h*numel(plot_samples);

    % total velocity + snowball formation and additional deceleration
    %helium_snowball_energy_loss = (randn(size(v_total))*100 + 50)/1000*eV; % in joule
    %     helium_snowball_energy_loss = randn(size(v_total))*0.1 + 0.4;
    %
    %     E_total = m* (v_total*100).^2/2;
    %     E_total = E_total.*(1 - helium_snowball_energy_loss).*b_inside;
    %     %b_negative = E_total<0;
    %
    %     v_total(b_inside) = sqrt(2*E_total(b_inside)/m)/100;
    
    % masses
    plot_samples = mass(b_select, t_id);
    xinterval = [min(edges_mass), max(edges_mass)];
    [h, sigma_h, centers_mass, ~, ~, ~] = bayes_hist(plot_samples/u, xinterval, true, 'r', edges_mass);
    h = h-min(h);

    histogram_data_mass(t_id, :) = h*numel(plot_samples);

    % total velocity
    %b_m = m_vectors(:, b_t)>127*u;
    %b_m = mass(:, t_id)== 127*u;
    %b_m = mass(:, t_id)== 127*u + 4*u*2;
    %b_m = abs( mass(:,end)/u - (127) )<1;
    %b_m = m_vectors(:,end)/u > 127 ;
    %b_m = abs( m_vectors(:,end)/u - (127) )<.1;  
    plot_samples = v_total;
    %b_v_total = v_total(:, 1/(time(2)-time(1));
    %plot_samples = v_total(b_m & b_both_outside);
    %plot_samples = v_total(b_inside);
    %plot_samples = plot_samples(~b_inside(b_m));
    plot_samples = plot_samples(b_select);

    xinterval = [min(edges_velocity), max(edges_velocity)];
    [h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
    h = h - min(h);

    histogram_data_total_velocity(t_id, :) = h*numel(plot_samples);
        h = h-min(h);


    % radius
    %plot_samples = L_droplet(~b_inside, b_t);

    plot_samples = r_vectors(:); 


    %plot_samples = depth(:,t_id);

    plot_samples = plot_samples(b_select);
    xinterval = [min(edges_radius), max(edges_radius)];
    [h, sigma_h, centers_radius, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_radius);
        h = h-min(h);
    histogram_data_radius(t_id, :) = h*numel(plot_samples);


    plot_samples = r_interatomic(:);
        xinterval = [min(edges_radius), max(edges_radius)];
    [h, sigma_h, centers_radius, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_interatomic_distance);
        h = h-min(h);
    histogram_data_interatomic_distance(t_id, :) = h*numel(plot_samples);

    % energy
    %plot_samples = m.*(100*v).^2/2/eV*1000;
    plot_samples = E_kin(:, t_id);
    plot_samples = plot_samples(b_select)*1000;

    xinterval = [min(edges_energy), max(edges_energy)];
    [h, sigma_h, centers_energy, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_energy);
        h = h-min(h);
    histogram_data_kinetic_energy(t_id, :) = h*numel(plot_samples);
    

    plot_samples = E_pot(b_select, t_id)*1000; % convert to meV
   

    xinterval = [min(edges_energy), max(edges_energy)];
    [h, sigma_h, centers_energy, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_energy);
        h = h-min(h);
    histogram_data_potential_energy(t_id, :) = h*numel(plot_samples);
    

    plot_samples = E_pot(:, t_id) + E_kin(:, t_id);
    plot_samples = plot_samples(b_select)*1000;

    xinterval = [min(edges_energy), max(edges_energy)];
    [h, sigma_h, centers_energy, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_energy);
        h = h-min(h);
    histogram_data_total_energy(t_id, :) = h*numel(plot_samples);

    % interatomic distance relative to droplet size
    plot_samples = diff_magnitude(b_select,t_id);
    xinterval = [min(edges_interatomic_distance), max(edges_interatomic_distance)];
    [h, sigma_h, centers_interatomic_distance, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_interatomic_distance);
        h = h-min(h);
    histogram_data_interatomic_distance(t_id,:) = h*numel(plot_samples);

%     if t == t_unique(1)
%         mean_E_start = mean(plot_samples);
% 
%     end
% 
%     if t== t_unique(end)
%         mean_E_end = mean(plot_samples);

%     end

end

save('post_process_checkpoint', 'L_droplet', 'vx_components',  'vy_components',  'vz_components', 'r_start', 'b_inside',...
    'v_total','vx_total', 'vy_total', 'vz_total', 'time','t_unique', 'E_kin', 'centers_radius', 'histogram_data_radius', 'centers_velocity', 'histogram_data_velocity',...
    'histogram_data_total_velocity','centers_mass','histogram_data_mass', 'number_inside', 'number_outside', 'depth', 'b_select','b_ion_outside', ...
    'centers_interatomic_distance', 'histogram_data_interatomic_distance', 'centers_energy', 'histogram_data_potential_energy', 'histogram_data_kinetic_energy');


end
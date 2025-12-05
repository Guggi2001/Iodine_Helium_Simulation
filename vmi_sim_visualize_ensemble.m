%% post processing
function vmi_sim_visualize_ensemble(input_crate)
desolve_struct(input_crate);
visualization_mode = 2
run physical_constants.m
% 

% 
% global single_pulse
% if single_pulse
%   load('single_pulse_simulation\neutral_propagation_checkpoint_single_pulse');
%  %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);
% 
%  close all
% 
% load('single_pulse_simulation\ion_propagation_checkpoint.mat');
% else
% 
%   load('neutral_propagation_checkpoint');
%  %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);
% 
%  close all
% 
% 
% load('ion_propagation_checkpoint.mat');
% end

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


t_id_interest = ceil(1/(time(2)-time(1)));
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
b_select = true(size(b_select));

%b_select = b_ion_outside(:,t_id) & r_vectors>droplet_radii;



close all
u = 1.66053907e-27; %kg

 v = VideoWriter('positions', 'MPEG-4');
 v.FrameRate = 60;
% 
 open(v)


if true
fprintf('generate density');




x_density = linspace(-80, 80, 400);
dx = x_density(2)- x_density(1);

y_density = x_density;

dy = y_density(2) - y_density(1);

density = zeros(length(x_density), length(y_density), length(time));
mass_density = zeros(size(density));

vx_density = linspace(-40,40, 400);
vy_density = linspace(-40,40,400);
dvx = vx_density(2) - vx_density(1);
dvy = vy_density(2) - vy_density(1);
v_density = zeros(length(vx_density), length(vy_density), length(time));
% 


molecule_select = b_select;



x_components = x_components(molecule_select, :);
y_components = y_components(molecule_select, :);

plot(mean(vx_components.^2 + vy_components.^2+ vz_components.^2, 2))

    %% method 1
    %     x = x_components(:,t==time);
    %     y = y_components(:,t==time);
    %
    %     for i=1:length(x_density)-1
    %         b_x = x> x_density(i) & x<x_density(i+1);
    %         for j=1:length(y_density)-1
    %             b_y = y>y_density(j) & y < y_density(j+1);
    %             density(i,j, t==time) = sum(b_x & b_y);
    %
    %         end
    %     end

    %% method 2
    for m = 1:length(time)

        for k = 1:size(x_components, 1)
            y = y_components(k,m);
            x = x_components(k,m);


            idx = floor((x -min(x_density))/dx)+1;
            idy = floor((y - min(y_density))/dy)+1;

            if idx>0 && idx<size(density,1) && idy>0 && idy<size(density,2)
                density(idx, idy, m) = density(idx, idy, m) + 1;
                mass_density(idx, idy, m)= mass_density(idx,idy, m) + (mass(k,m)/u - 127);
            end


            %vx = vx_components(k,m);
            %vy = vy_components(k,m);

            %vx = vx_total(k,m);
            %vy = vy_total(k,m);
            vx = vx_components(k,m);
            vy = vy_components(k,m);
            
            
            idvx = floor((vx -min(vx_density))/dvx)+1;
            idvy = floor((vy - min(vy_density))/dvy)+1;

            if idvx>0 && idvx<size(v_density,1) && idvy>0 && idvy<size(v_density,2)
                v_density(idvx, idvy, m) = v_density(idvx, idvy, m) + 1;
            end

        end
        disp(time(m));

    end




%save('vmi_sim_2d_output','-append', 'density', 'x_density', 'y_density');
end

switch visualization_mode
    case 1
%% visual particle velocity


id=1;


transformed_density = v_density;
%transformed_density = density;
transformed_density = transformed_density - min(transformed_density(:));
transformed_density = transformed_density / max(transformed_density(:));

video_time = 0;
dt = 0.02;

while video_time < time(end)

    id = get_closest_index(time, video_time);

    if id==1
        [xx, yy] = meshgrid(vx_density, vy_density);
        d = transformed_density(:,:, id);
        z = zeros(size(d));
        color_stack_0 = ones(size(d,1), size(d,2), 3);
% 
%         color_stack_0(:,:,1) = d;
%         color_stack_0(:,:,2) = z;
%         color_stack_0(:,:,3) = z;

        p1 = surf(xx, yy, zeros(size(xx)), color_stack_0, 'EdgeColor','none');
        shading interp
       % caxis([0,0.01*max(transformed_density(:))]);
        view(0,90);
        %colormap(viridis)
        hold on
        phi = 0:0.1:2*pi;
    r = mean(droplet_radii);

    xlabel(' v_x / Angström/ps');
    ylabel(' v_y / Angström/ps');

    x = r*cos(phi);
    y = r*sin(phi);
    p2 = plot(x,y, 'color', [0,0,0]);
    xlim([-25,25]);
    ylim([-25,25]);
    pbaspect([1,1,1]);
    set(gca, 'Fontsize', 15);
    end

    color_stack = color_stack_0;


    img = transformed_density(:,:,id);
    
    num_colors = 255;

%     cmap = circshift(colorcet('L03', 'N', num_colors), 0, 1);
%     red_channel = cmap( floor( img*(num_colors-1))+1,1);
%     red_channel = reshape(red_channel, size(img,1), size(img,2));
% 
%     green_channel = cmap( floor( img*(num_colors-1))+1,2);
%     green_channel = reshape(green_channel, size(img,1), size(img,2));
% 
%     blue_channel =  cmap( floor( img*(num_colors-1))+1,3);
%     blue_channel = reshape(blue_channel, size(img,1), size(img,2));
% 
%     b_c1 = img>0 & mass_density(:,:,id)<1;
%     color_stack(:,:,1) = color_stack(:,:,1) - b_c1;% + b_c1.*red_channel;
%     color_stack(:,:,2) = color_stack(:,:,2) ; %+ b_c1.*green_channel;
%     color_stack(:,:,3) = color_stack(:,:,3) ; %+ b_c1.*blue_channel;
%   

    cmap = winter(num_colors);
    red_channel = cmap( floor( img*(num_colors-1))+1,1);
    red_channel = reshape(red_channel, size(img,1), size(img,2));

    green_channel = cmap( floor( img*(num_colors-1))+1,2);
    green_channel = reshape(green_channel, size(img,1), size(img,2));

    blue_channel =  cmap( floor( img*(num_colors-1))+1,3);
    blue_channel = reshape(blue_channel, size(img,1), size(img,2));

    b_c2 = img>0;% & mass_density(:,:,id)>=1;
    color_stack(:,:,1) = color_stack(:,:,1) -b_c2; %+ b_c2.*red_channel;
    color_stack(:,:,2) = color_stack(:,:,2) - b_c2;% + b_c2.*green_channel;
    color_stack(:,:,3) = color_stack(:,:,3) - b_c2;% + b_c2.*blue_channel;


    p1.CData= color_stack;
    drawnow
    pause(0.005);
    
    

    f = getframe(gcf);

    writeVideo(v,f);


   
                %scatter(x, y);
    %pause((time(id+1) - time(id) )*0.1)
    title(sprintf('%.2f ps', video_time));

    video_time = video_time + dt;
        %video_time = video_time + time(id+1) - time(id);
end

 close(v);

    case 2
%% visualize atom positions
id=1;


transformed_density = density;
transformed_density = transformed_density - min(transformed_density(:));
transformed_density = transformed_density / max(transformed_density(:));

video_time = 0;
dt = 0.1;

while video_time < time(end)

    id = get_closest_index(time, video_time);

    if id==1
        [xx, yy] = meshgrid(x_density, y_density);
        d = transformed_density(:,:, id);
        z = zeros(size(d));
        color_stack_0 = ones(size(d,1), size(d,2), 3);
% 
%         color_stack_0(:,:,1) = d;
%         color_stack_0(:,:,2) = z;
%         color_stack_0(:,:,3) = z;

        p1 = surf(xx, yy, zeros(size(xx)), color_stack_0, 'EdgeColor','none');
        shading interp
       % caxis([0,0.01*max(transformed_density(:))]);
        view(0,90);
        %colormap(viridis)
        hold on
        phi = 0:0.1:2*pi;
    r = mean(droplet_radii);

    xlabel(' x / Angström');
    ylabel(' y / Angström');

    x = r*cos(phi);
    y = r*sin(phi);
    p2 = plot(x,y, 'color', [0,0,0]);
    xlim([-80,80]);
    ylim([-80,80]);
    pbaspect([1,1,1]);
    set(gca, 'Fontsize', 15);
    end

    color_stack = color_stack_0;


    img = transformed_density(:,:,id);
    
    num_colors = 255;
    cmap = circshift(colorcet('L03', 'N', num_colors), -10, 1);
    red_channel = cmap( floor( img*(num_colors-1))+1,1);
    red_channel = reshape(red_channel, size(img,1), size(img,2));

    green_channel = cmap( floor( img*(num_colors-1))+1,2);
    green_channel = reshape(green_channel, size(img,1), size(img,2));

    blue_channel =  cmap( floor( img*(num_colors-1))+1,3);
    blue_channel = reshape(blue_channel, size(img,1), size(img,2));

    b_c1 = img>0 & mass_density(:,:,id)<1;
    color_stack(:,:,1) = color_stack(:,:,1) - b_c1 + b_c1.*red_channel;
    color_stack(:,:,2) = color_stack(:,:,2) - b_c1 + b_c1.*green_channel;
    color_stack(:,:,3) = color_stack(:,:,3) - b_c1 + b_c1.*blue_channel;
  

    cmap = winter(num_colors);
    red_channel = cmap( floor( img*(num_colors-1))+1,1);
    red_channel = reshape(red_channel, size(img,1), size(img,2));

    green_channel = cmap( floor( img*(num_colors-1))+1,2);
    green_channel = reshape(green_channel, size(img,1), size(img,2));

    blue_channel =  cmap( floor( img*(num_colors-1))+1,3);
    blue_channel = reshape(blue_channel, size(img,1), size(img,2));

    b_c2 = img>0 & mass_density(:,:,id)>=1;
    color_stack(:,:,1) = color_stack(:,:,1) - b_c2 + b_c2.*red_channel;
    color_stack(:,:,2) = color_stack(:,:,2) - b_c2 + b_c2.*green_channel;
    color_stack(:,:,3) = color_stack(:,:,3) - b_c2 + b_c2.*blue_channel;


    p1.CData= color_stack;
    drawnow
    pause(0.005);
    


    f = getframe(gcf);
    writeVideo(v,f);
    
   
                %scatter(x, y);
    %pause((time(id+1) - time(id) )*0.1)
    title(sprintf('%.2f ps', video_time));

    %video_time = video_time + dt;
    video_time = video_time + time(id+1) - time(id);
end

 close(v);
end



%% velocity density
% id=1;
% time = t_vectors;
% 
% for t= time
% 
% 
%     if id==1
%         [xx, yy] = meshgrid(vx_density, vy_density);
%         p1 = surf(xx, yy, zeros(size(xx)), v_density(:,:, t==time), 'EdgeColor','none');
%         shading interp
%         %caxis([0,2]);
%         view(0,90);
% 
%     end
% 
%     p1.CData= v_density(:,:,t==time);
%     drawnow
%     pause(0.005);
% 
%     id = id+1;
% 
% end
% 












function vmi_sim_post_process(mode)

if nargin<1
    mode ='pumpprobe';
   mode = 'ion'
   % mode = 'neutral';

end

global effusive_dynamics
global use_single_droplet_size 
global single_pulse
global single_initial_position

if single_pulse
    neutral_data = load('single_pulse_simulation\neutral_propagation_checkpoint');
    %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);


       if use_single_droplet_size  && single_initial_position 
    ion_data = load('single_pulse_simulation\ion_propagation_checkpoint_hedft.mat');
       else
            ion_data = load('single_pulse_simulation\ion_propagation_checkpoint.mat');
       end
       
else

    neutral_data = load('neutral_propagation_checkpoint');
    %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_me3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);

    


        if use_single_droplet_size  && single_initial_position 
        fn = 'ion_propagation_checkpoint_hedft';
        ion_data = load(fn);
    else
        fn = 'ion_propagation_checkpoint';
        ion_data = load(fn);

            if effusive_dynamics
        ion_data =  load('ion_propagation_checkpoint_gas.mat');
    else
    ion_data =  load('ion_propagation_checkpoint.mat');
    end

    end

end

%% recollection data for switching between ion dynamics and neutral dynamics analysis
% this is not optimal, but necessary because i did not have the foresight
% to name the neutral and ion dynamics variables appropriately, so that
% post process can plot both
switch mode
    %% version for plotting neutral dynamics with ion simulation
    case 'pumpprobe'
        input_crate = struct();
        
        input_crate.b_invalid_spawn = neutral_data.b_invalid_spawn;
        
        % neutral dynamics data
        input_crate.mass =ion_data.mass_i_final;
        input_crate.num_molecules = neutral_data.num_molecules;

        input_crate.time = neutral_data.time;
        input_crate.droplet_radii = neutral_data.droplet_radii;

        input_crate.x_components = neutral_data.x_components;
        input_crate.y_components = neutral_data.y_components;
        input_crate.z_components = neutral_data.z_components;

        input_crate.vx_components = neutral_data.vx_components;
        input_crate.vy_components = neutral_data.vy_components;
        input_crate.vz_components = neutral_data.vz_components;

        input_crate.E_kin = neutral_data.E_kin;
        input_crate.E_pot = neutral_data.E_pot;
        input_crate.L_droplet = neutral_data.L_droplet;

        % ion dynamics data
        input_crate.vx_total = ion_data.vx_total;
        input_crate.vy_total = ion_data.vy_total;
        input_crate.vz_total = ion_data.vz_total;


        input_crate.x_ci  = ion_data.x_ci;
        input_crate.y_ci  = ion_data.y_ci;
        input_crate.z_ci  = ion_data.z_ci;


        input_crate.b_ion_outside = ion_data.b_ion_outside;
        

       % vmi_sim_visualize_particle_paths(input_crate);

        test = 1;
        post_process_function(input_crate);
        generate_post_process_figures;


        
        %vmi_sim_visualize_trajectory(input_crate);
        % debugging 
      %  figure
         %plot(sqrt(ion_data.vx_total.^2 + ion_data.vy_total.^2 + ion_data.vz_total.^2));

       % plot(sqrt(neutral_data.vx_components.^2 + neutral_data.vy_components.^2 + neutral_data.vz_components.^2)')
       % plot(sqrt(neutral_data.x_components(1:num_molecules,:)  ))

        %% version for plotting neutral dynamics with ion simulation
    case 'neutral'
        warning('visualizing neutral data only');
        
        input_crate = struct();


        % neutral dynamics data
        input_crate.b_invalid_spawn = neutral_data.b_invalid_spawn;
        input_crate.mass = neutral_data.mass;
        input_crate.num_molecules = neutral_data.num_molecules;

        input_crate.time = neutral_data.time;
        input_crate.droplet_radii = neutral_data.droplet_radii;

        input_crate.x_components = neutral_data.x_components;
        input_crate.y_components = neutral_data.y_components;
        input_crate.z_components = neutral_data.z_components;

        input_crate.vx_components = neutral_data.vx_components;
        input_crate.vy_components = neutral_data.vy_components;
        input_crate.vz_components = neutral_data.vz_components;

        input_crate.E_kin = neutral_data.E_kin;
        input_crate.E_pot = neutral_data.E_pot;
        input_crate.L_droplet = neutral_data.L_droplet;

        % ion dynamics data
        input_crate.vx_total = neutral_data.vx_components;
        input_crate.vy_total = neutral_data.vy_components;
        input_crate.vz_total = neutral_data.vz_components;

        input_crate.b_ion_outside = true(size(neutral_data.vx_components)); % all timesteps are considered outside ions
        
        input_crate.time_i = neutral_data.time;

        post_process_function(input_crate);

        generate_post_process_figures;
        
        %generate_post_process_figures;
        vmi_sim_visualize_particle_paths(input_crate);



        %% version for plotting ion dynamics
    case 'ion'
        input_crate = struct();
        input_crate.droplet_radii = neutral_data.droplet_radii;

        input_crate.b_invalid_spawn = neutral_data.b_invalid_spawn;

        input_crate.mass = ion_data.mass_i;
        input_crate.num_molecules = neutral_data.num_molecules;

        input_crate.time = ion_data.time_i;
        input_crate.droplet_radii = neutral_data.droplet_radii;

        input_crate.x_components = ion_data.x_ci;
        input_crate.y_components =  ion_data.y_ci;
        input_crate.z_components =  ion_data.z_ci;

        input_crate.vx_components = ion_data.vx_ci;
        input_crate.vy_components = ion_data.vy_ci;
        input_crate.vz_components = ion_data.vz_ci;

        input_crate.E_kin = ion_data.E_kin_ion;
        input_crate.E_pot = ion_data.E_pot_ion;
       %input_crate.L_droplet = ion_data.L_droplet_ion;

        input_crate.vx_total = ion_data.vx_ci*0;
        input_crate.vy_total = ion_data.vy_ci*0;
        input_crate.vz_total = ion_data.vz_ci*0;

        input_crate.b_ion_outside = ion_data.vx_ci*0 ==0; % spoof the ion ejection data, setting it so that the ion is always outside if ion dynamics are plotted!





end



%%

%% visualize selected neutral and ion trajectories
input_crate.x_ci = ion_data.x_ci;
input_crate.y_ci = ion_data.y_ci;
input_crate.z_ci = ion_data.z_ci;


%test = 1;
%return
%vmi_sim_visualize_ensemble(input_crate);

%vmi_sim_visualize_distributions(input_crate);

%return

vmi_sim_visualize_distributions(input_crate)



%post_process_function(input_crate);
%generate_post_process_figures;
%vmi_sim_visualize_particle_paths(input_crate);
%vmi_sim_visualize_trajectory(input_crate);


end
%% optional random seed
seed = 123;
%rng(seed);

global custom_DFT_start 



% Recreate the handle (use a robust form; exp(-2*log(v)) == v.^(-2) for v>0)
sigma_lookup = @(v) v.^(-2);
global DEBUG

run physical_constants.m

global geometric_scattering_crosssection_I
mean_free_path = 1/(density_droplet*  geometric_scattering_crosssection_I);

global neutral_scatter_angle_std

global scattering_probability

global hard_sphere_collision_mode

global single_initial_position

global deltaR0

global T_particles
%% simulation parameters
t_max = 200; % ps

%dt = 0.005;
dt = 0.005;
dt  =0.01;

if single_pulse
    t_max = dt*2;
end


num_timesteps = ceil(t_max/dt);

global num_molecules
m =127*u + zeros(num_molecules,1);


global partner_interaction;

%% choose between effusive and droplet dynamics
global effusive_dynamics

if effusive_dynamics
    relative_energy_loss = 0;
    he_direction_scattering =0;
    hard_sphere_collisions = 0;
    binding_energy_I_atom = 0;
    attach_he = 0;
else

    % simple energy loss relative to current atom energy
    relative_energy_loss = 0;
    he_direction_scattering =0;
    scatter_strength = 0.0*1E-2*he_direction_scattering;

    % hard sphere collisions
    hard_sphere_collisions =  1;

    %binding_energy_I_atom = 0.031; % binding energy in eV (this value was used a lot)
    binding_energy_I_atom = 318.43*k_B/eV; %value from the notes
    %binding_energy_I_atom = 0.015;
    %binding_energy_I_atom = 0.002;

    % attach heliums, this is currently not used
    attach_he = 0;
end


debug_plot = true;
if num_molecules>10
    debug_plot = false;
end




%% set initial kinetic energies based on pump laser wavelength and bandwidth
global lambda_pump
global E_diss



fwhm_lambda = 33; % nm
%fwhm_lambda = 0;

E_initial = hc/lambda_pump*eV;

% pump energy width
fwhm_E = hc/lambda_pump^2 * fwhm_lambda ; % in eV
fwhm_v = 1/2*fwhm_E*eV/sqrt(E_initial*127*u) / 100; % in Angström per picosecond




if partner_interaction
    mean_v = sqrt(E_initial./m) / 100;
else
    mean_v = sqrt((E_initial - E_diss)./m)/100;
end

% check in with Brauns simulation
% mean_v =mean_v*0 + 8;
% fwhm_v = fwhm_v*0.2;

if effusive_dynamics
    fwhm_v = 0;
end


%% parameters for cooling, only used with simple exponential energy loss
%velocity_decrease_factor = 0.0001;
%data_in = load("cooling_rate_parameters.mat"); % E_min, relative_energy_loss_per_ps
%E_min = data_in.E_min;
%relative_energy_loss_per_ps = data_in.relative_energy_loss_per_ps;
% these parameters were obtained from a fit of the energy loss of a xenon
% atom in the coppens2017 paper, DOI: 10.1039/C7CP03307A

global v_limit
E_min = (127*u)*v_limit^2/2/eV;

if E_min>binding_energy_I_atom
    warning('all neutrals will escape!');
end

%E_min = 1.645/1000; % from 50 m/s landau velocity of I atoms

%E_min = 6.58/ 1000; % from 100 m/s landau velocity of I atoms
%E_min = 5.33/ 1000; % from 90 m/s landau velocity of I atoms
%% generate droplet radii base on log normal distribution
global use_single_droplet_size
global single_droplet_size
if~use_single_droplet_size
    N = generate_droplet_sizes(num_molecules,false, true);
    %N = N/2;
else
    N = zeros(num_molecules,1) + single_droplet_size;
end

%% convert droplet sizes in num of helium atoms to radii in angström
droplet_radii = 2.22*N.^(1/3) ;

%rho_bulk = 0.0218;
%rho_droplet = rho_bulk*0.8;
%droplet_radii = (3/(4*pi*rho_droplet))^(1/3) * N.^(1/3); % this formula
%would make the droplets a bit bigger, 30 Angström instead of 28 Angström
%for a 2000 He droplet

fprintf('mean N: %.0f\n', mean(N));
fprintf('mean R: %.0f\n', mean(droplet_radii));


%% droplet solvation potential
%droplet_potential = @(beta, x) (erf((x-beta(3))/beta(1))*1+1)/2*beta(2);

%potential_steepness = 11.89;
potential_steepness = 14.2;
%potential_steepness = 1;

droplet_potential_atom = @(r) droplet_potential([potential_steepness, binding_energy_I_atom, 0],r);

figure
r = [0:0.1:80];
 plot(r, droplet_potential_atom(r - mean(droplet_radii) ))
 hold on
hline(E_min);

h = 0.000001;
droplet_force = @(x) (droplet_potential_atom(x+h) - droplet_potential_atom(x))/h;


%% random sampling of molecule initial positions and orientations
E_max = 200; % in meV
T_droplet = T_particles; % in K
%T_droplet = 4;
r_step = 0.01;
debug_plot_sample_generation = true;

% from ernesto dft result beta = [14.3324   26.9916   34.4431]
potential_steepness_molecule = 14.3324; % from fit of solvation potential DFT result
%binding_energy_molecule = 26.9; % in meV for I2 in He from fit of solvation potential DFT result
binding_energy_molecule =  573.3*k_B/eV*1000;

% Xenon as a stand in
%potential_steepness_molecule = 11.89;
%binding_energy_molecule = 18; % in meV for Xenon in HeN, assuming I2 behaves similarly

r_accepted = generate_radial_samples_3d(droplet_radii, ...
    @droplet_potential, potential_steepness_molecule, binding_energy_molecule, ...
    E_max, T_droplet, r_step, debug_plot_sample_generation);


%generate_radial_samples_3d_paper_fig(droplet_radii, ...
%  @droplet_potential, potential_steepness_molecule, binding_energy_molecule, ...
%    E_max, T_droplet, r_step, debug_plot_sample_generation);

% start molecules at positions in the potential corresponding to boltzmann
% energy distribution at 0.34 K
r0 = r_accepted; % molecule distance from center

if single_initial_position
    r0 = r0*0;
end


E_pot_molecule = droplet_potential([potential_steepness_molecule, binding_energy_molecule, 0], r_accepted - droplet_radii);


E_thermal_molecule = 3/2*k_B*T_droplet;
v_thermal_molecule = sqrt(2*E_thermal_molecule/(2*127*u))/100;

%scatter(r0 - droplet_radii, E_pot_molecule)


% the following does not work, produces more samples at the poles
% beta = rand(size(r0))*pi*2; % azimuth angle of molecule location inside droplet
% gamma = rand(size(r0))*pi*pi; %  polar angle

% try this version, fixes the random sampling of molecule positions
beta= rand(size(r0))*2*pi;
costheta = rand(size(r0))*2 -1;
gamma = acos(costheta);


% distance of atoms in the molecule at t0
global R0_GS
molecule_equilibrium_distance = R0_GS + randn(size(r0))*deltaR0;


if single_pulse 
    
    %anisotropic molecule angle
    alpha_accept = [];
    delta_accept = [];


    while length(alpha_accept)< num_molecules
        alpha = rand(size(r0))*2*pi;
        costheta = rand(size(r0))*2-1;
        delta = acos(costheta);

        cosphi = cos(alpha).*sin(delta); % angle of molecule to x-axis

        p = abs(cosphi.^2); % from: Molecular reorientation during dissociative multiphoton ionization, Physical Review A, 1993 - APS
        u_sample = rand(size(p));

        b_accept = u_sample<p;

        alpha_accept = [alpha_accept; alpha(b_accept)];
        delta_accept = [delta_accept; delta(b_accept)];

    end

    alpha = alpha_accept(1:num_molecules);
    delta = delta_accept(1:num_molecules);


else
    %isotropic molecule angle
    alpha = rand(size(r0))*2*pi;
    costheta = rand(size(r0))*2-1;
    delta = acos(costheta);
end


% all molecules with same angle and at edge of droplet
%
%  alpha = 0*alpha;
%  delta = 0*delta + pi/2;
% %
%  beta = beta*0+0;
%  gamma = gamma*0 + pi/2;
% r0 = r0*0 + droplet_radii(1);


%% initialize variables that hold information about the molecules



v0 = randn(num_molecules,1)*fwhm_v + mean_v;

if single_pulse
    v0 = v0*0;
end


x0 = r0.*cos(beta).*sin(gamma) + cos(alpha).*sin(delta).*molecule_equilibrium_distance/2;
y0 = r0.*sin(beta).*sin(gamma) + sin(alpha).*sin(delta).*molecule_equilibrium_distance/2;
z0 = r0.*cos(gamma) + cos(delta).*molecule_equilibrium_distance/2;


direction_molecule = [cos(beta).*sin(gamma), sin(beta).*sin(gamma), cos(gamma)];
direction_molecule = [direction_molecule; direction_molecule];

%scatter3(direction_molecule(:,1), direction_molecule(:,2), direction_molecule(:,3))

direction_atom1 = [cos(alpha).*sin(delta), sin(alpha).*sin(delta+pi), cos(delta)];
direction_atom2 = [cos(alpha).*sin(delta), sin(alpha).*sin(delta+pi), cos(delta+pi)];
%scatter3(direction_atom(:,1), direction_atom(:,2), direction_atom(:,3));
direction_atom = [direction_atom1; direction_atom2];

angle = sum( direction_atom.*direction_molecule,2);

figure
scatter3(x0,y0,z0, '.'); % check starting distribution

x0_twin = r0.*cos(beta).*sin(gamma) + cos(alpha).*sin(delta + pi).*molecule_equilibrium_distance/2;
y0_twin= r0.*sin(beta).*sin(gamma) + sin(alpha).*sin(delta + pi ).*molecule_equilibrium_distance/2;
z0_twin = r0.*cos(gamma) + cos(delta + pi).*molecule_equilibrium_distance/2;


vx0 = v0.*cos(alpha).*sin(delta);
vy0 = v0.*sin(alpha).*sin(delta);
vz0 = v0.*cos(delta);

vx0_twin = v0.*cos(alpha).*sin(delta + pi);
vy0_twin = v0.*sin(alpha).*sin(delta + pi);
vz0_twin = v0.*cos(delta + pi);

interatomic_axis_unit_vectors = zeros(num_molecules*2, 3);
interatomic_axis_unit_vectors(:,1) = [cos(alpha).*sin(delta);  cos(alpha).*sin(delta + pi)];
interatomic_axis_unit_vectors(:,2) = [sin(alpha).*sin(delta);  sin(alpha).*sin(delta + pi)];
interatomic_axis_unit_vectors(:,3) = [cos(delta);  cos(delta + pi)];

%scatter3(interatomic_axis_unit_vectors(:,1), interatomic_axis_unit_vectors(:,2), interatomic_axis_unit_vectors(:,3));


distance_travelled_since_last_collision = zeros(num_molecules*2,1);
collision_count = zeros(size(distance_travelled_since_last_collision));


% allocate memory for the positions and velocities at each timestep
% note that the the first num_molecules entries are for 1. atom of the molecule
% and the last num_molecules are for the 2. atom of the molecule

x_components = zeros(num_molecules*2, num_timesteps);
y_components = zeros(num_molecules*2, num_timesteps);
z_components = zeros(num_molecules*2, num_timesteps);

vx_components = zeros(num_molecules*2, num_timesteps);
vy_components = zeros(num_molecules*2, num_timesteps);
vz_components = zeros(num_molecules*2, num_timesteps);

ax_components = zeros(size(vx_components));
ay_components = zeros(size(vx_components));
az_components = zeros(size(vx_components));


x_components(:,1) = [x0; x0_twin ];
y_components(:,1) = [y0; y0_twin];
z_components(:,1) = [z0; z0_twin];


r_atoms = sqrt(x_components(:,1).^2 +y_components(:,1).^2 + z_components(:,1).^2);

sum(r_atoms>droplet_radii(1))/numel(r_atoms);

b_spawned_outside = r_atoms>[droplet_radii; droplet_radii];

b_invalid_spawn = b_spawned_outside(1:num_molecules) | b_spawned_outside(num_molecules+1:end);
b_invalid_spawn = [b_invalid_spawn; b_invalid_spawn];


vx_components(:,1) = [vx0; vx0_twin ];
vy_components(:,1) = [vy0; vy0_twin];
vz_components(:,1) = [vz0; vz0_twin];

ax_components(:,1) =[x0*0; x0*0];
ay_components(:,1) =[x0*0; x0*0];
ay_components(:,1) = [z0*0; z0*0];

time = zeros(1, num_timesteps);
E_dissip = zeros(size(vx_components));
E_kin = zeros(size(vx_components));
E_pot = zeros(size(vx_components));
%E_mass_attach_defect 23.01.25: no mass attachment in neutral propagation

droplet_radii = [droplet_radii; droplet_radii];

m = [m; m];
%m = [m; m*1.5];


mass = zeros(size(x_components));
mass(:,1) = m;


if custom_DFT_start
start_time = 3;

global Xdip_active
if Xdip_active
start_time = 0.8; 
end

%warning('start time at 0.9 ps!')

dft_fit = load('T:\github synchronized\Iodine_Helium_Simulation\HeDFT_MD_comparison_neutral\custom_start_interpolating_functions.mat');

% mimic the dynamics of TD-HeDFT
tau_id = 1;
tau = 0;
while tau < start_time


    % add R data from He DFT to random molecule offset position
x_components(1:num_molecules,tau_id) = interatomic_axis_unit_vectors(1:num_molecules,1)*dft_fit.R_interp(tau)/2 + r0.*cos(beta).*sin(gamma) ;
y_components(1:num_molecules,tau_id) = interatomic_axis_unit_vectors(1:num_molecules,2)*dft_fit.R_interp(tau)/2 + r0.*sin(beta).*sin(gamma) ;
z_components(1:num_molecules,tau_id) = interatomic_axis_unit_vectors(1:num_molecules,3)*dft_fit.R_interp(tau)/2 + r0.*cos(gamma);

x_components(num_molecules+1:end,tau_id) = interatomic_axis_unit_vectors(num_molecules+1:end,1)*dft_fit.R_interp(tau)/2 + r0.*cos(beta).*sin(gamma) ;
y_components(num_molecules+1:end,tau_id) = interatomic_axis_unit_vectors(num_molecules+1:end,2)*dft_fit.R_interp(tau)/2 + r0.*sin(beta).*sin(gamma) ;
z_components(num_molecules+1:end,tau_id) = interatomic_axis_unit_vectors(num_molecules+1:end,3)*dft_fit.R_interp(tau)/2 + r0.*cos(gamma);

% R_test = sqrt((x_components(1:num_molecules,tau_id) -x_components(num_molecules+1:end,tau_id) ).^2 + ...
%         	    (y_components(1:num_molecules,tau_id) -y_components(num_molecules+1:end,tau_id) ).^2 + ...
%                     (z_components(1:num_molecules,tau_id) -z_components(num_molecules+1:end,tau_id) ).^2);
% 
% 
% plot(interatomic_axis_unit_vectors(1:num_molecules,1));
% hold on
% plot(interatomic_axis_unit_vectors(num_molecules+1:end,1));
% 
% plot(x_components(1:num_molecules,tau_id) );
% hold on
% plot(x_components(num_molecules+1:end,tau_id));


vx_components(:,tau_id) = interatomic_axis_unit_vectors(:,1)*dft_fit.vz_interp(tau);
vy_components(:,tau_id) = interatomic_axis_unit_vectors(:,2)*dft_fit.vz_interp(tau);
vz_components(:,tau_id) = interatomic_axis_unit_vectors(:,3)*dft_fit.vz_interp(tau);

% sqrt(vx_components(:,tau_id).^2 + vy_components(:,tau_id) .^2 + vz_components(:,tau_id) .^2)   ;

%sqrt(interatomic_axis_unit_vectors(:,1).^2 + interatomic_axis_unit_vectors(:,2).^2 + interatomic_axis_unit_vectors(:,3).^2 );
plot3([0,vx_components(1,tau_id)], [0,vy_components(1,tau_id)], [0,vz_components(1,tau_id)])
hold on
plot3([0,vx_components(1+num_molecules,tau_id)], [0,vy_components(1+num_molecules,tau_id)], [0,vz_components(1+num_molecules,tau_id)])
tau = tau+ dt;
tau_id = tau_id+1;

time(:,tau_id) = tau;
mass(:,tau_id) = m;
end

t_id = tau_id-1;


else

    start_time = 0;
    time(:,1) = [start_time];
    t_id = 1;
end



% R_test = sqrt((x_components(1:num_molecules,1:tau_id) -x_components(num_molecules+1:end,1:tau_id) ).^2 + ...
%         	    (y_components(1:num_molecules,1:tau_id) -y_components(num_molecules+1:end,1:tau_id) ).^2 + ...
%                     (z_components(1:num_molecules,1:tau_id) -z_components(num_molecules+1:end,1:tau_id) ).^2);
% 
% 


% add observables
E_kin(:,1) = [m.*(vx_components(:,1).^2 + vy_components(:,1).^2 + vz_components(:,1).^2).^2/2/eV];

E_pot(:,1) = droplet_potential_atom( sqrt(x_components(:,1).^2+ y_components(:,1).^2 + z_components(:,1).^2) - droplet_radii);

L_droplet = zeros(size(x_components)); % pathlength traveled through droplet
% used for filtering later




%% debug plot for small number of moleculs
if num_molecules<10
    figure
    for i=1:size(vx_components,1)
        %disp(col(i,:))

        col = lines;


        [X,Y,Z] = sphere;
        X = X*mean(droplet_radii);
        Y = Y*mean(droplet_radii);
        Z = Z*mean(droplet_radii);

        s=  surface(X,Y,Z, 'FaceAlpha',0.1, 'FaceLighting','gouraud', 'EdgeColor',[0.1,0.3,0.9], 'EdgeAlpha',0.2, 'EdgeLighting','gouraud');
        s.CData = Z*0;
        colormap('winter')

        if i==1
            l = light;
            l.Color = [0.8,0.8,1];
        end

        hold on
        %p_center_of_mass = scatter3(r0(i,:)*cos(beta(i,:)),r0(i,:).*sin(beta(i,:)), 'MarkerEdgeColor',col(i,:));

        s1 = scatter3(x_components(i,t_id), y_components(i,t_id),z_components(i,t_id), 'MarkerEdgeColor',col(i+2,:));
        %scatter3(x0_twin(i,:), y0_twin(i,:),z0_twin(i,:), 'MarkerEdgeColor',col(i,:));

        %quiver3(x_components(i,t_id),y_components(i,t_id), z_components(i,t_id), vx_components(i,t_id),vy_components(i,t_id), vz_components(i,t_id),'Color',col(i,:));
        %quiver3(x0_twin(i,:),y0_twin(i,:),z0_twin(i,:),  vx0_twin(i,:),vy0_twin(i,:),vz0_twin(i,:), 'Color',col(i,:));

        %q1 = quiver3(x_components(i,t_id),y_components(i,t_id), z_components(i,t_id), ax_components(i,t_id),ay_components(i,t_id), az_components(i,t_id),'Color',col(i,:));

    end

end

% generate random sampled free paths for each particle
u_sample = rand(size(distance_travelled_since_last_collision));
free_path = -mean_free_path*log(1 - u_sample);
% free path distribution is approximately exponential: 10.1103/PhysRevE.77.041117


frog_step_crate = struct;
frog_step_crate.m = mass(:,1);
frog_step_crate.droplet_radii =droplet_radii;
frog_step_crate.droplet_force = droplet_force;
frog_step_crate.he_direction_scattering = he_direction_scattering;



%% normal propagation
while t_id<num_timesteps



    % get current masses
    m = mass(:,t_id);
    frog_step_crate.m = m;

    % leap frogging
    % get current position at t_id
    x0 = x_components(:,t_id);
    y0 = y_components(:,t_id);
    z0 = z_components(:,t_id);

    % get current velocities
    vx0 = vx_components(:,t_id);
    vy0 = vy_components(:,t_id);
    vz0 = vz_components(:, t_id);


    [x1,y1,z1, vx1, vy1, vz1, E_pot_partner] = frog_step_neutral(x0,y0,z0, vx0, vy0, vz0, frog_step_crate, dt);

    % calculate depth here!!
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
    depth = r1 - droplet_radii;


    % after collisions with droplet wall, implement energy loss
    v = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
    v_unit_x = vx1./v;
    v_unit_y = vy1./v;
    v_unit_z = vz1./v;




    E0 =(v*100).^2.*m/2/eV;
    E1 = E0;

    % in case of no energy loss, set new energy to old energy
    % this is used to keep track of the energy over time
    E1 = E0;

    if relative_energy_loss
        % relative energy loss
        dE = heaviside(-depth).*E0*relative_energy_loss_per_ps*dt.*(1 + randn(size(depth))*0.1);


        E1 = max( [E0 - dE, zeros(size(E0)) + min(E_min, E0) ],[],2); % energy after energy loss is limited: to E_min or E_0 if E_0 was already smaller than E_min
        v = sqrt(2*E1*eV./m)/100;

        vx1 = v.*v_unit_x;
        vy1 = v.*v_unit_y;
        vz1 = v.*v_unit_z;

    end


    if hard_sphere_collisions
        % this section is based on Andreas Braun's phd thesis, section D.2





        switch hard_sphere_collision_mode
            case 1
                trial_random_number = rand(size(E0));
                b_collision = trial_random_number < scattering_probability & (depth <0);

            case 2
                b_collision = (distance_travelled_since_last_collision> free_path) & (depth<0);
            
            case 3 
                    if t_id>1
                    distance_travelled_in_timestep =sqrt(   ( x_components(:, t_id) - x_components(:, t_id-1)).^2 + ...
                                                                     (y_components(:, t_id) - y_components(:, t_id-1)).^2 + ...
                                                                                  (z_components(:, t_id) - z_components(:, t_id-1)).^2   );

                    bulk_density_helium = 0.0219; % Angström ^-3
                    density_droplet = 0.8*bulk_density_helium;
                    trial_random_number = rand(size(E0));

                    if sigma_dependent_on_v
                        sigma = sigma_lookup(v)*0.5;
                    else
                        sigma = geometric_scattering_crosssection_I;
                    end


                    p_scatter = distance_travelled_in_timestep.*sigma*density_droplet ;
                    b_collision =(trial_random_number < p_scatter) & (depth<0);
  

                    test = 1;
                    % todo: implement v dependent crosssection here
                    else
                        b_collision = false(size(E0));
                    end


        end

        if DEBUG
            b_collision = true(size(b_collision));
        end

        % no collision if energy is already below minimum energy
        b_landau = E0 < E_min;
        if sum(b_landau) ~=0
            test = 1;
        end

        b_collision = b_collision & ~b_landau;



        %sum(b_collision)
        collision_count = collision_count + b_collision;

        % due to geometry of hard sphere collision in 3d, the probability
        % for a impact parameter is p(b) 2 pi b/ (pi R ^2) where R is the sum of
        % radii of the colliding spheres
        % the associated CDF is then CDF(b) = b^2/R^2
        % using inverse CDF sampling, given random numbers x between zero
        % and 1, the function ICDF(x) = sqrt(x) will yield random samples of b/R


        impact_parameter_norm = sqrt( rand(size(E0))); % b/R
        % sampling in this way guarantees that the fact that impact
        % parameters b/R closer to one have a higher probability than impact
        % parameters closer to zero (= grazing collisions are more likely
        % than head on collisions)

        COSTHETA = (2*impact_parameter_norm.^2 - 1) ; % costheta seems to be evenly distributed again
        SINTHETA = sqrt(1- COSTHETA.^2);


        COSTHETA(~b_collision ) = 1;
        SINTHETA(~b_collision ) = 0;


        % mass ratio
        global scatter_mass_neutral;
        RHO = (m/u)/scatter_mass_neutral ;


        % get new kinetic energy after scattering event
        E1 = E0.* (1 + 2*RHO.*COSTHETA + RHO.^2)./(1 + RHO).^2;

        % THETA is the angle in the center of mass frame
        % ==> transform to laboratory frame

        COStheta = (COSTHETA+ RHO)./sqrt(1 + 2*RHO.*COSTHETA+ RHO.^2);
        SINtheta= sqrt(1- COStheta.^2);
        % COStheta = COSTHETA;
        % SINtheta = SINTHETA;

        % additional angle variation
        theta = acos(COStheta(b_collision)) + randn(size(COStheta(b_collision)))*neutral_scatter_angle_std*pi/180;
        COStheta(b_collision)  = cos(theta);
        SINtheta(b_collision) = sqrt(1- COStheta(b_collision).^2);



        % calculate unit vectors of plane normal to old velocity vector
        % (velocity_normal_x and velocity_normal_y)
        velocity_unit_vectors = [v_unit_x, v_unit_y, v_unit_z];

        %reference_direction = [1,0,0]; % if reference direction is along one of the cartesian axis, there is a weird error where th scattering directions are not random but have a bias towards one direction


        % a = velocity_unit_vectors(:,1);
        %b = velocity_unit_vectors(:,2);
        %c = velocity_unit_vectors(:,3);
        % velocity_normal_1 = [b+c, c-a, -a-b]; % this also has a preferred direction..

        reference_direction = rand(size(velocity_unit_vectors))-0.5;
        reference_direction=  reference_direction./sqrt(sum( reference_direction.* reference_direction, 2));

        velocity_normal_1 = cross(velocity_unit_vectors, reference_direction, 2);

        velocity_normal_1 = velocity_normal_1./sqrt(sum(velocity_normal_1.*velocity_normal_1, 2));

        %velocity_normal_1 = cross(velocity_unit_vectors,repmat( reference_direction, size(v_unit_x,1), 1), 2);
        %velocity_normal_1= velocity_normal_1./sqrt( sum( velocity_normal_1.^2 ,2));

        velocity_normal_2 = cross(velocity_unit_vectors,velocity_normal_1, 2);

        % get random scattering direction in this plane

        COSBETA = (rand(size(E0))-0.5)*2; % beta angle selects the velocity change direction in the plane normal to the incident velocity vector
        SINBETA = sqrt(1 - COSBETA.^2);

        % calculate new velocities
        v_new = sqrt(2*E1*eV./m)/100;

        v_normal = v_new.*SINtheta;
        v_parallel = v_new.*COStheta;

        new_velocity_vectors = velocity_unit_vectors.*v_parallel  ...
            + velocity_normal_1.*COSBETA.*v_normal  + velocity_normal_2.*SINBETA.*v_normal;

        %completely random new velocity
        % new_velocity_vectors = reference_direction.*v_new;

        %new_velocity_vectors =[ reference_direction(:,1).*v_new, reference_direction(:,2).*v_new, reference_direction(:,3).*v_new];

        if sum(b_collision)>0
            test = 1;
        end


        % check out the three orthonormal vectors
        if DEBUG
            %
            % figure
            % plot_vector([0,0,0], velocity_unit_vectors(1,:));
            % hold on
            % plot_vector([0,0,0],  velocity_normal_1(1,:));
            % plot_vector([0,0,0],  velocity_normal_2(1,:));
            % pbaspect([1,1,1])
            %
            % dot(velocity_normal_1(1,:), velocity_normal_2(1,:))
            % dot(velocity_normal_1(1,:), velocity_unit_vectors(1,:))
            %
            % dot(velocity_unit_vectors(1,:), velocity_normal_2(1,:))
            %
            % dot(velocity_unit_vectors(1,:), velocity_normal_1(1,:))*180/pi
            % dot(velocity_unit_vectors(1,:), velocity_normal_2(1,:))*180/pi
            %
            % norm(velocity_unit_vectors(1,:))
            % norm(velocity_normal_2(1,:))
            % norm(velocity_normal_1(1,:))
            close all
            figure


            plot_vector([0,0,0], velocity_unit_vectors(1,:)*v(1));
            hold on
            plot_vector([0,0,0],  velocity_normal_1(1,:));
            plot_vector([0,0,0],  velocity_normal_2(1,:));

            pbaspect([1,1,1])
            xlim([-5,5])
            ylim([-5,5])
            zlim([-5,5]);

            plot_vector([0,0,0],velocity_unit_vectors(1,:).*v_parallel(1));
            hold on
            plot_vector([0,0,0],   + velocity_normal_1(1,:).*COSBETA(1).*v_normal(1)  + velocity_normal_2(1,:).*SINBETA(1).*v_normal(1))

            plot_vector([0,0,0],  new_velocity_vectors(1,:));

            legend({'initial v', 'vnx', 'vny', 'new v parallel', 'new v normal', 'new v'});

            id_collisions = 1:length(b_collision);

            id_collisions = id_collisions(b_collision);
            if length(id_collisions)>0
                % check if the angle between old and new vectors is correct
                fprintf('angle old vs new %.2f, acos(cos(theta)) %.2f, asin(sin(theta)) %.2f\n',dot(velocity_unit_vectors(id_collisions(1),:), new_velocity_vectors(id_collisions(1),:)') /norm(new_velocity_vectors(id_collisions(1),:))*180/pi,...
                    acos(COStheta(id_collisions(1)))*180/pi,...
                    asin(SINtheta(id_collisions(1)))*180/pi)
            end


            %return

        end


        vx1 = new_velocity_vectors(:,1);
        vy1 = new_velocity_vectors(:,2);
        vz1 = new_velocity_vectors(:,3);


    end





    actual_dE = E0 - E1;
    b_decelerated = actual_dE>0;
    %vel_angle = vel_angle + (rand(size(v))- 0.5)*2*2*pi*scattering_angle_width/360.*b_decelerated; %random direction change after deceleration

    if hard_sphere_collisions
        if sum(b_collision)>0


            scattering_angle = acos(mean(COSTHETA(b_collision)));
            percent_energy_loss = mean(actual_dE(b_collision))/mean(E0(b_collision));

            if hard_sphere_collision_mode==2
                % important step: set distance traveled to zero if collision happened!
                distance_travelled_since_last_collision(b_collision) = 0;

                % determine new free path for atoms that collided
                u_sample = rand(size(distance_travelled_since_last_collision(b_collision)));
                free_path(b_collision) = -mean_free_path*log(1 - u_sample);
            end


        end
    end


    %% attach heliums at droplet surface
    if attach_he
        [v, M] = attach_Helium_3d(v, m, depth, binding_energy_I_atom, ...
            droplet_potential_atom( sqrt(x1.^2+ y1.^2 + z1.^2) - droplet_radii)/empirical_potential_factor,...
            [v_unit_x, v_unit_y, v_unit_z], [r_unit_x1, r_unit_y1, r_unit_z1]);




        m1 = M;
    else
        m1 = m;
    end





    % end of timestep, assign new velocities, positions and accelerations
    % to arrays
    vx_components(:, t_id+1) = vx1;
    vy_components(:, t_id+1) = vy1;
    vz_components(:, t_id+1) = vz1;

    x_components(:, t_id+1) = x1;
    y_components(:, t_id+1) = y1;
    z_components(:, t_id+1) = z1;

    time(:, t_id+1) = time(:,t_id) + dt;

    %  ax_components(:,t_id + 1) = ax1;
    %   ay_components(:,t_id+1) = ay1;
    %    az_components(:,t_id+1) = az1;

    mass(:, t_id+1) = m1;

    v = sqrt(vx1.^2 + vy1.^2 + vz1.^2);

    E_kin(:, t_id + 1) = m1.*(v*100).^2/2/eV ;

    E_pot(:,t_id+1) = droplet_potential_atom( sqrt(x1.^2+ y1.^2 + z1.^2) - droplet_radii) + [E_pot_partner;E_pot_partner]/2;

    E_dissip(:,t_id+1) = E_dissip(:,t_id) + actual_dE;

    distance_traveled = sqrt( (x_components(:, t_id+1) - x_components(:, t_id)).^2 + ...
        (y_components(:, t_id+1) - y_components(:, t_id)).^2 + ...
        (z_components(:, t_id+1) - z_components(:, t_id)).^2);

    distance_travelled_since_last_collision = distance_travelled_since_last_collision + distance_traveled;

    % get distance travelled in this timestep
    lx = x1 - x0;
    ly = y1 - y0;
    lz = z1 - z0;

    % if still inside droplet (depth<0) add distance travel this timestep
    L_droplet(:, t_id+1) = L_droplet(:, t_id) + (depth<0).*sqrt(lx.^2 + ly.^2 + lz.^2);



    t_id = t_id+1;






end

%% debug plot that shows trajectories
if debug_plot & num_molecules<10

    for i=1:size(x_components,1)
        %r =sqrt(x_vectors(i,:).^2 + y_vectors(i,:).^2);
        %v = sqrt( vx_vectors(i,:).^2 + vy_vectors(i,:).^2);
        hold on
        %plot(x_vectors(i,:), y_vectors(i,:),  'color', [0.3,0,0.8,0.1]);
        % show path through droplet and acceleration due to droplet potential
        plot3(x_components(i,:), y_components(i,:), z_components(i,:), 'color', col(i+2,:));

        % hold on
        %quiver3(x_components(i,:), y_components(i,:),z_components(i,:),  ax_components(i,:), ay_components(i,:), az_components(i,:));


        % show difference vectors
        %quiver(x_vectors(i,end), y_vectors(i,end), -diff_vector_x(i,end), -diff_vector_y(i,end));
        %scatter(x_vectors(i,end)-diff_vector_x(i,end), y_vectors(i,end)-diff_vector_y(i,end))

        % show normal velocity vectors
        %quiver(x_vectors(i,end), y_vectors(i,end), vx_vectors(i,end), vy_vectors(i,end));

        % show coulomb velocity vectors
        %quiver(x_vectors(i,end), y_vectors(i,end), v_coulomb_x(i,end), v_coulomb_y(i,end));

        % show coulomb velocity vectors
        %quiver(x_vectors(i,end), y_vectors(i,end), vx_total(i,end), vy_total(i,end));

        % figure(2);
        % hold on
        % plot(r,m*v.^2/2/eV*1000)
        %
        % figure(3)
        % hold on
        % plot(r, v)
        pbaspect([1,1,1])
        xlim([-60,60])
        ylim([-60,60])
        zlim([-60,60])
        view(60,30)
        set(gcf, 'Position', [680 225 1015 653]);

    end
    return
end



close all
E_system = E_kin + E_pot+ E_dissip;
% check energy conservation
plot(sum(E_kin,1));
hold on
plot(sum(E_pot,1));
plot(sum(E_dissip,1));
plot(sum(E_system,1));

legend('E_{kin}', 'E_{pot}', 'E_{dissip}', 'E_{system}');
title('Energy balance neutral atoms');

savefig(gcf, 'T:\github synchronized\Iodine_Helium_Simulation\debug_images\neutral_energy');



test = 1;


% reduce size of saved quantities
global num_neutral_export_timestept

t_select = logspace(0, log10(t_max), num_neutral_export_timestept)-1;
%t_select = [logspace(0,log10(10), 20), 21:1:300];

%t_select = unique(time);
indices_select = [];

for t_of_interest = t_select
    indices_select(end+1) = find(abs(time - t_of_interest)== min(abs(time-t_of_interest)));

end


indices_select = unique(indices_select);
t_select_from_time = time(indices_select);

% figure
% plot(t_select, t_select);
% hold on
% scatter(t_select_from_time, t_select_from_time);




mass = mass(:, indices_select);


x_components = x_components(:, indices_select);
y_components = y_components(:,indices_select);
z_components = z_components(:,indices_select);
vx_components = vx_components(:,indices_select);
vy_components = vy_components(:,indices_select);
vz_components = vz_components(:,indices_select);
E_pot = E_pot(:,indices_select);
E_kin = E_kin(:,indices_select);
E_dissip = E_dissip(:, indices_select);

L_droplet = L_droplet(:, indices_select);

time = time(indices_select);

save('neutral_propagation_checkpoint', 'direction_atom', 'direction_molecule','angle',...
    'mass', 'vx_components', 'vy_components', 'vz_components',...
    'eV', 'x_components', 'y_components', 'z_components', ...
    'num_molecules',  ...
    'time', 'droplet_radii', 'u', 'E_kin', 'E_pot', 'L_droplet',...
    'r0', 'E_initial', 'E_dissip','m', 'effusive_dynamics', ...
    'binding_energy_I_atom', 'he_direction_scattering', 'E_min','b_invalid_spawn',...
    'hard_sphere_collision_mode', 'scattering_probability', 'mean_free_path',  'geometric_scattering_crosssection_I',...
    'mean_free_path');

function [x1,y1,z1, vx1, vy1, vz1, Epot] = frog_step_ion(x0, y0, z0, vx0, vy0, vz0, fixed_inputs, dt)
    
Epot = zeros(size(x0));

%obtain inputs from struct 
%desolve_struct(fixed_inputs);% this loads: m, droplet_radii, droplet_force, partner_interaction_ions


% explicit assignment (its faster)
m = fixed_inputs.m;
droplet_radii = fixed_inputs.droplet_radii;
droplet_force =     fixed_inputs.droplet_force;
partner_interaction_ions =     fixed_inputs.partner_interaction_ions;
charge_i =    fixed_inputs.charge_i;
he_direction_scattering =     fixed_inputs.he_direction_scattering ;
 additional_droplet_charges =     fixed_inputs.additional_droplet_charges;
 droplet_charge_model =    fixed_inputs.droplet_charge_model;
ion_charge_positions =     fixed_inputs.ion_charge_positions;


% leap frogging



    if sum(imag(vx0))>0
        test = 1;
    end
    



    % calculate droplet depth and acceleration
    %depth = sqrt(x0.^2 + y0.^2) - droplet_radius*ones(size(y_vectors(:,t_id)));
    r0 = sqrt(x0.^2 + y0.^2 + z0.^2);
    depth = r0 - droplet_radii;
    F = droplet_force(depth)*1.602E-9 ; % convert from eV/Angström to Newton

    a = -F./m; % force in m/s^2

    a = a*1E-14; % convert to Angström per picosecond^2


    % get unit vector of position
    r_unit_x0 =  x0./r0;
    r_unit_y0 = y0./r0;
    r_unit_z0 = z0./r0;
    
    % acceleration acts in opposite direction of radial unit vector
    ax0 = a.*r_unit_x0;
    ay0 = a.*r_unit_y0;
    az0 = a.*r_unit_z0;
    
    if partner_interaction_ions
    [a_potx, a_poty, a_potz, ~] = add_partner_interaction_ion(x0,y0,z0, m, charge_i);
    

    %%diff = sqrt((x0 - circshift(x0,num_molecules)).^2 + (y0 - circshift(y0, num_molecules)).^2 + (z0 - circshift(z0,num_molecules)).^2 );
    %scatter(diff, sqrt(a_potx.^2 + a_poty.^2 + a_potz.^2) )
    
    a0test = sqrt(ax0.^2 + ay0.^2 + az0.^2);

    apottest = sqrt(a_potx.^2 + a_poty.^2 + a_potz.^2);

    ax0 = ax0 + a_potx;
    ay0 = ay0 + a_poty;
    az0 = az0 + a_potz;
    end
    
    if additional_droplet_charges>0
            
        switch droplet_charge_model
            case 1 % repulsion from helium ion fixed at the droplet center
            [a_he_ionx, a_he_iony, a_he_ionz, ~] = add_helium_interaction_coulomb(x0,y0,z0, m, ion_charge_positions);
            test = 1;
            
            case 2   % charged droplet 
                   [a_he_ionx, a_he_iony, a_he_ionz, ~] = add_helium_interaction_charged_droplet(x0,y0,z0, m, droplet_radii);
        end

        ax0 = ax0+ a_he_ionx;
        ay0 = ay0+ a_he_iony;
        az0 = az0+ a_he_ionz;
    end


    v0 = sqrt(vx0.^2 + vy0.^2 + vz0.^2);
    v_unit_x = vx0./v0;
    v_unit_y = vy0./v0;
    v_unit_z = vz0./v0;



    if he_direction_scattering % change velocity unit vector if direction should be randomized
    direction_change = randn(size(v_unit_x,1),3)*scatter_strength;
    v_unit_x = v_unit_x + direction_change(:,1);
    v_unit_y = v_unit_y + direction_change(:,2);
    v_unit_y = v_unit_y + direction_change(:,3);

    v_renorm = sqrt(v_unit_x.^2 + v_unit_y.^2 + v_unit_z.^2);
    v_unit_x = v_unit_x./v_renorm;
    v_unit_y = v_unit_y./v_renorm;
    v_unit_z = v_unit_z./v_renorm;
    
    vx0 = v0.*v_unit_x;
    vy0 = v0.*v_unit_y;
    vz0 = v0.*v_unit_z;
    end



    % get new positions
    x1 = x0  + dt*vx0 + 1/2*ax0*dt^2;
    y1 = y0  + dt*vy0 + 1/2*ay0*dt^2;
    z1 = z0  + dt*vz0 + 1/2*az0*dt^2;
    
      % get next acceleration to calculate next velocity
    %depth = sqrt(x1.^2 + y1.^2) - droplet_radius*ones(size(y_vectors(:,t_id)));
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);

    depth = r1 - droplet_radii;
    
    try
    F = droplet_force(depth)*1.602E-9; % convert from eV/Angström to Newton
    catch ME
test = 2;
    end
    a = -F./m; % force in m/s^2
    a = a*1E-14; % convert to Angström per picosecond^2

    r_unit_x1 =  x1./r1;
    r_unit_y1 = y1./r1;
    r_unit_z1 = z1./r1;
    
    % acceleration acts in opposite direction of radial unit vector
    ax1 = a.*r_unit_x1;
    ay1 = a.*r_unit_y1;
    az1 = a.*r_unit_z1;
    
    if partner_interaction_ions
    [a_potx, a_poty, a_potz, Epot] = add_partner_interaction_ion(x1,y1,z1, m, charge_i);

    ax1 = ax1+ a_potx;
    ay1 = ay1+ a_poty;
    az1 = az1+ a_potz;
    end
    
    if  additional_droplet_charges>0
            
        
        
        switch droplet_charge_model
            case 1 % repulsion from helium ion fixed at the droplet center
            [a_he_ionx, a_he_iony, a_he_ionz,Epot_he_ions] = add_helium_interaction_coulomb(x1,y1,z1, m, ion_charge_positions);
            
            case 2   % charged droplet 
                   [a_he_ionx, a_he_iony, a_he_ionz,Epot_he_ions] = add_helium_interaction_charged_droplet(x1,y1,z1, m, droplet_radii);



        end

     
     
        ax1 = ax1+ a_he_ionx;
        ay1 = ay1+ a_he_iony;
        az1 = az1+ a_he_ionz;
    else
        Epot_he_ions = 0;
    end


    vx1 = vx0 + 1/2*(ax0 + ax1)*dt;
    vy1 = vy0 + 1/2*(ay0 + ay1)*dt;
    vz1 = vz0 + 1/2*(az0 + az1)*dt;

    v0test = sqrt(vx0.^2 + vy0.^2 + vz0.^2);
    v1test = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
    
        a0test = sqrt(ax0.^2 + ay0.^2 + az0.^2);
        a1test = sqrt(ax1.^2 + ay1.^2 + az1.^2);
    vdiff = v1test - v0test;

    Epot = Epot + Epot_he_ions;
end

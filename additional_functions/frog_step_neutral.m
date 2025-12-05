function [x1,y1,z1, vx1, vy1, vz1, E_pot] = frog_step_neutral(x0, y0, z0, vx0, vy0, vz0, fixed_inputs, dt)
    

% frog step according to kick-drift-kick form
% v(t+h/2) = x(t) + (h/2)
global partner_interaction

Epot = zeros(size(x0));

%obtain inputs from struct 
desolve_struct(fixed_inputs);% this loads: m, droplet_radii, droplet_force, partner_interaction_ions


    % calculate depth of atoms in droplet and force due to droplet
    % potential
    r0 = sqrt(x0.^2 + y0.^2 + z0.^2);
    depth = r0 - droplet_radii;
    F = droplet_force(depth)*1.602E-9; % convert from eV/Angström to Newton

    a = -F./m; % force in m/s^2

    a = a*1E-14; % convert to Angström per picosecond^2

    % get unit vector components of position
    r_unit_x0 =  x0./r0;
    r_unit_y0 = y0./r0;
    r_unit_z0 = z0./r0;

    % acceleration acts in opposite direction of radial unit vector
    ax0 = a.*r_unit_x0;
    ay0 = a.*r_unit_y0;
    az0 = a.*r_unit_z0;

    if partner_interaction
        [a_potx, a_poty, a_potz,  ~] = add_partner_interaction(x0,y0,z0, m);

        ax0 = ax0 + a_potx;
        ay0 = ay0 + a_poty;
        az0 = az0 + a_potz;
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


    % second leap frog step
    % get new positions
    x1 = x0  + dt*vx0 + 1/2*ax0*dt^2;
    y1 = y0  + dt*vy0 + 1/2*ay0*dt^2;
    z1 = z0  + dt*vz0 + 1/2*az0*dt^2;


        % get next acceleration to calculate next velocity
    %depth = sqrt(x1.^2 + y1.^2) - droplet_radius*ones(size(y_vectors(:,t_id)));
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);

    depth = r1 - droplet_radii;

    F = droplet_force(depth)*1.602E-9; % convert from eV/Angström to Newton
    a = -F./m; % force in m/s^2
    a = a*1E-14; % convert to Angström per picosecond^2

    r_unit_x1 =  x1./r1;
    r_unit_y1 = y1./r1;
    r_unit_z1 = z1./r1;

    % acceleration acts in opposite direction of radial unit vector
    ax1 = a.*r_unit_x1;
    ay1 = a.*r_unit_y1;
    az1 = a.*r_unit_z1;

    if partner_interaction
        [a_potx, a_poty, a_potz, E_pot] = add_partner_interaction(x1,y1,z1, m);

        ax1 = ax1+ a_potx;
        ay1 = ay1+ a_poty;
        az1 = az1+ a_potz;

        
    else
        E_pot = zeros(size(x1,1)/2,1);
    end
    

    vx1 = vx0 + 1/2*(ax0 + ax1)*dt;
    vy1 = vy0 + 1/2*(ay0 + ay1)*dt;
    vz1 = vz0 + 1/2*(az0 + az1)*dt;

    

end




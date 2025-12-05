function [ax,ay,az,Epot] = add_helium_interaction_coulomb(x,y,z, mass, ion_positions)

global additional_droplet_charges



num_particles = size(x,1);

r1 = [x(1:num_particles/2), y(1:num_particles/2), z(1:num_particles/2)]; % position of first iodine
r2 =[x(num_particles/2+1:end), y(num_particles/2+1:end), z(num_particles/2+1:end)]; %position of second iodine

u = 1.66053907e-27; %kg

h= 0.0001;
for i=1:additional_droplet_charges

    dr_vec = [r1; r2] - ion_positions(:,:,i); % vector between iodine and additional charge

    dr = sqrt(sum(dr_vec.^2,2));

    dr_unit_vec = dr_vec./dr;
    
    Epot = coulomb_interaction_potential(dr);

    if i==1
        F = (Epot - coulomb_interaction_potential(dr+h))/h .* dr_unit_vec;% force in eV/Angström
    else
        F = F + (Epot - coulomb_interaction_potential(dr+h))/h .* dr_unit_vec;
    end

end





a_vec = F./(mass/u); % in eV/(Angström u)

a_vec = a_vec*9648.53322; % in Angström / ps^2



ax = a_vec(:,1);
ay = a_vec(:,2);
az = a_vec(:,3);
%
% r = [2:0.1:20];
% plot(r, atom_interaction_potential(r))

end

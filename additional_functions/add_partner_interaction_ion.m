function [ax,ay,az, Epot] = add_partner_interaction_ion(x,y,z, mass, charge)

% problem: if starting with correct potential, the particles will not have
% enough energy to break the bond
% (555 m/s) = 0.4 eV, the X potential height is however 1.5 eV 
num_particles = size(x,1);

r1 = [x(1:num_particles/2), y(1:num_particles/2), z(1:num_particles/2)];
r2 =[x(num_particles/2+1:end), y(num_particles/2+1:end), z(num_particles/2+1:end)];
q1  = charge(1:num_particles/2);
q2 = charge(num_particles/2+1:end);

dr_vec = r1-r2;

dr = sqrt(sum(dr_vec.^2,2));

dr_unit_vec = dr_vec./dr;

h= 0.0001;

Epot = ion_interaction_potential(dr, q1,q2);
F = (Epot - ion_interaction_potential(dr+h, q1, q2))/h;% force in eV/Angström



u = 1.66053907e-27; %kg



F = [F; F];

a = F./(mass/u); % in eV/(Angström u)

a = a*9648.53322; % in Angström / ps^2
a = a/1;

% size of a should be size of 2*size(q2)
a_vec = a.*[dr_unit_vec; -dr_unit_vec];

ax = a_vec(:,1);
ay = a_vec(:,2);
az = a_vec(:,3);
% 
% r = [2:0.1:20];
% plot(r, atom_interaction_potential(r))


 Epot = [Epot; Epot]/2; % split potential energy between partners 
end

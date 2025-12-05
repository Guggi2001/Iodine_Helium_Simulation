function [ax,ay,az, Epot] = add_helium_interaction_charged_droplet(x,y,z, mass, R)


global additional_droplet_charges

num_particles = size(x,1);

r1 = [x(1:num_particles/2), y(1:num_particles/2), z(1:num_particles/2)];
r2 =[x(num_particles/2+1:end), y(num_particles/2+1:end), z(num_particles/2+1:end)];

r_vec = [r1; r2];

r = sqrt(sum(r_vec.^2,2));

r_unit_vec = r_vec./r;

% potential energy of test charge at distance r from homogeneously charged
% sphere with radius R
U = @(r, R)14.39964548 * (heaviside(r-R).*1./r  + heaviside(R-r) .* (3./R - r.^2./R.^3)*1/2 );

h= 0.0001;
Epot = U(r, R);
F = additional_droplet_charges*(Epot - U(r+h, R))/h;% force in eV/Angström


u = 1.66053907e-27; %kg



a = F./(mass/u); % in eV/(Angström u)

a = a*9648.53322; % in Angström / ps^2
a = a/1;

a_vec = a.*r_unit_vec;

ax = a_vec(:,1);
ay = a_vec(:,2);
az = a_vec(:,3);
% 
% r = [2:0.1:20];
% plot(r, atom_interaction_potential(r))

end

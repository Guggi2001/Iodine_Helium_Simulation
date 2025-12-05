function [vx,vy,vz] = add_single_ionization_vinit(x,y,z, vx0,vy0,vz0,mass, charge, state_select_id)

num_particles = size(x,1);
num_molecules = size(x,1)/2;

b_charge = charge(1:num_molecules) + charge(num_molecules+1:end)==1;

r1 = [x(1:num_particles/2), y(1:num_particles/2), z(1:num_particles/2)];
r2 =[x(num_particles/2+1:end), y(num_particles/2+1:end), z(num_particles/2+1:end)];
q1  = charge(1:num_particles/2);
q2 = charge(num_particles/2+1:end);

dr_vec = r1-r2;

dr = sqrt(sum(dr_vec.^2,2));

dr_unit_vec = dr_vec./dr;


test = 1;


% potential parameters of four lowerst curves of I_2^+
D_e = [2.7, 2.03,  1.26,   0.56]; %in eV
omega_e = [240,230,  141,   117];  % in cm^1
omega_e_x_e = [0.69,  0.29, 0.32,  0.38];% in cm^1
IP_rel = [0, 0.63, 1.68,  2.44];
R_e = [2.61, 2.61, 2.95, 2.95];
IP_0 = 9.36; % baseline ip of I2, lowest state



hf = 1.54980;
E_photon = hf*8;

r_array = dr;

U = arrayfun( @(i) ...
    morse_potential_I2plus(r_array(i), D_e(state_select_id(i)), omega_e(state_select_id(i)), omega_e_x_e(state_select_id(i)), R_e(state_select_id(i))) ... % the actual function
    - morse_potential_I2plus(2.666, D_e(state_select_id(i)), omega_e(state_select_id(i)), omega_e_x_e(state_select_id(i)), R_e(state_select_id(i))) ...  % subtract value at 2.666
    + IP_0 + IP_rel(state_select_id(i))...    % add reference IP and relative IP for each state
    ,  1:length(state_select_id) );



    E_excess = hf*(8+ (rand(size(U))'>0.5)) -  U';
    
u = 1.66053907e-27; %kg

eV = 1.602e-19; %joule

m1 = mass(1:num_molecules);
m2 = mass(num_molecules+1:end);

b_E = E_excess>0;
E_excess(~b_E) = 0; % dont add velocity if excess energy would be negative
E_excess(~b_charge) = 0; % dont add velocity if sum of charges of both partners is unequal to 1

v1 = sqrt(2*m2./m1 .* E_excess*eV./(m1+m2))/100;
v2 = sqrt(2*m1./m2 .* E_excess*eV./(m1+m2))/100;

dv1 = dr_unit_vec.*v1;
dv2 = -dr_unit_vec.*v2;

vx= vx0;
vy= vy0;
vz = vz0;


vx(1:num_molecules) = vx0(1:num_molecules) + dv1(:,1);
vy(1:num_molecules) = vy0(1:num_molecules) + dv1(:,2); % there was a terrible mistake here as well!!!!!!
vz(1:num_molecules) = vz0(1:num_molecules) + dv1(:,3);

vx(num_molecules+1:end) = vx0(num_molecules+1:end) + dv2(:,1);
vy(num_molecules+1:end) = vy0(num_molecules+1:end) + dv2(:,2); % there was a terrible mistake here!!!!!!
vz(num_molecules+1:end) = vz0(num_molecules+1:end) + dv2(:,3);

end

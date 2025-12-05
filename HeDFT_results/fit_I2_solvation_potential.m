load('I2_contraint_solvation.mat');
E = E - min(E);
k_B = 1.381E-23; % J/K

 r = [r; 100; ];
 E = [E; 574.12; ];



E = E*k_B; %energy in joule 

E = E/1.602E-19; % energy in eV
E = E*1000; % energy in meV

scatter(r, E);
hold on
beta = nlinfit(r, E, @droplet_potential, [11, 60, 20]);
r_new = 0:0.01:100;

plot(r_new, droplet_potential(beta, r_new));

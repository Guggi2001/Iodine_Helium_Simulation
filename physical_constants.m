%% definitions
e_charge = -1.602e-19; %C
epsilon_0 = 8.85418781e-12; %farads  / meter

eV_per_wavenumber = 1/8065.54429;

coulomb_energy = @(x)e_charge.^2./(4*pi*epsilon_0*x*1E-10); % distance has to be in angström

coulomb_velocity = @(x, m) sqrt(coulomb_energy(x)./m);

u = 1.66053907e-27; %kg

eV = 1.602e-19; %joule

hc = 1240; % eV nm


%& hard sphere collision parameters
% mean free path method
bulk_density_helium = 0.0219; % Angström ^-3
density_droplet = 0.8*bulk_density_helium;
%https://journals.aps.org/prb/pdf/10.1103/PhysRevB.58.3341


k_B = 1.380649E-23; %J / K

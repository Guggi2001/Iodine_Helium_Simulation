
eV_per_wavenumber = 1/8065.54429;

data_in = load('iodine_atom_binding_potential_andi_hauser_email');
data_in2 = load('iodine_atom_binding_potential_andi_hauser_email2.mat')

R1 = [data_in.radius; 100];
E1 = data_in.temperature;
E1 = [E1; E1(end) + 250];
E1 = E1 - min(E1);

R2 = [data_in2.R; 100];
E2 = data_in2.E;
E2 = [E2; E2(end)+250];
E2 = E2 - min(E2);

scatter(R1, E1*eV_per_wavenumber);
hold on
scatter(R2, E2*eV_per_wavenumber);

beta = nlinfit(R2, E2*eV_per_wavenumber, @droplet_potential, [20,5,0.03]);
plot(min(R2):0.1:max(R2), droplet_potential(beta, min(R2):0.1:max(R2)));
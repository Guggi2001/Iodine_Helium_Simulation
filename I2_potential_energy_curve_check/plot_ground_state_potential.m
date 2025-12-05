pot = get_morse_potential_X();

figure
r = 2.2:0.01:30;
plot(r, pot(r) - pot(r(end)));
xlabel('r / Angström')
ylabel('E / eV')



tl = tiledlayout(1,2)
nexttile
r = 2.2:0.01:30;
plot(r, pot(r) - pot(r(end)));
xlabel('r / Angström')
ylabel('E / eV')
nexttile

eV_to_kelvin = 11605;
plot(r,(pot(r) - pot(r(end)))*eV_to_kelvin);

xlabel('r / Angström')
ylabel('E / K')

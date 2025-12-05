% goal:plot the selected ionic pec for I2+
% from this paper 
%https://doi.org/10.1063/1.475194

% table V.
D_e = [2.7, 2.03,  1.26,   0.56]; %in eV
omega_e = [240,230,  141,   117];  % in cm^1
omega_e_x_e = [0.69,  0.29, 0.32,  0.38];% in cm^1

IP_rel = [0, 0.63, 1.68,  2.44];

% the R_e is obtained by comparing the assignment and main config. in table
% V to table IV, R_e
% i left out the state with Sigma, because they dont mention an R_e for
% this one

name = {'2Pi_g', '2Pi_g', '2Pi_u','2Pi_u'};

R_e = [2.61, 2.61, 2.95, 2.95];
r = [2.2:0.01:7];
for i = 1:4
    f = get_morse_potential_I2plus(D_e(i), omega_e(i), omega_e_x_e(i), R_e(i));


    plot(r, f(r) - f(2.666) + IP_rel(i) + 9.36);
    hold on
    

end


% plot dissociative photon excitation lines
hf = 1.54980;

for n = [8,9]
    hline(hf*n);
    hold on
end

% determine excess energies
for i=1:4

    f = get_morse_potential_I2plus(D_e(i), omega_e(i), omega_e_x_e(i), R_e(i));

    E_limit =  f(10) - f(2.666) + IP_rel(i) + 9.36;
    E_low = -E_limit+ hf*8
    E_high  = -E_limit+ hf*9

end


v_final = 650;

E_final = v_final.^2/(2)*(127*u) / eV;


v_final_target = 500;

dE_target = v_final_target.^2/(2)*(127*u)/eV

E_bind = E_final - dE_target;

%%

127*u*(40)^2/2/eV*1000


%%
q_ion = 1;
Q_droplet = 0.9408;

E_droplet_coulomb = 1/(4*pi*epsilon_0)*q_ion*Q_droplet*e_charge.^2/(mean(droplet_radii)*1E-10) /eV ;

v_landau = sqrt(2*E_min*eV/(127*u));

v_final = sqrt(v_landau.^2 +    2*eV*E_droplet_coulomb / (127*u)   );



%%
1/(4*pi*epsilon_0)*1/(9E-10)*e_charge.^2/eV


%% if atoms are ejected at v_landau, could droplet ionization at 15 ns cause the ions to have 500 m/s velocity
v_landau = 57;

x = v_landau*15*1E-9;% distance between atom and droplet in meters after 15 ns have passed
Q_droplet = 10;
E_droplet_coulomb = 1/(4*pi*epsilon_0)*1*Q_droplet*e_charge.^2/(x) /eV;


v_final = sqrt(v_landau.^2 +    2*eV*E_droplet_coulomb / (127*u)   );

%==> no, even a tenfold charged droplet would only lead to 170 m/s ion
%velocity at 15 ns


%% velocity expected for coulomb explosion from 9 Angström
E_c_9 = 1/(4*pi*epsilon_0)*1*1*e_charge.^2/(9E-10) /eV;

E_c_18 = 1/(4*pi*epsilon_0)*1*1*e_charge.^2/(18E-10) /eV;

E_c_35 = 1/(4*pi*epsilon_0)*1*1*e_charge.^2/(35E-10) /eV;


v_final = sqrt(eV*(E_c_9 - E_c_35) / (127*u));


E_kin_final = 2 * 127*u*v_final^2/2 / eV;

v_final_actual = 450;
E_kin_actual = 2*127*u*v_final_actual^2/2/eV;


energy_lost_to_droplet = (E_kin_final - E_kin_actual) /E_kin_final;

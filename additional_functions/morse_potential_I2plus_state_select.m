function U= morse_potential_I2plus_state_select(r_array)
%MORSE_POTENTIAL_I2PLUS_STATE_SELECT choose one of 4 ionization states,
%given the state_selet_id array with the same size as r and with entries
%ranging from 1 to 4
global I2plus_state_select_id

state_select_id = I2plus_state_select_id;
molecule_id = 1:length(state_select_id);

% otherwise, use morse potential of the lowest curves of I2+
D_e = [2.7, 2.03,  1.26,   0.56]; %in eV
omega_e = [240,230,  141,   117];  % in cm^1
omega_e_x_e = [0.69,  0.29, 0.32,  0.38];% in cm^1
IP_rel = [0, 0.63, 1.68,  2.44];

R_e = [2.61, 2.61, 2.95, 2.95];

 % randomly choose a I2+ pec for each iodine atom pair
% if this is done here in the potential function, 

IP_0 = 9.36; % baseline ip of I2, lowest state

%r_array = rand(200, 1)'*5 + 2;



% U = arrayfun( @(i) (morse_potential_I2plus(r_array(i), D_e(i), omega_e(i), omega_e_x_e(i), R_e(i)) - ...
%                     morse_potential_I2plus(2.666, D_e(i), omega_e(i), omega_e_x_e(i), R_e(i)) + IP_0 + IP_rel(i)), state_select_id);


U = arrayfun( @(i) ...
    morse_potential_I2plus(r_array(i), D_e(state_select_id(i)), omega_e(state_select_id(i)), omega_e_x_e(state_select_id(i)), R_e(state_select_id(i))) ... % the actual function
    - morse_potential_I2plus(2.666, D_e(state_select_id(i)), omega_e(state_select_id(i)), omega_e_x_e(state_select_id(i)), R_e(state_select_id(i))) ...  % subtract value at 2.666
    + IP_0 + IP_rel(state_select_id(i))...    % add reference IP and relative IP for each state
    ,  molecule_id);


end


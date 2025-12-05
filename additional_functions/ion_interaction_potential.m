function U = ion_interaction_potential(r, q1, q2)

global E_coulomb_scale

if nargin<2
    q1 = ones(size(r));
    q2 = ones(size(r));
end

% if both atoms have charge, use coulomb potential
U = E_coulomb_scale*q1.*q2.*14.39964548./r ; % in eV when r is in angström
% energy for test particle with charge 1!!!

b_single_ionization = q1 + q2 == 1;

global single_charge_ionization_allowed
if single_charge_ionization_allowed
    U_temp = U;

U = U_temp + (b_single_ionization.*morse_potential_I2plus_state_select(r)') ;
end


end


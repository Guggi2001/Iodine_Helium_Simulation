function U = coulomb_interaction_potential(r, q1, q2)


if nargin<2
    q1 = ones(size(r));
    q2 = ones(size(r));
end

% if both atoms have charge, use coulomb potential
U = q1.*q2.*14.39964548./r ; % in eV when r is in angström
% energy for test particle with charge 1!!!


global E_coulomb_scale

U = U*E_coulomb_scale;


end


function U = atom_interaction_potential(r)
%ATOM_INTERACTION_POTENTIAL 
% re = 2.665; % anström
% a = 2;
% De = 1.5; %eV
% 
% r = [1:0.01:6];
% U = De*( exp(-2*a*(r - re)) -2*exp(-a*(r - re)));

morse = get_morse_potential_X();



% plot(r, U);
% hold on
% plot(r, morse(r));

U = morse(r);

end


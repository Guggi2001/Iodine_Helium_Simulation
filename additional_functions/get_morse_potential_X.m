function function_handle = get_morse_potential_X()

global Xdip_active
angstrom_per_bohr_radius = 0.529177211;
eV_per_hartree = 27.2114; % eV per hartree
c = 299792458;
u =1.66054e-27;
h = 6.62607015*1E-34; %J s
eV_per_wavenumber = 1.23984/10000; % ev per cm^-1
joule_per_eV = 1.602E-19;
joule_per_wavenumber = joule_per_eV*eV_per_wavenumber;



% values from: J. Chem. Phys. 107, 9046 (1997); doi: 10.1063/1.475194 
omega_e = 214.5; %cm^-1, 
    omega_e_x_e= 0.65; % cm^-1 

    x_e = omega_e_x_e / omega_e;

    D_e = 1.556; 

    mu = 127/2*u;
    R_e = 2.666;

    a = omega_e*(c*1E2)*2*pi/sqrt(2* ( D_e*joule_per_eV) /mu)*1E-10;


    %
    
eV_to_kelvin = 11605;
    
    morse = @(r) D_e*(1 - exp(-a*(r - R_e))).^2      - 0.9*g(9, 0.3, r)*Xdip_active;


  %  warning('0.6 eV dip in X potential at 9 Angst')

   %plot(2:0.01:20, morse(2:0.01:20))

if false
   R = [1.8:0.01:20];
   plot(R, morse(R));
end

    function_handle = morse;
end

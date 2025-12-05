        % get new kinetic energy after scattering event
        M = 127;
        m = 4;

        M = 40;
        m = 4;

         RHO = M/m ;
        E0 = 1;

        theta = 180 /360*pi;

        COSTHETA = cos(theta);

        E1 = E0.* (1 + 2*RHO.*COSTHETA + RHO.^2)./(1 + RHO).^2;

        relative_energy_loss = (E1 - E0)/E0;


        relative_energy_loss2 = -4*m/M/(1 + m/M)^2;
        
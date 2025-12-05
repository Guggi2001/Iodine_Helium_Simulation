function [res, diagnostic] = f_fast(v, theta, lmax, V, mu)

debug = false;

% calculate crosssection according to 10.1109/NSSMIC.2004.1462242
rmin  =2;
rmax = 300;



u = 1.66053907e-27; %kg
hbar_SI = 1.05457182E-34;
eV = 1.602176487e-19; %joule;
fs = 1E-15;


hbar = 0.6582845318; % eV femtosecond

unitfactor1 = (2*3.8*u*1*eV/(hbar_SI^2) / (1/1E-10)^2 )  /  (2*3.8*1/hbar^2);% convert u/(fs^2 eV) into 1/Angström^2

%https://www.wolframalpha.com/input?i=u+%2F%28eV+fs%5E2%29+%2F+%28+1%2FAngstr%C3%B6m%5E2%29
% proof of unitfactor:
% direct: https://www.wolframalpha.com/input?i=2*3.8+u+*+1+eV+%2F+%28hbar%5E2+%29+%2F+%281%2FAngstr%C3%B6m%5E2%29
% using unitfactor1: https://www.wolframalpha.com/input?i=2*3.8+*+1+%2F+%280.6582%5E2%29+*+103.642


res = zeros(length(v),length(theta));


unitfactor2 = ((mu*u*1000/hbar_SI) /( 1/1E-10) )  / ( mu*10/hbar);%

k_array =  mu*v/hbar * unitfactor2; % in 1/Angström

%ktest = mu*u*1000/hbar_SI/(1/1E-10)
%ktest2 = mu*10/hbar*unitfactor2;

% control: (3.8779 u) * 30 Angström/ps /(planck constant/(2*pi)) /(1/Angström)
% k vector for 10 Angström /ps : mu*10/hbar*unitfactor2 = 1.5746 /Angström

% determine r0 for each k (closest approach)
% rsample = linspace(1,30, 1000);
% for i=1:length(k)
%     K = k(i);
%     lsample = 6;
%     plot(rsample, K^2 - 2*mu*V(rsample)/hbar^2*unitfactor1- (lsample+1/2)^2./rsample.^2);
%     yyaxis right
%     plot(rsample, K^2 - 2*mu*V(rsample)/hbar^2*unitfactor1- (lsample+1/2)^2./rsample.^2>0);
% end

r = linspace(rmin,rmax,4000);

ell = (0:lmax)';

r_array = repmat(r, size(ell,1), 1);

ell_array = repmat(ell,1, size(theta,2));
theta_array = repmat(theta, size(ell,1), 1);
fprintf('calculating legendre array \n');
legendre_array = legendreP(ell_array, cos(theta_array));
fprintf('done, continuing ..\n');

close all



for i=1:length(k_array)
    k = k_array(i);


    % check if k*R_pot >l
    %     r = 3:0.01:20;
    %     plot(r, abs(V(r))/ max(abs(V(r))));
    %     hold on
    %     b_small = abs(V(r))/max(abs(V(r))) < 0.01;
    %
    %     r_select = r(b_small);
    %
    %     plot(r_select, abs(V(r_select))/max(abs(V(r))));
    %     R_char = r_select(1);
    %     lmax = R_char*k;


    % find outer zeros of the temp expression below
    temp = k^2 - 2*mu*V(r)/hbar^2*unitfactor1- (ell+1/2).^2./r.^2;
    %
    %     temp(:,end) = -eps; % put an artificial one at the end
    %
    %     temp_flip = flip(temp,2);
    %     index_flip = flip(1:size(temp,2),2);
    %
    %
    %     indices = arrayfun( @(n) find(temp_flip(n,:)<0, 1,'first'), 1:size(temp,1));
    %     r0_old = r(index_flip(indices))';

    % try alternative way
    tempfun = @(R, l_id) k^2 - 2*mu*V(R)/hbar^2*unitfactor1- (ell(l_id)+1/2).^2./R.^2;

    fzero(@(R) tempfun(R, 1), max(r));

    zero_present = sum(temp>0,2)>0;
    ell_id = 1:size(ell);

    r0 = max(r)*ones(size(ell));

    % find zeros where they are present (this only occurs for some l
    % values)
    zero_found = arrayfun( @(i) fzero(@(R) tempfun(R, i), max(r)), ell_id(zero_present));

    r0(ell_id(zero_present)) = zero_found;

    if debug
        for ii = 1:length(ell)
            plot(r(r>r0(ii)), tempfun(r(r>r0(ii)), ii));
            hold on
        end
    end

    % assign r0
    %r0 = flip( r( -col + length(r) + 1)' ); % need to flip r0 because so that r0(1) corresponds to lowest ell value

    b_r = r>=r0;


    if debug
        for ii = 1:length(ell)
            plot(r(b_r(i,:)), tempfun(r(b_r(i,:)), i));
            hold on
        end
        set(gca,'YScale','log');

    end

    % debug plot in order to figure out if it is working correctly
    % col should have the indices pointing to the r value of the outer zero

    if debug

        figure
        [rr, ll] = meshgrid(r, ell);
        try
            surf(rr, ll, b_r*1, 'edgecolor', 'none');
        catch ME
            test = 0;
        end

        view(0,90);

    end





    V_array = repmat( V(r), size(ell,1), 1);

    eta_integrand = sqrt(  k.^2 - 2*mu*V_array/hbar^2 * unitfactor1 - (ell+1/2).^2./r_array.^2) - k  ;
    if sum(imag(eta_integrand(:).*b_r(:))>0 )
        warning('phase integrand has imaginary parts!');
        if debug
        figure
        for ii=1:size(eta_integrand,1)

            plot(r, real(eta_integrand(ii,:)));
            hold on

        end
        for ii=1:size(eta_integrand,1)

            plot(r, real(eta_integrand(ii,:)).*b_r(ii,:));
            hold on

        end
        figure
        for ii=1:size(eta_integrand,1)

            plot(r, imag(eta_integrand(ii,:)));
            hold on

        end
        figure
        for ii=1:size(eta_integrand,1)

            plot(r, imag(eta_integrand(ii,:)).*b_r(ii,:));
            hold on

        end
        end
    end

    eta_integrand(imag(eta_integrand)>0 ) = 0; % brute force this shit

    eta_integral = trapz(r, eta_integrand.*b_r , 2);



    % eta =  eta_part - eta_integral_prefactor .* eta_integral; % mistake here: there is no multiplicative prefactor!
    eta =  (ell + 1/2)*pi/2 - k*r0+ eta_integral; % this should be the correct expression now

    if debug
        % debug integration interval:
        plot(r, temp(1,:));
        hold on
        plot(r,imag(eta_integrand(1,:)))
        hold on
        plot(r, real(eta_integrand(1,:)));
        plot(r(b_r(1,:)), eta_integrand(1,b_r(1,:)));

        hold on
        vline(r0(1));
    end

    eta_array = repmat(eta, 1, size(theta,2));

    f_array = (2*ell_array+1).*(exp(2*1j*eta_array) - 1).*legendre_array;

    phase_term = (exp(2*1j*eta_array) - 1);

    fsum = sum(f_array, 1);


    disp(i/length(k_array));

    res(i, :) = fsum/(2*1j*k);


    % end
    %
    % diagnostic.r0 = r0;
    % diagnostic.temp = temp(7,:);
    % diagnostic.b_r = b_r(7,:);


    diagnostic = struct();

end


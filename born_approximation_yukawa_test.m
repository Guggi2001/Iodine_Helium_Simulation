%% goal: calculate I+ velocity dependent crosssection of I+ He collision
%% intermediate step: yukawa potential
addpath(fullfile(pwd, 'plot_utility'));
set_groot_properties % set graphics root properties so it looks nicer

% define physical constants
run physical_constants
hbar_SI = 1.05457182E-34;
fs = 1E-15;

m1 = 127; 
m2 = 4;
mu = m1*m2/(m1+m2); %reduced mass

%% define yukawa potential for testing
alpha = 1; %  range parameter in 1/Angström
V0 = 20; % scale of potential in eV

% function definition
V = @(x) V0*exp(-alpha*x)./x;


%% plot the potential
figure
x = 1.9:0.01:25;
plot(x, V(x)*1000);
xlabel('r / Angström'); ylabel('E / meV');


%% direct integral from book
% Quantum Mechanics, Volume 2 Angular Momentum, Spin, and Approximation Methods

rmin = 0; rmax = 50;

% arrays for the sigma calculation
v_array = (1:100:2200);
theta_array = linspace(-pi,pi, 100);

sigma = zeros(length(v_array), length(theta_array));
counter = 0;


%% direct computation of fourier transform for each value of velocity and angle
for i = 1:length(v_array)
    for j = 1:length(theta_array)
    	
        % fetch entries of arrays in each loop
        theta = theta_array(j);
        v = v_array(i);

        k0 = mu*u*v/hbar_SI; % wavevector of incoming wave, B-14, Chapter VIII

        % scattering wavevector
        K = 2*k0.*sin(theta/2); % complement C_Viii, equation (6)

        % prefactor for integral
        prefactor = -2*mu*u/hbar_SI.^2./K; % % complement C_Viii, equation (4)

        %% try different numerical integration methods

        %integrand = (r*1E-10).*sin((r*1E-10).*K).*V(r)*eV;
        %f = trapz(r, integrand)*prefactor;
        % ==> trapz for the radial integral is not good enough

        % better: use integral method which takes a function as the argument
        integrand = @(z) (z*1E-10).*sin((z*1E-10).*K).*V(z)*eV; % complement C_Viii, equation (4)
        f = integral(integrand, rmin, rmax, 'AbsTol',0.00000001)*prefactor;


        % put the obtained value for the crosssection into a matrix
        sigma(i,j) = abs(f).^2; % B-24

        % increase loop counter
        counter = counter+1;

        % display progress
        if mod(counter, ceil(numel(sigma)/20))==0
        fprintf('calculation progress: %.0f percent\n', counter/numel(sigma)*100);
        end

    end
end

% remove nans if there are any
if sum(isnan(sigma(:)))>0
    warning('nan in calculated sigma values, removing..')
    sigma(isnan(sigma)) = 0;
end


% determine total crosssection by integrating over the angles
sigma_total = trapz(theta_array, sigma,2);


%% fit v-dependence
figure

normalization = 1;

plot(v_array, sigma_total/normalization);
hold on 

% perform fit of the total crosssection using power law
fit_function = @(beta, logv) beta(1) + logv*beta(2);

beta = nlinfit(log(v_array)', log(sigma_total), fit_function, [log(max(sigma_total)), -1])
plot(v_array, exp(fit_function(beta, log(v_array))));
leg = ['fit, \sigma = ', replace( sprintf(' %.1e x v^{%.2f}', exp(beta(1)), beta(2) ) , 'x', '\times')];
    
set(gca, 'YScale', 'log');
xlabel('v / m/s')
ylabel('\sigma / m^{-2}')
legend(leg);

%% plot the numerical crosssection against the analytical result for the yukawa potential
figure
plot(v_array, sigma_total)
xlabel('v / m/s')
ylabel('\sigma / m^{-2}')
hold on 

% compare to analytial solution of yukawa potential
k = (mu*u)*v_array/hbar_SI; % need to redefine the k vector here

% calculate all sigma(v, theta) values at once, ensuring that sigma_yukawa
% has the same dimensions as the numerical sigma matrix (for plotting and
% comparing them later on)
sigma_yukawa =  4*(mu*u)^2*(V0*eV)^2/hbar_SI^4     *    1./       (    (alpha*1E10)^2   +      4*(k.^2)' .* sin(theta_array/2).^2      ).^2;
% note that the alpha value has to be multipled by 1E10 to ensure it has SI
% units of 1/m instead of 1/Angström

% calculate the total crosssection for all angles by integrating, as before
sigma_total_yukawa_numeric = trapz(theta_array, sigma_yukawa,2);

% plot the analytical solution
plot(v_array, sigma_total_yukawa_numeric , ':');

legend('numerical', 'analytical');


%% plot the sigma(v, theta) surfaces for the numerical sigma, the analytical sigma and the difference between the two

figure
colormap(colorcet('L07')); % colormaps, for available options see: https://colorcet.com/gallery.html#linear

%tl = tiledlayout(1,3);

% nexttile

surface(v_array', theta_array', sigma');

xlabel('v / m/s'); ylabel('\theta / radian');
cb = colorbar;
cb.Label.String = '\sigma / m^{-2}';

title('numerical solution')

% nexttile

surface(v_array', theta_array', sigma_yukawa' );

xlabel('v / m/s'); ylabel('\theta / radian');
cb = colorbar;
cb.Label.String = '\sigma / m^{-2}';

title('analytial solution')

% nexttile

surface(v_array', theta_array', sigma_yukawa' - sigma');

title('difference')

xlabel('v / m/s'); ylabel('\theta / radian');


f = gcf;
f.Position = [0.4130    0.5418    1.1520    0.4200]*1000;


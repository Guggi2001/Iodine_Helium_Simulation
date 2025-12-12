%% goal: calculate I+ velocity dependent crosssection of I+ He collision
%% intermediate step: yukawa potential

set_groot_properties % set graphics root properties so it looks nicer

% define physical constants
run physical_constants
hbar_SI = 1.05457182E-34;
fs = 1E-15;

m1 = 127; 
m2 = 4;
mu = m1*m2/(m1+m2); %reduced mass

%% define  potential
%I+ He potential from A. A. Buchachenko, T. V. Tscherbul, J. Klos, M. M.
%SzczÂ¸eÂ´sniak G. ChalasiÂ´nsk, R. Webb, and L.A. Viehland
%J. Chem. Phys. 122, 194311 (2005)

% type in the parameters from the EPAPS.txt repository file

alpha = 0.5132894e+1;
beta = -0.8349857e+1;
g = [0.8833500e+7, -0.1624297e+8,  0.1156310e+8,  -0.3647464e+7,  0.2590770e+6, 0.1324182e+6,  -0.2782933e+5 ,  -0.1278957e+2,   0.2773608e+3];

delta = 0.7159647e+3;
D4 =  0.1185215e+5;
D6 =  0.1117964e+6;
D8 = 0.7710727e+6;

eV_per_wavenumber = 1/8065.54429;

l = 0:8;
x = 1.9:0.01:25;

% short range part
VSR = @(R) sum( g(l+1).* (R').^l .* exp(-alpha*R' - beta), 2)';

% long range part
VLR = @(R) -D4*R.^(-4) - D6*R.^(-6) - D8*R.^(-8) ;

% switch function
sw = @(R) 1/2*(1+ tanh(1 + delta*R));

%% plot components of the potential
%  figure
%  plot(x, (VLR(x).*sw(x) + VSR(x)));
%  hold on
%  plot(x, VSR(x));
%  plot(x, VLR(x));  
%  plot(x, VLR(x)*1 + VSR(x), ':');
% xlabel(['R / ', char(0197)]);
% ylabel('E / eV');
% 
% legend('full version', 'short range part', 'long range part', 'full version, switch omitted');

% define the I+He potential, the switch function can be omitted because it
% does nothing
V_original= @(R) (VLR(R)*1 + VSR(R))* eV_per_wavenumber ;


% define a soft step function, which is 0 below R0 and 1 above R0
softstep = @(R, R0) (erf((R- R0)/0.01) + 1)/2;


% create modified potential which is well behaved below 2 AngstrÃ¶m
R_switch = 2;
V = @(R) V_original(R) .*softstep(R, R_switch) + V_original(R_switch)* (1-softstep(R, R_switch));



%% plot the potentials (modified and original)
figure
x = 0.01:0.01:25;
plot(x, V(x)*1000);
hold on
plot(x, V_original(x)*1000)
% ylim([min(V(x)*1000), max(V(x)*1000)])
ylim([0, 1e6])


%% Lennard-Jones Regularized Potentials 
% Parameters
eps = .01784;            % Depth of the potential well
sig = 3.25/2^(1/6);              % Characteristic length scale
regpar = 0.05 * sig;       % Soft-core regularization parameter
rc = 0.05 * sig;       % Cutoff radius for exponential damping
n = 8;                  % Damping exponent

% Distance range
r = linspace(0.5*sig, 3*sig, 1000);

% --- Lennard-Jones Potential
V_LJ = @(b, r) 4*b(1)*[(b(2)./r).^12 - (b(2)./r).^6];

% --- Soft-Core Regularization
V_soft = @(b, r) 4 * b(1) * ( ((b(2)^2 + b(3)^2) ./ (r.^2 + b(3)^2)).^6 - ((b(2)^2 + b(3)^2) ./ (r.^2 + b(3)^2)).^3 );

% --- Exponential Damping Regularization
damping = 1 - exp( - (r ./ rc).^n );
V_exp = @(b,r) V_LJ([b(1), b(2)], r) .* (1 - exp( - (r ./ b(3)).^b(4) ));

plot(x, V_LJ([eps, sig], x)*1000);
plot(x, V_soft([eps, sig, regpar], x)*1000);
plot(x, V_exp([eps, sig, rc, n], x)*1000);

xlim([0 5]);
ylim([-20 100]);
legend('clipped', 'original', 'Lennard-Jones','soft-core', 'exp cutoff')
xlabel('r / Angström'); ylabel('E / meV');

close gcf

%% fit Lennard Jones potential
% x_LJ = x(x>2.5);
% beta = nlinfit(x_LJ, V_original(x_LJ), V_LJ, [.018, 3]);

%% try out multiple R_switch values by looping over possible choices, and redefining the potential each time
figure
ax = axes;

% R_switch_values =0.8:0.2:2.8;
R_switch_values = 1;

leg = {}; % cell array to hold legend entries, which are filled during each loop
for id =1:length(R_switch_values)
R_switch = R_switch_values(id);

leg{end+1} = ['R_{switch} = ', sprintf('%.1f', R_switch)];

V = @(R) V_original(R) .*softstep(R, R_switch) + V_original(R_switch)* (1-softstep(R, R_switch));


%% direct integral from book
% Quantum Mechanics, Volume 2 Angular Momentum, Spin, and Approximation Methods

rmin = 0; rmax = 1e12;

% arrays for the sig calculation
v_array = (10:100:2200);
theta_array = linspace(-pi,pi, 100);
sigma = zeros(length(v_array), length(theta_array));
counter = 0;


%% direct computation of fourier transform for each value of velocity and
% angle
for i = 1:length(v_array)
    for j = 1:length(theta_array)
    	
        % fetch entries of arrays in each loop
        theta = theta_array(j);
        v = v_array(i);

        k0 = (mu*u)*v/hbar_SI; % wavevector of incoming wave, B-14, Chapter VIII

        % scattering wavevector
        K = 2*k0.*sin(theta/2); % complement C_Viii, equation (6)

        % prefactor for integral
        prefactor = -2*mu*u/hbar_SI.^2./K; % % complement C_Viii, equation (4)

        %% try different numberical integration methods

        %integrand = (r*1E-10).*sin((r*1E-10).*K).*V(r)*eV;
        %f = trapz(r, integrand)*prefactor;
        % ==> trapz for the radial integral is not good enough

        % better: use integral method which takes a function as the argument
%         integrand = @(z) (z*1E-10).*sin((z*1E-10).*K).*V(z)*eV; % complement C_Viii, equation (4)
         integrand = @(z) (z*1E-10).*sin((z*1E-10).*K).*V_exp([eps, sig, rc, n],z)*eV;
        f = integral(integrand, rmin, rmax, 'AbsTol',0.00000001)*prefactor;

        % V_LJ([eps, sig];
        % V_soft([eps, sig, regpar]
        % V_exp([eps, sig, rc, n]


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


%% plot the numerical crosssection against the analytical result for the yukawa potential
colors = colormap(colorcet('L08', 'N', length(R_switch_values)+1));

normalization = trapz(v_array, sigma_total);
normalization = 1;

plot(ax, v_array, sigma_total/normalization , 'color', colors(id,:));
hold on 

% perform fit of the total crosssection using power law
fit_function = @(beta, logv) beta(1) + logv*beta(2);

if numel(R_switch_values) ==1

%     beta = nlinfit(log(v_array)', log(sigma_total), fit_function, [log(max(sigma_total)), -1])
%     plot(v_array, exp(fit_function(beta, log(v_array))));
%     leg{end+1} = ['fit, \sigma = ', replace( sprintf(' %.1e x v^{%.2f}', exp(beta(1)), beta(2) ) , 'x', '\times')];
    
end


end

set(gca, 'YScale', 'log');

xlabel('v / m/s')
ylabel('\sigma / m^{-2}')

legend(leg);

%% plot the sigma(v, theta) surfaces for the numerical sigma

figure
colormap(colorcet('L07')); % colormaps, for available options see: https://colorcet.com/gallery.html#linear

% tl = tiledlayout(1,1);

% nexttile

surface(v_array', theta_array', sigma');

xlabel('v / m/s'); ylabel('\theta / radian');
cb = colorbar;
cb.Label.String = '\sigma / m^{-2}';
view(45,45);

close gcf

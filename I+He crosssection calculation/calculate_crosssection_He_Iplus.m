%% goal: calculate I+ velocity dependent crosssection of I+ He collision
% using the lippmann schwinger equation
run physical_constants

%I+ He potential from A. A. Buchachenko, T. V. Tscherbul, J. Klos, M. M.
%Szcz¸e´sniak G. Chalasi´nsk, R. Webb, and L.A. Viehland
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
x = 1.9:0.1:25;

VSR = @(R) sum( g(l+1).* (R').^l .* exp(-alpha*R' - beta), 2)';
test = VSR(x);

VLR = @(R) -D4*R.^(-4) - D6*R.^(-6) - D8*R.^(-8) ;


sw = @(R) 1/2*(1+ tanh(1 + delta*R));

figure
plot(x, (VLR(x).*sw(x) + VSR(x)));
hold on
plot(x, VSR(x));
plot(x, VLR(x));  
plot(x, VLR(x)*1 + VSR(x));


V = @(R) (VLR(R)*1 + VSR(R))* eV_per_wavenumber;
mu = 127*4/(127 + 4); % in u

% use Xe potential for debugging the code
% C4 = 29.122; % in eV Angström
% C6 = 198.43; % in eV Angström
% 
% C = [0,0,0,C4,0,C6];
% 
% b = 2.7526;
% A = 5037.96;
% B = 107.756;
% 
% f2 = @(n, R)( 1 - exp(-b*R').*sum( ((b*R').^[0:2*n])./factorial([0:2*n]),2))';
% 
% V_ref = @(R) A*exp(-b*R) - B*exp(-b*R/2) - f2(2,R)*C(4).*R.^(-4) -  f2(3,R)*C(6).*R.^(-6);
% 
% V = V_ref;
% 
% mu = 131*2/(131 + 131); % in u Xe Xe

figure
plot(x, V(x)*1000);
xlabel('r / Angström'); ylabel('E / meV');

xmin = x(find(V(x)==min(V(x)),1));

v = 1:1:30;

xchar = x(   find(   abs(V(x)) < 0.01*abs(V(xmin))  ,1 )   );
hbar_SI = 1.05457182E-34;
hbar = 0.6582845318; % eV femtosecond
unitfactor2 = ((mu*u*1000/hbar_SI) /( 1/1E-10) )  / ( mu*10/hbar);%
k_array =  mu*v/hbar * unitfactor2; % in 1/Angström

lmax_from_xchar = xchar*max(k_array);

%v = [100:50:3000] / 100; % velocity in angström / ps


E = mu*u * (v*100).^2 / 2/eV;

%E = logspace(-3,1, 50);

%v = sqrt(2*E*eV/(mu*u))/100;

%v = 4;
%v = [1];
%v = 1;
theta = linspace(0,pi,90);
%theta = [0.05];

leg = {};
res_cell = {};

% code fails at 190 due to dimansion mismatch in some boolean matrices related to b_r, i hope it converges before
% that
lmax_array = [1, 5, 10. 50, 100, 150, 180];
for lmax = lmax_array

[res, diagnostic_fast] = f_fast(v, theta,lmax, V, mu);
res_cell{end+1} = res;

leg{end+1} = sprintf('l_{max} = %.0f', lmax);


end

cmap = colorcet('L09', 'N',length(res_cell)*2);

%%

figure
tiledlayout(1,2);
ax1 = nexttile;
set(ax1,'Yscale', 'log');
ylabel(['Q_s / 1E-17 m^2'])
xlabel(['E / eV']);
hold on
set(gca, 'XScale', 'log');

ax2 = nexttile;
set(ax2,'Yscale', 'log');
xlabel('\Theta / radian');
ylabel(['Q_s / 1E-17 m^2']);
hold on


% [vv, thetatheta] = meshgrid(v, theta);
% surf(vv, thetatheta, log(Is'));
%surf(vv, thetatheta, log(abs((sin(theta).*Is)')));


for i=1:length(res_cell)
    res = res_cell{i};

Is = res.*conj(res)/4; % assume there is only one contribution (gerade or ungerade)



Qs = 2*pi*trapz(theta, sin(theta).*Is,2);

Qs = Qs * 1E-20 / (1E-17); 
scatter(ax1, E, Qs,  'markeredgecolor', cmap(i,:),'markerfacecolor', cmap(i,:))
scatter(ax2, theta, trapz(E, Is,1), 'markeredgecolor', cmap(i,:), 'markerfacecolor', cmap(i,:));
end
legend(ax1, leg);


figure
tiledlayout(1,2);
ax1 = nexttile;
set(ax1,'Yscale', 'log');
ylabel(['Q_S /', char(0197), '^2'])
xlabel(['v / ',char(0197), '/ps']);
hold on

ax2 = nexttile;
set(ax2,'Yscale', 'log');
xlabel('\Theta / radian');
ylabel(['Q_s /' char(0197), '^2']);
hold on


% [vv, thetatheta] = meshgrid(v, theta);
% surf(vv, thetatheta, log(Is'));
%surf(vv, thetatheta, log(abs((sin(theta).*Is)')));
add_letter_norm('a', 'topleft', 20, ax1);
add_letter_norm('b',  'topleft', 20, ax2);


fit_parameters = [];
fun = @(beta, logv) -beta(1)*logv + beta(2);




for i=1:length(res_cell)
    res = res_cell{i};

Is = res.*conj(res)/4; % assume there is only one contribution (gerade or ungerade)



Qs = 2*pi*trapz(theta, sin(theta).*Is,2);

scatter(ax1, v, Qs,  'markeredgecolor', cmap(i,:),'markerfacecolor', cmap(i,:))

beta = nlinfit(log(v), log(Qs)', fun, [100, 1]);

fit_parameters = [fit_parameters; beta];

plot(ax1, v, exp(fun(beta, log(v))),'color', cmap(i,:), 'HandleVisibility','off');


scatter(ax2, theta, trapz(v, Is,1), 'markeredgecolor', cmap(i,:), 'markerfacecolor', cmap(i,:));
end
legend(ax1, leg);


exportgraphics(gcf, 'crosssection_IplusHe.pdf', 'ContentType','vector');

figure
% plot fit parameters for different lmax
plot(lmax_array, fit_parameters(:,1))

xlabel('l_{max}');
ylabel('k / (unitless)');
title(sprintf('mean(k) : %.1f', mean(fit_parameters(:,1))));

yyaxis right
plot(lmax_array, fit_parameters(:,2));
ylabel('')


figure
% plot v, sigma(v) for the mean slop k for the function log(sigma) = -k log(v) + d 
sigma0 = 300;
v_test = 0:0.1:30;
sigma_test = exp( fun([2, 0], log(v_test)));
plot(v_test, sigma_test);

%sigma_lookup = @(v_test) interp1(v, Qs,v_test,"linear" );

mean_exponent = mean(fit_parameters(:,1));

sigma_lookup = @(v_test) exp(-2*log(v_test))

hold on
plot(v_test, sigma_lookup(v_test), '--');

save('crosssection_v_dependence', "sigma_lookup");

% PROBLEM: lmax determines the absolute value of sigma. how to determine
% lmax? 


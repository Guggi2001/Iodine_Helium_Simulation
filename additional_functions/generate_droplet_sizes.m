function N = generate_droplet_sizes(num_molecules,reduced_crosssection, debug_pickup_plot)


if nargin <1
    num_molecules = 500;
    debug_pickup_plot = true;
    reduced_crosssection = false;
end


% see lackner diss for formulas
global p_source
global T_source
p0 = p_source;
T0 = T_source;

%
%mu = (0.00184*p0 - 0.331)*T0 + (-0.0133*pi + 13.3747);


%N_median = 6250*2;
%N_mean =7700*2;

N_mean = get_dropletsize(p0, T0);
delta = 0.625; %use constant delta according to kornilo2009, see lackner diss

mu = log(N_mean)- delta^2/2;

%mu = log(N_median);

%delta = sqrt(2*(log(N_mean) - mu));

%FWHM_N = exp(mu - delta^2 + delta*sqrt(2*log(2))) - exp(mu - delta^2 - delta*sqrt(2*log(2)));

N = lognrnd(mu, delta, 1,num_molecules*10);


%samples_initial = N*2; %why the fuck is there a two here, remove this
samples_initial = N;






%% TD: simulate pickup of I2 from this distribution

num_colors =60;
cmap_pickup = colormap('jet');
cmap_pickup  = cmap_pickup(1:floor(length(cmap_pickup )/num_colors):length(cmap_pickup ), :);
bar_colors = cmap_pickup';


samples = samples_initial;

if size(samples,1)<size(samples,2)
    samples = samples';
end

% => get how size distribution of droplets with one I2 look like
% => decrease they size by E_thermal = 5/2 kT 
n_he = 2.18E28; % density of liquid helium in particles per m^3

n_droplet = n_he*0.8;

% droplet radius from number of hleium atoms in droplet
R_droplet = @(N)(3*N/(4*pi*n_droplet)).^(1/3);


% incorporate kinetic energy dissipation dependent crosssection
l = 4*1E-10; % mean free path
epsilon = 0.04; % relative energy loss per collision

E_solv = 14; % meV
E_kin_0 = 25.85 + 14; % meV of molecule at 300 K

%kappa = epsilon/l;

b_thresh =  @(N) real( sqrt(R_droplet(N).^2 - l^2 / (4*epsilon^2)*log(E_solv/E_kin_0)^2) );
%histogram(b_thresh(N));
%hold on
%histogram(R_droplet(N));
if reduced_crosssection
sigma_droplet = @(N)pi*(b_thresh(N)).^2; % crosssection of helium droplet that can pick up the dopant
else
sigma_droplet = @(N)pi*(R_droplet(N)).^2; % crosssection of helium droplet
end

mbar = 100; % mbar in pascal
eV_to_K = 11604.5250061657; 
kb = 1.381E-23; % boltzmann constant in J / K

T = 275;
pressure_puc = 0.5E-5;
n_gas = pressure_puc*mbar/(T*kb); 
a = 2E-2; % length of pickup region in m

p_pickup_event = @(sigma_droplet,a)(a.*n_gas.*sigma_droplet);




max_pickup = 1; % maximum number of pickup processes


%mu =@(N) -7.15398 + 11.53557*N.^(-1/3) + 1.95583*N.^(-2/3); % chemical potential in kelvin
mu = @(N) -7.21 + 17.71*N.^(-1/3) - 5.95*N.^(-2/3); % 10.1007/978-3-030-94896-2_1

num_evap = @(k, N) T./mu(N) + (k>1).*  (T*5/2)  ./mu(N);



% variables that change every loop
total_destroyed = 0; % number of droplets that have been destroyed
total_pickup = zeros(length(samples),1); % list that save how many atoms each droplet has picked up, where destroyed droplets and finished droplets are not included

samples_completed = []; % list that contains droplet sizes of droplets that have passed the pickup region 
total_pickup_completed = []; % list that contains how many atoms those droplets have picked up

id = 1; % running counter variable

%distance_traveled = zeros(length(samples),1);


leg = {};
while ~isempty(samples) && max(total_pickup)<max_pickup
    
    % !! the quantities that are used inside this loop are vectors with the
    % same length as samples
    
    
    % calculate pickup probability
    pickup_probabilities = p_pickup_event(a, sigma_droplet(samples));
    
   

    % boolean values that are true when pickup occurs
    pickup_bool =  rand(length(samples) ,1) < pickup_probabilities ;

    % increment total_pickup
    total_pickup = total_pickup + pickup_bool;
    
    % calculate how many helium atoms are evaporated for each droplet
    evap =  num_evap(total_pickup,samples).*pickup_bool;
    % calculate how many helium atoms are evaporated for each droplet
%     evap = zeros(length(samples),1);
%     for k=1:20
%         evap = evap + num_evap(k, samples).*(total_pickup==k);
%     end
%     

if debug_pickup_plot
     fprintf('mean pickup probability: %.2f \n', mean(pickup_probabilities));
    fprintf('mean evaporated He atoms: %.0f \n', mean(evap))
end

    % decrement droplet sizes in samples by the amount evaporated
    new_samples = samples + evap;
    
    % boolean vector that contains which droplets have been destroyed
    destroyed = new_samples<=0;
    
    % remove destroyed droplets from vectors
    new_samples(destroyed)= [];
    total_pickup(destroyed) = [];
   % distance_traveled(destroyed) = [];
    
    % increase total number of destroyed droplets
    total_destroyed = total_destroyed + sum(destroyed);
    
    % boolean list that is true for droplets that did not have pickup
    no_pickup = ( total_pickup == id - 1 );
    
    % add droplet sizes and number of atoms picked up to completed list
    total_pickup_completed = [total_pickup_completed; total_pickup(no_pickup)];
    samples_completed = [samples_completed; new_samples(no_pickup)];
    
    % remove entries of droplets that have passed the pickup cell (= did
    % not have pickup this round)
    new_samples(no_pickup) = [];
    total_pickup(no_pickup) = [];
    %distance_traveled(no_pickup) = [];
    
    % if droplet has picked atom up and is not destroyed, let it travel on
    % to a random location in the remaining cell
    % dont know if this is physical or not?
    %distance_traveled = distance_traveled + (a - distance_traveled).*rand(length(distance_traveled),1);
    
    
    % optional plotting of new_samples distribution
    if debug_pickup_plot
        hold on
        [h, sigma_h, centers, binwidth, barplot, errorplot] = bayes_hist(new_samples, [1,max(N) ], true, bar_colors(:,1+mod(id-1,length(bar_colors))));
        plot(centers, h, 'linewidth', 1.2, 'color', bar_colors(:,1+mod(id-1,length(bar_colors))));
    
    end
    leg{end+1} = sprintf(' %.0f pickup round', id);
    
    
    samples =new_samples;
    id = id+1;
    
end
test = 1;

samples_completed = [samples_completed; samples];
total_pickup_completed = [total_pickup_completed; total_pickup];


if max(total_pickup)<max_pickup
    disp('max_pickup reached');
end


if debug_pickup_plot
    legend(leg);
    xlabel('N');
    ylabel('p(N)');
    title('droplet size distributions after pickup events')
end



samples_1_pickup = samples_completed(total_pickup_completed==1);
figure
histogram(samples_initial);
hold on
histogram(samples_1_pickup);
legend('initial droplet sizes', 'droplet sizes with 1 molecule picked up');


%% after finishing the pickup simulation, pickup the correct amount of droplets
if length(samples_1_pickup)<num_molecules
    error('not enough droplets with 1 molecule, increase number of samples droplets!');
end

N = samples_1_pickup(1:num_molecules);

end
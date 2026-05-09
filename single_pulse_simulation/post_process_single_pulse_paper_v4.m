close all
run physical_constants.m
%% fat todo list:
% refactor and clean up this code, because it is horrible to work with
% use high s/n ratio measurements from this october 
% write code to view ion trajectories
% check if ion deceleration is similar to HeDFT result for I2 diss. in
% droplet


% use the small droplet power scan measurement from 16.10.24
%Probe only I+ He and gas # 43554 

%Probe only I+ gas 43555

%Probe only Ihe+  43556

%Probe power / mW 160

%frames 2500

%      



gas_phase_reference_measurement = 43555;%160 mW
I_plus_from_drop_reference_measurement = 43554; % 160 mW

gas_phase_reference_measurement = 43568;%600 mW
I_plus_from_drop_reference_measurement = 43567; % 600 mW

I_plus_He_from_drop_reference_measurement = 43556;

I_plus_He_from_drop_reference_measurement = 43563; % 300 mW

if effusive_dynamics


    global plot_processed_with_ROI
    plot_processed_with_ROI = false;

    % get gas phase I+ reference measurement
    res = plot_processed_VMI(gas_phase_reference_measurement  , 1, [482.9299  392.4866], 0);
    %res = plot_processed_VMI(gas_phase_reference_measurement  , 1, [477.7867  410.8572], 0);
    
    close all
    global velocity_factor


        vf = 8.6178; % for 7.95 kV repeller, new calibr.!
    
    leg_fig1 = {};

    fig1 = figure(1)
    tl = tiledlayout(2,1);
    tile1 = nexttile;

    plot(res.r*vf ,res.radial_distribution/max(res.radial_distribution))
    hold on
    leg_fig1{end+1} = 'm. 1: singlepulse, gas';

    % get radial distribution in limited angular range
   b_phi = (res.phi >0 & res.phi <pi/3) | (res.phi>5/6*2*pi & res.phi<2*pi);
    v =res.r*vf;
   y = mean(res.image_polar(b_phi,:),1);
   plot(v, y, ':')
    leg_fig1{end+1} = 'm. 1: singlepulse, gas, lim.angle';

    fun = @(beta, v) g(beta(1), beta(2), v).*beta(3) + g(beta(4), beta(5), v).*beta(6);
    beta = nlinfit(v, y, fun, [814, 1000, 1, 2000, 100, 1]);
    plot(v, fun(beta, v));
    leg_fig1{end+1} = 'double g fit';

    plot(v, beta(6)*g(beta(4), beta(5), v));
    leg_fig1{end+1} = 'high vel. g fit';

    b_r = res.r>250 & res.r<400;

    
    %leg_fig1 = { 'm. 3: p.p. at t<0','simulation'}


   tile2 = nexttile;
    y2 = mean(res.image_polar(:, b_r),2);
    y2 = movmean(y2,4);
    %y2 = y2 - min(y2);
    y2 = y2/max(y2);


    plot(res.phi,y2 );
    fun2 = @(beta, phi) [real(beta(1)*cos(phi).^(beta(2))); imag(real(beta(1)*cos(phi).^(beta(2))))] ;
    beta2 = nlinfit(res.phi', [y2; 0*y2], fun2, [1,2]);
    hold on
    fitresult = fun2(beta2, res.phi')
    plot(res.phi, fitresult(1:length(y2)));


    figure(3)
    hold on
    plot(v/100, g(beta(4), beta(5), v));

    figure(4)
    [rr, phiphi] = meshgrid(res.r, res.phi);
    image_polar_select = res.image_polar;
    image_polar_select(~b_phi,:) = 0;
    surf(rr, phiphi, image_polar_select);


else % comparison IHe+
    

        global plot_processed_with_ROI
    plot_processed_with_ROI = false;

    % from 26.8.24, 600 mW probe, 40 bar, 14 K
    %res = plot_processed_VMI(32910, 0, [528.4614  382.9408],0);
   % close all
    
   % load vmi from filenumbers
    res_Iplus_from_He = plot_processed_VMI(I_plus_from_drop_reference_measurement, true,[524.5297  380.8430], true);
    res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement , true, [524.5297  380.8430], true);

    res_Iplus_from_He = subtract_processed_data(res_Iplus_from_He, res_Iplus_gas);
    
    VM_center_Iplus_He = [ 509.3664  387.6409];

    res_Iplus_He = plot_processed_VMI(I_plus_He_from_drop_reference_measurement, true,VM_center_Iplus_He, true);
    
    % alternatively, load vmi from high snr stored mat files

    %data_in =  load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He');
    %res_Iplus_He = data_in.res_sum;

    global velocity_factor

    % load timescan probe only for comparison
    vf = velocity_factor;
    vf_timescan =  6.4149*0.9675;
    global GS_bleach_correction
    GS_bleach_correction = 1;
    res2 = mean_timescan_2d_VMI([296:297], false,[524.5297  380.8430], false);

    tiledlayout(1,1);
    tile1 = nexttile;

   % fig1 = figure(1)
    %plot(res.r*vf ,res.radial_distribution/max(res.radial_distribution))
     hold on

    % v =res.r*vf;
    % y = mean(res.image_polar(b_phi,:),1);
    % plot(v, y/max(y), ':')
    % 
    % fit 1 gaussian
    % fun = @(beta, v) g(beta(1), beta(2), v).*beta(3);% + g(beta(4), beta(5), v).*beta(6);
    % beta = nlinfit(v, y, fun, [814, 1000, 1]);
    % plot(v, fun(beta, v)/max(fun(beta, v)));
    % %plot(v, beta(3)*g(beta(1), beta(2), v));

    %b_r = res.r>150 & res.r<280;
    xlabel('v / m/s');
    ylabel('signal / arb. units');

    %y = sum(res2.data_probe_only(:,:),2);
    y =  res2.data(:, res2.t>150);
    y = y/max(y);

    %plot(res2.r*vf_timescan,y);

    %plot(res2.r*vf_timescan, y/max(y));
    vf_single =8.6178;

    plot(res_Iplus_He.r*vf_single, res_Iplus_He.radial_distribution / max(res_Iplus_He.radial_distribution ));
    %plot(res_Iplus_from_He.r*vf_single, res_Iplus_from_He.radial_distribution/max(res_Iplus_from_He.radial_distribution));

  leg_fig1 = { 'I_2:I^+He'};
    

 

end

linestyle = {'--', ':', '-.'};

%% plot
% load('T:\github synchronized\Iodine_Helium_Simulation\\neutral_propagation_checkpoint.mat');
% load('T:\github synchronized\Iodine_Helium_Simulation\ion_propagation_checkpoint.mat');

data_neutral = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');

if effusive_dynamics
data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint_gas.mat');
else
data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
end

if single_initial_position
    data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint_hedft.mat');
else
    data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
end

%for mass_select =127 + 4*[0:17]
line_id = 1;

for mass_select = 127+ 4*[0,1]
b_mass_select = round(data_ion.mass_i(:,end)/u)==mass_select;
%b_mass_select = true(size(b_mass_select));


fprintf('mass selection accounts for %.2f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
b_select = b_mass_select & data_ion.b_ion_outside;


vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);



% figure
% scatter3(vx_total(:,1), vy_total(:,1),vz_total(:,1));


v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);


axes(tile1);
edges_velocity = [0:0.05:35];
plot_samples = v_total_proj;


xinterval = [min(edges_velocity), max(edges_velocity)];
[h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
%histogram_data_velocity(1, :) = movmean(h,20);
vd_ion = movmean(h,40);
%vd_ion = vd_ion -min(vd_ion);
vd_ion = vd_ion - min(vd_ion);

hold on
plot(centers_velocity*100, vd_ion/max(vd_ion), 'linestyle', linestyle{line_id}, 'linewidth', 1.5, 'color', [0.1,0.1,0.1]);
leg_fig1{end+1} = sprintf('simulated, mass = %.0f u', mass_select);


line_id = line_id+1;





end
axes(tile1)
legend(leg_fig1)

f = gcf;
f.Position(3:4) = figsize*1;
f.Position(2) = f.Position(2) - figsize(2);
%add_letter_norm('a','topleft', 18, tile1, 0);
%add_letter_norm('b','topleft', 18, tile2, 0);


exportgraphics(f, 'compare_simulation_and_measurement_simpler.pdf', 'ContentType','vector');


figure

vr_total = sqrt(data_ion.vx_total.^2 + data_ion.vy_total.^2);
theta_total = atan2(data_ion.vx_total, data_ion.vy_total);

theta1 = theta_total(1:num_molecules) + pi;
theta2 = theta_total(num_molecules+1:end) + pi;

numbins = 90;
theta = linspace(0,2*pi, numbins);
dtheta = theta(2)- theta(1);

theta1 = theta1(b_select(1:num_molecules)  & b_select(num_molecules+1:end));
theta2 = theta2(b_select(1:num_molecules)  & b_select(num_molecules+1:end));

simulated_angular_covariance = zeros(length(theta), length(theta));

for i=1:length(theta1)
   
        simulated_angular_covariance( floor(theta1(i)/dtheta)+1 , floor(theta2(i)/dtheta)+1) = simulated_angular_covariance( floor(theta1(i)/dtheta)+1 , floor(theta2(i)/dtheta)+1) + 1;


end

surf(theta-pi, theta-pi, simulated_angular_covariance);
xlabel('\theta / radian'); ylabel('\theta / radian');
colormap(colorcet('L08')); cb = colorbar; cb.Label.String = 'counts';
view(0,90);
scatter(theta1, theta2, 'x');
pbaspect([1,1,1]);



figure
title('ion mass histrogram')
histogram(data_ion.mass_i(:,end)/u);
xlabel('m / u')

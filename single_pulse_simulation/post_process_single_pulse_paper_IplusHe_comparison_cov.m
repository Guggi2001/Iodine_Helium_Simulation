close all
run physical_constants.m

effusive_dynamics = false;
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

gas_phase_reference_measurement = 43562; % 300 mW

gas_phase_reference_measurement = 43562; % 300 mW

gas_phase_reference_measurement = 43632; %300 mW

%gas_phase_reference_measurement = 43568;%600 mW

I_plus_from_drop_reference_measurement = 43567; % 600 mW

I_plus_He_from_drop_reference_measurement = 43556;


I_plus_He_from_drop_reference_measurement = 43563; % 300 mW 16.10.24

I_plus_He_from_drop_reference_measurement = 43630; % 300 mW 17.10.24


I_plus_He_from_drop_reference_measurements = [45668, 45662, 45667, 45686]; % I+He measurements, low doping, 300 mW, 3.12.24

I_plus_He_from_drop_reference_measurements = [45668, 45662, 45667]; % I+He measurements, low doping, 300 mW, 3.12.24


%I_plus_He_from_drop_reference_measurement = 43569; % 600 mW

f1 = figure;

tl = tiledlayout(3,2, 'TileSpacing','tight', 'Padding','compact');
tile1 = nexttile(1, [1,1]);
tile2 = nexttile(2, [1,1]);
tile3 = nexttile(3, [1,2]);

tile4 = nexttile(5, [1,1]);
tile5 = nexttile(6, [1,1]);


colors = colormap(colorcet('L08', 'N', 5));

if effusive_dynamics


else % comparison IHe+
    

        global plot_processed_with_ROI
    plot_processed_with_ROI = false;

    % from 26.8.24, 600 mW probe, 40 bar, 14 K
    %res = plot_processed_VMI(32910, 0, [528.4614  382.9408],0);
   % close all
    
   % load vmi from filenumbers
    res_Iplus_from_He = plot_processed_VMI(I_plus_from_drop_reference_measurement, true,[524.5297  380.8430], true);

    %res_Iplus_gas= plot_processed_VMI(gas_phase_reference_measurement  , 1, [482.9299  392.4866],true);
    
    %data_in_hsnr_gas = load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2I^+.mat');
    %res_Iplus_gas = data_in_hsnr_gas.res_sum;
     res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement  , 1, [482.9299  392.4866], true);


    %res_Iplus_from_He = subtract_processed_data(res_Iplus_from_He, res_Iplus_gas);

   % res_Iplus_He = plot_processed_VMI(I_plus_He_from_drop_reference_measurement, true,[524.5297  380.8430], true);
   
   % average multiple I+He measurements
   for k=1:length(I_plus_He_from_drop_reference_measurements)
       fn = I_plus_He_from_drop_reference_measurements(k);
       if k==1
             res_Iplus_He = plot_processed_VMI(fn, true,[524.5297  380.8430], true);

       else
            res_temp = plot_processed_VMI(fn, true,[524.5297  380.8430], true);
            res_Iplus_He = add_processed_data(res_Iplus_He, res_temp);
       end
   end
   res_Iplus_He = multiply_processed_data(res_Iplus_He,1/length(I_plus_He_from_drop_reference_measurements));

    % alternatively, load vmi from high snr stored mat files

    %data_in_hsnr =  load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He');
    %res_Iplus_He = data_in_hsnr.res_sum;

    global velocity_factor

    % load timescan probe only for comparison
    vf = 8.6178;
    vf_timescan = 5.636;
    global GS_bleach_correction
    GS_bleach_correction = 1;
    res2 = mean_timescan_2d_VMI([296:297], false,[524.5297  380.8430], false);


    axes(tile3);

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

    xlabel(['v / ', char(0197), '/ps']);

    ylabel('signal / arb. units');

    %y = sum(res2.data_probe_only(:,:),2);
    y =  res2.data(:, res2.t>150);
    y = y/max(y);

    %plot(res2.r*vf_timescan,y);

    %plot(res2.r*vf_timescan, y/max(y));
    vf_single =8.6178;
    y= res_Iplus_He.radial_distribution;
    y  = movmean(y, 1);
    

    mass_correction_factor= sqrt((127) /(131));

        hold on
    plot(res_Iplus_gas.r*vf_single/100, res_Iplus_gas.radial_distribution/max( res_Iplus_gas.radial_distribution), 'color',colors(1,:));

    plot(res_Iplus_He.r*vf_single/100*mass_correction_factor, y / max(y ),  ':', 'color',colors(2,:), 'LineWidth',2);
    %plot(res_Iplus_from_He.r*vf_single, res_Iplus_from_He.radial_distribution/max(res_Iplus_from_He.radial_distribution));
    
    v_exp = res_Iplus_He.r*vf_single/100;
    fprintf('contributions at v > 20 A/ps: %.0f \n', trapz(v_exp(v_exp>20),y(v_exp>20) ) /trapz(v_exp,y )*100 );

    fprintf('contributions at v < 20 A/ps: %.0f \n', trapz(v_exp(v_exp<20),y(v_exp<20) ) /trapz(v_exp,y )*100 );



  leg_fig1 = {'I_2:I^+','I_2He_N:I^+He'};
    

 
    figure
    tl2 = tiledlayout(1,1);
    ax2 = nexttile;

    % choose between timescan probe only and high s/n measurement
    % using timescan probe only to compare to because it has better s/n
    % ratio
   % res = res2.data_2d_sum_probe_only;
    res = res_Iplus_He;

    %plot(res.r*vf_timescan, res.radial_distribution);
    VMIN_ANGULAR_DISTR =0;
        b_r = res.r*vf_single> VMIN_ANGULAR_DISTR;% & res.r*vf_single<3000;
            b_phi = (res.phi >0 & res.phi <2*pi);
    y2 = mean(res.image_polar(:, b_r),2);
    y2 = y2/max(y2);

    
    plot(res.phi,y2 );
    hold on
    % fun2 = @(beta, phi) [real(beta(1)*cos(phi).^(beta(2)) + beta(3)); imag(real(beta(1)*cos(phi).^(beta(2))) + beta(3))] ;
    % beta2 = nlinfit(res.phi', [y2; 0*y2], fun2, [1,2,1]);
    % hold on
    % fitresult = fun2(beta2, res.phi')
    % plot(res.phi, fitresult(1:length(y2)));
    

    xlabel('phi / radian');
    ylabel('signal / arb. units');


    % figure(4)
    % [rr, phiphi] = meshgrid(res.r, res.phi);
    % image_polar_select = res.image_polar;
    % image_polar_select(~b_phi,:) = 0;
    % surf(rr, phiphi, image_polar_select);
    % xlabel('phi / radian');
    % ylabel('r / pixel');

end


%% plot simulation result of MD
axes(tile3);

data_neutral = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');
data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');


mass_select = 127+ 4*1;
b_mass_select = round(data_ion.mass_i(:,end)/u)==mass_select;
%b_mass_select = true(size(b_mass_select));


fprintf('mass selection accounts for %.2f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
b_select = b_mass_select & data_ion.b_ion_outside;


vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);


edges_phi = [0:0.05:2*pi];

phi_sim = atan2(vy_total, vx_total) + pi;

plot_samples = phi_sim;
xinterval = [min(edges_phi), max(edges_phi)];
[h, sigma_h, centers_phi, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_phi);
h = movmean(h,15);
% 


axes(ax2)
plot(centers_phi, h/max(h));
legend("I_2:I^+He", "simulation n=0", "simulation n=1", "simulation n=2")
xlabel('phi / radian');




axes(tile2);

%velocity_bins = [-flip(edges_velocity), edges_velocity(2:end)];
velocity_bins = [-35:0.5:35];

velocity_map = zeros(size(velocity_bins,2),size(velocity_bins,2));

for i=1:size(vx_total,1)
    vx_id = get_closest_index(vx_total(i,1), velocity_bins);
    vy_id = get_closest_index(vy_total(i,1), velocity_bins);

    try
        velocity_map (vx_id, vy_id) = velocity_map (vx_id, vy_id) +1;
    end

end


[vxvx, vyvy] = meshgrid(velocity_bins, velocity_bins);

s1 = surf(vxvx,vyvy, velocity_map');
s1.CData = s1.ZData;
s1.ZData = s1.ZData*0;
hold on
% circle1 = draw_circle_on_surface(tile2, [ 1212]/100, [0, 0]);
% circle1.LineWidth = 2;
% circle2 = draw_circle_on_surface(tile2, 7.8, [0, 0]);
% circle2.LineWidth = 2;
% circle3= draw_circle_on_surface(tile2, 18, [0, 0]);
% circle3.LineWidth = 2;

  cmap = colorcet( 'L08');
    cmap = rotate_colors(cmap, 0/180*pi);
   % colormap( clamp(cmap*1.0,0.3));

cb = colorbar;
cb.Label.String = 'signal / counts';

clim([min(s1.CData(:)),max(s1.CData(:))*0.3]);

view(0,90)
xlabel(['v_x / ', char(0197), '/ps']);
ylabel(['v_y / ', char(0197), '/ps']);
pbaspect([1,1,1]);
xlim([-35,35]);
ylim([-35,35]);



axes(tile1)
VX= (res.X- res.image_center_x)*vf_single/100*mass_correction_factor;
VY = (res.Y - res.image_center_y)*vf_single/100*mass_correction_factor;

s2 = surf( VX, VY, res.image);
s2.CData = s2.ZData;
s2.ZData = s2.ZData*0;
hold on
% circle4 = draw_circle_on_surface(tile1, [  1212]/100, [0, 0]);
% circle4.LineWidth = 2;
% circle5 = draw_circle_on_surface(tile1, 7.8, [0, 0]);
% circle5.LineWidth = 2;
% circle6 = draw_circle_on_surface(tile1, 18, [0, 0]);
% circle6.LineWidth = 2;
view(0,90)
xlabel(['v_x / ', char(0197), '/ps']);
ylabel(['v_y / ', char(0197), '/ps']);
pbaspect([1,1,1]);
xlim([-1,1]* max(VX(:)));
ylim([-1,1]* max(VY(:)));
clim([0,max(s2.CData(:))*0.4]);
colormap(colorcet('L08')); %purple yellow)
cb = colorbar
cb.Label.String = 'signal / cps';




letter = add_letter_norm('a','topleft', 20, tile1);
letter.Color = [1,1,1];
letter=add_letter_norm('b','topleft', 20, tile2);
letter.Color = [1,1,1];

letter = add_letter_norm('c', 'topleft', 20, tile3);

%letter_c = add_letter_norm('c','topleft', 18, tile3, 0);
%letter_c.Color = [1,1,1];
%%letter_d = add_letter_norm('d','topleft', 18, tile4, 0);
%letter_d.Color = [1,1,1];


vmin = 4; vmax = 22;

% add covariance data
%result = generate_VMI_covariance_matrices(data_in_hsnr.filenumbers,[vmin, vmax]*100/velocity_factor);

%result = generate_VMI_covariance_matrices(I_plus_He_from_drop_reference_measurement,[vmin, vmax]*100/velocity_factor);

center = autocenter_from_extended_data(I_plus_He_from_drop_reference_measurements)

apply_angular_filter = false;
event_filter = true;
theta_target = pi;
theta_range = 40/180*pi;
result = generate_VMI_covariance_matrices(I_plus_He_from_drop_reference_measurements,[0,600], center, [90,90], apply_angular_filter,  event_filter, theta_target, theta_range);

axes(tile4);


cov_angular = result.cov_angular - diag(diag(result.cov_angular));

surf(result.theta,result.theta,cov_angular , 'EdgeColor','none');
view(0,90);
clim([0,max(cov_angular(:))*0.7])
xlabel('angle / radian');
ylabel('angle / radian');
xlim([-pi, pi]);
ylim([-pi,pi]);
colormap(colorcet('L08'))
cb = colorbar;
cb.Label.String = 'covariance /arb. units';
pbaspect([1,1,1]);
letter = add_letter_norm('d', 'topleft', 20, tile4);
letter.Color = [1,1,1];

axes(tile5)

cov_radial = result.cov_radial - diag(diag(result.cov_radial));


cov_radial = movmean(cov_radial, 2, 1);
cov_radial = movmean(cov_radial, 2, 2);

surfdata =surf(result.r*velocity_factor/100*mass_correction_factor, result.r*velocity_factor/100*mass_correction_factor, cov_radial, 'EdgeColor','none');
surfdata.CData = surfdata.ZData;
surfdata.ZData = surfdata.ZData*0;

view(0,90);

xlabel(['v / ', char(0197), '/ps']);
ylabel(['v / ', char(0197), '/ps']);
colormap(colorcet('L08'))
cb = colorbar;
cb.Label.String = 'covariance /arb. units';
clim([0,max(cov_radial(:))*0.7]);
pbaspect([1,1,1]);

letter = add_letter_norm('e', 'topleft', 20, tile5);
letter.Color = [1,1,1];

xlim([0,max(result.r*velocity_factor/100)]);
ylim([0,max(result.r*velocity_factor/100)]);


axes(tile5)
hold on
%plot([0,35], [vmin, vmin], 'linestyle', '--', 'color', [1,1,1],'HandleVisibility','off');
%plot([0, 35], [vmax, vmax], 'linestyle', '--',  'color',[1,1,1],'HandleVisibility','off');


axes(tile3)
v_radial_corr = result.r*velocity_factor/100*mass_correction_factor;

b_v = v_radial_corr>vmin & v_radial_corr<vmax;

vd_radial_corr =  sum(cov_radial(b_v,:),1)/2; % divide by two because symmetric matrix



%vd_radial_corr = movmean(vd_radial_corr, 5);

scatter(v_radial_corr,vd_radial_corr / max(vd_radial_corr));
leg_fig1{end+1} = "v-cov. trace";
legend(leg_fig1)


f1.Position = [191.4000   84.2000  779.2000  953.6000];

%f1.Position = [191.4000   84.2000  779.2000  953.6000*2/3];

exportgraphics(f1, 'I+He_image.pdf', 'ContentType','vector');


savefig(f1, 'I+He_image');

% compare experimental angular covariance to simulated angular covariance
%%
figure

vr_total = sqrt(data_ion.vx_total.^2 + data_ion.vy_total.^2);
theta_total = atan2(data_ion.vx_total, data_ion.vy_total);

theta1 = theta_total(1:num_molecules) + pi;
theta2 = theta_total(num_molecules+1:end) + pi;

vr1 = vr_total(1:num_molecules); vr2 = vr_total(num_molecules+1:end);
vr1 = vr1(b_select(1:num_molecules)  & b_select(num_molecules+1:end));
vr2 = vr2(b_select(1:num_molecules)  & b_select(num_molecules+1:end));

numbins = 90;

vr = linspace(0,max(result.r),numbins)*velocity_factor/100;
dvr = vr(2)-vr(1);


theta = linspace(0,2*pi, numbins);
dtheta = theta(2)- theta(1);

theta1 = theta1( b_select(1:num_molecules)  & b_select(num_molecules+1:end) );
theta2 = theta2( b_select(1:num_molecules)  & b_select(num_molecules+1:end) );

simulated_angular_covariance = zeros(length(theta), length(theta));
simulated_radial_covariance = zeros(length(vr), length(vr));

for i=1:length(theta1)
   
        simulated_angular_covariance( floor(theta1(i)/dtheta)+1 , floor(theta2(i)/dtheta)+1) = simulated_angular_covariance( floor(theta1(i)/dtheta)+1 , floor(theta2(i)/dtheta)+1) + 1;


end



for i=1:length(vr1)
    simulated_radial_covariance(floor(vr1(i)/dvr)+1, floor(vr2(i)/dvr+1) )= simulated_radial_covariance(floor(vr1(i)/dvr)+1, floor(vr2(i)/dvr+1))  + 1;
end


tl = tiledlayout(1,3);

nexttile
surf(theta-pi, theta-pi, simulated_angular_covariance);
xlabel('\theta / radian'); ylabel('\theta / radian');
colormap(colorcet('L08')); cb = colorbar; cb.Label.String = 'counts';
view(0,90);
%scatter(theta1, theta2, 'x');
pbaspect([1,1,1]);

nexttile

line = plot(theta - pi, sum(simulated_angular_covariance,1));
line.YData = movmean(line.YData,3);
line.YData = line.YData - min(line.YData);
line.YData = line.YData / max(line.YData);

hold on

line = plot(result.theta, sum(cov_angular,1), '--');
line.YData = movmean(line.YData,3);
line.YData = line.YData - min(line.YData);
line.YData = line.YData / max(line.YData);

nexttile


line = plot(vr, sum(simulated_radial_covariance,1));
hold on
line.YData = movmean(line.YData,3);
line.YData = line.YData - min(line.YData);
line.YData = line.YData / max(line.YData);

line = plot(result.r*velocity_factor/100, sum(cov_radial,1), '--');
line.YData = movmean(line.YData,3);
line.YData = line.YData - min(line.YData);
line.YData = line.YData / max(line.YData);

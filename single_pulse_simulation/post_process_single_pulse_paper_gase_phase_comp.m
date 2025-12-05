close all
run physical_constants.m

effusive_dynamics = true;

E_coulomb_scale = 0.8;
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

gas_phase_reference_measurement = 43632; %300 mW


%gas_phase_reference_measurement = 43568;%600 mW

I_plus_from_drop_reference_measurement = 43567; % 600 mW

I_plus_He_from_drop_reference_measurement = 43556;

f1 = figure
tl = tiledlayout(3,2, 'TileSpacing','tight', 'Padding','compact');


tile1 = nexttile(1, [1,1]);
tile2 = nexttile(2, [1,1]);
hold on
tile3 = nexttile(3, [1,2]);
tile4 = nexttile(5, [1,1]);
tile5 = nexttile(6, [1,1]);



if effusive_dynamics


    global plot_processed_with_ROI
    plot_processed_with_ROI = false;

    % get gas phase I+ reference measurement
    res = plot_processed_VMI(gas_phase_reference_measurement  , 1, [482.9299  392.4866], true);

    %data_in_hsnr = load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2I^+.mat');
   % res = data_in_hsnr.res_sum;

    %res = plot_processed_VMI(gas_phase_reference_measurement  , 1, [477.7867  410.8572], 0);

    global velocity_factor

    vf = 8.6178; % for 7.95 kV repeller!
    
    leg_fig1 = {};
    axes(tile3);
 
    plot(res.r*vf/100 ,res.radial_distribution/max(res.radial_distribution))
    hold on
    leg_fig1{end+1} = 'I_2:I^+';
        xlabel(['v / ', char(0197), '/ps']);
    
    ylabel('signal / arb. units');

    % get radial distribution in limited angular range
   % b_phi = (res.phi >0 & res.phi <pi/3) | (res.phi>5/6*2*pi & res.phi<2*pi);
   %  v =res.r*vf;
   % y = mean(res.image_polar(b_phi,:),1);
   % plot(v, y, ':')
   %  leg_fig1{end+1} = 'm. 1: singlepulse, gas, lim.angle';
   % 
   %  fun = @(beta, v) g(beta(1), beta(2), v).*beta(3) + g(beta(4), beta(5), v).*beta(6);
   %  beta = nlinfit(v, y, fun, [814, 1000, 1, 2000, 100, 1]);
   %  plot(v, fun(beta, v));
   %  leg_fig1{end+1} = 'double g fit';
   % 
   %  plot(v, beta(6)*g(beta(4), beta(5), v));
   %  leg_fig1{end+1} = 'high vel. g fit';
   % 
   %  b_r = res.r>250 & res.r<400;

    
    %leg_fig1 = { 'm. 3: p.p. at t<0','simulation'}


    figure(2)
    y2 = mean(res.image_polar(:, :),2);
    y2 = movmean(y2,4);
    %y2 = y2 - min(y2);
    y2 = y2/max(y2);


    plot(res.phi,y2 );
    fun2 = @(beta, phi) [real(beta(1)*cos(phi).^(beta(2))); imag(real(beta(1)*cos(phi).^(beta(2))))] ;
    beta2 = nlinfit(res.phi', [y2; 0*y2], fun2, [1,2]);
    hold on
    fitresult = fun2(beta2, res.phi')
    plot(res.phi, fitresult(1:length(y2)));




else % comparison IHe+
    



end


%% plot
% load('T:\github synchronized\I2HeN_velocity_simulation\\neutral_propagation_checkpoint.mat');
% load('T:\github synchronized\I2HeN_velocity_simulation\ion_propagation_checkpoint.mat');

data_neutral = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');
data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint_gas.mat');



%for mass_select =127 + 4*[0:17]
for mass_select = 127
b_mass_select = round(data_ion.mass_i(:,end)/u)==mass_select;
%b_mass_select = true(size(b_mass_select));

fprintf('mass selection accounts for %.0f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
b_select = b_mass_select;

vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);



% figure
% scatter3(vx_total(:,1), vy_total(:,1),vz_total(:,1));


v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);



edges_velocity = [0:0.05:35];
plot_samples = v_total_proj;


xinterval = [min(edges_velocity), max(edges_velocity)];
[h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
%histogram_data_velocity(1, :) = movmean(h,20);
vd_ion = movmean(h,20);
vd_ion = vd_ion -min(vd_ion);

axes(tile3);
 
plot(centers_velocity, vd_ion/max(vd_ion), '--');
leg_fig1{end+1} = 'simulation';




test = 1;

figure

edges_phi = [0:0.05:2*pi];

phi_sim = atan2(vy_total, vx_total) + pi;

plot_samples = phi_sim;
xinterval = [min(edges_phi), max(edges_phi)];
[h, sigma_h, centers_phi, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_phi);
h = movmean(h,15);

legend('fit measurement', 'simulation');

figure(2);
plot(centers_phi, h/max(h));
legend("AD TS I+He", "fit", "simulation")
xlabel('phi / radian');





figure

axes(tile2);

%velocity_bins = [-flip(edges_velocity), edges_velocity(2:end)];
velocity_bins = [-35:0.4:35];

velocity_map = zeros(size(velocity_bins,2),size(velocity_bins,2));

for i=1:size(vx_total,1)
    vx_id = get_closest_index(vx_total(i,1), velocity_bins);
    vy_id = get_closest_index(vy_total(i,1), velocity_bins);

    try
        velocity_map (vx_id, vy_id) = velocity_map (vx_id, vy_id) +1;
    end

end


[vxvx, vyvy] = meshgrid(velocity_bins, velocity_bins);

s2 = surf(vxvx,vyvy, velocity_map');



view(0,90)
xlabel(['v_x / ', char(0197), '/ps']);
ylabel(['v_y / ', char(0197), '/ps']);
pbaspect([1,1,1]);
xlim([-35,35]);
ylim([-35,35]);
cb = colorbar;
cb.Label.String = 'signal / counts';
clim([0,max(s2.CData(:))*0.4]);
letter = add_letter_norm('b', 'topleft', 20, tile2);
letter.Color = [1,1,1];

axes(tile1);


VX = (res.X- res.image_center_x)*vf/100;
VY =  (res.Y - res.image_center_y)*vf/100;

s1 = surf( VX,VY , res.image);
view(0,90)
xlabel(['v_x / ', char(0197), '/ps']);
ylabel(['v_y / ', char(0197), '/ps']);
pbaspect([1,1,1]);
xlim([-1,1]*max(VX(:)));
ylim([-1,1]*max(VY(:)));


colormap(colorcet('L08')); %purple yellow)
cb = colorbar;
cb.Label.String = 'signal / cps';
clim([0,max(s1.CData(:))*0.1]);

letter = add_letter_norm('a', 'topleft', 20, tile1);
letter.Color = [1,1,1];

letter = add_letter_norm('c', 'topleft', 20, tile3);


vmin = 16; vmax = 24;
%vmin = 0; vmax = 10;
%vmin = 24; vmax = 30;
% add covariance data
%result = generate_VMI_covariance_matrices(data_in_hsnr.filenumbers,[vmin, vmax]*100/velocity_factor);



center = autocenter_from_extended_data(gas_phase_reference_measurement);
result = generate_VMI_covariance_matrices(gas_phase_reference_measurement,[0,600], center, [90,90], false,false, pi, 40/180*pi);


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

s = surf(result.r*velocity_factor/100, result.r*velocity_factor/100, cov_radial, 'EdgeColor','none');
view(0,90);
s.CData = s.ZData; s.ZData = s.ZData*0;

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



axes(tile3)

v_radial_corr = result.r*velocity_factor/100;
b_v = v_radial_corr>vmin & v_radial_corr<vmax;

vd_radial_corr =  sum(cov_radial(b_v,:),1)/2; % divide by two because symmetric matrix




%vd_radial_corr = movmean(vd_radial_corr, 5);

scatter(result.r*velocity_factor/100,vd_radial_corr / max(vd_radial_corr));
leg_fig1{end+1} = "v-cov. trace";
l = legend(tile3, leg_fig1)

%axes(tile5)
%hold on
%plot([0,35], [vmin, vmin], 'linestyle', '--', 'color', [1,1,1],'HandleVisibility','off');
%plot([0, 35], [vmax, vmax], 'linestyle', '--',  'color',[1,1,1],'HandleVisibility','off');


end


f1.Position = [191.4000  0.2000  779.2000  953.6000];


exportgraphics(f1, 'T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\gas_phase_results\single_pulse_gas_phase.pdf', 'ContentType','vector');


savefig(f1, 'T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\gas_phase_results\single_pulse_gas_phase');
% if 

if E_coulomb_scale<1
f = figure

newAx = copyobj( tile3, f);

% Set the new axes as the current axes in the new figure
set(f, 'CurrentAxes', newAx);

% Optionally, adjust the position of the new axes
newAx.OuterPosition = [0 0 1 1]; % Full figure size

newAx.Children(2).String = '';
f.Position(3:4) = figsize;


set(newAx, 'fontsize', 15);

legend(newAx, leg_fig1);

%exportgraphics(f, 'T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\gas_phase_results\single_pulse_gas_phase_Ec_80percent.pdf', 'ContentType','vector');

end

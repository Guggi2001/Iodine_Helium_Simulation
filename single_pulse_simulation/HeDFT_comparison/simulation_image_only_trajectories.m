addpath 'T:\github synchronized\Iodine_Helium_Simulation'
% if running without prior simulation, do this:
%run inputfiles_dft_comparison\single_pulse_droplet_distribution.m 

run physical_constants.m

global abel_inv_post
abel_inv_post = true;

if exist('diff','var')
    clear diff
end

%% load marti data

global R0_GS
switch R0_GS
    case 9
        addpath('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\HeDFT_comparison\9Angström\')
        [t, v] = importfile_v2('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\HeDFT_comparison\9Angström\data_vabs2.csv');
        v = abs(v);
        
        [tR, R] = importfile_R1_R2('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\HeDFT_comparison\9Angström\R1-R2.csv');


    case 18
        addpath('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\HeDFT_comparison\18Angström');
        importpath = 'T:\Cloud\MATLAB iodine\MartiPiNotes\impurities_dynamics\impurities_dynamics\';

        R1 = importfile_marti([importpath, 'rimp.I_1']);
        R2 = importfile_marti([importpath,'rimp.I_2']);
        tR = R1.t;
        R = sqrt(sum((R1.vec - R2.vec).^2,2));

        V1 = importfile_marti([importpath,'vimp.I_1']);
        V2 = importfile_marti([importpath,'vimp.I_2']);

        t = V2.t;

        v = sqrt(sum(V2.vec.^2,2));
end

%% post processing
close all
global single_pulse
if single_pulse
  data_neutral = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\neutral_propagation_checkpoint');
 %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);

end

data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint_hedft.mat');
%%
% run 'T:\github synchronized\Iodine_Helium_Simulation\inputfiles_dft_comparison\single_pulse_N2000.m'
% for i = 1:3
%     global E_coulomb_scale
%     E_coulomb_scale = 1.1-i*0.1;
%     sigma_dependent_on_v = false;
%     vmi_sim_3d_neutral_propa_HeDFT_mimic;
%     sigma_dependent_on_v = true;
%     vmi_sim_3d_ion_propa;
%     data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint_hedft.mat');
%     disp('-----------------------')
%     fprintf('Quantify the mismatch across the whole overlap interval for a Coulomb Energy of %3f \n', E_coulomb_scale)
%     
%     Nmol = size(data_ion.x_ci,1)/2;
%     dx = data_ion.x_ci(1:Nmol,:) - data_ion.x_ci(1+Nmol:end,:);
%     dy = data_ion.y_ci(1:Nmol,:) - data_ion.y_ci(1+Nmol:end,:);
%     dz = data_ion.z_ci(1:Nmol,:) - data_ion.z_ci(1+Nmol:end,:);
%     dR_mean = mean(sqrt(dx.^2+dy.^2+dz.^2), 1);
%     
%     tmax = min(max(tR), max(data_ion.time_i));
%     mask_md = data_ion.time_i <= tmax;
%     maskR   = tR <= tmax;
%     
%     t_md   = data_ion.time_i(mask_md).';
%     d_md   = dR_mean(mask_md).';
%     tR_use = tR(maskR).';
%     R_use  = R(maskR).';
%     
%     dR_mean_on_tR = interp1(t_md, d_md, tR_use, 'linear');
%     
%     good  = isfinite(dR_mean_on_tR) & isfinite(R_use);
%     rmse  = sqrt(mean((dR_mean_on_tR(good) - R_use(good)).^2));
%     ratio = mean(dR_mean_on_tR(good) ./ R_use(good));
%     
%     fprintf('RMSE = %.3f Å, mean ratio = %.4f\n', rmse, ratio);
% end


%%


disp('-----------------------')
disp('Quantify the mismatch across the whole overlap interval')

Nmol = size(data_ion.x_ci,1)/2;
dx = data_ion.x_ci(1:Nmol,:) - data_ion.x_ci(1+Nmol:end,:);
dy = data_ion.y_ci(1:Nmol,:) - data_ion.y_ci(1+Nmol:end,:);
dz = data_ion.z_ci(1:Nmol,:) - data_ion.z_ci(1+Nmol:end,:);
dR_mean = mean(sqrt(dx.^2+dy.^2+dz.^2), 1);

tmax = min(max(tR), max(data_ion.time_i));
mask_md = data_ion.time_i <= tmax;
maskR   = tR <= tmax;

t_md   = data_ion.time_i(mask_md).';
d_md   = dR_mean(mask_md).';
tR_use = tR(maskR).';
R_use  = R(maskR).';

dR_mean_on_tR = interp1(t_md, d_md, tR_use, 'linear');

good  = isfinite(dR_mean_on_tR) & isfinite(R_use);
rmse  = sqrt(mean((dR_mean_on_tR(good) - R_use(good)).^2));
ratio = mean(dR_mean_on_tR(good) ./ R_use(good));

fprintf('RMSE = %.3f Å, mean ratio = %.4f\n', rmse, ratio);



%%
figure



plot(tR, R);

hold on

num_molecules = size(data_ion.x_ci,1)/2;

for particle_id = 1:num_molecules/10

pos1 = [data_ion.x_ci(particle_id,:); data_ion.y_ci(particle_id,:); data_ion.z_ci(particle_id,:)];



pos2 = [data_ion.x_ci(particle_id+num_molecules,:); data_ion.y_ci(particle_id+num_molecules,:); data_ion.z_ci(particle_id+num_molecules,:)];

dR_vector = pos1 - pos2;

dR = sqrt(sum( dR_vector.^2, 1));

if particle_id>1
plot(data_ion.time_i, dR , 'color', [1,0.2,0.6, 0.1], 'HandleVisibility','off');
else
plot(data_ion.time_i, dR , 'color', [1,0.2,0.6,0.1], 'HandleVisibility','on');
end

hold on
end

xlim([0,6]);
ylim([8, 40]);

ylabel('R_1 - R_2 / A');
xlabel('t / ps');


legend({ 'HeDFT', 'MD trajectories'});

f1 = figure
tl = tiledlayout(2,2, 'TileSpacing','tight','Padding','compact');



tile3 = nexttile(1,[1,2]);

tile4 = nexttile(3,[1,2]);

axes(tile3);

plot(t, v, 'LineStyle',':', 'LineWidth',2);
hold on

v_total = sqrt(data_ion.vx_ci(:,:).^2+ data_ion.vy_ci(:,:).^2+ data_ion.vz_ci(:,:).^2);


for particle_id = 1:round(num_molecules/15):num_molecules
vel1 = [data_ion.vx_ci(particle_id,:); data_ion.vy_ci(particle_id,:); data_ion.vz_ci(particle_id,:)];


vel2 = [data_ion.vx_ci(particle_id+num_molecules,:); data_ion.vy_ci(particle_id+num_molecules,:); data_ion.vz_ci(particle_id+num_molecules,:)];

v1 = sqrt(sum( vel1.^2, 1));
v2 = sqrt(sum( vel2.^2, 1));

if particle_id>1
plot(data_ion.time_i,v1,  'color',[0.2,0.2,0.6, 0.1], 'HandleVisibility','off');
hold on
plot(data_ion.time_i,v2,  'color',[0.2,0.2,0.6, 0.1], 'HandleVisibility','off');
else
plot(data_ion.time_i,v1,  'color',[0.2,0.2,0.6, 0.1], 'HandleVisibility','on');
hold on
plot(data_ion.time_i,v2,  'color',[0.2,0.2,0.6, 0.1], 'HandleVisibility','off');

end

end

plot(data_ion.time_i, mean(v_total,1), '--','color', [0,0,0]);

xlim([0,12]);
%ylim([0,6]);
xlabel('t / ps');
ylabel('v / A/ps');



legend({'HeDFT','MD velocity', 'mean MD velocity'} );


return

axes(tile4)
% plot I+ gas and I+He hsnr data
%data_in_hsnr =  load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He');
%res_Iplus_He = data_in_hsnr.res_sum;

I_plus_He_from_drop_reference_measurement = 43630; % 300 mW 17.10.24
I_plus_He_from_drop_reference_measurements = [45668, 45662, 45667, 45686]; % I+He measurements, low doping, 300 mW, 3.12.24

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


res_Iplus_He.image(res_Iplus_He.image<0) = 0;
res_Iplus_He.image = movmean(res_Iplus_He.image, 3,1);
res_Iplus_He.image = movmean(res_Iplus_He.image, 3,2);

res_Iplus_He = abel_invert_processed_VMI(res_Iplus_He);

%res_Iplus_He = plot_processed_VMI(I_plus_He_from_drop_reference_measurement, true,[524.5297  380.8430], true);
    
    
%data_in_hsnr_gas = load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2I^+.mat');
%res_Iplus_gas = data_in_hsnr_gas.res_sum;
gas_phase_reference_measurement = 43632; %300 mW
res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement  , 1, [482.9299  392.4866], true);
res_Iplus_gas = abel_invert_processed_VMI(res_Iplus_gas);

colors = colormap(colorcet('L08', 'N', 5));


xlabel(['v / ', char(0197), '/ps']);

ylabel('signal / arb. units');


vf_single =8.6178;

y= res_Iplus_He.radial_distribution;
y  = movmean(y, 1);

hold on
b_v = res_Iplus_gas.r*vf_single/100>4;

plot(res_Iplus_gas.r*vf_single/100, res_Iplus_gas.radial_distribution/max( res_Iplus_gas.radial_distribution(b_v)), 'color',colors(1,:));

mass_correction_factor = sqrt(127/131);

plot(res_Iplus_He.r*vf_single/100*mass_correction_factor, y / max(y ),  ':', 'color',colors(2,:));





% % plot simulation result for I+He and I+He2
% data_neutral = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');
% data_ion = load('T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
% 
% 
% 
% for mass_select = 127+ 4*[1,2]
% b_mass_select = round(data_ion.mass_i(:,end)/u)==mass_select;
% 
% 
% fprintf('mass selection accounts for %.2f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
% b_select = b_mass_select & data_ion.b_ion_outside;
% 
% 
% vx_total = data_ion.vx_total(b_select, :);
% vy_total = data_ion.vy_total(b_select, :);
% vz_total = data_ion.vz_total(b_select, :);
% 
% v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);
% v_total = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2 + vz_total(:,1).^2);
% 
% edges_velocity = [0:0.04:26];
% 
% if abel_inv_post
%     plot_samples = v_total;
% else
% plot_samples = v_total_proj;
% end
% 
% xinterval = [min(edges_velocity), max(edges_velocity)];
% [h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
% vd_ion = movmean(h,15);
% 
% xlim([0,28]);
% 
% ylim([0,1.1]);
% 
% hold on
% 
% 
% vd_ion =vd_ion - min(vd_ion);
% 
% if mass_select==131
% plot(centers_velocity, vd_ion/max(vd_ion), '--',  'color',colors(3,:));
% else
%     plot(centers_velocity, vd_ion/max(vd_ion), '-.',  'color',colors(4,:));
% end
% 
% 
% end
% 
% 
% 
% %f1.Position(3:4) = figsize2.*[1,2];
% 
%   leg = {'I_2:I^+','I_2He_N:I^+He',  'simulation I^+He',  'simulation I^+He_2'};
%   legend(leg);
% 
% f1.Position = [    1.0306    0.3746    0.8216    0.5776]*1E3;
% 
% add_letter_norm('a', 'topleft', 20, tile3)
% add_letter_norm('b', 'topleft', 20, tile4)
% 
% 
% 
% 
% 
% if abel_inv_post
% exportgraphics(gcf, 'simulation_image_inv.pdf', 'ContentType', 'vector')
% else
%     exportgraphics(gcf, 'simulation_image.pdf', 'ContentType', 'vector')
% end
% 
% 
% 
% if E_coulomb_scale <1
% f = figure
% 
% newAx = copyobj( tile4, f);
% 
% % Set the new axes as the current axes in the new figure
% set(f, 'CurrentAxes', newAx);
% 
% % Optionally, adjust the position of the new axes
% newAx.OuterPosition = [0 0 1 1]; % Full figure size
% 
% newAx.Children(1).String = '';
% f.Position(3:4) = figsize;
% 
% set(newAx, 'fontsize', 15);
% 
% legend(newAx, leg);
% 
%   exportgraphics(gcf, 'simulation_image_Ec_80percent.pdf', 'ContentType', 'vector')
% end



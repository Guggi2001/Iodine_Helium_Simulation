

global abel_inv_post; 
abel_inv_post = true;
reeval = 0;

gas_phase_reference_measurement = 43632; %300 mW

I_plus_He_from_drop_reference_measurements = [45668, 45662, 45667]; % I+He measurements, low doping, 300 mW, 3.12.24

f1  = figure
ax1 = axes;
hold on
xlabel(['v / ', char(0197), '/ps']);
ylabel('signal / arb. units');
f1.Position = [432   591   922   335];

figure 
ax2  = axes;
hold on
xlabel('angle / °');
ylabel('signal / arb. units');
set(ax1, 'fontsize', 18);
set(ax2, 'fontsize', 18);

%% gas phase 
res = plot_processed_VMI(gas_phase_reference_measurement  , reeval, [482.9299  392.4866], true);
res = abel_invert_processed_VMI(res);

VX = (res.X- res.image_center_x)*vf/100;
VY =  (res.Y - res.image_center_y)*vf/100;

figure
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


%letter = add_letter_norm('a', 'topleft', 20, gca);
%letter.Color = [1,1,1];

clim([0,20])

plot(ax1, res.r*velocity_factor/100, res.radial_distribution/max(res.radial_distribution));
plot(ax2, res.phi*180/pi, res.angular_distribution/max(res.angular_distribution));

%% I+He


   % average multiple I+He measurements
   for k=1:length(I_plus_He_from_drop_reference_measurements)
       fn = I_plus_He_from_drop_reference_measurements(k);
       if k==1
             res_Iplus_He = plot_processed_VMI(fn, reeval,[524.5297  380.8430], true);

       else
            res_temp = plot_processed_VMI(fn, reeval,[524.5297  380.8430], true);
            res_Iplus_He = add_processed_data(res_Iplus_He, res_temp);
       end
   end
   res = multiply_processed_data(res_Iplus_He,1/length(I_plus_He_from_drop_reference_measurements));

res = abel_invert_processed_VMI(res);

VX = (res.X- res.image_center_x)*vf/100;
VY =  (res.Y - res.image_center_y)*vf/100;

figure
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


%letter = add_letter_norm('a', 'topleft', 20, gca);
%letter.Color = [1,1,1];

clim([0,2])

v_smooth_size = 200;
movmean_size = ceil(v_smooth_size/((res.r(2) - res.r(1))*velocity_factor));
plot(ax1, res.r*velocity_factor/100,movmean( res.radial_distribution/max(res.radial_distribution), movmean_size));
plot(ax2, res.phi*180/pi, res.angular_distribution / max(res.angular_distribution));


legend(ax1,'I^+ from isolated I_2', ...
    'I^+He from inside droplet');



%% I+He simulation
figure
if effusive_dynamics
data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint_gas.mat');
else
%data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
if single_initial_position
    data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint_hedft.mat');
else
    data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
end
end


mass_select = 127+ 4*1;
b_mass_select = round(data_ion.mass_i(:,end)/u)==mass_select;
%b_mass_select = round(data_ion.mass_i(:,end)/u)>140 &  round(data_ion.mass_i(:,end)/u)>127 <  127+4*10;

%b_mass_select = true(size(b_mass_select)); 



fprintf('mass selection accounts for %.2f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
b_select = b_mass_select & data_ion.b_ion_outside;


vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);


%velocity_bins = [-flip(edges_velocity), edges_velocity(2:end)];
velocity_bins = [-35:0.5:35];

velocity_map = zeros(size(velocity_bins,2),size(velocity_bins,2));

v_total = sqrt(vx_total.^2 + vy_total.^2 + vz_total.^2);

ux = vx_total./sqrt(vx_total.^2 + vy_total.^2);
uy = vy_total./sqrt(vy_total.^2 + vx_total.^2);

for i=1:size(vx_total,1)

   %     vx_id = get_closest_index(vx_total(i,1), velocity_bins);
  %  vy_id = get_closest_index(vy_total(i,1), velocity_bins);

    vx_id = get_closest_index(ux(i,1)*v_total(i,1), velocity_bins);
    vy_id = get_closest_index(uy(i,1)*v_total(i,1), velocity_bins);

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
colormap(colorcet('L08'));

cb = colorbar;
cb.Label.String = 'signal / counts';

%clim([min(s1.CData(:)),max(s1.CData(:))*0.3]);

view(0,90)
xlabel(['v_x / ', char(0197), '/ps']);
ylabel(['v_y / ', char(0197), '/ps']);
pbaspect([1,1,1]);
xlim([-35,35]);
ylim([-35,35]);


%% 

global velocity_factor


filenumbers = gas_phase_reference_measurement;
mass_correction_factor = 1;

filenumbers = I_plus_He_from_drop_reference_measurements;
    mass_correction_factor= sqrt((127) /(131));

center = autocenter_from_extended_data(filenumbers(1) );
center = [505, 350];

apply_angular_filter = false;
event_filter = false;
theta_target = pi;
theta_range = 30/180*pi;
result = generate_VMI_covariance_matrices(filenumbers ,[0,600], center, [150,150], apply_angular_filter,  event_filter, theta_target, theta_range);



figure
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



figure

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


xlim([0,max(result.r*velocity_factor/100)]);
ylim([0,max(result.r*velocity_factor/100)]);


figure
 plot(result.r*velocity_factor/100*mass_correction_factor, sum(cov_radial,1))
xlabel(['v / ', char(0197), '/ps']);
ylabel('total covariance signal');

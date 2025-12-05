

run physical_constants.m

%% post processing
close all
global single_pulse
if single_pulse
  data_neutral = load('single_pulse_simulation\neutral_propagation_checkpoint');
 %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);



data_ion = load('single_pulse_simulation\ion_propagation_checkpoint.mat');
else

data_neutral = load('neutral_propagation_checkpoint');
 %  [vx_total, vy_total, vz_total] = add_cei_vel_to_vel_3d(x_components, y_components,z_components,  vx_components, vy_components , vz_components, mass);


if effusive_dynamics

data_ion = load('ion_propagation_checkpoint_gas.mat');
else
data_ion = load('ion_propagation_checkpoint.mat');
end


end

E_initial = data_neutral.E_initial;
E_d = data_neutral.E_initial;
m = data_neutral.m;
mass = data_neutral.mass;



load('post_process_checkpoint')


figure
tiledlayout(1,3);
nexttile
num_molecules = size(data_ion.x_ci,1)/2;

for particle_id = 1:ceil(num_molecules/100):num_molecules

pos1 = [data_ion.x_ci(particle_id,:); data_ion.y_ci(particle_id,:); data_ion.z_ci(particle_id,:)];

pos2 = [data_ion.x_ci(particle_id+num_molecules,:); data_ion.y_ci(particle_id+num_molecules,:); data_ion.z_ci(particle_id+num_molecules,:)];

dR_vector = pos1 - pos2;

dR = sqrt(sum( dR_vector.^2, 1));

plot(data_ion.time_i, dR );
hold on
end




ylabel('R_1 - R_2 / A');
xlabel('t / ps');
title('i+ interatomic distance')


nexttile
num_molecules = size(data_ion.x_ci,1)/2;

for particle_id = 1:ceil(num_molecules/100):num_molecules

pos = [data_ion.x_ci(particle_id,:); data_ion.y_ci(particle_id,:); data_ion.z_ci(particle_id,:)];



dR_vector = pos;

dR = sqrt(sum( dR_vector.^2, 1));

plot(data_ion.time_i, dR );
hold on
end



ylabel('R / A');
xlabel('t / ps');
title('i+ distance from origin')


nexttile

for particle_id = 1:ceil(num_molecules/100):num_molecules
vel1 = [data_ion.vx_ci(particle_id,:); data_ion.vy_ci(particle_id,:); data_ion.vz_ci(particle_id,:)];


vel2 = [data_ion.vx_ci(particle_id+num_molecules,:); data_ion.vy_ci(particle_id+num_molecules,:); data_ion.vz_ci(particle_id+num_molecules,:)];

v1 = sqrt(sum( vel1.^2, 1));
v2 = sqrt(sum( vel2.^2, 1));

line = plot(data_ion.time_i,v1);
hold on
plot(data_ion.time_i,-v2, 'color', line.Color);

%ylim([0,6]);
xlabel('t / ps');
ylabel('v / A/ps');


end


title('i+ velocities')

f = gcf;
f.Position =[0.3458    0.5010    1.2656    0.4200]*1E3;

test = 1;


figure
tl = tiledlayout(1,3);

nexttile
scatter(L_droplet(:,end), sqrt( vx_components(:,end).^2 + vy_components(:,end).^2) )
xlabel('distance traveled in droplet');
ylabel('final velocity');

starting_radius = r_start;
nexttile
scatter(starting_radius(~b_inside),  sqrt( vx_components(~b_inside,end).^2 + vy_components(~b_inside,end).^2));
xlabel('starting distance from center');
ylabel('final velocity');

nexttile

scatter(starting_radius,  v_total);
xlabel('starting distance from center');
ylabel('total ion velocity');

% figure
% scatter(m_vectors(~b_inside,end)/u, v_total(~b_inside));



f = gcf;
f.Position = 1.0e+03 *[1.0600    0.4720    0.7810    0.4445];

try
figure

colormap(viridis);

tl = tiledlayout(3,4, 'Padding','tight', 'TileSpacing','tight')

nexttile

data = movmean(movmean(histogram_data_interatomic_distance,1,1),5,2);
%data = movmean(movmean(histogram_data_radius,1,1),5,2);


[tt,xx] = meshgrid(centers_radius, t_unique);
surf(xx, tt, data, 'EdgeColor','none')
shading interp
view(0,90)
%clim([0,max(data(end,:))*0.8]);



title('radius')
xlabel('t / ps');
ylabel('r / Angström');
set(gca,'XScale', 'log');
xlim([min(t_unique), max(t_unique)])
colorbar


nexttile
[tt,xx] = meshgrid(centers_velocity, t_unique);

data = movmean(movmean(histogram_data_velocity,1,1),5,2);
surf(xx, tt, data, 'EdgeColor','none')
shading interp
view(0,90)
try
clim([0,max(data(end,:))*0.8]);
catch ME
clim([0,max(data(:))*0.8]);
end
colorbar


try
%clim([min(histogram_data_velocity(end,:)), max(histogram_data_velocity(end,:))])
catch  ME
end

 title('velocity')
 xlabel('t / ps');

 ylabel('v / A/ps');
 
set(gca,'XScale', 'log');
xlim([min(t_unique), max(t_unique)])

nexttile

%%

t_local = t_unique';
v_local = centers_velocity;
smooth_fun = @(mu, sig, x) 1/sqrt(2 *pi*sig^2)*exp( -(x - mu).^2 / (2*sig^2));

v_distr_local = histogram_data_total_velocity;
v_distr_smoothed = 0*v_distr_local;

normalization = v_distr_smoothed;

sigma_t_smooth = 0.1; %0.1 ps
sigma_v_smooth = 0.1; % 0.1 A/ps

% time smoothing
for k=1:size(v_distr_local,2)
    for l=1:size(v_distr_local,1)
        t_center = t_local(l);
        v_distr_smoothed(l,k) = trapz(t_local, smooth_fun(t_center, sigma_t_smooth, t_local).*v_distr_local(:,k)   )./trapz(t_local, smooth_fun(t_center,sigma_t_smooth, t_local));
        normalization(l,k) = trapz(t_local, smooth_fun(t_center, sigma_t_smooth, t_local));
    end
end

v_distr_local = v_distr_smoothed;

for l=1:size(v_distr_local,1)
    for k=1:size(v_distr_local,2)
        v_center = v_local(k);
        v_distr_smoothed(l,k) = trapz(v_local, smooth_fun(v_center, sigma_v_smooth , v_local).*v_distr_local(l,:))./trapz(v_local, smooth_fun( v_center, sigma_v_smooth , v_local));

    end
end




%data = movmean(movmean(histogram_data_total_velocity,5,1),5,2);
data = v_distr_smoothed;
[xx,tt] = meshgrid(centers_velocity, t_unique);
surf(tt, xx, data, 'EdgeColor','none')
colorbar

try
clim([0,max(data(end,:))*0.8]);
catch ME
clim([0,max(data(:))*0.8]);
end


colorbar

shading interp
view(0,90)
try
%caxis([min(histogram_data_total_velocity(end,:)), max(histogram_data_total_velocity(end,:))])
catch ME
    warning(ME.message)
end
title('velocity total')
xlabel('t / ps');
ylabel('v / A/ps');
hold on
r_landau = 2.665 + 2*0.50*time; % in angström
v_ref =sqrt((E_initial - E_d)./m(1))/100;

r_ref = 2.665 + 2*v_ref*time;
% check if this is correct!
plot3(time, sqrt( (coulomb_energy(r_ref))/mass(1)  + (v_ref*100)^2)/100, zeros(size(time)) + max(histogram_data_total_velocity(:)) ,'color', [1,1,1], 'LineWidth',0.8);
hold on
plot3(time, sqrt( (2*coulomb_energy(r_ref))/mass(1)  + (v_ref*100)^2)/100, zeros(size(time)) + max(histogram_data_total_velocity(:)) ,'color', [1,1,1],'LineWidth',0.8, 'LineStyle','--');
plot3(time, sqrt( (4*coulomb_energy(r_ref))/mass(1)  + (v_ref*100)^2)/100, zeros(size(time)) + max(histogram_data_total_velocity(:)) ,'color', [1,1,1],'LineWidth',0.8,'LineStyle',':');
set(gca,'XScale', 'log');
xlim([min(t_unique), max(t_unique)])
ylim([0,22])




nexttile
plot(time, number_inside);
hold on
plot(time, number_outside);

xlabel('t / ps')
ylabel('number of atoms');
legend('inside', 'outside');
set(gca,'XScale', 'log');
xlim([min(t_unique), max(t_unique)])
nexttile
plot(time, mean(E_kin,1))

fun = @(beta, t) beta(1)*exp(-t/beta(2)) + beta(3);
try
beta = nlinfit(time, mean(E_kin,1), fun, [50, 1,0.1]);
catch ME
    beta = [nan,nan,nan];
end

xlabel('t  / ps');
ylabel('mean kinetic energy neutral');

hold on
plot(time, fun(beta,time));
fprintf('relative energy loss: %.2f', 1/beta(2)*100)

Deltat = 1;
set(gca,'XScale', 'log');
%legend(sprintf('relative energy loss: %.2f', 1/beta(2)*100));
legend([sprintf('relative energy loss per ps: %.2f ', (exp(-Deltat/beta(2)) - 1)/Deltat*100), '%']);
xlim([min(t_unique), max(t_unique)])



nexttile


plot(time, sum(b_ion_outside,1)/size(b_ion_outside,1));
hold on
%total_velocity = sqrt(vx_total.^2 + vy_total.^2 + vz_total.^2);
%plot(time, sum(b_ion_outside & total_velocity<7.50,1)/size(b_ion_outside,1));

xlabel('t / ps');
ylabel('fraction escaped ions');

set(gca,'XScale', 'log');
%legend('all ion velocities', 'v_{total}<750 m/s');
xlim([min(t_unique), max(t_unique)])

f = gcf;
f.Position  =  [75 391 1143 567.5000] ;








nexttile
[tt,xx] = meshgrid(centers_energy, t_unique);
data = histogram_data_potential_energy;
data = movmean(data,5,2);

colorbar


surf(xx, tt, data, 'EdgeColor','none')
shading interp
view(0,90)
ylabel('E_{pot} / meV');
xlabel('t / ps');
colorbar

try
%clim([0,max(data(:))*0.1]);
catch ME
end


nexttile
[tt,xx] = meshgrid(centers_mass, t_unique);

data =histogram_data_mass;
data = movmean(data,2,2);
surf(xx, tt, data, 'EdgeColor','none')
shading interp
view(0,90)
ylabel('mass');
xlabel('t / ps');
try
clim([0,max(data(end,:))*0.8]);
catch ME
clim([0,max(data(:))*0.8]);
end
colorbar


catch ME

end



function custom_tick_settings(ax)

xticks = ax.XAxis.TickValues;
ax.XAxis.TickLabels = {};
for xtk = xticks
    if xtk<1
    ax.XAxis.TickLabels{end+1} = sprintf('%.1f', xtk);
    else
    ax.XAxis.TickLabels{end+1} = sprintf('%.0f', xtk);
    end

end

end

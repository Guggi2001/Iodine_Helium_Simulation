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


if effusive_dynamics
    path = 'T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures\timescan figures\effusive example\';
    data_in = load([path, 'ts_figure_evaluated_data.mat']);
    res2 = data_in.ts_inverted;
       vf_timescan =  5.636*sqrt(5/3.4);


else % comparison IHe+
    

    % load timescan probe only for comparison
    vf = velocity_factor;
    vf_timescan = 5.636;
    

    path = 'T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures\timescan figures\IplusHe_7.95kV\';
    data_in = load([path, 'ts_figure_evaluated_data.mat']);
    res2 = data_in.ts_inverted;
    %res2 = data_in.ts_res_no_GSB_correction;
    vf_timescan =  8.9838*sqrt(127/131);


end


%% plot

figure
[tt,xx] = meshgrid(res2.t, res2.r*vf_timescan/100);
data_in = load('T:\github synchronized\I2HeN_velocity_simulation\post_process_checkpoint.mat')
s = surf(tt,xx,res2.data, 'HandleVisibility','off');
s.CData = s.ZData;
s.ZData = s.ZData*0;

view(0,90);
xlim([min(data_in.time), max(data_in.time)]);
ylim([min(res2.r*vf_timescan), max(res2.r*vf_timescan)]/100)
set(gca,'XScale','log')

hold on
% plot mean total velocity of simulation
support = data_in.centers_velocity;
distribution = data_in.histogram_data_total_velocity;
distribution = distribution - min(distribution(:,:),[],2);
plot(data_in.time, trapz(support,distribution.*support,2)./ trapz(support,distribution,2), 'color', [1,1,1])
set(gca,'XScale','log')

xlabel('t / ps');
ylabel(['v / ',char(0197), '/ps'])
colormap(viridis);



hold on
% plot mean total velocity of simulation
support = res2.r*vf_timescan/100;
distribution =res2.data';
distribution(distribution<0) = 0;

%distribution = distribution - min(distribution(:,:),[],2);
plot(res2.t, trapz(support,distribution.*support,2)./ trapz(support,distribution,2), '--','color', [1,1,1])
set(gca,'XScale','log')
cb = colorbar;
cb.Label.String = 'signal / counts';


ids = [];

% find maximum point of distribution
% for i=1:length(res2.t)
%     y = movmean(res2.data(:,i), 50);
%         peak_ids =  findLocalMaxima(y,31,90,60);
% 
% 
%     % plot(res2.r, res2.data(:,i));
%     % hold on
%     % plot(res2.r, y);
%     % 
%     % 
%     % hold on
%     % scatter(res2.r(peak_ids), res2.data(peak_ids,i))
%     % test = 1;
% 
%     ids = [ids, peak_ids(end-2)];
% 
% 
% end

%plot(tile1, res2.t, res2.r(ids)*vf_timescan/100, ':',  'Color',[1,1,1]);
f= gcf;
f.Position(3:4) = figsize2;


ylim([0,22]);

xlim([0.2,200]);

ax1 = gca;

figure
[tt,xx] = meshgrid(data_in.time, data_in.centers_velocity);
s = surf(tt,xx,data_in.histogram_data_total_velocity'-min(data_in.histogram_data_total_velocity',[],1), 'HandleVisibility','off' );
s.CData = s.ZData;
s.ZData = s.ZData*0;
ax2 = gca;

view(0,90);
xlim([min(data_in.time), max(data_in.time)]);
%ylim([min(res2.r*vf_timescan), max(res2.r*vf_timescan)]/100)
ylim([0,2200]);

set(gca,'XScale','log')
hold on
% % plot mean total velocity
% plot(data_in.time, trapz(data_in.centers_velocity,data_in.histogram_data_total_velocity.*data_in.centers_velocity,2)./ trapz(data_in.centers_velocity,data_in.histogram_data_total_velocity,2) )
% set(gca,'XScale','log')

% plot mean total velocity of simulation
support = data_in.centers_velocity;
distribution = data_in.histogram_data_total_velocity;

distribution = distribution - min(distribution(:,:),[],2);
plot(data_in.time, trapz(support,distribution.*support,2)./ trapz(support,distribution,2) , 'color', [1,1,1])
colormap(viridis);

xlabel('t / ps');
ylabel(['v / ',char(0197), '/ps'])

ylim([0,22]);

xlim([0.2,200]);

cb = colorbar;
cb.Label.String = 'signal / counts';


f2= gcf;
f2.Position(3:4) = figsize2;




axes(ax1);

l = legend( 'mean(v) simulation','mean(v) experiment');
l.Box = 'off';
l.TextColor = [1,1,1];

let = add_letter_norm('a','topleft', 18, ax1); let.Color = [1,1,1];
let = add_letter_norm('b','topleft', 18, ax2); let.Color = [1,1,1];



exportgraphics(f, 'T:\github synchronized\I2HeN_velocity_simulation\pumpprobe_simulation_images\exp.pdf');

exportgraphics(f2, 'T:\github synchronized\I2HeN_velocity_simulation\pumpprobe_simulation_images\sim.pdf');


figure
support = data_in.centers_velocity;
distribution = data_in.histogram_data_velocity;
distribution = distribution - min(distribution(:,:),[],2);

plot(data_in.time, trapz(support,distribution.*support,2)./ trapz(support,distribution,2) )

set(gca,'XScale','log')
title('mean velocity');



figure
support = data_in.centers_radius;

distribution = data_in.histogram_data_radius;
distribution = distribution - min(distribution(:,:),[],2);


plot(data_in.time, trapz(support,distribution.*support,2)./ trapz(support,distribution,2) )

title('mean interatomic distance');
set(gca,'XScale','log')


test = 1;



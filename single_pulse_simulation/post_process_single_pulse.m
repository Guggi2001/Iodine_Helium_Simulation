


if effusive_dynamics
    global plot_processed_with_ROI
    plot_processed_with_ROI = false;
    res = plot_processed_VMI(28735, 1, [489.4394  377.3878], 1);
    close all
    global velocity_factor

    vf = velocity_factor;
    vf =  6.4149*0.9675;
    figure(1)
    plot(res.r*vf ,res.radial_distribution/max(res.radial_distribution))
    hold on
    b_phi = (res.phi >0 & res.phi <pi/3) | (res.phi>5/6*2*pi & res.phi<2*pi);
    v =res.r*vf;
    y = mean(res.image_polar(b_phi,:),1);
    plot(v, y, ':')

    fun = @(beta, v) g(beta(1), beta(2), v).*beta(3) + g(beta(4), beta(5), v).*beta(6);
    beta = nlinfit(v, y, fun, [814, 1000, 1, 2000, 100, 1]);
    plot(v, fun(beta, v));
    plot(v, beta(6)*g(beta(4), beta(5), v));

    b_r = res.r>250 & res.r<300;

    
    leg_fig1 = {'m. 1: singlepulse','fit', 'm. 2: probe only from ts', 'm. 3: p.p. at t<0','simulation'}


    figure(2)
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
    res = plot_processed_VMI(32910, 0, [528.4614  382.9408],0);
    close all
    global velocity_factor

    %vf = velocity_factor;
    %vf =  6.4149*0.9675;
    global GS_bleach_correction
    GS_bleach_correction = 1;
    res2 = mean_timescan_2d_VMI([296:297], false,[524.5297  380.8430], false)

    vf = velocity_factor;

    figure(1)
    %plot(res.r*vf ,res.radial_distribution/max(res.radial_distribution))
    hold on
    b_phi = (res.phi >0 & res.phi <2*pi);
    v =res.r*vf;
    y = mean(res.image_polar(b_phi,:),1);
    plot(v, y/max(y), ':')
    
    % fit 1 gaussian
    fun = @(beta, v) g(beta(1), beta(2), v).*beta(3);% + g(beta(4), beta(5), v).*beta(6);
    beta = nlinfit(v, y, fun, [814, 1000, 1]);
    plot(v, fun(beta, v)/max(fun(beta, v)));
    %plot(v, beta(3)*g(beta(1), beta(2), v));

    b_r = res.r>150 & res.r<280;
    xlabel('v / m/s');
    y = sum(res2.data_probe_only(:,:),2);
    plot(res2.r*velocity_factor, y/max(y));
    
    y = res2.data(:,1);
    hold on
    plot(res2.r*velocity_factor, y/max(y))
  leg_fig1 = {'m. 1: singlepulse','fit', 'm. 2: probe only from ts', 'm. 3: p.p. at t<0','simulation'}
    
  %      figure(1)
 %   hold on
%plot(v, g(beta(1), beta(2), v));
 
    figure(2)
    y2 = mean(res.image_polar(:, b_r),2);
    y2 = movmean(y2,4);
    %y2 = y2 - min(y2);
    y2 = y2/max(y2);

    
    plot(res.phi,y2 );
    fun2 = @(beta, phi) [real(beta(1)*cos(phi).^(beta(2)) + beta(3)); imag(real(beta(1)*cos(phi).^(beta(2))) + beta(3))] ;
    beta2 = nlinfit(res.phi', [y2; 0*y2], fun2, [1,2,1]);
    hold on
    fitresult = fun2(beta2, res.phi')
    plot(res.phi, fitresult(1:length(y2)));

    xlabel('phi / radian');


    figure(4)
    [rr, phiphi] = meshgrid(res.r, res.phi);
    image_polar_select = res.image_polar;
    image_polar_select(~b_phi,:) = 0;
    surf(rr, phiphi, image_polar_select);
    xlabel('phi / radian');
    ylabel('r / pixel');

end

% load('T:\github synchronized\I2HeN_velocity_simulation\\neutral_propagation_checkpoint.mat');
% load('T:\github synchronized\I2HeN_velocity_simulation\ion_propagation_checkpoint.mat');

load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');
load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');
figure
scatter3(vx_total(:,1), vy_total(:,1),vz_total(:,1));

figure(3)
v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);


figure(1)
edges_velocity = [0:0.05:22];
plot_samples = v_total_proj;


xinterval = [min(edges_velocity), max(edges_velocity)];
[h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
%histogram_data_velocity(1, :) = movmean(h,20);
vd_ion = movmean(h,20);
vd_ion = vd_ion -min(vd_ion);
plot(centers_velocity*100, vd_ion/max(vd_ion));


legend(leg_fig1)
figure

edges_phi = [0:0.1:2*pi];

phi_sim = atan2(vy_total, vx_total) + pi;

plot_samples = phi_sim;
xinterval = [min(edges_phi), max(edges_phi)];
[h, sigma_h, centers_phi, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_phi);
histogram_data_phi(1, :) = h;

legend('fit measurement', 'simulation');

figure(2);
plot(centers_phi, histogram_data_phi/max(histogram_data_phi));
legend("measurement", "theory", "simulation")
xlabel('phi / radian');



figure
%velocity_bins = [-flip(edges_velocity), edges_velocity(2:end)];
velocity_bins = [-22:0.5:22];

velocity_map = zeros(size(velocity_bins,2),size(velocity_bins,2));

for i=1:size(vx_total,1)
    vx_id = get_closest_index(vx_total(i,1), velocity_bins);
    vy_id = get_closest_index(vy_total(i,1), velocity_bins);

    try
        velocity_map (vx_id, vy_id) = velocity_map (vx_id, vy_id) +1;
    end

end


[vxvx, vyvy] = meshgrid(velocity_bins, velocity_bins);

surf(vxvx,vyvy, velocity_map);


view(90,90)
xlabel('vx / A/ps');
ylabel('vy / A/ps');





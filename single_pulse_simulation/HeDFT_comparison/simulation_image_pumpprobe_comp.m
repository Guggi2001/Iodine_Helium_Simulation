    % load timescan probe only for comparison

    abel_inv_post = true;

    vf = velocity_factor;
    vf_timescan = 5.636;
    

    path = 'T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures\timescan figures\IplusHe_7.95kV\';
    data_in = load([path, 'ts_figure_evaluated_data.mat']);
    res2 = data_in.ts_inverted;
    vf_timescan =  8.9838*sqrt(127/131);


% plot mean total velocity of simulation
support = res2.r*vf_timescan/100;
distribution =res2.data';
distribution(distribution<0) = 0;
b_t = res2.t>100;

%distribution = distribution - min(distribution(:,:),[],2);

distr = trapz(res2.t(b_t),distribution(b_t,:),1);
plot(support, distr/max(distr), 'color', [0,0,0])





% plot simulation result for I+He and I+He2
data_neutral = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\neutral_propagation_checkpoint.mat');
data_ion = load('T:\github synchronized\I2HeN_velocity_simulation\single_pulse_simulation\ion_propagation_checkpoint.mat');

droplet_radii = data_neutral.droplet_radii;
b_invalid_spawn = data_neutral.b_invalid_spawn;

R_ion_final = sqrt(data_ion.x_ci_final.^2 + data_ion.y_ci_final.^2 + data_ion.z_ci_final.^2);
b_ion_outside = R_ion_final > [droplet_radii];

for mass_select = 127+ 4*[0,1]
b_mass_select = (round(data_ion.mass_i(:,end)/u)==mass_select )  & (~data_neutral.b_invalid_spawn) & data_ion.b_ion_outside;



fprintf('mass selection accounts for %.2f percent of total ions\n', 100*sum(b_mass_select)/numel(b_mass_select));
b_select = b_mass_select & data_ion.b_ion_outside;

%b_select = true(size(b_select));

vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);

v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);
v_total = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2 + vz_total(:,1).^2);

edges_velocity = [0:0.04:26];

if abel_inv_post
    plot_samples = v_total;
else
plot_samples = v_total_proj;
end

xinterval = [min(edges_velocity), max(edges_velocity)];
[h, sigma_h, centers_velocity, ~, ~, ~] = bayes_hist(plot_samples, xinterval, true, 'r', edges_velocity);
vd_ion = movmean(h,15);

xlim([0,28]);

ylim([0,1.1]);

hold on


vd_ion =vd_ion - min(vd_ion);

vd_ion = movmean(vd_ion,15);

if mass_select==131
plot(centers_velocity, vd_ion/max(vd_ion), '--',  'color',[0.1, 0.7, 0.1]);
else
    plot(centers_velocity, vd_ion/max(vd_ion), '-.',  'color',[0.0 0.3, 0.8]);
end

xlabel(['v / ', char(0197), '/ps'])
ylabel('counts / arb. units');







end

  leg = {'I_2He_N:I^+He',  'simulation I^+', 'simulation I^+He'};
  legend(leg);

  f = gcf;
f.Position(3:4) = figsize;

  exportgraphics(f,[ sprintf('comparison_pumpprobe_%.2fK', T_particles), '.pdf'])





  %% add a plot: initial depth and orientation vs final velocity
figure
  test = 1;

num_molecules = data_neutral.num_molecules;
atom1 = 1:num_molecules;
atom2 = num_molecules+1:2*num_molecules;

  X0_1 = [data_neutral.x_components(atom1,1), data_neutral.y_components(atom1,1), data_neutral.z_components(atom1, 1)];

  X0_2 =  [data_neutral.x_components(atom2,1), data_neutral.y_components(atom2,1), data_neutral.z_components(atom2, 1)];

  X0 = [X0_1; X0_2];

r0_atoms = sqrt(sum( ([X0_1; X0_2]).^2,2));

  com = (X0_1 + X0_2)/2;



r0_test = sqrt(sum(com.^2, 2));

r0 = data_neutral.r0;

depth = droplet_radii - r0_atoms;

plot(depth(~data_neutral.b_invalid_spawn))


vec_from_com = X0 - [com;com];

scalar_prod = sum( [com;com].*vec_from_com, 2)./ (sqrt(sum( [com;com].^2, 2)).*sqrt(sum( vec_from_com.^2, 2)));

angle = acos(scalar_prod); 

%scatter([r0;r0], angle)

%hold on
%scatter([data_neutral.r0; data_neutral.r0], acos(data_neutral.angle))


% something about the angle calculation later is off. 
% plot directly with initial angle
mass_select = 131;
b_mass_select = (round(data_ion.mass_i(:,end)/u)>=mass_select );

b_select = ~b_invalid_spawn & b_ion_outside; % & b_mass_select;

vx_total = data_ion.vx_total(b_select, :);
vy_total = data_ion.vy_total(b_select, :);
vz_total = data_ion.vz_total(b_select, :);

v_total_proj = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2);
v_total = sqrt(vx_total(:,1).^2 + vy_total(:,1).^2 + vz_total(:,1).^2);


%scatter(depth(b_select), acos(data_neutral.angle(b_select)))

scatter3(depth(b_select), acos(data_neutral.angle(b_select)), v_total, 'markeredgecolor', [0, 0.8, 0.8]);

scatter(acos(data_neutral.angle(b_select)), v_total)

figure
histogram(acos(data_neutral.angle(b_select)));


figure
histogram(v_total);

figure

% form probability distribution for depth
 [h, sigma_h, centers, binwidth, barplot, errorplot] = bayes_hist(depth,[min(depth),max(depth)], true);



p_depth = h*binwidth;


% form probability distribution for depth, given escape
 [h, sigma_h, centers2, binwidth, barplot, errorplot] = bayes_hist(depth(b_select),[min(depth),max(depth)], true);
p_depth_given_escaped = h*binwidth;


plot(centers, p_depth);
hold on
plot(centers2, p_depth_given_escaped);

p_escape = sum(b_select)/numel(b_select);


% bayes
p_escape_given_depth = p_depth_given_escaped;


figure
histogram2(v_total, depth(b_select), 'NumBins',[30,30]);
xlabel(['escaped ion velocity / ', char(0197), '/ps']);
ylabel(['depth / ', char(0197)]);
zlabel('# of ions')

exportgraphics(gcf, 'histogram2d_depth_ion_velocity.pdf');

figure
histogram2(acos(data_neutral.angle(b_select))*180/pi, depth(b_select), 'NumBins',[30,30]);
xlabel('angle to surface normal');
ylabel(['depth / ', char(0197)]);
zlabel('# of ions')

exportgraphics(gcf, 'histogram2d_depth_angle.pdf');







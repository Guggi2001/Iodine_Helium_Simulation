function vmi_sim_visualize_distributions(input_crate)

desolve_struct(input_crate);

%% velocity
figure
    v_total = sqrt(vx_components(:, :).^2 + vy_components(:,:).^2 + vz_components(:,:).^2);


data = v_total;
edges_in= linspace(0,25,60);

leg = {};
for i=40:50:length(time)
%for i=1:5:40
    [n, edges, bin] = histcounts(data(:,i), edges_in);

    binwidth = edges(2)-edges(1);
    centers = edges(1:end-1) + binwidth/2;
    plot(centers, n);
    hold on

    leg{end+1} = sprintf('%.2f ps', time(i));
end

recolor_lines
legend(leg);
xlabel(['velocity / ', char(0197),'/ps']);
ylabel('number of ions');


%% depth
figure
r = sqrt(x_components(:, :).^2 + y_components(:,:).^2 + z_components(:,:).^2);
depth = mean(droplet_radii) - r;

data = depth;
edges_in= linspace(-40,mean(droplet_radii), 60);

leg = {};
%for i=40:50:length(time)
for i=1:20:160
    [n, edges, bin] = histcounts(data(:,i), edges_in);

    binwidth = edges(2)-edges(1);
    centers = edges(1:end-1) + binwidth/2;
    plot(centers, n);
    hold on

    leg{end+1} = sprintf('%.2f ps', time(i));
end

recolor_lines
legend(leg);
xlabel(['depth / ', char(0197)]);
ylabel('number of ions');


%% energy
figure


data = E_kin;
edges_in= linspace(0,2.2, 60);

E_mean_array = [];
leg = {};
%for i=40:50:length(time)
time_indices = 20:10:ceil(length(time)/2);

for i=time_indices
    [n, edges, bin] = histcounts(data(:,i), edges_in);

    E_mean = mean(data(:,i));
    E_mean_array = [E_mean_array, E_mean];

    binwidth = edges(2)-edges(1);
    centers = edges(1:end-1) + binwidth/2;
    plot(centers, n);
    hold on

    leg{end+1} = sprintf('%.2f ps', time(i));
end

recolor_lines
legend(leg);
xlabel(['E_{kin} / eV']);
ylabel('number of ions');

figure
plot(time(time_indices), E_mean_array);
%weights = ones(size(time_indices))';

%weights = (time(time_indices)>0)'+eps;

%fun = @(beta, t) beta(1)*exp(-t/beta(2)) + beta(3);
%[Ypred, prediction_error, beta, parameter_error, ErrorModelInfo] = nonlinear_fit_wrapper_2(time(time_indices), E_mean_array,...
%  weights, 0.05, fun, [1,1,1]);
% hold on
%plot(time(time_indices), fun(beta, time(time_indices)));


ylabel('<E_{kin}> / eV');
xlabel('t / ps');



figure
local_time =time(time_indices);

dt = local_time(2:end) - local_time(1:end-1);
difftime = local_time(1:end-1) + dt/2;
plot(difftime, diff(E_mean_array)./dt);
ylabel('energy loss per ps');
xlabel('time / ps');


end

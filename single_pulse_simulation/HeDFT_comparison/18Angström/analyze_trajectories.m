%%
R1 = importfile_marti('rimp.I_1');
R2 = importfile_marti('rimp.I_2');

X = R1.vec - R2.vec;

%plot(R1.t, sqrt(sum(X.^2,2)));

V1 = importfile_marti('vimp.I_1');


%% calculate coulomb force at each time step
addpath('T:\github synchronized\I2HeN_velocity_simulation\additional_functions\');
global E_coulomb_scale
E_coulomb_scale = 1;
global single_charge_ionization_allowed
single_charge_ionization_allowed  = false;


dr = sqrt(sum(X.^2,2));
dr_unit_vec = X./dr;

Epot = ion_interaction_potential(dr, 1,1);
h= 0.0001;

F = (Epot - ion_interaction_potential(dr+h, 1, 1))/h;% force in eV/Angström


u = 1.66053907e-27; %kg
mass = 127*u;


a = F./(mass/u); % in eV/(Angström u)

a = a*9648.53322; % in Angström / ps^2


% size of a should be size of 2*size(q2)
a1_calc_coulomb = a.*dr_unit_vec;
a2_calc_coulomb = -a.*dr_unit_vec;



%%
A1 = importfile_marti('aimp.I_1');
A2 = importfile_marti('aimp.I_2');

V1_calc = cumtrapz(A1.t, A1.vec,1);
R1_calc = cumtrapz(A1.t, V1_calc,1) + R1.vec(1,:);


%%
plot(A1.t, sqrt(sum(R1_calc.^2,2)));
hold on
plot(A1.t, sqrt(sum(R1.vec.^2,2)), ':');




%% plot initial positions in space, and acceleration vectors at these positions
addpath('T:\github synchronized\I2HeN_velocity_simulation\plot_utility\');
figure
scatter3( R1.vec(1,1),  R1.vec(1,2),  R1.vec(1,3));
hold on

%plot_vector(R1.vec(1,:)', A1.vec(1,:)', 1E0);
%plot_vector(R1.vec(1,:)', a1_calc_coulomb(1,:)', 1E0);
plot_vector(R1.vec(1,:)', A1.vec(1,:)'- a1_calc_coulomb(1,:)', 1E0);

scatter3( R2.vec(1,1),  R2.vec(1,2),  R2.vec(1,3));
hold on

%plot_vector(R2.vec(1,:)', A2.vec(1,:)', 1E0);
plot_vector(R2.vec(1,:)', A2.vec(1,:)' - a2_calc_coulomb(1,:)', 1E0);

xlim([-50,50]);
ylim([-50,50]);
zlim([-50,50]);
pbaspect([1,1,1]);



%% quiver plot of decelleration (does not really work)
a1_decel = A1.vec(:,:)- a1_calc_coulomb(:,:);
a2_decel = A2.vec(:,:)- a2_calc_coulomb(:,:);


figure
plotstep = 50;
id = 1:plotstep:size(R1.vec,1);

quiver3(R1.vec(id ,1), R1.vec(id ,2), R1.vec(id ,3), a1_decel(id ,1), a1_decel(id ,2), a1_decel(id ,3), 100);

xlim([-50,50]);
ylim([-50,50]);
zlim([-50,50]);
pbaspect([1,1,1]);



%% plot coulomb and decel, vs full acceleration
plot(A1.t, sqrt(sum(A1.vec.^2,2)));
hold on
plot(A1.t, sqrt(sum(a1_calc_coulomb.^2,2)));

plot(A1.t, sqrt(sum(a1_decel.^2,2)));

%% plot coulomb and decel, vs full acceleration (only z component
component_id = 3;

figure
plot(A1.t, A1.vec(:,component_id));
hold on
plot(A1.t, a1_calc_coulomb(:,component_id));

plot(A1.t, a1_decel(:,component_id));

plot(A1.t,a1_decel_average);

a1_decel_average = movmean(a1_decel(:,component_id),50);


legend('a_z', 'a_z (C)', 'a_z (fric)', 'movmean(a_z (fric))')


%% plot avarage deceleration as function of position and velocity

plot3(R1.vec(:,3), V1.vec(1:end-1,3), a1_decel_average)

xlabel('z position'); ylabel('v_z'); zlabel('a_{decel}, average');


%% do it again, but with norm


figure
plot(A1.t, sqrt(sum(A1.vec.^2,2)));
hold on
plot(A1.t, sqrt(sum(a1_calc_coulomb.^2,2)));

plot(A1.t, sqrt(sum( a1_decel.^2, 2)));


a1_decel_average = movmean(sqrt(sum(a1_decel.^2,2)),50000);

plot(A1.t,a1_decel_average);



legend('a', 'a (C)', 'a (fric)', 'movmean(a (fric))')


% plot avarage deceleration as function of position and velocity

plot3(sqrt(sum(R1.vec.^2,2)), sqrt(sum( V1.vec(1:end-1,:).^2,2)), a1_decel_average);

xlabel('z position'); ylabel('v'); zlabel('a_{decel}, average');



% %% plot coulomb and decel, vs full acceleration (only z component
% figure
% plot(A2.t, A2.vec(:,component_id));
% hold on
% plot(A2.t, a2_calc_coulomb(:,component_id));
% 
% plot(A2.t, a2_decel(:,component_id));
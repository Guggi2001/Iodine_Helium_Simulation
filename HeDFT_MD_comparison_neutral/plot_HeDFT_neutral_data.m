% fig 14 of notes from 11.11.24
tiledlayout(2,1);

data_in = load('vz_neutral.mat');
nexttile
plot(data_in.t, data_in.vz);

xlabel('t / ps');
ylabel('vz / Angstroem/ps');

vz_interp = @(tnew) interp1(data_in.t,data_in.vz, tnew);




hold on
t=[0:0.001:10];



scatter(t,vz_interp(t));

data_in = load('vz_longer');
plot(data_in.t,abs(data_in.vz));
hline(0.4);

R0 = 2.666;

% fig 15 of notes from 11.11.24
data_in = load('R1_R2_neutral.mat');
nexttile
plot(data_in.t, data_in.R/data_in.R(end) * (data_in.R(end)-R0) + R0 );
xlabel('t / ps');
ylabel('R / Ansgroem')


hold on



%R_interp = cumtrapz(t,vz_interp*2 );

R_interp =  @(tnew) interp1(data_in.t,data_in.R/data_in.R(end) * (data_in.R(end)-R0) + R0, tnew);
scatter(t, R_interp(t));


save('custom_start_interpolating_functions', 'vz_interp', 'R_interp');





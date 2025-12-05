figure


data_neutral = load('T:\github synchronized\I2HeN_velocity_simulation\neutral_propagation_checkpoint.mat');


HeDFT_data = load('T:\github synchronized\I2HeN_velocity_simulation\HeDFT_MD_comparison_neutral\vz_neutral.mat');

v = sqrt(data_neutral.vx_components.^2 + data_neutral.vy_components.^2 + data_neutral.vz_components.^2);


plot(data_neutral.time, v(1:5:end,:)', 'color', [0.1,0.5,0.7,0.1], 'HandleVisibility','off');
hold on
plot(data_neutral.time, mean(v,1), 'color', [0.1,0.1,0.5,1]);

b_t = HeDFT_data.t<=3; 
line = plot(HeDFT_data.t(b_t), HeDFT_data.vz(b_t),'--')


HeDFT_data_longer = load('T:\github synchronized\I2HeN_velocity_simulation\HeDFT_MD_comparison_neutral\vz_longer.mat');

b_t = HeDFT_data_longer.t>3;
plot(HeDFT_data_longer.t(b_t), abs(HeDFT_data_longer.vz(b_t)), '--', 'Color',line.Color, 'HandleVisibility', 'off')


legend('MD simulation', 'TDHeDFT');
set(gca,'XScale', 'log');
 xlabel('t / ps')
ylabel(['|v| / ', char(0197), '/ps'])
grid on
xlim([ 0.0450  199.0000]);

exportgraphics(gcf, 'neutral_dynamics_comparison.pdf');

%%
figure
 

vz = HeDFT_data.vz;
 

 

r_from_integral = 2*cumtrapz(HeDFT_data.t,  HeDFT_data.vz)+ 2.666;
 
E_kin =127*u*(HeDFT_data.vz*100).^2/2/eV *2 ;
 

total_energy_available =  E_kin + e_charge^2./(4*pi*epsilon_0*r_from_integral*1E-10)  /eV;
 

 


 

plot(HeDFT_data.t, total_energy_available/2);
 hold on
plot(HeDFT_data.t, E_kin/2*10, ':')
a = gca;
 


 

hold on
 

plot([min(HeDFT_data.t), max(HeDFT_data.t)], [0.3,0.3], '--')
 


 



ylabel('E / eV');

 


yyaxis right
 

a = gca;

 


line = plot(HeDFT_data.t, r_from_integral, '-o' );
 
line.Color = [0.5,0.1,0.4];

 a.YAxis(2).Color = line.Color;

ylabel(['R(t) / ', char(0197)]);
 
grid on
xlabel('t / ps');
 
legend('E_{C}/2 + E_{kin}^{neutral}' ,'E_{kin}^{neutral} \times 10', 'E_{solv} I^+', 'R(t)');
exportgraphics(gcf, 'TDHeDFT_max_energy_gain.pdf');

%%

r = sqrt((data_neutral.x_components(1:data_neutral.num_molecules,:) - data_neutral.x_components(data_neutral.num_molecules+1:end,:) ).^2 + ...
    (data_neutral.y_components(1:data_neutral.num_molecules,:) - data_neutral.y_components(data_neutral.num_molecules+1:end,:)).^2 + ...
    (data_neutral.z_components(1:data_neutral.num_molecules,:) - data_neutral.z_components(data_neutral.num_molecules+1:end,:)).^2);

HeDFT_data = load('T:\github synchronized\I2HeN_velocity_simulation\HeDFT_MD_comparison_neutral\R1_R2_neutral.mat');

figure



plot(data_neutral.time, r'/2,'Color',[0.1,0.1,0.5,0.01]);
hold on
plot(HeDFT_data.t, (HeDFT_data.R+2.666)/2, 'color',[1,0,0]);


